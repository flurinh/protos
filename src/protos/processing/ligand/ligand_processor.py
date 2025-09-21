"""LigandProcessor: registry-aligned manager for ligand data."""

from __future__ import annotations

import json
import math
import pickle
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.io.ingest.utils.ligand_utils import (
    calculate_molecular_properties,
    sanitize_smiles_filename,
    smiles_to_inchi,
    validate_smiles,
    create_sdf_from_smiles,
)
from protos.processing.property import PropertyProcessor

try:  # Optional RDKit dependency
    from rdkit import Chem  # pragma: no cover - optional
except ImportError:  # pragma: no cover - optional
    Chem = None


class LigandProcessor(BaseProcessor):
    """Manage ligand entities stored as SMILES strings or extracted structures."""

    processor_type = "ligand"

    def __init__(
        self,
        name: str = "ligand_processor",
        *,
        structure_processor: Optional[Any] = None,
    ) -> None:
        super().__init__(name=name)

        self._smiles_cache: Dict[str, Dict[str, Any]] = {}
        self._structure_processor = structure_processor

        self.smiles_dir = Path(self.get_subdirectory_path("sdf_dir"))
        self.cache_dir = Path(self.get_subdirectory_path("cache_dir"))
        self.datasets_dir = Path(self.get_subdirectory_path("datasets_dir"))
        self.smiles_dir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.datasets_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Dependency helpers
    # ------------------------------------------------------------------
    @property
    def dependencies_available(self) -> bool:
        return Chem is not None

    # ------------------------------------------------------------------
    # BaseProcessor abstract implementations
    # ------------------------------------------------------------------
    def load_entity(self, name: str) -> Optional[Dict[str, Any]]:
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info is None:
            return None

        metadata = dict(entity_info.metadata or {})
        kind = metadata.get("kind", "smiles")

        payload: Dict[str, Any] = {
            "name": entity_info.original_id,
            "kind": kind,
            "metadata": metadata,
        }

        if kind == "structure" and entity_info.file_path:
            file_path = Path(self.paths.data_root) / entity_info.file_path
            if file_path.exists():
                with open(file_path, "rb") as handle:
                    payload["structure"] = pickle.load(handle)
        else:
            payload["smiles"] = entity_info.original_id
        return payload

    def save_entity(
        self,
        name: str,
        data: Union[str, Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        if isinstance(data, str):
            ligand_record = {"smiles": data}
        else:
            ligand_record = dict(data)

        ligand_record.setdefault("input_name", name)
        stored_name = self._persist_smiles_entity(ligand_record, metadata or {})
        return stored_name

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def register_smiles_dataset(
        self,
        smiles_map: Dict[str, str],
        *,
        dataset_name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Tuple[str, List[str]]:
        entity_names: List[str] = []
        for label, smiles in smiles_map.items():
            record = {
                "smiles": smiles,
                "label": label,
            }
            entity_name = self._persist_smiles_entity(record, metadata or {})
            entity_names.append(entity_name)

        dataset_id = dataset_name or self._sanitize_filename("smiles_dataset")
        dataset_meta = {
            "kind": "smiles",
            "entity_count": len(entity_names),
        }
        if metadata:
            dataset_meta.update(metadata)

        if self.dataset_manager.dataset_exists(dataset_id):
            self.dataset_manager.delete_dataset(dataset_id)

        self.create_dataset(dataset_id, entity_names, dataset_meta)
        return dataset_id, entity_names

    def list_ligands(self, dataset: Optional[str] = None) -> List[str]:
        if dataset:
            return self.dataset_manager.get_dataset_entities(dataset)
        return self.entity_registry.list_entities(self.processor_type)

    # ------------------------------------------------------------------
    # Structure integration
    # ------------------------------------------------------------------
    def extract_ligands_from_structure(
        self,
        structure_id: str,
        *,
        structure_processor: Optional[Any] = None,
        dataset_name: Optional[str] = None,
        min_atoms: int = 3,
    ) -> Tuple[str, List[str]]:
        struct_proc = structure_processor or self._structure_processor
        if struct_proc is None:
            from protos.processing.structure import StructureProcessor

            struct_proc = StructureProcessor()

        structure_df = struct_proc.load_entity(structure_id)
        if structure_df is None:
            raise ValueError(f"Structure '{structure_id}' could not be loaded")

        hetero = structure_df[structure_df["group"].str.upper() == "HETATM"]
        if hetero.empty:
            return dataset_name or f"{structure_id}__ligands", []

        grouping_cols = [
            col
            for col in ("auth_asym_id", "auth_seq_id", "auth_comp_id", "pdbx_PDB_ins_code")
            if col in hetero.columns
        ] or ["auth_seq_id"]

        entity_names: List[str] = []
        for group_values, ligand_df in hetero.groupby(grouping_cols, dropna=False):
            if len(ligand_df) < min_atoms:
                continue

            descriptor = self._format_ligand_descriptor(group_values)
            entity_name = self._persist_structure_ligand(
                structure_id,
                descriptor,
                ligand_df,
            )
            entity_names.append(entity_name)

        if not entity_names:
            return dataset_name or f"{structure_id}__ligands", []

        dataset_id = dataset_name or self._sanitize_filename(f"{structure_id}__ligands")
        metadata = {
            "structure_id": structure_id,
            "entity_count": len(entity_names),
            "kind": "structure",
        }

        if self.dataset_manager.dataset_exists(dataset_id):
            self.dataset_manager.delete_dataset(dataset_id)

        self.create_dataset(dataset_id, entity_names, metadata)
        return dataset_id, entity_names

    def compute_interactions(
        self,
        structure_id: str,
        *,
        structure_processor: Optional[Any] = None,
        ligands: Optional[Iterable[str]] = None,
        distance_cutoff: float = 4.0,
    ) -> pd.DataFrame:
        struct_proc = structure_processor or self._structure_processor
        if struct_proc is None:
            from protos.processing.structure import StructureProcessor

            struct_proc = StructureProcessor()

        structure_df = struct_proc.load_entity(structure_id)
        if structure_df is None:
            raise ValueError(f"Structure '{structure_id}' could not be loaded")

        coord_cols = self._coordinate_columns(structure_df)

        protein = structure_df[structure_df["group"].str.upper() == "ATOM"].copy()
        lig_atoms = structure_df[structure_df["group"].str.upper() == "HETATM"].copy()
        if lig_atoms.empty or protein.empty:
            return pd.DataFrame()

        protein_coords = protein[coord_cols].to_numpy(float)
        chain_col, resseq_col, resname_col = self._structure_metadata_columns(protein)
        protein_meta = protein[[chain_col, resseq_col, resname_col]].to_numpy()

        records: List[Dict[str, Any]] = []
        grouping_cols = [
            col
            for col in ("auth_asym_id", "auth_seq_id", "auth_comp_id")
            if col in lig_atoms.columns
        ]
        if not grouping_cols:
            grouping_cols = [col for col in ("chain_id", "label_seq_id", "res_name") if col in lig_atoms.columns]
        ligand_groups = lig_atoms.groupby(grouping_cols, dropna=False)

        ligand_entity_lookup = {name: name for name in (ligands or [])}
        for group_values, ligand_df in ligand_groups:
            descriptor = self._format_ligand_descriptor(group_values)
            if ligands is not None and descriptor not in ligand_entity_lookup:
                continue

            ligand_coords = ligand_df[coord_cols].to_numpy(float)
            if ligand_coords.size == 0:
                continue

            distances = self._pairwise_distances(ligand_coords, protein_coords)
            contact_indices = np.where(distances <= distance_cutoff)
            if contact_indices[0].size == 0:
                continue

            for lig_idx, prot_idx in zip(*contact_indices):
                residue_info = protein_meta[prot_idx]
                records.append(
                    {
                        "structure_id": structure_id,
                        "ligand_descriptor": descriptor,
                        "protein_chain": residue_info[0],
                        "protein_resseq": residue_info[1],
                        "protein_resname": residue_info[2],
                        "distance": float(distances[lig_idx, prot_idx]),
                        "interaction_type": "contact",
                    }
                )

        df = pd.DataFrame(records)
        if df.empty:
            return df

        df["ligand_entity"] = df["ligand_descriptor"].map(
            lambda d: self._sanitize_filename(d) if d in ligand_entity_lookup else None
        )

        scope_entries: List[List[Dict[str, str]]] = []
        for _, row in df.iterrows():
            scope = [{"format": "structure", "name": structure_id}]
            if row["ligand_entity"]:
                scope.append({"format": self.processor_type, "name": row["ligand_entity"]})
            scope_entries.append(scope)

        df["scope"] = scope_entries
        return df

    def record_interactions(
        self,
        table_name: str,
        interactions: pd.DataFrame,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        allow_create: bool = True,
    ) -> pd.DataFrame:
        prop_proc = PropertyProcessor()
        return prop_proc.record_properties(
            table_name,
            interactions,
            metadata=metadata,
            allow_create=allow_create,
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _persist_smiles_entity(
        self,
        record: Dict[str, Any],
        metadata: Dict[str, Any],
    ) -> str:
        smiles = record.get("smiles")
        if not smiles:
            raise ValueError("Ligand record missing 'smiles'")

        is_valid, canonical_smiles = validate_smiles(smiles)
        if not is_valid or not canonical_smiles:
            raise ValueError(f"Invalid SMILES string: {smiles}")

        entity_name = canonical_smiles
        safe_name = sanitize_smiles_filename(canonical_smiles)
        file_path = self.smiles_dir / f"{safe_name}.sdf"

        if self.dependencies_available:
            sdf_content = create_sdf_from_smiles(canonical_smiles, record.get("properties"))
            if sdf_content:
                file_path.parent.mkdir(parents=True, exist_ok=True)
                with open(file_path, "w") as handle:
                    handle.write(sdf_content)
                relative_path = str(file_path.relative_to(self.paths.data_root))
            else:
                relative_path = ""
        else:
            relative_path = ""

        properties = record.get("properties")
        if properties is None:
            properties = calculate_molecular_properties(canonical_smiles)

        standard_meta = {
            "kind": "smiles",
            "input_name": record.get("input_name", canonical_smiles),
            "properties": properties or {},
        }

        inchi_data = smiles_to_inchi(canonical_smiles)
        if inchi_data:
            standard_meta.update(inchi_data)

        standard_meta.update(metadata)

        self.entity_registry.register_entity(
            name=entity_name,
            format_type=self.processor_type,
            file_path=relative_path,
            metadata=standard_meta,
        )

        self._smiles_cache[entity_name] = {
            "smiles": canonical_smiles,
            "metadata": standard_meta,
        }
        return entity_name

    def _persist_structure_ligand(
        self,
        structure_id: str,
        descriptor: str,
        ligand_df: pd.DataFrame,
    ) -> str:
        entity_name = self._sanitize_filename(f"{structure_id}__{descriptor}")
        file_path = self.cache_dir / f"{entity_name}.pkl"
        file_path.parent.mkdir(parents=True, exist_ok=True)
        ligand_df.to_pickle(file_path)

        metadata = {
            "kind": "structure",
            "structure_id": structure_id,
            "descriptor": descriptor,
            "atom_count": len(ligand_df),
        }

        relative_path = str(file_path.relative_to(self.paths.data_root))
        self.entity_registry.register_entity(
            name=entity_name,
            format_type=self.processor_type,
            file_path=relative_path,
            metadata=metadata,
        )

        try:
            self.entity_registry.add_relationship(
                source_name=entity_name,
                target_name=structure_id,
                rel_type="derived_from",
                metadata={"descriptor": descriptor},
            )
        except ValueError:
            self.logger.warning(
                "Structure '%s' not registered; skipping relationship for ligand '%s'",
                structure_id,
                entity_name,
            )

        return entity_name

    def _pairwise_distances(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        diff = a[:, None, :] - b[None, :, :]
        return np.linalg.norm(diff, axis=-1)

    def _format_ligand_descriptor(self, values: Union[str, Tuple[Any, ...]]) -> str:
        if isinstance(values, tuple):
            parts = [str(v) for v in values if v is not None and str(v) != "nan"]
            return "_".join(parts)
        return str(values)

    def _coordinate_columns(self, df: pd.DataFrame) -> List[str]:
        for candidate in ("x_coord", "y_coord", "z_coord"), ("x", "y", "z"):
            if all(col in df.columns for col in candidate):
                return list(candidate)
        raise ValueError("Structure dataframe lacks coordinate columns ('x_coord'/'y_coord'/'z_coord')")

    def _structure_metadata_columns(self, df: pd.DataFrame) -> Tuple[str, str, str]:
        chain_candidates = ["auth_asym_id", "auth_chain_id", "chain_id", "asym_id"]
        resseq_candidates = ["auth_seq_id", "label_seq_id", "seq_id", "residue_number"]
        resname_candidates = ["auth_comp_id", "label_comp_id", "res_name", "comp_id"]

        def pick(candidates: List[str]) -> str:
            for col in candidates:
                if col in df.columns:
                    return col
            raise ValueError(f"None of {candidates} present in structure dataframe")

        return pick(chain_candidates), pick(resseq_candidates), pick(resname_candidates)
