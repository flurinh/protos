
"""GraphProcessor: manages structure-derived graph entities."""

from __future__ import annotations

import json
import math
import pickle
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.processing.structure import StructureProcessor

try:  # Optional PyTorch / PyG support
    import torch
except ImportError:  # pragma: no cover - optional dependency
    torch = None

try:  # pragma: no cover - optional dependency
    from torch_geometric.data import Data as PyGData
except ImportError:  # pragma: no cover - optional dependency
    PyGData = None


_ELEMENT_TO_Z: Dict[str, int] = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'P': 15,
    'S': 16,
    'SE': 34,
}

# Column name mappings for normalization
# NOTE: These are used for _get_column (finding columns) and _normalize_columns (renaming)
# Be careful not to create duplicate columns - each column should only appear in ONE alias list
_COLUMN_ALIASES = {
    'chain_id': ['chain_id', 'auth_chain_id', 'auth_asym_id', 'label_asym_id'],
    'residue_number': ['residue_number', 'auth_seq_id', 'label_seq_id', 'res_id'],
    'residue_name': ['residue_name', 'res_name3l', 'auth_comp_id', 'label_comp_id', 'res_name'],
    'atom_name': ['atom_name', 'label_atom_id', 'auth_atom_id'],
    'element': ['element', 'type_symbol'],  # 'atom_name' removed - it's a different column
    'x': ['x', 'x_coord', 'Cartn_x'],
    'y': ['y', 'y_coord', 'Cartn_y'],
    'z': ['z', 'z_coord', 'Cartn_z'],
    'group': ['group', 'record_type'],
    'grn': ['grn', 'grn_position', 'generic_number'],
}


class GraphProcessor(BaseProcessor):
    """Processor for PyTorch Geometric compatible graphs derived from structures."""

    processor_type = "graph"

    def __init__(
        self,
        name: str = "graph_processor",
        *,
        default_cutoff: float = 5.0,
        default_level: str = "atom",
    ) -> None:
        super().__init__(name=name)

        self.graphs_dir = self.get_subdirectory_path('graphs_dir')
        self.datasets_dir = self.get_subdirectory_path('datasets_dir')
        self.graphs_dir.mkdir(parents=True, exist_ok=True)
        self.datasets_dir.mkdir(parents=True, exist_ok=True)

        self.default_cutoff = float(default_cutoff)
        self.default_level = default_level

        # Cached StructureProcessor instance
        self._struct_proc: Optional[StructureProcessor] = None

        self._check_dependencies()

    # ------------------------------------------------------------------
    # Dependency helpers
    # ------------------------------------------------------------------
    def _check_dependencies(self) -> None:
        if torch is None or PyGData is None:
            self.logger.warning(
                "PyTorch Geometric is not available. Graphs will be saved as serialised "
                "NumPy dictionaries. Install with: pip install torch torch_geometric"
            )

    @property
    def dependencies_available(self) -> bool:
        return torch is not None and PyGData is not None

    @property
    def struct_proc(self) -> StructureProcessor:
        """Lazy-loaded StructureProcessor instance."""
        if self._struct_proc is None:
            self._struct_proc = StructureProcessor()
        return self._struct_proc

    # ------------------------------------------------------------------
    # Column normalization
    # ------------------------------------------------------------------
    def _get_column(self, df: pd.DataFrame, canonical_name: str) -> Optional[str]:
        """Find the actual column name for a canonical name."""
        aliases = _COLUMN_ALIASES.get(canonical_name, [canonical_name])
        for alias in aliases:
            if alias in df.columns:
                return alias
        return None

    def _normalize_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize column names to canonical form.

        Only renames if the canonical name doesn't already exist to avoid duplicates.
        """
        df = df.copy()
        rename_map = {}
        for canonical, aliases in _COLUMN_ALIASES.items():
            # Skip if canonical column already exists
            if canonical in df.columns:
                continue
            # Find first alias that exists and rename it
            for alias in aliases:
                if alias in df.columns and alias != canonical:
                    rename_map[alias] = canonical
                    break
        if rename_map:
            df = df.rename(columns=rename_map)
        return df

    # ------------------------------------------------------------------
    # Structure filtering
    # ------------------------------------------------------------------
    def _filter_structure(
        self,
        df: pd.DataFrame,
        *,
        chain: Optional[str] = None,
        grn_positions: Optional[List[str]] = None,
        near_hetatm: Optional[Tuple[str, float]] = None,
        protein_only: bool = True,
    ) -> pd.DataFrame:
        """Filter structure DataFrame based on criteria.

        Args:
            df: Structure DataFrame
            chain: Filter to specific chain ID
            grn_positions: Filter to residues with these GRN positions
            near_hetatm: Tuple of (hetatm_name, distance) to filter residues near a ligand
            protein_only: If True, exclude HETATM records (waters, ligands, etc.)

        Returns:
            Filtered DataFrame
        """
        df = self._normalize_columns(df)

        # Filter to protein atoms only
        if protein_only:
            group_col = self._get_column(df, 'group')
            if group_col and group_col in df.columns:
                df = df[df[group_col] == 'ATOM']

        # Filter by chain
        if chain:
            chain_col = self._get_column(df, 'chain_id')
            if chain_col and chain_col in df.columns:
                df = df[df[chain_col] == chain]

        # Filter by GRN positions
        if grn_positions:
            grn_col = self._get_column(df, 'grn')
            if grn_col and grn_col in df.columns:
                df = df[df[grn_col].isin(grn_positions)]
            else:
                self.logger.warning("GRN column not found, cannot filter by GRN positions")

        # Filter by proximity to HETATM
        if near_hetatm:
            hetatm_name, distance = near_hetatm
            df = self._filter_near_hetatm(df, hetatm_name, distance)

        return df

    def _filter_near_hetatm(
        self,
        df: pd.DataFrame,
        hetatm_name: str,
        distance: float,
    ) -> pd.DataFrame:
        """Filter to residues within distance of a HETATM (ligand)."""
        # Need the original unfiltered data to find the HETATM
        # Get coordinates of the ligand
        group_col = self._get_column(df, 'group')
        resname_col = self._get_column(df, 'residue_name')

        # We need access to the full structure to find the ligand
        # This is a limitation - we need the full DataFrame passed in
        # For now, assume the ligand atoms are in the DataFrame if protein_only=False was used earlier

        x_col = self._get_column(df, 'x')
        y_col = self._get_column(df, 'y')
        z_col = self._get_column(df, 'z')

        if not all([x_col, y_col, z_col]):
            self.logger.warning("Coordinate columns not found, cannot filter by distance")
            return df

        # This method expects the full structure with HETATM records
        # Find ligand atoms
        if group_col and resname_col:
            ligand_mask = (df[resname_col] == hetatm_name)
            if group_col in df.columns:
                ligand_mask = ligand_mask | ((df[group_col] == 'HETATM') & (df[resname_col] == hetatm_name))
        else:
            self.logger.warning("Cannot identify ligand atoms")
            return df

        ligand_atoms = df[ligand_mask]
        if ligand_atoms.empty:
            self.logger.warning(f"No atoms found for HETATM '{hetatm_name}'")
            return df

        ligand_coords = ligand_atoms[[x_col, y_col, z_col]].values

        # Get protein atoms
        if group_col and group_col in df.columns:
            protein_df = df[df[group_col] == 'ATOM'].copy()
        else:
            protein_df = df.copy()

        if protein_df.empty:
            return protein_df

        protein_coords = protein_df[[x_col, y_col, z_col]].values

        # Compute minimum distance from each protein atom to any ligand atom
        # Using broadcasting for efficiency
        min_distances = np.min(
            np.linalg.norm(
                protein_coords[:, None, :] - ligand_coords[None, :, :],
                axis=-1
            ),
            axis=1
        )

        # Get residue numbers of atoms within distance
        resnum_col = self._get_column(protein_df, 'residue_number')
        chain_col = self._get_column(protein_df, 'chain_id')

        close_mask = min_distances <= distance

        if resnum_col:
            # Get unique residues that have at least one atom within distance
            close_residues = protein_df.loc[close_mask, [resnum_col]].drop_duplicates()
            if chain_col and chain_col in protein_df.columns:
                close_residues = protein_df.loc[close_mask, [chain_col, resnum_col]].drop_duplicates()
                result = protein_df.merge(close_residues, on=[chain_col, resnum_col])
            else:
                result = protein_df.merge(close_residues, on=[resnum_col])
            return result
        else:
            # Fall back to atom-level filtering
            return protein_df[close_mask]

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def generate_graph(
        self,
        structure_id: str,
        *,
        structure: Optional[pd.DataFrame] = None,
        # Filtering options
        chain: Optional[str] = None,
        grn_positions: Optional[List[str]] = None,
        near_hetatm: Optional[Tuple[str, float]] = None,
        protein_only: bool = True,
        # Graph construction options
        level: str = "atom",
        cutoff: Optional[float] = None,
        include_hydrogens: bool = False,
        # Output options
        dataset_name: Optional[str] = None,
        entity_name: Optional[str] = None,
    ) -> str:
        """Generate and persist a graph for a single structure.

        Args:
            structure_id: Entity name of the structure to load
            structure: Optional pre-loaded DataFrame (if None, loads from entity registry)
            chain: Filter to specific chain ID
            grn_positions: Filter to residues with these GRN positions
            near_hetatm: Tuple of (hetatm_name, distance) to filter residues near a ligand
            protein_only: If True, exclude HETATM records
            level: 'atom' or 'residue' level graph
            cutoff: Distance cutoff for edges (default: self.default_cutoff)
            include_hydrogens: Include hydrogen atoms
            dataset_name: Optional dataset to associate with
            entity_name: Optional custom entity name for the graph

        Returns:
            Entity name of the generated graph
        """
        cutoff = float(cutoff or self.default_cutoff)
        level = level or self.default_level

        # Load structure if not provided
        if structure is None:
            structure = self.struct_proc.load_entity(structure_id)
        if structure is None:
            raise ValueError(f"Structure '{structure_id}' could not be loaded")

        # For near_hetatm filtering, we need to keep HETATM temporarily
        if near_hetatm:
            # Filter with protein_only=False first to keep ligand for distance calc
            filtered = self._filter_structure(
                structure,
                chain=chain,
                grn_positions=grn_positions,
                near_hetatm=near_hetatm,
                protein_only=False,  # Keep HETATM for distance calculation
            )
            # The _filter_near_hetatm already returns only protein atoms
        else:
            filtered = self._filter_structure(
                structure,
                chain=chain,
                grn_positions=grn_positions,
                protein_only=protein_only,
            )

        if filtered.empty:
            raise ValueError(f"No atoms remaining after filtering structure '{structure_id}'")

        # Build graph
        graph_payload = self._structure_to_graph(
            structure_id,
            filtered,
            level=level,
            cutoff=cutoff,
            include_hydrogens=include_hydrogens,
        )

        # Add filter metadata
        graph_payload["graph_metadata"]["filters"] = {
            "chain": chain,
            "grn_positions": grn_positions,
            "near_hetatm": near_hetatm,
            "protein_only": protein_only,
        }

        # Persist
        result_name = self._persist_graph(
            structure_id,
            graph_payload,
            dataset_name=dataset_name,
            entity_name=entity_name,
        )

        return result_name

    def generate_graphs_from_dataset(
        self,
        structure_dataset: str,
        *,
        dataset_name: Optional[str] = None,
        # Filtering options
        chain: Optional[Union[str, Dict[str, str]]] = None,
        grn_positions: Optional[List[str]] = None,
        near_hetatm: Optional[Union[Tuple[str, float], Dict[str, Tuple[str, float]]]] = None,
        protein_only: bool = True,
        # Graph construction options
        level: str = "atom",
        cutoff: Optional[float] = None,
        include_hydrogens: bool = False,
    ) -> Tuple[str, List[str]]:
        """Generate graphs for every structure in a dataset.

        Args:
            structure_dataset: Name of the structure dataset
            dataset_name: Name for the output graph dataset
            chain: Chain filter (str for all, or dict mapping structure_id -> chain)
            grn_positions: GRN positions to include
            near_hetatm: Ligand proximity filter (tuple for all, or dict mapping structure_id -> tuple)
            protein_only: Exclude HETATM records
            level: 'atom' or 'residue' level
            cutoff: Edge distance cutoff
            include_hydrogens: Include hydrogen atoms

        Returns:
            Tuple of (dataset_name, list of graph entity names)
        """
        cutoff = float(cutoff or self.default_cutoff)
        level = level or self.default_level

        structures = self.struct_proc.load_dataset(
            structure_dataset,
            return_format="dict",
        )
        if not isinstance(structures, dict) or len(structures) == 0:
            raise ValueError(f"Structure dataset '{structure_dataset}' is empty")

        graph_names: List[str] = []
        for structure_id, structure_df in structures.items():
            try:
                # Resolve per-structure filters
                struct_chain = chain[structure_id] if isinstance(chain, dict) else chain
                struct_near = near_hetatm[structure_id] if isinstance(near_hetatm, dict) else near_hetatm

                graph_name = self.generate_graph(
                    structure_id,
                    structure=structure_df,
                    chain=struct_chain,
                    grn_positions=grn_positions,
                    near_hetatm=struct_near,
                    protein_only=protein_only,
                    level=level,
                    cutoff=cutoff,
                    include_hydrogens=include_hydrogens,
                    dataset_name=None,  # postpone dataset registration until the end
                )
                graph_names.append(graph_name)
            except Exception as exc:
                self.logger.warning(
                    "Failed to generate graph for %s: %s", structure_id, exc
                )

        if not graph_names:
            raise RuntimeError("No graphs generated from the structure dataset")

        graph_dataset_name = dataset_name or self._sanitize_filename(
            f"{structure_dataset}__{level}__{cutoff:.1f}A"
        )

        metadata = {
            "source_structure_dataset": structure_dataset,
            "level": level,
            "cutoff": cutoff,
            "include_hydrogens": include_hydrogens,
            "entity_count": len(graph_names),
            "filters": {
                "grn_positions": grn_positions,
                "protein_only": protein_only,
            },
        }

        if self.dataset_manager.dataset_exists(graph_dataset_name):
            self.dataset_manager.delete_dataset(graph_dataset_name)

        self.create_dataset(graph_dataset_name, graph_names, metadata)

        return graph_dataset_name, graph_names

    # ------------------------------------------------------------------
    # BaseProcessor abstract implementations
    # ------------------------------------------------------------------
    def load_entity(self, name: str) -> Optional[Dict[str, Any]]:
        try:
            payload = self.load_graph(name)
            return payload
        except Exception:
            return None

    def save_entity(
        self,
        name: str,
        data: Dict[str, Any],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        graph_meta = dict(data.get("graph_metadata", {}))
        if metadata:
            graph_meta.update(metadata)

        structure_id = graph_meta.get("structure_id", name)
        dataset_name = graph_meta.get("dataset")
        entity_name = self._sanitize_filename(name)

        return self._persist_graph(
            structure_id,
            data,
            dataset_name=dataset_name,
            entity_name=entity_name,
        )

    def load_graph(self, entity_name: str) -> Dict[str, Any]:
        """Load a persisted graph entity."""

        entity_info = self.entity_registry.find_entity(entity_name, self.processor_type)
        if entity_info is None:
            raise ValueError(f"Graph entity '{entity_name}' not found")

        file_path = Path(self.paths.data_root) / entity_info.file_path
        if not file_path.exists():
            raise FileNotFoundError(f"Graph file not found: {file_path}")

        with open(file_path, "rb") as handle:
            return pickle.load(handle)

    def list_graphs(self, dataset: Optional[str] = None) -> List[str]:
        """List registered graph entities."""

        if dataset:
            return self.dataset_manager.get_dataset_entities(dataset)
        return self.entity_registry.list_entities(self.processor_type)

    def to_pyg(self, graph_dict: Dict[str, Any]) -> "PyGData":
        """Convert a stored graph dictionary to a PyG Data object."""

        if not self.dependencies_available:
            raise RuntimeError("PyTorch Geometric is not installed")

        pos = torch.tensor(graph_dict["node_positions"], dtype=torch.float32)
        edge_index = torch.tensor(graph_dict["edge_index"], dtype=torch.long)

        edge_weight = graph_dict.get("edge_weight")
        if edge_weight is not None:
            edge_attr = torch.tensor(edge_weight, dtype=torch.float32).unsqueeze(1)
        else:
            edge_attr = None

        x_data = graph_dict.get("node_features")
        if x_data is not None:
            x = torch.tensor(x_data, dtype=torch.float32)
        else:
            x = None

        pyg_data = PyGData(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_attr,
            pos=pos,
        )
        return pyg_data

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _structure_to_graph(
        self,
        structure_id: str,
        structure_df: pd.DataFrame,
        *,
        level: str,
        cutoff: float,
        include_hydrogens: bool,
    ) -> Dict[str, Any]:
        """Convert a structure DataFrame into a graph payload."""

        df = self._normalize_columns(structure_df.copy())

        # Find coordinate columns
        x_col = self._get_column(df, 'x')
        y_col = self._get_column(df, 'y')
        z_col = self._get_column(df, 'z')

        if not all([x_col, y_col, z_col]):
            raise ValueError("Structure dataframe missing coordinate columns")
        coord_cols = [x_col, y_col, z_col]

        # Find element column
        element_col = self._get_column(df, 'element')
        if element_col is None:
            element_col = "element"
            df[element_col] = "C"  # default fallback

        # Filter hydrogens
        if not include_hydrogens:
            df = df[df[element_col].astype(str).str.upper().str[0] != "H"]

        if df.empty:
            raise ValueError("Structure has no atoms after filtering")

        # Get column names for grouping
        chain_col = self._get_column(df, 'chain_id') or 'chain_id'
        resnum_col = self._get_column(df, 'residue_number') or 'residue_number'
        resname_col = self._get_column(df, 'residue_name') or 'residue_name'
        atomname_col = self._get_column(df, 'atom_name') or 'atom_name'
        grn_col = self._get_column(df, 'grn')

        if level == "residue":
            # Group by residue
            group_cols = [col for col in [chain_col, resnum_col, resname_col] if col in df.columns]
            if not group_cols:
                group_cols = [resnum_col] if resnum_col in df.columns else df.columns[:1].tolist()

            grouped = df.groupby(group_cols, dropna=False)
            pos = grouped[coord_cols].mean().to_numpy(dtype=np.float32)

            elements = grouped[element_col].agg(
                lambda series: Counter(series).most_common(1)[0][0]
            ).tolist()

            node_metadata = []
            for key in grouped.indices.keys():
                if isinstance(key, tuple):
                    meta = {
                        "chain": key[0] if len(key) > 0 else None,
                        "residue_number": key[1] if len(key) > 1 else None,
                        "residue_name": key[2] if len(key) > 2 else None,
                    }
                else:
                    meta = {"residue_number": key}

                # Add GRN if available
                if grn_col and grn_col in df.columns:
                    grp_df = grouped.get_group(key)
                    grn_vals = grp_df[grn_col].dropna().unique()
                    meta["grn"] = grn_vals[0] if len(grn_vals) > 0 else None

                node_metadata.append(meta)
        else:
            # Atom level
            pos = df[coord_cols].to_numpy(dtype=np.float32)
            elements = df[element_col].tolist()

            node_metadata = []
            for i in range(len(df)):
                row = df.iloc[i]
                meta = {
                    "atom_name": row.get(atomname_col) if atomname_col in df.columns else None,
                    "residue_name": row.get(resname_col) if resname_col in df.columns else None,
                    "chain": row.get(chain_col) if chain_col in df.columns else None,
                    "residue_number": row.get(resnum_col) if resnum_col in df.columns else None,
                }
                if grn_col and grn_col in df.columns:
                    meta["grn"] = row.get(grn_col)
                node_metadata.append(meta)

        node_count = len(pos)
        if node_count < 1:
            raise ValueError("No nodes available for graph construction")

        atomic_numbers = np.array(
            [self._atomic_number(element) for element in elements], dtype=np.float32
        ).reshape(-1, 1)

        edge_index, edge_weight = self._build_edges(pos, cutoff=cutoff)

        graph_metadata = {
            "structure_id": structure_id,
            "level": level,
            "cutoff": cutoff,
            "include_hydrogens": include_hydrogens,
            "node_count": node_count,
            "edge_count": int(edge_index.shape[1]) if edge_index.size else 0,
        }

        payload = {
            "node_positions": pos,
            "node_features": atomic_numbers,
            "node_metadata": node_metadata,
            "edge_index": edge_index,
            "edge_weight": edge_weight,
            "graph_metadata": graph_metadata,
        }

        return payload

    def _persist_graph(
        self,
        structure_id: str,
        graph_payload: Dict[str, Any],
        *,
        dataset_name: Optional[str],
        entity_name: Optional[str] = None,
    ) -> str:
        """Persist a graph dictionary and register it in the entity registry."""

        graph_meta = graph_payload["graph_metadata"]
        if entity_name is None:
            name_parts = [
                structure_id,
                graph_meta["level"],
                f"{graph_meta['cutoff']:.1f}A",
            ]
            if dataset_name:
                name_parts.append(dataset_name)
            entity_name = self._sanitize_filename("__".join(map(str, name_parts)))
        else:
            entity_name = self._sanitize_filename(entity_name)

        file_path = self.graphs_dir / f"{entity_name}.pkl"
        with open(file_path, "wb") as handle:
            pickle.dump(graph_payload, handle)

        relative_path = str(file_path.relative_to(self.paths.data_root))

        metadata = {
            **graph_meta,
            "dependencies_available": self.dependencies_available,
        }
        if dataset_name:
            metadata["dataset"] = dataset_name

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
                metadata={
                    "level": graph_meta["level"],
                    "cutoff": graph_meta["cutoff"],
                },
            )
        except ValueError:
            self.logger.debug(
                "Structure '%s' is not registered; relationship skipped", structure_id
            )

        return entity_name

    def _build_edges(
        self,
        positions: np.ndarray,
        *,
        cutoff: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute edges based on Euclidean distance cutoff."""

        if positions.shape[0] <= 1:
            return np.empty((2, 0), dtype=np.int64), np.empty((0,), dtype=np.float32)

        diff = positions[:, None, :] - positions[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)

        mask = (dist <= cutoff) & (dist > 0)
        src, dst = np.where(mask)
        edge_index = np.vstack([src, dst]).astype(np.int64)
        edge_weight = dist[src, dst].astype(np.float32)

        return edge_index, edge_weight

    def _atomic_number(self, element: Any) -> float:
        if isinstance(element, str):
            key = element.strip().upper()[:2]  # Handle 2-letter elements
            if key in _ELEMENT_TO_Z:
                return float(_ELEMENT_TO_Z[key])
            # Try first letter only
            return float(_ELEMENT_TO_Z.get(key[0], 0))
        return 0.0
