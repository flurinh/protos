"""Exporter for structure-derived ligand files."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from protos.io.formats.molecule_utils import canonical_dataframe_to_mol, require_rdkit
from protos.processing.molecule import MoleculeProcessor

try:  # Optional RDKit dependency
    from rdkit import Chem
except ImportError:  # pragma: no cover - optional
    Chem = None


class LigandExporter:
    """Convert canonical structure frames into SDF representations."""

    def __init__(
        self,
        structure_processor,
        *,
        molecule_processor: Optional[MoleculeProcessor] = None,
    ) -> None:
        self.structure_processor = structure_processor
        self.molecule_processor = molecule_processor or MoleculeProcessor(name=f"{structure_processor.name}_molecules")

    def export_structure_to_sdf(
        self,
        structure_id: str,
        *,
        out_path: Optional[Path] = None,
        register: bool = True,
        overwrite: bool = True,
    ) -> Path:
        require_rdkit()
        df = self.structure_processor.load_entity(structure_id)
        if df is None:
            raise ValueError(f"Structure '{structure_id}' could not be loaded")

        entity_info = self.structure_processor.entity_registry.find_entity(structure_id, format_type='structure')
        metadata = dict(entity_info.metadata or {}) if entity_info else {}

        mol = canonical_dataframe_to_mol(df, metadata=metadata)
        out_path = self._resolve_output_path(structure_id, out_path)
        if out_path.exists() and not overwrite:
            raise FileExistsError(f"File {out_path} already exists")

        out_path.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(out_path))
        writer.write(mol)
        writer.close()

        if register:
            self._register_export(structure_id, out_path, metadata)

        return out_path

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _resolve_output_path(self, structure_id: str, out_path: Optional[Path]) -> Path:
        if out_path is not None:
            return Path(out_path)

        sdf_dir = self.structure_processor.path_sdf_dir
        sdf_dir.mkdir(parents=True, exist_ok=True)
        return sdf_dir / f"{structure_id}.sdf"

    def _register_export(self, structure_id: str, out_path: Path, metadata: dict) -> None:
        rel_path = out_path.relative_to(self.structure_processor.paths.data_root)
        record = {
            'kind': 'sdf_export',
            'structure_id': structure_id,
            'file_path': str(rel_path),
            'smiles': metadata.get('canonical_smiles'),
        }
        export_name = self.molecule_processor.save_entity(out_path.stem, record, metadata={'kind': 'sdf'})
        try:
            self.structure_processor.entity_registry.add_relationship(
                source_name=export_name,
                target_name=structure_id,
                rel_type='exported_structure',
                metadata={'format': 'sdf'},
            )
        except ValueError:
            pass

    def export_dataset_to_sdf(
        self,
        dataset_name: str,
        *,
        out_path: Optional[Path] = None,
        overwrite: bool = True,
        register: bool = False,
    ) -> Path:
        """Export all structures in a dataset into a single multi-mol SDF.

        Writes to structure/sdf/{dataset_name}.sdf by default.
        """
        require_rdkit()
        entities = self.structure_processor.get_dataset_entities(dataset_name)
        if not entities:
            raise ValueError(f"Dataset '{dataset_name}' has no entities")

        sdf_dir = self.structure_processor.path_sdf_dir
        target = Path(out_path) if out_path is not None else (sdf_dir / f"{dataset_name}.sdf")
        if target.exists() and not overwrite:
            raise FileExistsError(f"Export already exists: {target}")

        target.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(target))
        try:
            for name in entities:
                df = self.structure_processor.load_entity(name)
                if df is None:
                    continue
                entity_info = self.structure_processor.entity_registry.find_entity(name, format_type='structure')
                meta = dict(entity_info.metadata or {}) if entity_info else {}
                mol = canonical_dataframe_to_mol(df, metadata=meta)
                if not mol.HasProp('Name'):
                    mol.SetProp('Name', name)
                writer.write(mol)
        finally:
            writer.close()

        if register:
            rel_path = target.relative_to(self.structure_processor.paths.data_root)
            dataset_record = {
                'kind': 'sdf_dataset_export',
                'dataset_name': dataset_name,
                'file_path': str(rel_path),
                'entity_count': len(entities),
            }
            self.molecule_processor.save_entity(f"{dataset_name}_sdf", dataset_record, metadata={'kind': 'sdf'})

        return target


__all__ = ["LigandExporter"]
