"""Loader utilities for ligand ingestion."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from protos.io.core.base_loader import BaseLoader
from protos.processing.structure import StructureProcessor
from protos.processing.molecule import MoleculeProcessor
from protos.io.paths import sanitize_storage_name

try:  # Optional RDKit dependency
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover - optional
    Chem = None

from protos.io.formats.molecule_utils import mol_to_canonical_dataframe, sdf_to_molecules, require_rdkit
import shutil


class LigandLoader(BaseLoader):
    """Loader for ligand data (SMILES, SDF, etc.) into structure storage."""

    loader_type = "structure"

    def __init__(
        self,
        name: str = "ligand_loader",
        *,
        structure_processor: Optional[StructureProcessor] = None,
        ligand_processor: Optional[MoleculeProcessor] = None,
    ) -> None:
        super().__init__(name=name)
        self._structure_processor = structure_processor
        self._ligand_processor = ligand_processor

    def _get_structure_processor(self) -> StructureProcessor:
        if self._structure_processor is None:
            self._structure_processor = StructureProcessor(name=f"{self.name}_structure")
        return self._structure_processor

    def _get_ligand_processor(self) -> MoleculeProcessor:
        if self._ligand_processor is None:
            self._ligand_processor = MoleculeProcessor(name=f"{self.name}_molecule")
        return self._ligand_processor

    # BaseLoader abstract requirements -------------------------------------------------

    def parse_identifier(self, identifier: str) -> Dict[str, Any]:
        path = Path(identifier)
        if not path.exists():
            raise FileNotFoundError(identifier)
        return {
            'id': sanitize_storage_name(path.stem),
            'source': 'local',
            'type': 'sdf',
            'filename': path.name,
            'original_id': identifier,
        }

    def fetch_entity(self, identifier: str, **kwargs) -> Optional[Path]:
        """Copy the SDF into structure-managed storage and return its path."""

        info = self.parse_identifier(identifier)
        source_path = Path(info['original_id'])
        structure_proc = self._get_structure_processor()
        target_dir = structure_proc.path_sdf_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        target_path = target_dir / info['filename']
        if source_path.resolve() != target_path.resolve():
            shutil.copy2(source_path, target_path)
        return target_path

    def import_sdf(
        self,
        file_path: str,
        *,
        dataset_name: Optional[str] = None,
        chain_id: str = 'L',
    ) -> Tuple[str, List[str]]:
        """Import ligands from an SDF file and register them as structure entities."""

        require_rdkit()
        stored_path = self.fetch_entity(file_path)
        path = Path(stored_path) if stored_path is not None else Path(file_path)
        if not path.exists():
            raise FileNotFoundError(file_path)

        structure_proc = self._get_structure_processor()
        ligand_proc = self._get_ligand_processor()

        # Register raw SDF artifact under structure domain
        try:
            rel_path = path.relative_to(self.paths.data_root)
        except ValueError:
            rel_path = path
        try:
            self.entity_registry.register_entity(
                name=sanitize_storage_name(path.stem),
                format_type='structure',
                file_path=str(rel_path),
                metadata={
                    'source_file': path.name,
                    'source_format': 'sdf',
                    'loader': self.__class__.__name__,
                },
            )
        except Exception:
            # Ignore if raw artifact already registered
            pass

        dataset_entities: List[str] = []
        molecule_entities: List[str] = []

        for index, mol in enumerate(sdf_to_molecules(path), start=1):
            if mol.GetNumAtoms() == 0:
                continue

            base_name = mol.GetProp('_Name') if mol.HasProp('_Name') and mol.GetProp('_Name') else path.stem
            structure_id = sanitize_storage_name(f"{base_name}_{index}")
            df = mol_to_canonical_dataframe(mol, structure_id=structure_id, chain_id=chain_id)

            smiles = Chem.MolToSmiles(mol)
            metadata = {
                'entity_type': 'ligand',
                'source_format': 'sdf',
                'source_file': str(path.name),
                'canonical_smiles': smiles,
                'mol_block': Chem.MolToMolBlock(mol),
            }
            structure_proc.save_entity(structure_id, df, metadata=metadata)
            dataset_entities.append(structure_id)

            molecule_name = ligand_proc.save_entity(
                structure_id,
                {
                    'smiles': smiles,
                    'input_name': structure_id,
                    'kind': 'structure_record',
                },
                metadata={
                    'structure_entity': structure_id,
                    'source_file': path.name,
                },
            )
            molecule_entities.append(molecule_name)
            try:
                self.entity_registry.add_relationship(
                    source_name=molecule_name,
                    target_name=structure_id,
                    rel_type='has_structure',
                )
            except ValueError:
                pass

        if not dataset_entities:
            raise RuntimeError(f"No valid molecules found in {file_path}")

        dataset_id = dataset_name or sanitize_storage_name(path.stem)
        dataset_meta = {
            'source_file': path.name,
            'entity_count': len(dataset_entities),
            'kind': 'structure_ligand',
        }

        if structure_proc.dataset_manager.dataset_exists(dataset_id):
            structure_proc.dataset_manager.delete_dataset(dataset_id)

        structure_proc.create_dataset(dataset_id, dataset_entities, dataset_meta)
        molecule_dataset_id = f"{dataset_id}_molecules"
        molecule_metadata = {
            'source_file': path.name,
            'entity_count': len(molecule_entities),
            'kind': 'smiles',
        }

        if ligand_proc.dataset_manager.dataset_exists(molecule_dataset_id):
            ligand_proc.dataset_manager.delete_dataset(molecule_dataset_id)
        ligand_proc.create_dataset(molecule_dataset_id, molecule_entities, molecule_metadata)

        return dataset_id, dataset_entities

    def import_smiles(
        self,
        smiles_map: Dict[str, str],
        *,
        dataset_name: Optional[str] = None,
        chain_id: str = 'L',
        generate_3d: bool = True,
    ) -> Tuple[str, List[str]]:
        """Import ligands from SMILES strings and register as structure entities.

        All path management is handled by processors; no filesystem writes occur here.
        """
        require_rdkit()

        structure_proc = self._get_structure_processor()
        ligand_proc = self._get_ligand_processor()

        dataset_entities: List[str] = []
        molecule_entities: List[str] = []

        for idx, (name, smiles) in enumerate(smiles_map.items(), start=1):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            if generate_3d:
                mol = Chem.AddHs(mol)
                try:
                    AllChem.EmbedMolecule(mol, randomSeed=0xF00D)
                    AllChem.UFFOptimizeMolecule(mol)
                except Exception:
                    # Fall back to 2D coords if embedding fails
                    pass
            if not mol.HasProp('_Name'):
                mol.SetProp('_Name', name)

            structure_id = sanitize_storage_name(f"{name}_{idx}")
            df = mol_to_canonical_dataframe(mol, structure_id=structure_id, chain_id=chain_id)
            metadata = {
                'entity_type': 'ligand',
                'source_format': 'smiles',
                'canonical_smiles': smiles,
            }
            structure_proc.save_entity(structure_id, df, metadata=metadata)
            dataset_entities.append(structure_id)

            molecule_name = ligand_proc.save_entity(
                structure_id,
                {
                    'smiles': smiles,
                    'input_name': structure_id,
                    'kind': 'smiles_record',
                },
                metadata={'structure_entity': structure_id},
            )
            molecule_entities.append(molecule_name)
            try:
                self.entity_registry.add_relationship(
                    source_name=molecule_name,
                    target_name=structure_id,
                    rel_type='has_structure',
                )
            except ValueError:
                pass

        if not dataset_entities:
            raise RuntimeError("No valid SMILES provided")

        ds_name = dataset_name or self._sanitize_filename("smiles_ligands")
        ds_meta = {'entity_count': len(dataset_entities), 'kind': 'structure_ligand'}
        if structure_proc.dataset_manager.dataset_exists(ds_name):
            structure_proc.dataset_manager.delete_dataset(ds_name)
        structure_proc.create_dataset(ds_name, dataset_entities, ds_meta)

        mol_ds_name = f"{ds_name}_molecules"
        mol_meta = {'entity_count': len(molecule_entities), 'kind': 'smiles'}
        if ligand_proc.dataset_manager.dataset_exists(mol_ds_name):
            ligand_proc.dataset_manager.delete_dataset(mol_ds_name)
        ligand_proc.create_dataset(mol_ds_name, molecule_entities, mol_meta)

        return ds_name, dataset_entities

__all__ = ['LigandLoader']
