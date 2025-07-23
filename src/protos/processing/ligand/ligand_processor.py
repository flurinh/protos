"""
Ligand processor for Protos.

This processor handles ligand data including:
- Small molecule structures and properties
- Bioactivity data
- Target-ligand interactions
- QSAR modeling
- Virtual screening

The processor follows Protos' BaseProcessor architecture and integrates
with the entity registry and dataset management system.
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Union, Any, Tuple
from datetime import datetime

from protos.core.base_processor import BaseProcessor
from protos.io.paths.path_config import ProtosPaths
from protos.io.entity_registry import EntityRegistry
from protos.io.dataset_manager import DatasetManager
# from protos.loaders.chembl_loader import ChEMBLDL  # Avoid circular import
from protos.loaders.ligand_utils import (
    sanitize_smiles_filename, validate_smiles, smiles_to_inchi,
    calculate_molecular_properties, is_drug_like, parse_activity_value
)

logger = logging.getLogger(__name__)

# Try to import optional dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    from rdkit.Chem.Fingerprints import FingerprintMols
    from rdkit.DataStructs import TanimotoSimilarity
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some ligand functionality will be limited.")
    HAS_RDKIT = False


class LigandProcessor(BaseProcessor):
    """
    Processor for handling ligand data in Protos.
    
    This processor manages:
    - Ligand structures (SMILES, SDF, MOL)
    - Bioactivity data from various sources
    - Molecular properties and descriptors
    - Target-ligand interactions
    - QSAR models and predictions
    """
    
    def __init__(self, name: str = "ligand_processor", 
                 paths: Optional[ProtosPaths] = None,
                 cif_processor=None):
        """
        Initialize the LigandProcessor.
        
        Args:
            name: Name for this processor instance
            paths: ProtosPaths instance (optional)
            cif_processor: CifBaseProcessor instance for structure integration
        """
        super().__init__(name=name, paths=paths)
        
        # Set processor type
        self.processor_type = 'ligand'
        
        # Initialize data paths
        self.data_path = Path(self.paths.get_processor_path(self.processor_type))
        self.sdf_dir = self.data_path / 'sdf'
        self.tables_dir = self.data_path / 'tables'
        self.cache_dir = self.data_path / 'cache'
        self.models_dir = self.data_path / 'models'
        self.pockets_dir = self.data_path / 'pockets'
        
        # Ensure directories exist
        for dir_path in [self.sdf_dir, self.tables_dir, self.cache_dir, 
                        self.models_dir, self.pockets_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize ChEMBL loader (lazy import to avoid circular dependency)
        self.chembl_loader = None
        try:
            from protos.loaders.chembl_loader import ChEMBLDL
            self.chembl_loader = ChEMBLDL(data_root=self.paths.data_root)
        except ImportError:
            logger.warning("ChEMBL loader not available")
        
        # Store reference to CifBaseProcessor if provided
        self.cif_processor = cif_processor
        
        # Cache for loaded data
        self._ligand_cache = {}
        self._fingerprint_cache = {}
    
    def _get_processor_type(self) -> str:
        """Return the processor type identifier."""
        return 'ligand'
    
    # ===== Entity Management Methods =====
    
    def load_entity(self, name: str) -> Optional[Dict]:
        """
        Load a ligand entity by name (SMILES or alias).
        
        Args:
            name: Entity name (SMILES, ChEMBL ID, InChI Key, etc.)
            
        Returns:
            Dictionary with ligand data, or None if not found
        """
        # Check cache first
        if name in self._ligand_cache:
            return self._ligand_cache[name]
        
        # Try to find entity in registry
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        
        if entity_info:
            # Load from registered path
            ligand_data = self._load_ligand_data(entity_info)
            if ligand_data:
                self._ligand_cache[name] = ligand_data
                return ligand_data
        
        # Try to interpret as SMILES directly
        is_valid, canonical_smiles = validate_smiles(name)
        if is_valid:
            # Create minimal ligand data
            ligand_data = {
                'smiles': canonical_smiles,
                'original_input': name,
                'properties': calculate_molecular_properties(canonical_smiles)
            }
            self._ligand_cache[name] = ligand_data
            return ligand_data
        
        return None
    
    def save_entity(self, name: str, data: Dict, metadata: Optional[Dict] = None):
        """
        Save a ligand entity.
        
        Args:
            name: Entity name (SMILES)
            data: Ligand data dictionary
            metadata: Optional metadata
        """
        # Validate SMILES
        smiles = data.get('smiles', name)
        is_valid, canonical_smiles = validate_smiles(smiles)
        if not is_valid:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Save SDF file
        safe_filename = sanitize_smiles_filename(canonical_smiles)
        sdf_path = self.sdf_dir / f"{safe_filename}.sdf"
        
        # Create SDF content
        from protos.loaders.ligand_utils import create_sdf_from_smiles
        properties = data.get('properties', {})
        if 'chembl_id' in data:
            properties['CHEMBL_ID'] = data['chembl_id']
        
        sdf_content = create_sdf_from_smiles(canonical_smiles, properties)
        if sdf_content:
            with open(sdf_path, 'w') as f:
                f.write(sdf_content)
        
        # Prepare metadata
        if metadata is None:
            metadata = {}
        
        # Add standard metadata
        metadata.update({
            'molecular_properties': data.get('properties', {}),
            'targets': data.get('targets', []),
            'activities': data.get('activities', [])
        })
        
        # Add aliases
        aliases = []
        if 'chembl_id' in data:
            aliases.append(data['chembl_id'])
        if 'inchi_key' in data:
            aliases.append(data['inchi_key'])
        metadata['aliases'] = aliases
        
        # Register entity
        self.entity_registry.register_entity(
            name=canonical_smiles,
            format_type=self.processor_type,
            file_path=str(sdf_path.relative_to(Path(self.paths.data_root))),
            metadata=metadata
        )
        
        # Update cache
        self._ligand_cache[canonical_smiles] = data
    
    def list_entities(self) -> List[str]:
        """
        List all ligand entities.
        
        Returns:
            List of entity names (SMILES)
        """
        # Get all entities of type 'ligand'
        all_entities = self.entity_registry.list_entities()
        ligand_entities = []
        
        for entity_id in all_entities:
            # Check if this entity has ligand format
            # find_entity returns hash_id, we need to check the internal registry
            hash_id = self.entity_registry._resolve_to_hash(entity_id)
            if hash_id and hash_id in self.entity_registry._registry:
                entity_data = self.entity_registry._registry[hash_id]
                if 'ligand' in entity_data.get('formats', {}):
                    # The entity ID returned by list_entities is already the human-readable name
                    ligand_entities.append(entity_id)
        
        return ligand_entities
    
    def entity_exists(self, name: str) -> bool:
        """
        Check if a ligand entity exists.
        
        Args:
            name: Entity name (SMILES or alias)
            
        Returns:
            True if entity exists
        """
        return self.entity_registry.find_entity(name, self.processor_type) is not None
    
    # ===== Dataset Management Methods =====
    
    def create_dataset(self, dataset_name: str, entity_names: List[str], 
                      metadata: Optional[Dict] = None):
        """
        Create a ligand dataset.
        
        Args:
            dataset_name: Name for the dataset
            entity_names: List of ligand names (SMILES or aliases)
            metadata: Optional metadata for the dataset
        """
        # Validate entities exist
        valid_entities = []
        for name in entity_names:
            if self.entity_exists(name) or validate_smiles(name)[0]:
                valid_entities.append(name)
            else:
                logger.warning(f"Entity not found: {name}")
        
        if not valid_entities:
            raise ValueError("No valid entities found for dataset")
        
        # Create dataset
        super().create_dataset(dataset_name, valid_entities, metadata)
    
    def load_dataset(self, dataset_name: str) -> Dict[str, Dict]:
        """
        Load a ligand dataset.
        
        Args:
            dataset_name: Name of the dataset
            
        Returns:
            Dictionary mapping entity names to ligand data
        """
        # Get dataset info
        dataset_info = self.dataset_manager.load_dataset(dataset_name)
        if not dataset_info:
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        # Check if this is a table-based dataset
        if 'metadata' in dataset_info and 'data_file' in dataset_info['metadata']:
            # Load from CSV table
            table_path = self.data_path / dataset_info['metadata']['data_file'].lstrip('../')
            if table_path.exists():
                df = pd.read_csv(table_path)
                
                # Convert to dictionary format
                result = {}
                for _, row in df.iterrows():
                    smiles = row['smiles']
                    ligand_data = {
                        'smiles': smiles,
                        'chembl_id': row.get('chembl_id', ''),
                        'activity_type': row.get('activity_type', ''),
                        'value_nm': row.get('value_nm', np.nan),
                        'properties': {
                            'mw': row.get('mw', np.nan),
                            'logp': row.get('logp', np.nan),
                            'hba': row.get('hba', np.nan),
                            'hbd': row.get('hbd', np.nan),
                            'tpsa': row.get('tpsa', np.nan)
                        }
                    }
                    result[smiles] = ligand_data
                
                return result
        
        # Load individual entities
        result = {}
        for entity_name in dataset_info.get('entities', []):
            try:
                ligand_data = self.load_entity(entity_name)
                if ligand_data:
                    result[entity_name] = ligand_data
            except Exception as e:
                logger.error(f"Failed to load entity {entity_name}: {e}")
        
        return result
    
    # ===== ChEMBL Integration Methods =====
    
    def get_protein_ligands(self, protein_id: str, reload: bool = False,
                           activity_types: Optional[List[str]] = None,
                           min_pchembl: float = 5.0) -> List[Dict]:
        """
        Get ligands for a protein from ChEMBL.
        
        Args:
            protein_id: Protein identifier (UniProt, PDB, gene name)
            reload: Force reload from ChEMBL
            activity_types: Activity types to retrieve
            min_pchembl: Minimum pChEMBL value
            
        Returns:
            List of ligand dictionaries
        """
        # Set reload flag on loader
        self.chembl_loader.reload = reload
        
        # Download ligands
        results = self.chembl_loader.download_protein_ligands(
            protein_id, 
            activity_types=activity_types,
            min_pchembl=min_pchembl
        )
        
        return results.get(protein_id, [])
    
    def search_similar_ligands(self, smiles: str, similarity: float = 0.7,
                             dataset: Optional[str] = None) -> List[Tuple[str, float]]:
        """
        Search for similar ligands in the database.
        
        Args:
            smiles: Query SMILES string
            similarity: Similarity threshold (0-1)
            dataset: Optional dataset to search within
            
        Returns:
            List of (smiles, similarity_score) tuples
        """
        if not HAS_RDKIT:
            logger.error("RDKit required for similarity search")
            return []
        
        # Validate query SMILES
        is_valid, query_smiles = validate_smiles(smiles)
        if not is_valid:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Get query fingerprint
        query_fp = self._get_fingerprint(query_smiles)
        if query_fp is None:
            return []
        
        # Get ligands to search
        if dataset:
            ligands = self.load_dataset(dataset)
        else:
            # Search all ligands
            ligands = {}
            for entity_name in self.list_entities():
                ligand_data = self.load_entity(entity_name)
                if ligand_data:
                    ligands[entity_name] = ligand_data
        
        # Calculate similarities
        results = []
        for ligand_smiles, ligand_data in ligands.items():
            if ligand_smiles == query_smiles:
                continue
            
            ligand_fp = self._get_fingerprint(ligand_data.get('smiles', ligand_smiles))
            if ligand_fp is not None:
                sim = TanimotoSimilarity(query_fp, ligand_fp)
                if sim >= similarity:
                    results.append((ligand_smiles, sim))
        
        # Sort by similarity
        results.sort(key=lambda x: x[1], reverse=True)
        
        return results
    
    # ===== Property Calculation Methods =====
    
    def calculate_properties(self, smiles: str) -> Optional[Dict]:
        """
        Calculate molecular properties for a SMILES string.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Dictionary of properties, or None if failed
        """
        return calculate_molecular_properties(smiles)
    
    def filter_drug_like(self, entity_names: List[str], strict: bool = False) -> List[str]:
        """
        Filter ligands by drug-like properties.
        
        Args:
            entity_names: List of ligand names
            strict: Apply stricter criteria
            
        Returns:
            List of drug-like ligand names
        """
        drug_like = []
        
        for name in entity_names:
            ligand_data = self.load_entity(name)
            if ligand_data:
                smiles = ligand_data.get('smiles', name)
                if is_drug_like(smiles, strict=strict):
                    drug_like.append(name)
        
        return drug_like
    
    # ===== Structure Integration Methods =====
    
    def find_ligand_in_structures(self, ligand_name: str) -> List[str]:
        """
        Find structures containing a specific ligand.
        
        Args:
            ligand_name: Three-letter ligand code (e.g., 'ATP', 'NAD')
            
        Returns:
            List of PDB IDs containing the ligand
        """
        if not self.cif_processor:
            raise ValueError("No CifBaseProcessor attached")
        
        has_ligand = []
        for pdb_id in self.cif_processor.pdb_ids:
            # Filter for HETATM records matching the ligand name
            ligand_data = self.cif_processor.data[
                (self.cif_processor.data['pdb_id'] == pdb_id) &
                (self.cif_processor.data['group'] == 'HETATM') &
                (self.cif_processor.data['res_name3l'] == ligand_name)
            ]
            
            if not ligand_data.empty:
                has_ligand.append(pdb_id)
        
        return has_ligand
    
    def get_ligand_from_structure(self, pdb_id: str, ligand_name: str) -> pd.DataFrame:
        """
        Extract ligand coordinates from a structure.
        
        Args:
            pdb_id: PDB identifier
            ligand_name: Three-letter ligand code
            
        Returns:
            DataFrame with ligand atom coordinates
        """
        if not self.cif_processor:
            raise ValueError("No CifBaseProcessor attached")
        
        structure_data = self.cif_processor.get_structure_by_idx(pdb_id)
        ligand_data = structure_data[
            (structure_data['group'] == 'HETATM') &
            (structure_data['res_name3l'] == ligand_name)
        ]
        
        return ligand_data
    
    # ===== Private Helper Methods =====
    
    def _load_ligand_data(self, entity_info) -> Optional[Dict]:
        """Load ligand data from entity info."""
        file_path = self.paths.data_root / entity_info.file_path
        
        if file_path.suffix == '.sdf':
            # Read SDF file
            if HAS_RDKIT:
                try:
                    suppl = Chem.SDMolSupplier(str(file_path))
                    for mol in suppl:
                        if mol:
                            smiles = Chem.MolToSmiles(mol)
                            props = dict(mol.GetPropsAsDict())
                            return {
                                'smiles': smiles,
                                'properties': calculate_molecular_properties(smiles),
                                'sdf_properties': props
                            }
                except Exception as e:
                    logger.error(f"Failed to read SDF file {file_path}: {e}")
        
        # Try to reconstruct from metadata
        metadata = entity_info.metadata
        return {
            'smiles': entity_info.name,
            'chembl_id': metadata.get('chembl_id', ''),
            'properties': metadata.get('molecular_properties', {}),
            'targets': metadata.get('targets', [])
        }
    
    def _get_fingerprint(self, smiles: str):
        """Get molecular fingerprint for a SMILES string."""
        if not HAS_RDKIT:
            return None
        
        # Check cache
        if smiles in self._fingerprint_cache:
            return self._fingerprint_cache[smiles]
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = FingerprintMols.FingerprintMol(mol)
                self._fingerprint_cache[smiles] = fp
                return fp
        except Exception as e:
            logger.debug(f"Failed to generate fingerprint for {smiles}: {e}")
        
        return None

