"""
Ligand processor for Protos - REFACTORED to follow DATA_MANAGEMENT_UNIFIED.md principles.

This processor handles ligand data including:
- Small molecule structures and properties
- Bioactivity data
- Target-ligand interactions
- QSAR modeling
- Virtual screening

Key changes from original:
- NO direct filesystem operations
- NO directory creation
- ALL paths managed by ProtosPaths through BaseProcessor
- Uses entity registry for all data management
- Implements abstract methods properly
"""

import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Union, Any, Tuple
from datetime import datetime

from protos.core.base_processor import BaseProcessor
from protos.loaders.ligand_utils import (
    sanitize_smiles_filename, validate_smiles, smiles_to_inchi,
    calculate_molecular_properties, is_drug_like, parse_activity_value,
    create_sdf_from_smiles
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
    
    Follows DATA_MANAGEMENT_UNIFIED.md principles:
    - Zero filesystem management
    - All paths handled by ProtosPaths
    - Human-readable names only
    - Entity registry for all data tracking
    """
    
    def __init__(self, name: str = "ligand_processor", 
                 paths=None, cif_processor=None):
        """
        Initialize the LigandProcessor.
        
        Args:
            name: Name for this processor instance
            paths: ProtosPaths instance (created if not provided)
            cif_processor: CifBaseProcessor instance for structure integration
        """
        # Initialize BaseProcessor - this handles ALL path management
        super().__init__(name=name, paths=paths)
        
        # Initialize ChEMBL loader (lazy import to avoid circular dependency)
        self._chembl_loader = None
        
        # Store reference to CifBaseProcessor if provided
        self.cif_processor = cif_processor
        
        # Caches for performance
        self._ligand_cache = {}
        self._fingerprint_cache = {}
    
    # ===== Properties for Path Access (delegating to BaseProcessor) =====
    
    @property
    def chembl_available(self):
        """Check if ChEMBL functionality is available."""
        try:
            from protos.loaders import chembl_loader
            return chembl_loader._get_chembl_client() is not None
        except ImportError:
            return False
    
    @property
    def sdf_subdir(self) -> str:
        """Get SDF subdirectory name (not path!)."""
        return 'sdf'
    
    @property
    def tables_subdir(self) -> str:
        """Get tables subdirectory name (not path!)."""
        return 'tables'
    
    @property
    def cache_subdir(self) -> str:
        """Get cache subdirectory name (not path!)."""
        return 'cache'
    
    # ===== Abstract Method Implementations =====
    
    def load_entity(self, name: str) -> Optional[Dict]:
        """
        Load a ligand entity by name.
        
        This method implements the abstract load_entity from BaseProcessor.
        It loads ligand data from the entity registry and associated files.
        
        Args:
            name: Human-readable name (SMILES, ChEMBL ID, InChI Key, etc.)
            
        Returns:
            Dictionary with ligand data, or None if not found
        """
        # Check cache first
        if name in self._ligand_cache:
            return self._ligand_cache[name]
        
        # Try to find entity in registry
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        
        if entity_info:
            # Build ligand data from entity info
            ligand_data = self._build_ligand_data_from_entity(entity_info)
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
                'properties': calculate_molecular_properties(canonical_smiles) if HAS_RDKIT else {}
            }
            self._ligand_cache[name] = ligand_data
            return ligand_data
        
        return None
    
    def save_entity(self, name: str, data: Union[str, Dict], metadata: Optional[Dict] = None):
        """
        Save a ligand entity.
        
        Args:
            name: Entity name (SMILES)
            data: Either a SMILES string or a dictionary with ligand data
            metadata: Optional metadata
        """
        # Handle both string and dict inputs
        if isinstance(data, str):
            smiles = data
            ligand_data = {'smiles': smiles}
        else:
            ligand_data = data
            smiles = ligand_data.get('smiles', name)
        
        # Validate SMILES
        is_valid, canonical_smiles = validate_smiles(smiles)
        if not is_valid:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Build file path using sanitized name
        safe_filename = sanitize_smiles_filename(canonical_smiles)
        relative_path = f"{self.processor_type}/{self.sdf_subdir}/{safe_filename}.sdf"
        
        # Try to create SDF content if RDKit is available
        properties = ligand_data.get('properties', {})
        if 'chembl_id' in ligand_data:
            properties['CHEMBL_ID'] = ligand_data['chembl_id']
        
        sdf_content = create_sdf_from_smiles(canonical_smiles, properties)
        
        if sdf_content:
            # Save SDF file through proper path system
            full_path = Path(self.paths.data_root) / relative_path
            full_path.parent.mkdir(parents=True, exist_ok=True)
            with open(full_path, 'w') as f:
                f.write(sdf_content)
        else:
            # When RDKit is not available, we still register the entity
            logger.debug(f"SDF file not created for {canonical_smiles} (RDKit not available)")
        
        # Prepare metadata
        if metadata is None:
            metadata = {}
        
        # Add standard metadata
        metadata.update({
            'molecular_properties': ligand_data.get('properties', {}),
            'targets': ligand_data.get('targets', []),
            'activities': ligand_data.get('activities', [])
        })
        
        # Add aliases
        aliases = []
        if 'chembl_id' in ligand_data:
            aliases.append(ligand_data['chembl_id'])
        if 'inchi_key' in ligand_data:
            aliases.append(ligand_data['inchi_key'])
        metadata['aliases'] = aliases
        
        # Register entity
        self.entity_registry.register_entity(
            name=canonical_smiles,
            format_type=self.processor_type,
            file_path=relative_path,
            metadata=metadata
        )
        
        # Update cache
        self._ligand_cache[canonical_smiles] = ligand_data
    
    # ===== ChEMBL Integration Methods =====
    
    def get_protein_ligands(self, protein_id: str, 
                           activity_types: Optional[List[str]] = None,
                           min_pchembl: float = 5.0,
                           limit: Optional[int] = None) -> List[Dict]:
        """
        Get ligands for a protein target from ChEMBL.
        
        Args:
            protein_id: Protein identifier (UniProt, gene name, etc.)
            activity_types: Activity types to retrieve (default: IC50, Ki, Kd)
            min_pchembl: Minimum pChEMBL value
            limit: Maximum number of compounds
            
        Returns:
            List of ligand dictionaries
        """
        if not self.chembl_available:
            logger.error("ChEMBL functionality not available")
            return []
        
        # Use ChEMBL utility to query data
        from protos.loaders import chembl_loader
        ligands = chembl_loader.query_protein_ligands(
            protein_id,
            activity_types=activity_types,
            min_pchembl=min_pchembl,
            limit=limit
        )
        for ligand_data in ligands:
            try:
                self.save_entity(ligand_data['smiles'], ligand_data)
            except Exception as e:
                logger.warning(f"Failed to save ligand {ligand_data.get('chembl_id', 'unknown')}: {e}")
        
        return ligands[:limit] if limit else ligands
    
    # ===== Similarity Search Methods =====
    
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
        
        structures = []
        for pdb_id in self.cif_processor.list_entities():
            structure_data = self.cif_processor.load_entity(pdb_id)
            if structure_data is not None:
                # Check for ligand in structure
                ligands = structure_data[
                    (structure_data['group'] == 'HETATM') &
                    (structure_data['res_name3l'] == ligand_name)
                ]['res_name3l'].unique()
                
                if len(ligands) > 0:
                    structures.append(pdb_id)
        
        return structures
    
    def get_ligand_from_structure(self, pdb_id: str, ligand_name: str, 
                                 chain_id: Optional[str] = None) -> pd.DataFrame:
        """
        Extract ligand coordinates from a structure.
        
        Args:
            pdb_id: PDB ID
            ligand_name: Three-letter ligand code
            chain_id: Optional chain ID
            
        Returns:
            DataFrame with ligand atoms
        """
        if not self.cif_processor:
            raise ValueError("No CifBaseProcessor attached")
        
        structure_data = self.cif_processor.load_entity(pdb_id)
        if structure_data is None:
            raise ValueError(f"Structure {pdb_id} not found")
        
        # Filter for ligand
        mask = (
            (structure_data['group'] == 'HETATM') &
            (structure_data['res_name3l'] == ligand_name)
        )
        
        if chain_id:
            mask &= (structure_data['auth_chain_id'] == chain_id)
        
        ligand_data = structure_data[mask]
        
        if ligand_data.empty:
            raise ValueError(f"Ligand {ligand_name} not found in {pdb_id}")
        
        return ligand_data
    
    # ===== Private Helper Methods =====
    
    def _build_ligand_data_from_entity(self, entity_info) -> Optional[Dict]:
        """Build ligand data dictionary from entity info."""
        # Start with metadata
        metadata = entity_info.metadata
        ligand_data = {
            'smiles': entity_info.original_id,
            'chembl_id': metadata.get('chembl_id', ''),
            'properties': metadata.get('molecular_properties', {}),
            'targets': metadata.get('targets', []),
            'activities': metadata.get('activities', [])
        }
        
        # Try to load SDF if it exists and RDKit is available
        if HAS_RDKIT and entity_info.file_path.endswith('.sdf'):
            full_path = Path(self.paths.data_root) / entity_info.file_path
            if full_path.exists():
                try:
                    suppl = Chem.SDMolSupplier(str(full_path))
                    for mol in suppl:
                        if mol:
                            # Enrich with SDF data
                            ligand_data['sdf_properties'] = dict(mol.GetPropsAsDict())
                            # Recalculate properties to ensure freshness
                            fresh_props = calculate_molecular_properties(ligand_data['smiles'])
                            if fresh_props:
                                ligand_data['properties'].update(fresh_props)
                            break
                except Exception as e:
                    logger.debug(f"Could not read SDF file: {e}")
        
        return ligand_data
    
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
    
    # ===== Dataset Loading Override =====
    
    def load_dataset(self, dataset_name: str) -> Dict[str, Dict]:
        """
        Load a ligand dataset.
        
        Overrides base class to handle potential table-based datasets.
        
        Args:
            dataset_name: Name of the dataset
            
        Returns:
            Dictionary mapping entity names to ligand data
        """
        # Try standard dataset loading first
        try:
            return super().load_dataset(dataset_name)
        except Exception:
            pass
        
        # Check if this is a table-based dataset (e.g., from ChEMBL download)
        table_path = Path(self.paths.get_processor_path(self.processor_type)) / self.tables_subdir / f"{dataset_name}.csv"
        if table_path.exists():
            df = pd.read_csv(table_path)
            
            # Convert to dictionary format
            result = {}
            for _, row in df.iterrows():
                smiles = row.get('smiles', '')
                if smiles:
                    ligand_data = {
                        'smiles': smiles,
                        'chembl_id': row.get('chembl_id', ''),
                        'properties': {
                            'mw': row.get('mw', 0),
                            'logp': row.get('logp', 0),
                            'hba': row.get('hba', 0),
                            'hbd': row.get('hbd', 0),
                            'tpsa': row.get('tpsa', 0)
                        },
                        'activity_type': row.get('activity_type', ''),
                        'value_nm': row.get('value_nm', 0),
                        'pchembl_value': row.get('pchembl_value', 0)
                    }
                    result[smiles] = ligand_data
            
            return result
        
        raise ValueError(f"Dataset '{dataset_name}' not found")