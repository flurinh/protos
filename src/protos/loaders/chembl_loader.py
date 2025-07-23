"""
ChEMBL data loader for Protos.

This loader downloads ligand and bioactivity data from ChEMBL database,
following Protos' standardized loader patterns.
"""

import os
import json
import logging
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Union, Tuple
from datetime import datetime

from protos.io.paths.path_config import ProtosPaths
from protos.io.entity_registry import EntityRegistry
from protos.loaders.ligand_utils import (
    sanitize_smiles_filename, validate_smiles, smiles_to_inchi,
    calculate_molecular_properties, extract_protein_mapping,
    parse_activity_value, create_sdf_from_smiles
)

logger = logging.getLogger(__name__)

# Try to import ChEMBL client
try:
    from chembl_webresource_client.new_client import new_client as chembl_client
    HAS_CHEMBL = True
except ImportError:
    logger.warning("chembl_webresource_client not available. ChEMBL functionality will be limited.")
    HAS_CHEMBL = False
    chembl_client = None


class ChEMBLDL:
    """
    Downloads ligand data from ChEMBL following Protos loader patterns.
    
    This loader handles:
    - Protein target to ligand mapping
    - Activity data retrieval
    - Compound structure downloads
    - Automatic entity registration
    """
    
    def __init__(self, data_root: Optional[str] = None, reload: bool = False, 
                 limit: Optional[int] = None):
        """
        Initialize ChEMBL loader with ProtosPaths.
        
        Args:
            data_root: Root directory for data. If None, uses default.
            reload: Whether to reload data even if cached.
            limit: Maximum number of compounds to download per target.
        """
        # Initialize ProtosPaths
        self.paths = ProtosPaths(data_root=data_root)
        
        # Get standard directories
        self.ligand_dir = Path(self.paths.get_processor_path('ligand'))
        self.sdf_dir = self.ligand_dir / 'sdf'
        self.chembl_dir = self.ligand_dir / 'chembl'
        self.cache_dir = self.ligand_dir / 'cache'
        self.tables_dir = self.ligand_dir / 'tables'
        self.datasets_dir = self.ligand_dir / 'datasets'
        
        # Create directories
        for dir_path in [self.sdf_dir, self.chembl_dir, self.cache_dir, 
                        self.tables_dir, self.datasets_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize entity registry
        self.entity_registry = EntityRegistry(self.paths)
        
        # Settings
        self.reload = reload
        self.limit = limit
        
        # Initialize ChEMBL client if available
        if HAS_CHEMBL:
            self.chembl = chembl_client
        else:
            self.chembl = None
            
        # Cache for protein ID mappings
        self.protein_mapping_cache = {}
        self._load_protein_mapping_cache()
    
    def _load_protein_mapping_cache(self):
        """Load cached protein ID mappings."""
        cache_file = self.cache_dir / 'protein_mappings.json'
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    self.protein_mapping_cache = json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load protein mapping cache: {e}")
                self.protein_mapping_cache = {}
    
    def _save_protein_mapping_cache(self):
        """Save protein ID mappings to cache."""
        cache_file = self.cache_dir / 'protein_mappings.json'
        try:
            with open(cache_file, 'w') as f:
                json.dump(self.protein_mapping_cache, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save protein mapping cache: {e}")
    
    def map_protein_to_chembl_target(self, protein_id: str) -> Optional[str]:
        """
        Map protein identifier to ChEMBL target ID.
        
        Args:
            protein_id: Protein identifier (UniProt, PDB, gene name)
            
        Returns:
            ChEMBL target ID, or None if not found
        """
        if not self.chembl:
            logger.error("ChEMBL client not available")
            return None
        
        # Check cache first
        if protein_id in self.protein_mapping_cache:
            return self.protein_mapping_cache[protein_id]
        
        # Extract identifier type
        id_info = extract_protein_mapping(protein_id)
        
        try:
            # Search based on identifier type
            if id_info['type'] == 'uniprot_acc':
                targets = self.chembl.target.filter(
                    target_components__accession=id_info['uniprot_id']
                )
            elif id_info['type'] == 'uniprot_name':
                # Try gene name first
                gene = id_info.get('gene_name', '')
                targets = self.chembl.target.filter(pref_name__iexact=gene)
            elif id_info['type'] == 'pdb':
                # Search via cross-references
                targets = self.chembl.target.filter(
                    target_components__xrefs__xref_id=id_info['pdb_id']
                )
            elif id_info['type'] == 'chembl':
                # Already a ChEMBL ID
                return id_info['chembl_id']
            else:
                # Try as gene name
                targets = self.chembl.target.filter(pref_name__iexact=protein_id)
            
            # Get first valid target
            for target in targets:
                if target['target_chembl_id']:
                    chembl_id = target['target_chembl_id']
                    # Cache the mapping
                    self.protein_mapping_cache[protein_id] = chembl_id
                    self._save_protein_mapping_cache()
                    return chembl_id
                    
        except Exception as e:
            logger.error(f"Failed to map {protein_id} to ChEMBL: {e}")
        
        return None
    
    def download_protein_ligands(self, protein_ids: Union[str, List[str]], 
                               activity_types: Optional[List[str]] = None,
                               min_pchembl: float = 5.0,
                               max_value_nm: float = 10000,
                               save_sdf: bool = True) -> Dict[str, List[Dict]]:
        """
        Download ligands for given proteins from ChEMBL.
        
        Args:
            protein_ids: Protein identifier(s) (UniProt/PDB/gene name)
            activity_types: Types of bioactivities to retrieve (default: IC50, Ki, Kd)
            min_pchembl: Minimum pChEMBL value (negative log of activity)
            max_value_nm: Maximum activity value in nM
            save_sdf: Whether to save SDF files
            
        Returns:
            Dict mapping protein_id to list of ligand data
        """
        if not self.chembl:
            logger.error("ChEMBL client not available")
            return {}
        
        # Handle single protein
        if isinstance(protein_ids, str):
            protein_ids = [protein_ids]
        
        # Default activity types
        if activity_types is None:
            activity_types = ['IC50', 'Ki', 'Kd', 'EC50']
        
        results = {}
        
        for protein_id in protein_ids:
            logger.info(f"Downloading ligands for {protein_id}")
            
            # Check cache first
            cache_file = self.chembl_dir / f"{protein_id}_activities.json"
            if cache_file.exists() and not self.reload:
                logger.info(f"Loading cached data for {protein_id}")
                with open(cache_file, 'r') as f:
                    ligand_data = json.load(f)
                results[protein_id] = ligand_data
                continue
            
            # Map to ChEMBL target
            chembl_target = self.map_protein_to_chembl_target(protein_id)
            if not chembl_target:
                logger.warning(f"Could not map {protein_id} to ChEMBL target")
                continue
            
            try:
                # Query activities
                activities = self.chembl.activity.filter(
                    target_chembl_id=chembl_target,
                    type__in=activity_types,
                    pchembl_value__gte=min_pchembl,
                    assay_type='B'  # Binding assays
                )
                
                ligand_data = []
                count = 0
                
                for activity in activities:
                    # Apply limit if set
                    if self.limit and count >= self.limit:
                        break
                    
                    # Extract compound info
                    compound_id = activity.get('molecule_chembl_id')
                    if not compound_id:
                        continue
                    
                    # Get compound details
                    try:
                        compound = self.chembl.molecule.get(compound_id)
                        smiles = compound.get('molecule_structures', {}).get('canonical_smiles')
                        
                        if not smiles:
                            continue
                        
                        # Validate SMILES
                        is_valid, canonical_smiles = validate_smiles(smiles)
                        if not is_valid:
                            continue
                        
                        # Parse activity value
                        value = activity.get('value')
                        units = activity.get('units', 'nM')
                        value_nm = parse_activity_value(value, units)
                        
                        if value_nm and value_nm <= max_value_nm:
                            ligand_info = {
                                'smiles': canonical_smiles,
                                'chembl_id': compound_id,
                                'activity_type': activity.get('type'),
                                'value': value,
                                'units': units,
                                'value_nm': value_nm,
                                'pchembl_value': activity.get('pchembl_value'),
                                'assay_id': activity.get('assay_chembl_id'),
                                'target_id': chembl_target,
                                'protein_id': protein_id
                            }
                            
                            # Add molecular properties
                            props = calculate_molecular_properties(canonical_smiles)
                            if props:
                                ligand_info['properties'] = props
                            
                            # Add InChI identifiers
                            inchi_data = smiles_to_inchi(canonical_smiles)
                            if inchi_data:
                                ligand_info.update(inchi_data)
                            
                            ligand_data.append(ligand_info)
                            count += 1
                            
                            # Register ligand entity
                            self._register_ligand(ligand_info)
                            
                            # Save SDF if requested
                            if save_sdf:
                                self._save_ligand_sdf(ligand_info)
                                
                    except Exception as e:
                        logger.debug(f"Failed to process compound {compound_id}: {e}")
                        continue
                
                # Cache the results
                with open(cache_file, 'w') as f:
                    json.dump(ligand_data, f, indent=2)
                
                # Save activity table
                if ligand_data:
                    self.save_ligand_table(protein_id, ligand_data)
                
                results[protein_id] = ligand_data
                logger.info(f"Downloaded {len(ligand_data)} ligands for {protein_id}")
                
            except Exception as e:
                logger.error(f"Failed to download ligands for {protein_id}: {e}")
                results[protein_id] = []
        
        return results
    
    def save_ligand_table(self, protein_id: str, ligand_data: List[Dict], 
                         force_overwrite: bool = True) -> Optional[Path]:
        """
        Save ligand activity data to CSV table.
        
        Args:
            protein_id: Protein identifier
            ligand_data: List of ligand dictionaries
            force_overwrite: Whether to overwrite existing file
            
        Returns:
            Path to saved file, or None if no data
        """
        if not ligand_data:
            return None
        
        # Prepare DataFrame
        df_data = []
        for ligand in ligand_data:
            row = {
                'smiles': ligand['smiles'],
                'chembl_id': ligand['chembl_id'],
                'activity_type': ligand['activity_type'],
                'value': ligand['value'],
                'units': ligand['units'],
                'value_nm': ligand['value_nm'],
                'pchembl_value': ligand.get('pchembl_value', ''),
                'mw': ligand.get('properties', {}).get('mw', ''),
                'logp': ligand.get('properties', {}).get('logp', ''),
                'hba': ligand.get('properties', {}).get('hba', ''),
                'hbd': ligand.get('properties', {}).get('hbd', ''),
                'tpsa': ligand.get('properties', {}).get('tpsa', ''),
                'inchi_key': ligand.get('inchi_key', '')
            }
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        # Save to CSV
        table_path = self.tables_dir / f"{protein_id}_chembl_activities.csv"
        if table_path.exists() and not force_overwrite:
            logger.warning(f"Table already exists: {table_path}")
            return None
        
        df.to_csv(table_path, index=False)
        logger.info(f"Saved activity table: {table_path}")
        
        # Create dataset entry
        self._create_ligand_dataset(protein_id, ligand_data)
        
        return table_path
    
    def _register_ligand(self, ligand_info: Dict):
        """
        Register ligand in entity system with SMILES as primary ID.
        
        Args:
            ligand_info: Dictionary with ligand information
        """
        smiles = ligand_info['smiles']
        safe_filename = sanitize_smiles_filename(smiles)
        
        # Prepare metadata
        metadata = {
            'chembl_id': ligand_info['chembl_id'],
            'targets': [ligand_info['protein_id']],
            'activity_types': [ligand_info['activity_type']],
            'min_value_nm': ligand_info['value_nm']
        }
        
        # Add aliases
        aliases = [ligand_info['chembl_id']]
        if 'inchi_key' in ligand_info:
            aliases.append(ligand_info['inchi_key'])
        metadata['aliases'] = aliases
        
        # Add properties if available
        if 'properties' in ligand_info:
            metadata['molecular_properties'] = ligand_info['properties']
        
        # Register entity
        try:
            self.entity_registry.register_entity(
                name=smiles,
                format_type='ligand',
                file_path=f"ligand/sdf/{safe_filename}.sdf",
                metadata=metadata
            )
        except Exception as e:
            logger.debug(f"Failed to register ligand {smiles}: {e}")
    
    def _save_ligand_sdf(self, ligand_info: Dict) -> Optional[Path]:
        """
        Save ligand structure as SDF file.
        
        Args:
            ligand_info: Dictionary with ligand information
            
        Returns:
            Path to saved file, or None if failed
        """
        smiles = ligand_info['smiles']
        safe_filename = sanitize_smiles_filename(smiles)
        sdf_path = self.sdf_dir / f"{safe_filename}.sdf"
        
        # Skip if already exists
        if sdf_path.exists() and not self.reload:
            return sdf_path
        
        # Create SDF with properties
        properties = {
            'CHEMBL_ID': ligand_info['chembl_id'],
            'ACTIVITY_TYPE': ligand_info['activity_type'],
            'VALUE_NM': str(ligand_info['value_nm']),
            'TARGET': ligand_info['protein_id']
        }
        
        sdf_content = create_sdf_from_smiles(smiles, properties)
        if sdf_content:
            try:
                with open(sdf_path, 'w') as f:
                    f.write(sdf_content)
                return sdf_path
            except Exception as e:
                logger.error(f"Failed to save SDF for {smiles}: {e}")
        
        return None
    
    def _create_ligand_dataset(self, protein_id: str, ligand_data: List[Dict]):
        """
        Create dataset entry for protein's ligands.
        
        Args:
            protein_id: Protein identifier
            ligand_data: List of ligand dictionaries
        """
        dataset_name = f"{protein_id}_chembl_ligands"
        dataset_path = self.datasets_dir / f"{dataset_name}.json"
        
        # Extract entity names (SMILES)
        entity_names = [ligand['smiles'] for ligand in ligand_data]
        
        # Create dataset metadata
        dataset_info = {
            'name': dataset_name,
            'description': f"ChEMBL ligands for {protein_id}",
            'entities': entity_names,
            'entity_count': len(entity_names),
            'data_file': f"../tables/{protein_id}_chembl_activities.csv",
            'metadata': {
                'protein_id': protein_id,
                'source': 'ChEMBL',
                'activity_types': list(set(l['activity_type'] for l in ligand_data)),
                'created': datetime.now().isoformat(),
                'min_pchembl': min(l.get('pchembl_value', 0) for l in ligand_data if l.get('pchembl_value')),
                'max_pchembl': max(l.get('pchembl_value', 0) for l in ligand_data if l.get('pchembl_value'))
            }
        }
        
        # Save dataset JSON
        with open(dataset_path, 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        logger.info(f"Created dataset: {dataset_name}")
    
    def download_compound_structures(self, chembl_ids: Union[str, List[str]], 
                                   format: str = 'sdf') -> Dict[str, Path]:
        """
        Download molecular structures from ChEMBL.
        
        Args:
            chembl_ids: ChEMBL compound ID(s)
            format: Output format ('sdf' or 'mol')
            
        Returns:
            Dictionary mapping ChEMBL ID to file path
        """
        if not self.chembl:
            logger.error("ChEMBL client not available")
            return {}
        
        # Handle single ID
        if isinstance(chembl_ids, str):
            chembl_ids = [chembl_ids]
        
        results = {}
        
        for chembl_id in chembl_ids:
            try:
                # Get compound
                compound = self.chembl.molecule.get(chembl_id)
                smiles = compound.get('molecule_structures', {}).get('canonical_smiles')
                
                if not smiles:
                    logger.warning(f"No structure found for {chembl_id}")
                    continue
                
                # Create ligand info
                ligand_info = {
                    'smiles': smiles,
                    'chembl_id': chembl_id,
                    'activity_type': 'unknown',
                    'value_nm': 0,
                    'protein_id': 'unknown'
                }
                
                # Calculate properties
                props = calculate_molecular_properties(smiles)
                if props:
                    ligand_info['properties'] = props
                
                # Register and save
                self._register_ligand(ligand_info)
                sdf_path = self._save_ligand_sdf(ligand_info)
                
                if sdf_path:
                    results[chembl_id] = sdf_path
                    
            except Exception as e:
                logger.error(f"Failed to download structure for {chembl_id}: {e}")
        
        return results
    
    def search_similar_compounds(self, smiles: str, similarity: float = 0.8, 
                               limit: int = 100) -> List[Dict]:
        """
        Search for similar compounds in ChEMBL.
        
        Args:
            smiles: Query SMILES string
            similarity: Similarity threshold (0-1)
            limit: Maximum number of results
            
        Returns:
            List of similar compound dictionaries
        """
        if not self.chembl:
            logger.error("ChEMBL client not available")
            return []
        
        try:
            # Validate SMILES
            is_valid, canonical_smiles = validate_smiles(smiles)
            if not is_valid:
                logger.error(f"Invalid SMILES: {smiles}")
                return []
            
            # Search similar molecules
            similar = self.chembl.similarity.filter(
                smiles=canonical_smiles,
                similarity=int(similarity * 100)  # ChEMBL uses percentage
            ).only(['molecule_chembl_id', 'similarity', 'molecule_structures'])[:limit]
            
            results = []
            for mol in similar:
                mol_smiles = mol.get('molecule_structures', {}).get('canonical_smiles')
                if mol_smiles:
                    result = {
                        'chembl_id': mol['molecule_chembl_id'],
                        'smiles': mol_smiles,
                        'similarity': mol['similarity'] / 100.0
                    }
                    
                    # Calculate properties
                    props = calculate_molecular_properties(mol_smiles)
                    if props:
                        result['properties'] = props
                    
                    results.append(result)
            
            return results
            
        except Exception as e:
            logger.error(f"Similarity search failed: {e}")
            return []