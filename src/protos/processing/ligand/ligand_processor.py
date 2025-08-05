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
    
    @property
    def databases_subdir(self) -> str:
        """Get databases subdirectory name (not path!)."""
        return 'databases'
    
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
    
    # ===== SDF File Operations =====
    
    def load_sdf_file(self, sdf_name: str, as_entities: bool = True) -> Union[List[Dict], pd.DataFrame]:
        """
        Load an SDF file.
        
        Args:
            sdf_name: Name of the SDF file (without .sdf extension)
            as_entities: Register molecules as entities (default: True)
            
        Returns:
            List of molecule dictionaries or DataFrame
        """
        # Build path to SDF file
        sdf_path = Path(self.paths.get_processor_path(self.processor_type)) / self.sdf_subdir / f"{sdf_name}.sdf"
        
        if not sdf_path.exists():
            raise FileNotFoundError(f"SDF file not found: {sdf_path}")
        
        # Use format registry to read
        from protos.io.formats import FormatRegistry
        registry = FormatRegistry()
        
        # Read as list of dictionaries
        molecules = registry.read(str(sdf_path), format_type='sdf')
        
        # Register as entities if requested
        if as_entities:
            for mol_data in molecules:
                if 'smiles' in mol_data:
                    try:
                        self.save_entity(mol_data['smiles'], mol_data)
                    except Exception as e:
                        logger.warning(f"Failed to register molecule: {e}")
        
        return molecules
    
    def save_sdf_file(self, sdf_name: str, data: Union[List[str], List[Dict], pd.DataFrame],
                      include_properties: bool = True) -> str:
        """
        Save molecules to an SDF file.
        
        Args:
            sdf_name: Name for the SDF file (without extension)
            data: Can be:
                - List of SMILES strings
                - List of entity names to load
                - List of molecule dictionaries
                - DataFrame with 'smiles' column
            include_properties: Include molecular properties
            
        Returns:
            Path to the created SDF file
        """
        # Convert data to appropriate format
        if isinstance(data, list) and data and isinstance(data[0], str):
            # List of SMILES or entity names
            molecules = []
            for item in data:
                # Try to load as entity first
                entity_data = self.load_entity(item)
                if entity_data:
                    molecules.append(entity_data)
                else:
                    # Treat as SMILES
                    molecules.append({'smiles': item})
        elif isinstance(data, pd.DataFrame):
            # DataFrame - will be handled by SDFHandler
            molecules = data
        else:
            # Assume list of dictionaries
            molecules = data
        
        # Build output path
        sdf_path = Path(self.paths.get_processor_path(self.processor_type)) / self.sdf_subdir / f"{sdf_name}.sdf"
        sdf_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Use format registry to write
        from protos.io.formats import FormatRegistry
        registry = FormatRegistry()
        
        registry.write(str(sdf_path), molecules, format_type='sdf')
        
        logger.info(f"Saved {len(molecules)} molecules to {sdf_path}")
        return str(sdf_path)
    
    def convert_to_structure_format(self, ligand_name: str, 
                                   output_format: str = 'cif',
                                   chain_id: str = 'L',
                                   res_name: Optional[str] = None) -> Optional[str]:
        """
        Convert ligand to structure format for CifBaseProcessor compatibility.
        
        Args:
            ligand_name: Ligand entity name or SMILES
            output_format: Output format ('cif', 'pdb', or 'mol2') - default is 'cif'
            chain_id: Chain identifier (default: 'L')
            res_name: Residue name (optional)
            
        Returns:
            Path to the converted file, or None if conversion failed
        """
        if not HAS_RDKIT:
            logger.error("RDKit required for structure conversion")
            return None
        
        # Load ligand data
        ligand_data = self.load_entity(ligand_name)
        if not ligand_data:
            logger.error(f"Ligand {ligand_name} not found")
            return None
        
        smiles = ligand_data.get('smiles')
        if not smiles:
            logger.error("No SMILES found for ligand")
            return None
        
        try:
            # Create 3D structure
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Failed to parse SMILES: {smiles}")
                return None
            
            # Add hydrogens and generate 3D coordinates
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Save to file
            safe_name = sanitize_smiles_filename(smiles)
            if output_format == 'cif':
                output_path = Path(self.paths.get_processor_path(self.processor_type)) / 'structures' / f"{safe_name}.cif"
                output_path.parent.mkdir(parents=True, exist_ok=True)
                # Use our CIF conversion
                return self.save_as_cif(ligand_name, str(output_path), chain_id=chain_id, res_name=res_name)
            elif output_format == 'pdb':
                output_path = Path(self.paths.get_processor_path(self.processor_type)) / 'structures' / f"{safe_name}.pdb"
                output_path.parent.mkdir(parents=True, exist_ok=True)
                Chem.MolToPDBFile(mol, str(output_path))
            elif output_format == 'mol2':
                output_path = Path(self.paths.get_processor_path(self.processor_type)) / 'structures' / f"{safe_name}.mol2"
                output_path.parent.mkdir(parents=True, exist_ok=True)
                Chem.MolToMol2File(mol, str(output_path))
            else:
                logger.error(f"Unsupported output format: {output_format}")
                return None
            
            logger.info(f"Converted ligand to {output_format}: {output_path}")
            return str(output_path)
            
        except Exception as e:
            logger.error(f"Failed to convert ligand: {e}")
            return None
    
    def convert_to_cif_dataframe(self, ligand_name: str, 
                                 chain_id: str = 'L',
                                 res_name: Optional[str] = None,
                                 res_id: int = 1) -> Optional[pd.DataFrame]:
        """
        Convert ligand to CIF DataFrame format compatible with StructureProcessor.
        
        This creates a DataFrame with the same columns as used by CifBaseProcessor,
        allowing seamless integration between ligand and structure data.
        
        Args:
            ligand_name: Ligand entity name or SMILES
            chain_id: Chain identifier for the ligand (default: 'L')
            res_name: Residue name (default: 'LIG' or from ChEMBL ID)
            res_id: Residue sequence number (default: 1)
            
        Returns:
            DataFrame in CIF format, or None if conversion failed
        """
        if not HAS_RDKIT:
            logger.error("RDKit required for structure conversion")
            return None
        
        # Load ligand data
        ligand_data = self.load_entity(ligand_name)
        if not ligand_data:
            logger.error(f"Ligand {ligand_name} not found")
            return None
        
        smiles = ligand_data.get('smiles')
        if not smiles:
            logger.error("No SMILES found for ligand")
            return None
        
        try:
            # Create 3D structure
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Failed to parse SMILES: {smiles}")
                return None
            
            # Add hydrogens and generate 3D coordinates
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Determine residue name
            if res_name is None:
                if 'chembl_id' in ligand_data and ligand_data['chembl_id']:
                    # Use last 3 chars of ChEMBL ID (e.g., CHEMBL25 -> L25)
                    res_name = ligand_data['chembl_id'][-3:]
                else:
                    res_name = 'LIG'
            
            # Build CIF DataFrame
            rows = []
            conf = mol.GetConformer()
            
            for atom_idx in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(atom_idx)
                pos = conf.GetAtomPosition(atom_idx)
                
                # Create row matching CIF format
                row = {
                    'group': 'HETATM',
                    'atom_id': atom_idx + 1,
                    'atom_name': f"{atom.GetSymbol()}{atom_idx + 1}",
                    'res_name': res_name,
                    'auth_chain_id': chain_id,
                    'auth_seq_id': res_id,
                    'x': round(pos.x, 3),
                    'y': round(pos.y, 3),
                    'z': round(pos.z, 3),
                    'element': atom.GetSymbol(),
                    'occupancy': 1.00,
                    'b_factor': 30.00,
                    'entity_id': 2,  # Typically 2 for ligands (1 is protein)
                    'model_num': 1,
                    'res_name3l': res_name,
                    'res_id': f"{res_name}_{res_id}_{chain_id}",
                    'label_chain_id': chain_id,
                    'label_seq_id': 1,
                    'alt_id': '.',
                    'charge': atom.GetFormalCharge() if atom.GetFormalCharge() != 0 else None
                }
                
                rows.append(row)
            
            # Create DataFrame
            df = pd.DataFrame(rows)
            
            # Add any additional metadata
            if 'pdb_id' in ligand_data:
                df['pdb_id'] = ligand_data['pdb_id']
            
            logger.info(f"Converted ligand to CIF DataFrame: {len(df)} atoms")
            return df
            
        except Exception as e:
            logger.error(f"Failed to convert ligand to CIF format: {e}")
            return None
    
    def save_as_cif(self, ligand_name: str, output_path: str,
                    chain_id: str = 'L', res_name: Optional[str] = None) -> Optional[str]:
        """
        Save ligand as a CIF file using the standard CIF handler.
        
        Args:
            ligand_name: Ligand entity name or SMILES
            output_path: Path for output CIF file
            chain_id: Chain identifier (default: 'L')
            res_name: Residue name (optional)
            
        Returns:
            Path to saved file, or None if failed
        """
        # Convert to CIF DataFrame
        df = self.convert_to_cif_dataframe(ligand_name, chain_id=chain_id, res_name=res_name)
        if df is None:
            return None
        
        try:
            # Use CIF handler to write
            from protos.io.cif_handler import CifHandler
            handler = CifHandler()
            
            # Write the file with force_overwrite
            written_path = handler.write_with_versioning(
                output_path, df, 
                versioned=False, 
                force_overwrite=True
            )
            logger.info(f"Saved ligand as CIF: {written_path}")
            return written_path
            
        except Exception as e:
            logger.error(f"Failed to save CIF file: {e}")
            return None
    
    # ===== Database Integration Methods =====
    
    def get_ccd_ligand(self, ccd_id: str, download_if_missing: bool = True) -> Optional[Dict]:
        """
        Get ligand from PDB Chemical Component Dictionary.
        
        This accesses the local CCD database, downloading it first if needed.
        The CCD contains all small molecules found in PDB structures.
        
        Args:
            ccd_id: Three-letter CCD code (e.g., 'ATP', 'NAD', 'HEM')
            download_if_missing: Auto-download CCD if not present
            
        Returns:
            Dictionary with ligand data, or None if not found
            
        Example:
            >>> atp = processor.get_ccd_ligand('ATP')
            >>> print(atp['smiles'])
            >>> print(atp['name'])
        """
        from protos.loaders import ccd_loader
        
        # Get database path using ProtosPaths
        db_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'ccd'
        
        # Check if CCD is downloaded
        if not ccd_loader.is_ccd_downloaded(db_path):
            if download_if_missing:
                logger.info("CCD not found. Downloading PDB Chemical Component Dictionary...")
                db_path.mkdir(parents=True, exist_ok=True)
                ccd_loader.download_ccd_components(db_path)
            else:
                logger.error("CCD database not found. Set download_if_missing=True to download.")
                return None
        
        # Load component
        component = ccd_loader.load_ccd_component(db_path, ccd_id)
        if not component:
            return None
        
        # Get SMILES
        smiles = ccd_loader.get_ccd_smiles(db_path, ccd_id)
        if not smiles:
            logger.warning(f"No SMILES found for CCD component {ccd_id}")
            return None
        
        # Build ligand data
        ligand_data = {
            'smiles': smiles,
            'ccd_id': ccd_id,
            'name': component.get('name', ccd_id),
            'formula': component.get('formula', ''),
            'type': component.get('type', 'non-polymer'),
            'source': 'CCD',
            'properties': calculate_molecular_properties(smiles) if HAS_RDKIT else {}
        }
        
        # Add InChI if available
        if 'inchi' in component:
            ligand_data['inchi'] = component['inchi']
        if 'inchi_key' in component:
            ligand_data['inchi_key'] = component['inchi_key']
        
        return ligand_data
    
    def search_qm9_by_properties(self, property_filters: Dict[str, Tuple[float, float]], 
                                 limit: Optional[int] = None,
                                 download_if_missing: bool = True) -> List[Dict]:
        """
        Search QM9 dataset by quantum chemical properties.
        
        QM9 contains ~134k small organic molecules with DFT-calculated properties.
        
        Args:
            property_filters: Dict mapping property names to (min, max) ranges
                Available properties: 'gap', 'homo', 'lumo', 'dipole', 'alpha', 
                                    'zpve', 'cv', 'u0', 'u298', 'h298', 'g298'
            limit: Maximum number of results
            download_if_missing: Auto-download QM9 if not present
            
        Returns:
            List of ligand dictionaries matching criteria
            
        Example:
            >>> # Find molecules with small HOMO-LUMO gap
            >>> molecules = processor.search_qm9_by_properties({
            ...     'gap': (0.1, 0.3),
            ...     'dipole': (0, 2.0)
            ... }, limit=100)
        """
        from protos.loaders import qm9_loader
        
        # Get database path using ProtosPaths
        db_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'qm9'
        
        # Check if QM9 is downloaded
        if not qm9_loader.is_qm9_downloaded(db_path):
            if download_if_missing:
                logger.info("QM9 not found. Downloading QM9 quantum chemistry dataset...")
                db_path.mkdir(parents=True, exist_ok=True)
                qm9_loader.download_qm9_dataset(db_path)
            else:
                logger.error("QM9 database not found. Set download_if_missing=True to download.")
                return []
        
        # Search by each property
        results = None
        for prop_name, (min_val, max_val) in property_filters.items():
            matches = qm9_loader.search_qm9_by_property(db_path, prop_name, min_val, max_val)
            
            if results is None:
                results = matches
            else:
                # Intersection of results
                result_ids = {m['id'] for m in results}
                match_ids = {m['id'] for m in matches}
                keep_ids = result_ids & match_ids
                results = [m for m in results if m['id'] in keep_ids]
        
        if results is None:
            results = []
        
        # Apply limit
        if limit and len(results) > limit:
            results = results[:limit]
        
        # Convert to ligand format
        ligands = []
        for mol in results:
            ligand_data = {
                'smiles': mol['smiles'],
                'qm9_id': mol['id'],
                'source': 'QM9',
                'quantum_properties': mol['properties'],
                'xyz_coordinates': mol.get('xyz', None)
            }
            
            # Add RDKit properties if available
            if HAS_RDKIT:
                ligand_data['properties'] = calculate_molecular_properties(mol['smiles'])
            
            ligands.append(ligand_data)
        
        return ligands
    
    def get_qm9_molecule(self, qm9_id: int, download_if_missing: bool = True) -> Optional[Dict]:
        """
        Get specific molecule from QM9 by ID.
        
        Args:
            qm9_id: QM9 molecule ID (1-133885)
            download_if_missing: Auto-download QM9 if not present
            
        Returns:
            Ligand dictionary with quantum properties, or None if not found
        """
        from protos.loaders import qm9_loader
        
        # Get database path using ProtosPaths
        db_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'qm9'
        
        # Check if QM9 is downloaded
        if not qm9_loader.is_qm9_downloaded(db_path):
            if download_if_missing:
                logger.info("QM9 not found. Downloading QM9 quantum chemistry dataset...")
                db_path.mkdir(parents=True, exist_ok=True)
                qm9_loader.download_qm9_dataset(db_path)
            else:
                logger.error("QM9 database not found. Set download_if_missing=True to download.")
                return None
        
        # Load molecule
        mol_data = qm9_loader.load_qm9_molecule(db_path, qm9_id)
        if not mol_data:
            return None
        
        # Convert to ligand format
        ligand_data = {
            'smiles': mol_data['smiles'],
            'qm9_id': qm9_id,
            'source': 'QM9',
            'quantum_properties': mol_data['properties'],
            'xyz_coordinates': mol_data.get('xyz', None)
        }
        
        # Add RDKit properties if available
        if HAS_RDKIT:
            ligand_data['properties'] = calculate_molecular_properties(mol_data['smiles'])
        
        return ligand_data
    
    def search_enamine_by_similarity(self, query_smiles: str, 
                                    dataset: str = 'diversity_1k',
                                    similarity: float = 0.7,
                                    limit: Optional[int] = 100,
                                    download_if_missing: bool = True) -> List[Tuple[str, float, Dict]]:
        """
        Search Enamine REAL database by structural similarity.
        
        Note: Enamine is a commercial database. Credentials must be set via
        environment variables (ENAMINE_USERNAME, ENAMINE_PASSWORD).
        
        Args:
            query_smiles: Query molecule SMILES
            dataset: Dataset name (e.g., 'diversity_1k', 'hit2lead_10k')
                    See enamine_loader.list_available_datasets() for options
            similarity: Tanimoto similarity threshold (0-1)
            limit: Maximum number of results
            download_if_missing: Auto-download dataset if not present
            
        Returns:
            List of (smiles, similarity_score, properties) tuples
            
        Example:
            >>> # Search for molecules similar to aspirin
            >>> similar = processor.search_enamine_by_similarity(
            ...     "CC(=O)Oc1ccccc1C(=O)O",
            ...     dataset='diversity_1k',  # Small test dataset
            ...     similarity=0.7
            ... )
        """
        if not HAS_RDKIT:
            logger.error("RDKit required for similarity search")
            return []
        
        from protos.loaders import enamine_loader
        
        # Get database path using ProtosPaths
        db_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'enamine'
        
        # Check if dataset is downloaded
        if not enamine_loader.is_dataset_downloaded(db_path, dataset):
            if download_if_missing:
                logger.info(f"Downloading Enamine dataset '{dataset}'...")
                db_path.mkdir(parents=True, exist_ok=True)
                if not enamine_loader.download_enamine_dataset(db_path, dataset):
                    logger.error(f"Failed to download dataset '{dataset}'")
                    return []
            else:
                logger.error(f"Enamine dataset '{dataset}' not found. Set download_if_missing=True to download.")
                return []
        
        # Get query fingerprint
        query_fp = self._get_fingerprint(query_smiles)
        if query_fp is None:
            return []
        
        # Search database
        results = []
        for compound in enamine_loader.stream_enamine_compounds(db_path, dataset):
            comp_fp = self._get_fingerprint(compound['smiles'])
            if comp_fp is not None:
                sim = TanimotoSimilarity(query_fp, comp_fp)
                if sim >= similarity:
                    results.append((compound['smiles'], sim, compound))
                    if limit and len(results) >= limit:
                        break
        
        # Sort by similarity
        results.sort(key=lambda x: x[1], reverse=True)
        
        return results
    
    def create_ccd_dataset(self, dataset_name: str, ccd_ids: List[str],
                          download_if_missing: bool = True) -> List[str]:
        """
        Create a dataset from CCD components.
        
        This is useful for working with known PDB ligands.
        
        Args:
            dataset_name: Name for the dataset
            ccd_ids: List of CCD three-letter codes
            download_if_missing: Auto-download CCD if not present
            
        Returns:
            List of successfully added entity names
            
        Example:
            >>> # Create dataset of common cofactors
            >>> processor.create_ccd_dataset(
            ...     "cofactors",
            ...     ["ATP", "NAD", "FAD", "COA", "HEM", "PLP"]
            ... )
        """
        successful = []
        
        for ccd_id in ccd_ids:
            ligand_data = self.get_ccd_ligand(ccd_id, download_if_missing)
            if ligand_data:
                try:
                    # Save as entity
                    self.save_entity(ligand_data['smiles'], ligand_data)
                    successful.append(ligand_data['smiles'])
                except Exception as e:
                    logger.warning(f"Failed to save CCD component {ccd_id}: {e}")
        
        if successful:
            # Create dataset
            self.create_dataset(dataset_name, successful)
            logger.info(f"Created dataset '{dataset_name}' with {len(successful)} CCD components")
        
        return successful
    
    def list_enamine_datasets(self) -> Dict[str, Dict]:
        """
        List all available Enamine datasets.
        
        Returns:
            Dictionary with dataset names and information
            
        Example:
            >>> datasets = processor.list_enamine_datasets()
            >>> for name, info in datasets.items():
            ...     print(f"{name}: {info['description']} ({info['size']} compounds)")
        """
        from protos.loaders import enamine_loader
        return enamine_loader.list_available_datasets()
    
    def get_enamine_dataset_info(self, dataset_name: str) -> Optional[Dict]:
        """
        Get detailed information about an Enamine dataset.
        
        Args:
            dataset_name: Name of the dataset
            
        Returns:
            Dataset information including download status
        """
        from protos.loaders import enamine_loader
        
        info = enamine_loader.get_dataset_info(dataset_name)
        if not info:
            return None
        
        # Add download status
        db_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'enamine'
        info['downloaded'] = enamine_loader.is_dataset_downloaded(db_path, dataset_name)
        
        # Add metadata if downloaded
        if info['downloaded']:
            metadata = enamine_loader.load_dataset_metadata(db_path, dataset_name)
            if metadata:
                info['download_date'] = metadata.get('downloaded')
        
        return info
    
    def get_database_statistics(self) -> Dict[str, Dict]:
        """
        Get statistics about available databases.
        
        Returns:
            Dictionary with database names and their statistics
        """
        from protos.loaders import qm9_loader, ccd_loader, enamine_loader
        
        stats = {}
        
        # Check CCD
        ccd_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'ccd'
        if ccd_loader.is_ccd_downloaded(ccd_path):
            stats['CCD'] = {
                'downloaded': True,
                'path': str(ccd_path),
                'component_count': len(list(ccd_path.glob('*.cif'))),
                'description': 'PDB Chemical Component Dictionary'
            }
        else:
            stats['CCD'] = {'downloaded': False, 'description': 'PDB Chemical Component Dictionary'}
        
        # Check QM9
        qm9_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'qm9'
        if qm9_loader.is_qm9_downloaded(qm9_path):
            stats['QM9'] = {
                'downloaded': True,
                'path': str(qm9_path),
                'molecule_count': 133885,
                'properties': ['gap', 'homo', 'lumo', 'dipole', 'alpha', 'zpve', 
                              'cv', 'u0', 'u298', 'h298', 'g298'],
                'description': 'Quantum chemistry dataset of small organic molecules'
            }
        else:
            stats['QM9'] = {'downloaded': False, 'description': 'Quantum chemistry dataset'}
        
        # Check Enamine
        enamine_path = Path(self.paths.get_processor_path(self.processor_type)) / self.databases_subdir / 'enamine'
        
        # Get info about all available datasets
        available_datasets = enamine_loader.list_available_datasets()
        downloaded_datasets = []
        
        if enamine_path.exists():
            for dataset_name in available_datasets:
                if enamine_loader.is_dataset_downloaded(enamine_path, dataset_name):
                    downloaded_datasets.append(dataset_name)
        
        stats['Enamine'] = {
            'downloaded': len(downloaded_datasets) > 0,
            'path': str(enamine_path) if enamine_path.exists() else None,
            'available_datasets': len(available_datasets),
            'downloaded_datasets': downloaded_datasets,
            'description': 'Enamine REAL commercial compound library',
            'credentials_set': bool(enamine_loader.get_enamine_credentials()[0])
        }
        
        return stats