#!/usr/bin/env python3
"""
Build a fast CCD index using gemmi for high-throughput access.
"""

import json
import gzip
import logging
from pathlib import Path
from typing import Dict, Optional
import pickle

try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    HAS_GEMMI = False
    logging.warning("gemmi not available. CCD indexing will be slower.")

logger = logging.getLogger(__name__)


def build_ccd_index_with_gemmi(ccd_gz_path: Path, output_dir: Path) -> bool:
    """
    Build CCD index using gemmi for fast parsing.
    Creates both JSON index and pickled data for ultra-fast access.
    
    Args:
        ccd_gz_path: Path to components.cif.gz
        output_dir: Directory to save index files
        
    Returns:
        True if successful
    """
    if not HAS_GEMMI:
        logger.error("gemmi is required for fast CCD indexing")
        return False
        
    if not ccd_gz_path.exists():
        logger.error(f"CCD file not found: {ccd_gz_path}")
        return False
    
    logger.info(f"Building CCD index from {ccd_gz_path}")
    
    try:
        # Parse the entire CCD file with gemmi
        # First extract to a temporary file since gemmi needs a filename
        import tempfile
        import shutil
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
            with gzip.open(ccd_gz_path, 'rt') as f:
                shutil.copyfileobj(f, tmp)
            tmp_path = tmp.name
        
        # Now read with gemmi
        doc = gemmi.cif.read(tmp_path)
        
        # Clean up temp file
        Path(tmp_path).unlink()
        
        # Build comprehensive index
        index = {
            'components': {},
            'by_type': {},
            'by_formula': {},
            'with_smiles': [],
            'with_inchi': []
        }
        
        # Also build a ligand data cache
        ligand_data = {}
        
        logger.info(f"Processing {len(doc)} components...")
        
        for block in doc:
            code = block.name
            
            # Extract key information
            comp_data = {
                'id': code,
                'name': None,
                'formula': None,
                'type': None,
                'smiles': None,
                'smiles_canonical': None,
                'inchi': None,
                'inchi_key': None
            }
            
            # Get basic properties
            try:
                comp_group = block.find_loop('_chem_comp.id')
                if comp_group:
                    row = comp_group.find_row(code)
                    if row:
                        comp_data['name'] = comp_group.find_value('_chem_comp.name', row)
                        comp_data['formula'] = comp_group.find_value('_chem_comp.formula', row)
                        comp_data['type'] = comp_group.find_value('_chem_comp.type', row)
            except:
                # Try direct access
                comp_data['name'] = block.find_value('_chem_comp.name')
                comp_data['formula'] = block.find_value('_chem_comp.formula')
                comp_data['type'] = block.find_value('_chem_comp.type')
            
            # Get SMILES from descriptor table
            try:
                desc_loop = block.find_loop('_pdbx_chem_comp_descriptor.comp_id')
                if desc_loop:
                    for row in desc_loop:
                        desc_type = desc_loop.find_value('_pdbx_chem_comp_descriptor.type', row)
                        if desc_type == 'SMILES':
                            comp_data['smiles'] = desc_loop.find_value('_pdbx_chem_comp_descriptor.descriptor', row)
                        elif desc_type == 'SMILES_CANONICAL':
                            comp_data['smiles_canonical'] = desc_loop.find_value('_pdbx_chem_comp_descriptor.descriptor', row)
                        elif desc_type == 'InChI':
                            comp_data['inchi'] = desc_loop.find_value('_pdbx_chem_comp_descriptor.descriptor', row)
                        elif desc_type == 'InChIKey':
                            comp_data['inchi_key'] = desc_loop.find_value('_pdbx_chem_comp_descriptor.descriptor', row)
            except:
                pass
            
            # Store in index
            index['components'][code] = {
                'name': comp_data['name'],
                'formula': comp_data['formula'],
                'type': comp_data['type'],
                'has_smiles': bool(comp_data['smiles'] or comp_data['smiles_canonical']),
                'has_inchi': bool(comp_data['inchi'])
            }
            
            # Index by type
            if comp_data['type']:
                comp_type = comp_data['type']
                if comp_type not in index['by_type']:
                    index['by_type'][comp_type] = []
                index['by_type'][comp_type].append(code)
            
            # Index by formula
            if comp_data['formula'] and comp_data['formula'] != '?':
                formula = comp_data['formula']
                if formula not in index['by_formula']:
                    index['by_formula'][formula] = []
                index['by_formula'][formula].append(code)
            
            # Track components with chemical data
            if comp_data['smiles'] or comp_data['smiles_canonical']:
                index['with_smiles'].append(code)
            if comp_data['inchi']:
                index['with_inchi'].append(code)
            
            # Store full data for fast access
            ligand_data[code] = comp_data
        
        # Save index as JSON
        index_path = output_dir / "ccd_index.json"
        with open(index_path, 'w') as f:
            json.dump(index, f, indent=2)
        logger.info(f"Saved index to {index_path}")
        
        # Save full data as pickle for ultra-fast loading
        data_path = output_dir / "ccd_ligands.pkl"
        with open(data_path, 'wb') as f:
            pickle.dump(ligand_data, f)
        logger.info(f"Saved ligand data to {data_path}")
        
        # Print statistics
        logger.info(f"\nCCD Index Statistics:")
        logger.info(f"  Total components: {len(index['components'])}")
        logger.info(f"  With SMILES: {len(index['with_smiles'])}")
        logger.info(f"  With InChI: {len(index['with_inchi'])}")
        logger.info(f"  Unique types: {len(index['by_type'])}")
        logger.info(f"  Unique formulas: {len(index['by_formula'])}")
        
        return True
        
    except Exception as e:
        logger.error(f"Failed to build CCD index: {e}")
        import traceback
        traceback.print_exc()
        return False


def load_ccd_ligand_fast(cache_dir: Path, comp_id: str) -> Optional[Dict]:
    """
    Load CCD component from pre-built index (ultra-fast).
    
    Args:
        cache_dir: Directory containing index files
        comp_id: Component ID (e.g., "ATP")
        
    Returns:
        Component data dictionary
    """
    cache_dir = Path(cache_dir)
    data_path = cache_dir / "ccd_ligands.pkl"
    
    if not data_path.exists():
        logger.error(f"CCD ligand data not found. Run build_ccd_index_with_gemmi first.")
        return None
    
    try:
        # Load pickled data (very fast)
        with open(data_path, 'rb') as f:
            ligand_data = pickle.load(f)
        
        return ligand_data.get(comp_id.upper())
        
    except Exception as e:
        logger.error(f"Failed to load CCD data: {e}")
        return None


def search_ccd_by_property_fast(cache_dir: Path, property_name: str, value: str) -> list:
    """
    Fast search using pre-built index.
    
    Args:
        cache_dir: Directory containing index files
        property_name: Property to search (type, formula, etc.)
        value: Value to search for
        
    Returns:
        List of matching component IDs
    """
    cache_dir = Path(cache_dir)
    index_path = cache_dir / "ccd_index.json"
    
    if not index_path.exists():
        logger.error("CCD index not found. Run build_ccd_index_with_gemmi first.")
        return []
    
    try:
        with open(index_path, 'r') as f:
            index = json.load(f)
        
        if property_name == 'type':
            return index['by_type'].get(value, [])
        elif property_name == 'formula':
            return index['by_formula'].get(value, [])
        elif property_name == 'has_smiles':
            return index['with_smiles'] if value else []
        elif property_name == 'has_inchi':
            return index['with_inchi'] if value else []
        else:
            # Search in components
            matches = []
            for code, data in index['components'].items():
                if data.get(property_name) == value:
                    matches.append(code)
            return matches
            
    except Exception as e:
        logger.error(f"Failed to search index: {e}")
        return []


if __name__ == "__main__":
    # Test the builder
    import sys
    if len(sys.argv) > 1:
        ccd_file = Path(sys.argv[1])
        output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else ccd_file.parent
        
        success = build_ccd_index_with_gemmi(ccd_file, output_dir)
        
        if success:
            # Test fast loading
            print("\nTesting fast loader...")
            for code in ['ATP', 'NAD', 'HEM']:
                ligand = load_ccd_ligand_fast(output_dir, code)
                if ligand:
                    print(f"{code}: {ligand.get('name')} - {ligand.get('formula')}")
                    if ligand.get('smiles'):
                        print(f"  SMILES: {ligand['smiles'][:50]}...")
        else:
            sys.exit(1)
    else:
        print("Usage: python ccd_index_builder.py <components.cif.gz> [output_dir]")