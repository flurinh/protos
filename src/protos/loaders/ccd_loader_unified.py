"""
Unified CCD loader for Protos.

The CCD contains all chemical components (ligands, modified residues, etc.) found in 
PDB structures. This loader provides download, indexing, and fast access functionality.

This loader follows Protos principles:
- NO path management (paths come from outside)
- NO directory creation
- Pure utility functions for data operations
"""

import os
import json
import gzip
import logging
import pickle
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Set
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Try to import RDKit and gemmi for enhanced functionality
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some CCD functionality will be limited.")
    HAS_RDKIT = False

try:
    import gemmi
    HAS_GEMMI = True
except ImportError:
    logger.warning("gemmi not available. CCD indexing will use fallback parser.")
    HAS_GEMMI = False

# CCD URLs
CCD_BASE_URL = "https://files.wwpdb.org/pub/pdb/data/monomers"
CCD_COMPONENTS_URL = f"{CCD_BASE_URL}/components.cif.gz"


def download_ccd(cache_dir: Path, force: bool = False) -> bool:
    """
    Download CCD components file.
    
    Args:
        cache_dir: Directory to store the CCD data
        force: Force re-download even if exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if components_path.exists() and not force:
        logger.info(f"CCD already exists at {components_path}")
        return True
    
    logger.info(f"Downloading CCD to {components_path}")
    
    try:
        response = urllib.request.urlopen(CCD_COMPONENTS_URL)
        total_size = int(response.headers.get('Content-Length', 0))
        
        # Download with progress bar
        with open(components_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, 
                     desc="Downloading CCD") as pbar:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))
        
        logger.info("CCD downloaded successfully")
        return True
        
    except Exception as e:
        logger.error(f"Failed to download CCD: {e}")
        if components_path.exists():
            components_path.unlink()
        return False


def build_ccd_index(cache_dir: Path, force: bool = False) -> bool:
    """
    Build index for fast CCD access.
    
    Creates:
    - ccd_index.json: Component metadata and lookup indices
    - ccd_ligands.pkl: Full component data for fast loading
    
    Args:
        cache_dir: Directory containing CCD data
        force: Force rebuild even if index exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    index_path = cache_dir / "ccd_index.json"
    data_path = cache_dir / "ccd_ligands.pkl"
    
    # Check if already indexed
    if not force and index_path.exists() and data_path.exists():
        logger.info("CCD index already exists")
        return True
    
    if not components_path.exists():
        logger.error(f"CCD not found at {components_path}. Download first.")
        return False
    
    logger.info("Building CCD index...")
    
    # Try gemmi first for fast parsing
    if HAS_GEMMI:
        try:
            return _build_index_with_gemmi(components_path, index_path, data_path)
        except Exception as e:
            logger.warning(f"Gemmi indexing failed: {e}. Falling back to manual parser.")
    
    # Fallback to manual parsing
    return _build_index_manual(components_path, index_path, data_path)


def _build_index_with_gemmi(ccd_path: Path, index_path: Path, data_path: Path) -> bool:
    """Build index using gemmi (fast)."""
    try:
        # Gemmi can read .gz files directly for smaller files
        # For large files, we need to extract first
        import tempfile
        
        logger.info("Extracting CCD for gemmi parsing...")
        with tempfile.NamedTemporaryFile(suffix='.cif', delete=False) as tmp:
            with gzip.open(ccd_path, 'rb') as f_in:
                # Extract in chunks with progress
                total_size = 0
                with tqdm(desc="Extracting CCD", unit='B', unit_scale=True) as pbar:
                    while True:
                        chunk = f_in.read(1024 * 1024)  # 1MB chunks
                        if not chunk:
                            break
                        tmp.write(chunk)
                        total_size += len(chunk)
                        pbar.update(len(chunk))
            tmp_path = Path(tmp.name)
        
        logger.info(f"Parsing {total_size / 1024 / 1024:.1f} MB CCD file with gemmi...")
        doc = gemmi.cif.read(str(tmp_path))
        
        # Build index
        index = {
            'components': {},
            'by_type': {},
            'by_formula': {},
            'with_smiles': [],
            'with_inchi': []
        }
        ligand_data = {}
        
        logger.info(f"Processing {len(doc)} components...")
        for i, block in enumerate(tqdm(doc, desc="Indexing components")):
            code = block.name
            comp_data = _extract_component_data_gemmi(block)
            
            if comp_data:
                # Add to index
                index['components'][code] = {
                    'name': comp_data['name'],
                    'formula': comp_data['formula'],
                    'type': comp_data['type'],
                    'has_smiles': bool(comp_data.get('smiles') or comp_data.get('smiles_canonical')),
                    'has_inchi': bool(comp_data.get('inchi'))
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
                if comp_data.get('smiles') or comp_data.get('smiles_canonical'):
                    index['with_smiles'].append(code)
                if comp_data.get('inchi'):
                    index['with_inchi'].append(code)
                
                # Store full data
                ligand_data[code] = comp_data
        
        # Clean up temp file
        tmp_path.unlink()
        
        # Save index and data
        with open(index_path, 'w') as f:
            json.dump(index, f, indent=2)
        
        with open(data_path, 'wb') as f:
            pickle.dump(ligand_data, f)
        
        logger.info(f"Indexed {len(ligand_data)} components")
        logger.info(f"  With SMILES: {len(index['with_smiles'])}")
        logger.info(f"  With InChI: {len(index['with_inchi'])}")
        
        return True
        
    except Exception as e:
        logger.error(f"Gemmi indexing failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def _extract_component_data_gemmi(block) -> Dict:
    """Extract component data from gemmi block."""
    comp_data = {
        'id': block.name,
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
        comp_data['name'] = block.find_value('_chem_comp.name')
        comp_data['formula'] = block.find_value('_chem_comp.formula')
        comp_data['type'] = block.find_value('_chem_comp.type')
    except:
        pass
    
    # Get SMILES from descriptor table
    try:
        desc_loop = block.find_loop('_pdbx_chem_comp_descriptor.comp_id')
        if desc_loop:
            for row in desc_loop:
                desc_type = desc_loop.find_value('_pdbx_chem_comp_descriptor.type', row)
                descriptor = desc_loop.find_value('_pdbx_chem_comp_descriptor.descriptor', row)
                
                if desc_type == 'SMILES':
                    comp_data['smiles'] = descriptor
                elif desc_type == 'SMILES_CANONICAL':
                    comp_data['smiles_canonical'] = descriptor
                elif desc_type == 'InChI':
                    comp_data['inchi'] = descriptor
                elif desc_type == 'InChIKey':
                    comp_data['inchi_key'] = descriptor
    except:
        pass
    
    return comp_data


def _build_index_manual(ccd_path: Path, index_path: Path, data_path: Path) -> bool:
    """Build index using manual parsing (slower but more robust)."""
    try:
        index = {
            'components': {},
            'by_type': {},
            'by_formula': {},
            'with_smiles': [],
            'with_inchi': []
        }
        ligand_data = {}
        
        with gzip.open(ccd_path, 'rt') as f:
            current_block = []
            current_id = None
            
            logger.info("Parsing CCD file (this may take a few minutes)...")
            for line in tqdm(f, desc="Reading CCD"):
                if line.startswith('data_'):
                    # Process previous block if exists
                    if current_block and current_id:
                        comp_data = _parse_component_block(current_block)
                        if comp_data:
                            # Add to index
                            index['components'][current_id] = {
                                'name': comp_data.get('name'),
                                'formula': comp_data.get('formula'),
                                'type': comp_data.get('type'),
                                'has_smiles': bool(comp_data.get('smiles') or comp_data.get('smiles_canonical')),
                                'has_inchi': bool(comp_data.get('inchi'))
                            }
                            
                            # Index by type
                            if comp_data.get('type'):
                                comp_type = comp_data['type']
                                if comp_type not in index['by_type']:
                                    index['by_type'][comp_type] = []
                                index['by_type'][comp_type].append(current_id)
                            
                            # Store full data
                            ligand_data[current_id] = comp_data
                    
                    # Start new block
                    current_id = line[5:].strip()
                    current_block = [line]
                else:
                    current_block.append(line)
            
            # Process last block
            if current_block and current_id:
                comp_data = _parse_component_block(current_block)
                if comp_data:
                    ligand_data[current_id] = comp_data
        
        # Build remaining indices
        logger.info("Building indices...")
        for code, comp_data in ligand_data.items():
            # Index by formula
            if comp_data.get('formula') and comp_data['formula'] != '?':
                formula = comp_data['formula']
                if formula not in index['by_formula']:
                    index['by_formula'][formula] = []
                index['by_formula'][formula].append(code)
            
            # Track components with chemical data
            if comp_data.get('smiles') or comp_data.get('smiles_canonical'):
                index['with_smiles'].append(code)
            if comp_data.get('inchi'):
                index['with_inchi'].append(code)
        
        # Save index and data
        with open(index_path, 'w') as f:
            json.dump(index, f, indent=2)
        
        with open(data_path, 'wb') as f:
            pickle.dump(ligand_data, f)
        
        logger.info(f"Indexed {len(ligand_data)} components")
        return True
        
    except Exception as e:
        logger.error(f"Manual indexing failed: {e}")
        return False


def _parse_component_block(lines: List[str]) -> Dict:
    """Parse a single component block (simplified)."""
    data = {'id': None}
    
    for line in lines:
        line = line.strip()
        
        # Extract component ID
        if line.startswith('data_'):
            data['id'] = line[5:]
        
        # Extract basic fields
        elif line.startswith('_chem_comp.'):
            parts = line.split(None, 1)
            if len(parts) == 2:
                field = parts[0].split('.')[1]
                value = parts[1].strip(' "\'\n')
                data[field] = value
        
        # Extract SMILES
        elif data.get('id') and 'SMILES' in line and line.startswith(data['id']):
            if 'SMILES_CANONICAL' in line:
                # Extract canonical SMILES
                parts = line.split('"')
                if len(parts) >= 2:
                    data['smiles_canonical'] = parts[-2]
            elif 'SMILES' in line:
                # Extract regular SMILES
                parts = line.split('"')
                if len(parts) >= 2:
                    data['smiles'] = parts[-2]
    
    return data


def is_ccd_ready(cache_dir: Path) -> bool:
    """
    Check if CCD is downloaded and indexed.
    
    Args:
        cache_dir: Directory containing CCD data
        
    Returns:
        True if ready to use, False otherwise
    """
    cache_dir = Path(cache_dir)
    return (cache_dir / "ccd_ligands.pkl").exists()


def ensure_ccd_ready(cache_dir: Path) -> bool:
    """
    Ensure CCD is downloaded and indexed.
    
    Args:
        cache_dir: Directory for CCD data
        
    Returns:
        True if ready, False otherwise
    """
    cache_dir = Path(cache_dir)
    
    # Check if already ready
    if is_ccd_ready(cache_dir):
        return True
    
    # Download if needed
    if not (cache_dir / "components.cif.gz").exists():
        logger.info("CCD not found, downloading...")
        if not download_ccd(cache_dir):
            return False
    
    # Build index
    logger.info("Building CCD index...")
    return build_ccd_index(cache_dir)


def get_ccd_component(cache_dir: Path, comp_id: str) -> Optional[Dict]:
    """
    Get CCD component data (fast indexed access).
    
    Args:
        cache_dir: Directory containing CCD data
        comp_id: Component ID (e.g., "ATP", "HEM")
        
    Returns:
        Component data dictionary, or None if not found
    """
    cache_dir = Path(cache_dir)
    data_path = cache_dir / "ccd_ligands.pkl"
    
    if not data_path.exists():
        logger.error("CCD not indexed. Run ensure_ccd_ready() first.")
        return None
    
    try:
        with open(data_path, 'rb') as f:
            ligand_data = pickle.load(f)
        return ligand_data.get(comp_id.upper())
    except Exception as e:
        logger.error(f"Failed to load CCD data: {e}")
        return None


def search_ccd(cache_dir: Path, query_type: str, query_value: str) -> List[str]:
    """
    Search CCD by various criteria.
    
    Args:
        cache_dir: Directory containing CCD data
        query_type: Type of search ('type', 'formula', 'has_smiles')
        query_value: Value to search for
        
    Returns:
        List of matching component IDs
    """
    cache_dir = Path(cache_dir)
    index_path = cache_dir / "ccd_index.json"
    
    if not index_path.exists():
        logger.error("CCD index not found. Run ensure_ccd_ready() first.")
        return []
    
    try:
        with open(index_path, 'r') as f:
            index = json.load(f)
        
        if query_type == 'type':
            return index['by_type'].get(query_value, [])
        elif query_type == 'formula':
            return index['by_formula'].get(query_value, [])
        elif query_type == 'has_smiles' and query_value:
            return index['with_smiles']
        elif query_type == 'has_inchi' and query_value:
            return index['with_inchi']
        else:
            return []
            
    except Exception as e:
        logger.error(f"Failed to search CCD: {e}")
        return []


def get_ccd_statistics(cache_dir: Path) -> Dict:
    """
    Get statistics about the CCD database.
    
    Args:
        cache_dir: Directory containing CCD data
        
    Returns:
        Dictionary with statistics
    """
    cache_dir = Path(cache_dir)
    index_path = cache_dir / "ccd_index.json"
    
    stats = {
        'total_components': 0,
        'with_smiles': 0,
        'with_inchi': 0,
        'types': {},
        'indexed': index_path.exists()
    }
    
    if not index_path.exists():
        return stats
    
    try:
        with open(index_path, 'r') as f:
            index = json.load(f)
        
        stats['total_components'] = len(index['components'])
        stats['with_smiles'] = len(index['with_smiles'])
        stats['with_inchi'] = len(index['with_inchi'])
        
        # Count by type
        for comp_type, codes in index['by_type'].items():
            stats['types'][comp_type] = len(codes)
            
    except Exception as e:
        logger.error(f"Failed to load statistics: {e}")
    
    return stats


# Common ligand IDs for quick reference
COMMON_LIGANDS = [
    'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',  # Nucleotides
    'NAD', 'NAP', 'FAD', 'FMN',  # Cofactors
    'HEM', 'HEC', 'HEA', 'HEB',  # Heme variants
    'SAM', 'SAH',  # SAM/SAH
    'COA', 'ACO',  # Coenzyme A
    'PLP', 'PMP',  # Vitamin B6
]