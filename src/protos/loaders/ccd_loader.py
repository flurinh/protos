"""
PDB Chemical Component Dictionary (CCD) loader for Protos.

The CCD contains all chemical components (ligands, modified residues, etc.) found in 
PDB structures. This loader provides utilities for downloading and accessing CCD data.

References:
    https://www.wwpdb.org/data/ccd

This loader follows Protos principles:
- NO path management
- NO directory creation
- NO entity registry management  
- Pure utility functions only
"""

import os
import json
import gzip
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Try to import RDKit for molecular operations
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some CCD functionality will be limited.")
    HAS_RDKIT = False

# CCD URLs
CCD_BASE_URL = "https://files.wwpdb.org/pub/pdb/data/monomers"
CCD_COMPONENTS_URL = f"{CCD_BASE_URL}/components.cif.gz"
CCD_AA_VARIANTS_URL = f"{CCD_BASE_URL}/aa-variants-v1.cif.gz"

# Common CCD categories
CCD_CATEGORIES = [
    'peptide',
    'rna',
    'dna', 
    'saccharide',
    'non-polymer',
    'water'
]


def is_ccd_downloaded(cache_dir: Path) -> bool:
    """
    Check if CCD components are downloaded.
    
    Args:
        cache_dir: Directory where CCD data is stored
        
    Returns:
        True if downloaded, False otherwise
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    extracted_path = cache_dir / "components.cif"
    
    # Check if either compressed or extracted version exists
    return components_path.exists() or extracted_path.exists()


def download_ccd_components(cache_dir: Path, force_download: bool = False) -> bool:
    """
    Download CCD components file.
    
    Args:
        cache_dir: Directory to store the CCD data
        force_download: Force re-download even if exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if components_path.exists() and not force_download:
        logger.info(f"CCD components already exist at {components_path}")
        return True
    
    logger.info(f"Downloading CCD components to {components_path}")
    
    try:
        # Download with progress bar
        response = urllib.request.urlopen(CCD_COMPONENTS_URL)
        total_size = int(response.headers.get('Content-Length', 0))
        
        with open(components_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, 
                     desc="Downloading CCD components") as pbar:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))
        
        logger.info("CCD components downloaded successfully")
        return True
        
    except Exception as e:
        logger.error(f"Failed to download CCD components: {e}")
        if components_path.exists():
            components_path.unlink()
        return False


def download_ccd_aa_variants(cache_dir: Path, force_download: bool = False) -> bool:
    """
    Download CCD amino acid variants file.
    
    Args:
        cache_dir: Directory to store the CCD data
        force_download: Force re-download even if exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    variants_path = cache_dir / "aa-variants-v1.cif.gz"
    
    if variants_path.exists() and not force_download:
        logger.info(f"CCD AA variants already exist at {variants_path}")
        return True
    
    logger.info(f"Downloading CCD AA variants to {variants_path}")
    
    try:
        # Download with progress bar
        response = urllib.request.urlopen(CCD_AA_VARIANTS_URL)
        total_size = int(response.headers.get('Content-Length', 0))
        
        with open(variants_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, 
                     desc="Downloading AA variants") as pbar:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))
        
        logger.info("CCD AA variants downloaded successfully")
        return True
        
    except Exception as e:
        logger.error(f"Failed to download CCD AA variants: {e}")
        if variants_path.exists():
            variants_path.unlink()
        return False


def parse_ccd_cif_block(block_lines: List[str]) -> Optional[Dict]:
    """
    Parse a single CCD component block from CIF format.
    
    Args:
        block_lines: Lines containing one component block
        
    Returns:
        Dictionary with component data, or None if parsing fails
    """
    data = {}
    
    # Extract component ID from data_ line
    for line in block_lines:
        if line.startswith('data_'):
            data['id'] = line[5:].strip()
            break
    
    if 'id' not in data:
        return None
    
    # Parse key-value pairs
    i = 0
    while i < len(block_lines):
        line = block_lines[i].strip()
        
        if line.startswith('_chem_comp.'):
            # Split into field name and value
            parts = line.split(None, 1)
            if len(parts) < 2:
                i += 1
                continue
                
            field_name = parts[0]
            value_part = parts[1] if len(parts) > 1 else ''
            key = field_name.split('.')[1]
            
            # Handle quoted values
            if value_part.startswith('"') and value_part.endswith('"'):
                value = value_part[1:-1]
            elif value_part.startswith("'") and value_part.endswith("'"):
                value = value_part[1:-1]
            elif value_part.startswith('"') or value_part.startswith("'"):
                # Multi-line quoted value
                quote_char = value_part[0]
                value_lines = [value_part[1:]]
                i += 1
                while i < len(block_lines) and not block_lines[i].strip().endswith(quote_char):
                    value_lines.append(block_lines[i].rstrip())
                    i += 1
                if i < len(block_lines):
                    value_lines.append(block_lines[i].rstrip().rstrip(quote_char))
                value = '\n'.join(value_lines)
            else:
                # Simple unquoted value
                value = value_part
            
            data[key] = value
        
        # Parse descriptor section for SMILES
        elif data.get('id') and line.strip().startswith(data['id']) and 'SMILES' in line:
            # Handle format: ATP SMILES_CANONICAL "OpenEye OEToolkits" 1.5.0 "smiles..."
            parts = line.strip().split(None, 2)  # Split into 3 parts max
            if len(parts) >= 3:
                comp_id = parts[0]
                smiles_type = parts[1]  # SMILES or SMILES_CANONICAL
                rest = parts[2]
                
                # Extract SMILES value (last quoted string)
                # Find the last quoted string which contains the SMILES
                if '"' in rest:
                    # Find all quoted sections
                    quote_start = rest.rfind('"')
                    # Find the matching opening quote by going backwards
                    quote_count = 0
                    for i in range(quote_start - 1, -1, -1):
                        if rest[i] == '"':
                            quote_count += 1
                            if quote_count % 2 == 1:  # Found opening quote
                                smiles_value = rest[i+1:quote_start]
                                break
                    else:
                        # Fallback: take everything after last space
                        smiles_value = rest.split()[-1].strip('"')
                else:
                    # No quotes, take last token
                    smiles_value = rest.split()[-1]
                
                # Store different SMILES types
                if smiles_type == "SMILES":
                    data['smiles'] = smiles_value
                elif smiles_type == "SMILES_CANONICAL":
                    data['smiles_canonical'] = smiles_value
        
        # Parse InChI fields
        elif data.get('id') and line.strip().startswith(data['id']) and 'InChI' in line:
            parts = line.strip().split(None, 3)
            if len(parts) >= 4:
                comp_id = parts[0]
                inchi_type = parts[1]
                program = parts[2]
                value = parts[3].strip('"')
                
                if inchi_type == "InChI":
                    data['inchi'] = value
                elif inchi_type == "InChIKey":
                    data['inchi_key'] = value
        
        i += 1
    
    return data


def load_ccd_component(cache_dir: Path, comp_id: str) -> Optional[Dict]:
    """
    Load a single CCD component by ID.
    
    Args:
        cache_dir: Directory containing CCD data
        comp_id: Component ID (e.g., "ATP", "HEM")
        
    Returns:
        Dictionary with component data, or None if not found
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if not components_path.exists():
        logger.error(f"CCD components not found at {components_path}")
        return None
    
    comp_id = comp_id.upper()
    
    try:
        with gzip.open(components_path, 'rt') as f:
            in_block = False
            block_lines = []
            
            for line in f:
                if line.startswith(f'data_{comp_id}'):
                    in_block = True
                    block_lines = [line]
                elif in_block:
                    if line.startswith('data_') and line.strip() != f'data_{comp_id}':
                        # End of block
                        break
                    block_lines.append(line)
            
            if block_lines:
                return parse_ccd_cif_block(block_lines)
                
    except Exception as e:
        logger.error(f"Failed to load component {comp_id}: {e}")
    
    return None


def search_ccd_by_name(cache_dir: Path, name_pattern: str, limit: int = 100) -> List[Dict]:
    """
    Search CCD components by name pattern.
    
    Args:
        cache_dir: Directory containing CCD data
        name_pattern: Name pattern to search (case-insensitive)
        limit: Maximum number of results
        
    Returns:
        List of matching components
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if not components_path.exists():
        logger.error(f"CCD components not found at {components_path}")
        return []
    
    name_pattern = name_pattern.lower()
    results = []
    
    try:
        with gzip.open(components_path, 'rt') as f:
            current_block = []
            
            for line in f:
                if line.startswith('data_'):
                    # Process previous block if exists
                    if current_block:
                        comp_data = parse_ccd_cif_block(current_block)
                        if comp_data:
                            name = comp_data.get('name', '').lower()
                            if name_pattern in name:
                                results.append(comp_data)
                                if len(results) >= limit:
                                    break
                    
                    current_block = [line]
                else:
                    current_block.append(line)
            
            # Process last block
            if current_block and len(results) < limit:
                comp_data = parse_ccd_cif_block(current_block)
                if comp_data:
                    name = comp_data.get('name', '').lower()
                    if name_pattern in name:
                        results.append(comp_data)
                        
    except Exception as e:
        logger.error(f"Failed to search CCD components: {e}")
    
    logger.info(f"Found {len(results)} components matching '{name_pattern}'")
    return results


def search_ccd_by_formula(cache_dir: Path, formula: str, limit: int = 100) -> List[Dict]:
    """
    Search CCD components by chemical formula.
    
    Args:
        cache_dir: Directory containing CCD data
        formula: Chemical formula (e.g., "C10H14N2")
        limit: Maximum number of results
        
    Returns:
        List of matching components
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if not components_path.exists():
        logger.error(f"CCD components not found at {components_path}")
        return []
    
    results = []
    
    try:
        with gzip.open(components_path, 'rt') as f:
            current_block = []
            
            for line in f:
                if line.startswith('data_'):
                    # Process previous block if exists
                    if current_block:
                        comp_data = parse_ccd_cif_block(current_block)
                        if comp_data and comp_data.get('formula') == formula:
                            results.append(comp_data)
                            if len(results) >= limit:
                                break
                    
                    current_block = [line]
                else:
                    current_block.append(line)
            
            # Process last block
            if current_block and len(results) < limit:
                comp_data = parse_ccd_cif_block(current_block)
                if comp_data and comp_data.get('formula') == formula:
                    results.append(comp_data)
                    
    except Exception as e:
        logger.error(f"Failed to search CCD by formula: {e}")
    
    logger.info(f"Found {len(results)} components with formula '{formula}'")
    return results


def get_ccd_categories(cache_dir: Path) -> Dict[str, List[str]]:
    """
    Get CCD components organized by category.
    
    Args:
        cache_dir: Directory containing CCD data
        
    Returns:
        Dictionary mapping categories to component IDs
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if not components_path.exists():
        logger.error(f"CCD components not found at {components_path}")
        return {}
    
    categories = {cat: [] for cat in CCD_CATEGORIES}
    categories['other'] = []
    
    try:
        with gzip.open(components_path, 'rt') as f:
            current_block = []
            
            for line in f:
                if line.startswith('data_'):
                    # Process previous block if exists
                    if current_block:
                        comp_data = parse_ccd_cif_block(current_block)
                        if comp_data:
                            comp_type = comp_data.get('type', '').lower()
                            comp_id = comp_data.get('id', '')
                            
                            # Categorize
                            categorized = False
                            for cat in CCD_CATEGORIES:
                                if cat in comp_type:
                                    categories[cat].append(comp_id)
                                    categorized = True
                                    break
                            
                            if not categorized:
                                categories['other'].append(comp_id)
                    
                    current_block = [line]
                else:
                    current_block.append(line)
            
            # Process last block
            if current_block:
                comp_data = parse_ccd_cif_block(current_block)
                if comp_data:
                    comp_type = comp_data.get('type', '').lower()
                    comp_id = comp_data.get('id', '')
                    
                    # Categorize
                    categorized = False
                    for cat in CCD_CATEGORIES:
                        if cat in comp_type:
                            categories[cat].append(comp_id)
                            categorized = True
                            break
                    
                    if not categorized:
                        categories['other'].append(comp_id)
                        
    except Exception as e:
        logger.error(f"Failed to categorize CCD components: {e}")
    
    # Log summary
    for cat, ids in categories.items():
        logger.info(f"Category '{cat}': {len(ids)} components")
    
    return categories


def get_ccd_smiles(cache_dir: Path, comp_id: str) -> Optional[str]:
    """
    Get SMILES string for a CCD component.
    
    Args:
        cache_dir: Directory containing CCD data
        comp_id: Component ID
        
    Returns:
        SMILES string, or None if not found
    """
    comp_data = load_ccd_component(cache_dir, comp_id)
    if not comp_data:
        return None
    
    # Try different SMILES fields
    smiles_fields = [
        'pdbx_smiles_canonical',
        'pdbx_smiles',
        'smiles_canonical',
        'smiles'
    ]
    
    for field in smiles_fields:
        smiles = comp_data.get(field)
        if smiles and smiles != '?':
            return smiles
    
    return None


def create_ccd_index(cache_dir: Path, force_rebuild: bool = False) -> bool:
    """
    Create an index file for fast CCD lookups.
    
    Args:
        cache_dir: Directory containing CCD data
        force_rebuild: Force rebuild even if index exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    index_path = cache_dir / "ccd_index.json"
    
    if index_path.exists() and not force_rebuild:
        logger.info(f"CCD index already exists at {index_path}")
        return True
    
    if not components_path.exists():
        logger.error(f"CCD components not found at {components_path}")
        return False
    
    logger.info("Building CCD index...")
    
    index = {
        'components': {},
        'by_type': {},
        'by_formula': {},
        'has_smiles': []
    }
    
    try:
        with gzip.open(components_path, 'rt') as f:
            current_block = []
            
            for line in f:
                if line.startswith('data_'):
                    # Process previous block if exists
                    if current_block:
                        comp_data = parse_ccd_cif_block(current_block)
                        if comp_data:
                            comp_id = comp_data.get('id', '')
                            comp_type = comp_data.get('type', '')
                            formula = comp_data.get('formula', '')
                            
                            # Add to main index
                            index['components'][comp_id] = {
                                'name': comp_data.get('name', ''),
                                'type': comp_type,
                                'formula': formula
                            }
                            
                            # Index by type
                            if comp_type:
                                if comp_type not in index['by_type']:
                                    index['by_type'][comp_type] = []
                                index['by_type'][comp_type].append(comp_id)
                            
                            # Index by formula
                            if formula and formula != '?':
                                if formula not in index['by_formula']:
                                    index['by_formula'][formula] = []
                                index['by_formula'][formula].append(comp_id)
                            
                            # Track components with SMILES
                            smiles = get_ccd_smiles(cache_dir, comp_id)
                            if smiles:
                                index['has_smiles'].append(comp_id)
                    
                    current_block = [line]
                else:
                    current_block.append(line)
        
        # Save index
        with open(index_path, 'w') as f:
            json.dump(index, f, indent=2)
        
        logger.info(f"CCD index created with {len(index['components'])} components")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create CCD index: {e}")
        return False


def load_ccd_index(cache_dir: Path) -> Optional[Dict]:
    """
    Load CCD index for fast lookups.
    
    Args:
        cache_dir: Directory containing CCD index
        
    Returns:
        Index dictionary, or None if not found
    """
    cache_dir = Path(cache_dir)
    index_path = cache_dir / "ccd_index.json"
    
    if not index_path.exists():
        logger.warning(f"CCD index not found at {index_path}")
        return None
    
    try:
        with open(index_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load CCD index: {e}")
        return None


def get_common_ligands() -> List[str]:
    """
    Get list of common ligand IDs found in PDB structures.
    
    Returns:
        List of common ligand component IDs
    """
    return [
        'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP',  # Nucleotides
        'NAD', 'NAP', 'FAD', 'FMN',  # Cofactors
        'HEM', 'HEC', 'HEA', 'HEB',  # Heme variants
        'SO4', 'PO4', 'GOL', 'EDO',  # Common ions/solvents
        'SAM', 'SAH',  # SAM/SAH
        'COA', 'ACO',  # Coenzyme A
        'PLP', 'PMP',  # Vitamin B6
        'TPP', 'TDP',  # Thiamine
        'B12', 'CNC',  # Vitamin B12
        'F3S', 'F4S', 'SF4',  # Iron-sulfur clusters
        'ZN', 'MG', 'CA', 'MN', 'FE', 'CU', 'NI',  # Metal ions
        'CL', 'NA', 'K',  # Common ions
    ]


def ensure_ccd_ready(cache_dir: Path) -> bool:
    """
    Ensure CCD is downloaded and ready to use.
    
    CCD is read directly from gzip, so no extraction needed.
    
    Args:
        cache_dir: Directory for CCD data
        
    Returns:
        True if CCD is ready to use, False otherwise
    """
    cache_dir = Path(cache_dir)
    components_path = cache_dir / "components.cif.gz"
    
    if components_path.exists():
        return True
    
    logger.warning("CCD not found. Call download_ccd_components() first.")
    return False


def get_ccd_ligand_safe(cache_dir: Path, comp_id: str) -> Optional[Dict]:
    """
    Safely load a CCD component, ensuring CCD is ready.
    Uses fast indexed access if available.
    
    Args:
        cache_dir: Directory containing CCD data
        comp_id: Component ID (e.g., "ATP", "HEM")
        
    Returns:
        Dictionary with component data including SMILES, or None
    """
    if not ensure_ccd_ready(cache_dir):
        return None
    
    # Try fast indexed access first
    try:
        from .ccd_index_builder import load_ccd_ligand_fast
        ligand = load_ccd_ligand_fast(cache_dir, comp_id)
        if ligand and (ligand.get('smiles') or ligand.get('smiles_canonical')):
            # Return in expected format
            return {
                'id': comp_id.upper(),
                'name': ligand.get('name', comp_id),
                'formula': ligand.get('formula', ''),
                'type': ligand.get('type', 'non-polymer'),
                'smiles': ligand.get('smiles') or ligand.get('smiles_canonical'),
                'inchi': ligand.get('inchi'),
                'inchi_key': ligand.get('inchi_key')
            }
    except:
        logger.debug("Fast CCD access not available, falling back to slow method")
    
    # Fall back to slow method
    comp_data = load_ccd_component(cache_dir, comp_id)
    if not comp_data:
        return None
    
    # Get SMILES
    smiles = get_ccd_smiles(cache_dir, comp_id)
    if not smiles:
        logger.warning(f"No SMILES found for {comp_id}")
        return None
    
    # Build complete ligand data
    result = {
        'id': comp_id,
        'name': comp_data.get('name', comp_id),
        'formula': comp_data.get('formula', ''),
        'type': comp_data.get('type', 'non-polymer'),
        'smiles': smiles
    }
    
    # Add InChI if available
    if 'inchi' in comp_data:
        result['inchi'] = comp_data['inchi']
    if 'inchi_key' in comp_data:
        result['inchi_key'] = comp_data['inchi_key']
    
    return result