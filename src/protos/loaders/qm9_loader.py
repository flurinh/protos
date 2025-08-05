"""
QM9 dataset loader for Protos.

The QM9 dataset contains ~134k molecules with 19 properties calculated at DFT level.
This loader provides utility functions for downloading and accessing QM9 data locally.

References:
    Ramakrishnan et al. "Quantum chemistry structures and properties of 134 kilo 
    molecules." Scientific Data 1 (2014): 140022.

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
import tarfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import urllib.request
import urllib.error
from tqdm import tqdm

logger = logging.getLogger(__name__)

# QM9 dataset URL
QM9_URL = "https://figshare.com/ndownloader/files/3195389"
QM9_FILENAME = "dsgdb9nsd.xyz.tar.bz2"

# Property names in QM9
QM9_PROPERTIES = [
    'rcA', 'rcB', 'rcC',  # Rotational constants
    'mu',  # Dipole moment
    'alpha',  # Isotropic polarizability
    'homo',  # HOMO energy
    'lumo',  # LUMO energy
    'gap',  # HOMO-LUMO gap
    'r2',  # Electronic spatial extent
    'zpve',  # Zero point vibrational energy
    'U0',  # Internal energy at 0K
    'U',  # Internal energy at 298.15K
    'H',  # Enthalpy at 298.15K
    'G',  # Free energy at 298.15K
    'Cv',  # Heat capacity at 298.15K
    'omega1',  # Highest vibrational frequency
]


def is_qm9_downloaded(cache_dir: Path) -> bool:
    """
    Check if QM9 dataset is downloaded.
    
    Args:
        cache_dir: Directory where QM9 data is stored
        
    Returns:
        True if downloaded, False otherwise
    """
    cache_dir = Path(cache_dir)
    dataset_path = cache_dir / QM9_FILENAME
    extracted_dir = cache_dir / "dsgdb9nsd"
    
    # Check if either archive or extracted directory exists
    return dataset_path.exists() or extracted_dir.exists()


def download_qm9_dataset(cache_dir: Path, force_download: bool = False) -> bool:
    """
    Download QM9 dataset to cache directory.
    
    Args:
        cache_dir: Directory to store the dataset
        force_download: Force re-download even if exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    dataset_path = cache_dir / QM9_FILENAME
    
    if dataset_path.exists() and not force_download:
        logger.info(f"QM9 dataset already exists at {dataset_path}")
        return True
    
    logger.info(f"Downloading QM9 dataset to {dataset_path}")
    
    try:
        # Download with progress bar
        response = urllib.request.urlopen(QM9_URL)
        total_size = int(response.headers.get('Content-Length', 0))
        
        with open(dataset_path, 'wb') as f:
            with tqdm(total=total_size, unit='B', unit_scale=True, 
                     desc="Downloading QM9 dataset") as pbar:
                while True:
                    chunk = response.read(8192)
                    if not chunk:
                        break
                    f.write(chunk)
                    pbar.update(len(chunk))
        
        logger.info("QM9 dataset downloaded successfully")
        return True
        
    except Exception as e:
        logger.error(f"Failed to download QM9 dataset: {e}")
        if dataset_path.exists():
            dataset_path.unlink()
        return False


def extract_qm9_dataset(cache_dir: Path) -> bool:
    """
    Extract QM9 dataset from tar.bz2 archive.
    
    Args:
        cache_dir: Directory containing the dataset
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    dataset_path = cache_dir / QM9_FILENAME
    extract_dir = cache_dir / "qm9_molecules"
    
    if not dataset_path.exists():
        logger.error(f"QM9 dataset not found at {dataset_path}")
        return False
    
    if extract_dir.exists() and any(extract_dir.iterdir()):
        logger.info(f"QM9 dataset already extracted at {extract_dir}")
        return True
    
    logger.info(f"Extracting QM9 dataset to {extract_dir}")
    
    try:
        extract_dir.mkdir(exist_ok=True)
        
        with tarfile.open(dataset_path, "r:bz2") as tar:
            members = tar.getmembers()
            with tqdm(total=len(members), desc="Extracting QM9 dataset") as pbar:
                for member in members:
                    tar.extract(member, extract_dir)
                    pbar.update(1)
        
        logger.info("QM9 dataset extracted successfully")
        return True
        
    except Exception as e:
        logger.error(f"Failed to extract QM9 dataset: {e}")
        return False


def parse_qm9_xyz(xyz_path: Path) -> Optional[Dict]:
    """
    Parse a single QM9 XYZ file.
    
    QM9 XYZ format:
    Line 1: Number of atoms
    Line 2: Properties (separated by spaces)
    Line 3+: Element x y z mulliken_charge
    Last lines: Vibrational frequencies
    Last line: SMILES InChI
    
    Args:
        xyz_path: Path to XYZ file
        
    Returns:
        Dictionary with molecule data, or None if parsing fails
    """
    try:
        with open(xyz_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 3:
            return None
        
        # Parse number of atoms
        n_atoms = int(lines[0].strip())
        
        # Parse properties
        props = lines[1].strip().split()
        if len(props) < len(QM9_PROPERTIES):
            logger.warning(f"Incomplete properties in {xyz_path}")
            return None
        
        properties = {name: float(value) for name, value in zip(QM9_PROPERTIES, props)}
        
        # Parse atoms
        atoms = []
        for i in range(2, 2 + n_atoms):
            parts = lines[i].strip().split()
            if len(parts) >= 5:
                atom = {
                    'element': parts[0],
                    'x': float(parts[1]),
                    'y': float(parts[2]),
                    'z': float(parts[3]),
                    'mulliken_charge': float(parts[4])
                }
                atoms.append(atom)
        
        # Parse SMILES and InChI from last line
        last_line = lines[-1].strip().split()
        smiles = last_line[0] if len(last_line) > 0 else None
        inchi = last_line[1] if len(last_line) > 1 else None
        
        # Get molecule ID from filename
        mol_id = xyz_path.stem
        
        return {
            'id': mol_id,
            'n_atoms': n_atoms,
            'atoms': atoms,
            'properties': properties,
            'smiles': smiles,
            'inchi': inchi
        }
        
    except Exception as e:
        logger.error(f"Failed to parse {xyz_path}: {e}")
        return None


def load_qm9_molecule(cache_dir: Path, mol_id: str) -> Optional[Dict]:
    """
    Load a single QM9 molecule by ID.
    
    Args:
        cache_dir: Directory containing extracted QM9 data
        mol_id: Molecule ID (e.g., "dsgdb9nsd_000001")
        
    Returns:
        Dictionary with molecule data, or None if not found
    """
    cache_dir = Path(cache_dir)
    mol_path = cache_dir / "qm9_molecules" / f"{mol_id}.xyz"
    
    if not mol_path.exists():
        logger.error(f"Molecule {mol_id} not found at {mol_path}")
        return None
    
    return parse_qm9_xyz(mol_path)


def search_qm9_by_property(cache_dir: Path, 
                          property_name: str,
                          min_value: Optional[float] = None,
                          max_value: Optional[float] = None,
                          limit: Optional[int] = None) -> List[Dict]:
    """
    Search QM9 molecules by property range.
    
    Args:
        cache_dir: Directory containing extracted QM9 data
        property_name: Name of property to search
        min_value: Minimum property value
        max_value: Maximum property value
        limit: Maximum number of results
        
    Returns:
        List of molecules matching criteria
    """
    if property_name not in QM9_PROPERTIES:
        logger.error(f"Invalid property name: {property_name}")
        return []
    
    cache_dir = Path(cache_dir)
    mol_dir = cache_dir / "qm9_molecules"
    
    if not mol_dir.exists():
        logger.error(f"QM9 molecules not found at {mol_dir}")
        return []
    
    results = []
    count = 0
    
    # Iterate through molecule files
    for xyz_file in sorted(mol_dir.glob("*.xyz")):
        if limit and count >= limit:
            break
        
        mol_data = parse_qm9_xyz(xyz_file)
        if not mol_data:
            continue
        
        prop_value = mol_data['properties'].get(property_name)
        if prop_value is None:
            continue
        
        # Check property range
        if min_value is not None and prop_value < min_value:
            continue
        if max_value is not None and prop_value > max_value:
            continue
        
        results.append(mol_data)
        count += 1
    
    logger.info(f"Found {len(results)} molecules with {property_name} in range")
    return results


def get_qm9_statistics(cache_dir: Path) -> Optional[Dict]:
    """
    Calculate statistics for QM9 dataset.
    
    Args:
        cache_dir: Directory containing extracted QM9 data
        
    Returns:
        Dictionary with dataset statistics
    """
    cache_dir = Path(cache_dir)
    mol_dir = cache_dir / "qm9_molecules"
    
    if not mol_dir.exists():
        logger.error(f"QM9 molecules not found at {mol_dir}")
        return None
    
    # Count molecules
    xyz_files = list(mol_dir.glob("*.xyz"))
    n_molecules = len(xyz_files)
    
    if n_molecules == 0:
        return {'n_molecules': 0}
    
    # Initialize property statistics
    prop_stats = {prop: {'min': float('inf'), 'max': float('-inf'), 'sum': 0, 'count': 0} 
                  for prop in QM9_PROPERTIES}
    
    # Calculate statistics
    for xyz_file in xyz_files[:1000]:  # Sample first 1000 for speed
        mol_data = parse_qm9_xyz(xyz_file)
        if not mol_data:
            continue
        
        for prop, value in mol_data['properties'].items():
            if prop in prop_stats:
                prop_stats[prop]['min'] = min(prop_stats[prop]['min'], value)
                prop_stats[prop]['max'] = max(prop_stats[prop]['max'], value)
                prop_stats[prop]['sum'] += value
                prop_stats[prop]['count'] += 1
    
    # Calculate means
    for prop in prop_stats:
        if prop_stats[prop]['count'] > 0:
            prop_stats[prop]['mean'] = prop_stats[prop]['sum'] / prop_stats[prop]['count']
        else:
            prop_stats[prop]['mean'] = None
    
    return {
        'n_molecules': n_molecules,
        'properties': prop_stats,
        'sampled': min(1000, n_molecules)
    }


def create_qm9_index(cache_dir: Path, force_rebuild: bool = False) -> bool:
    """
    Create an index file for fast molecule lookup.
    
    Args:
        cache_dir: Directory containing extracted QM9 data
        force_rebuild: Force rebuild even if index exists
        
    Returns:
        True if successful, False otherwise
    """
    cache_dir = Path(cache_dir)
    mol_dir = cache_dir / "qm9_molecules"
    index_path = cache_dir / "qm9_index.json"
    
    if index_path.exists() and not force_rebuild:
        logger.info(f"QM9 index already exists at {index_path}")
        return True
    
    if not mol_dir.exists():
        logger.error(f"QM9 molecules not found at {mol_dir}")
        return False
    
    logger.info("Building QM9 index...")
    
    index = {
        'molecules': {},
        'by_smiles': {},
        'by_n_atoms': {}
    }
    
    for xyz_file in mol_dir.glob("*.xyz"):
        mol_data = parse_qm9_xyz(xyz_file)
        if not mol_data:
            continue
        
        mol_id = mol_data['id']
        
        # Add to main index
        index['molecules'][mol_id] = {
            'smiles': mol_data['smiles'],
            'n_atoms': mol_data['n_atoms'],
            'properties': mol_data['properties']
        }
        
        # Index by SMILES
        if mol_data['smiles']:
            if mol_data['smiles'] not in index['by_smiles']:
                index['by_smiles'][mol_data['smiles']] = []
            index['by_smiles'][mol_data['smiles']].append(mol_id)
        
        # Index by number of atoms
        n_atoms = str(mol_data['n_atoms'])
        if n_atoms not in index['by_n_atoms']:
            index['by_n_atoms'][n_atoms] = []
        index['by_n_atoms'][n_atoms].append(mol_id)
    
    # Save index
    try:
        with open(index_path, 'w') as f:
            json.dump(index, f, indent=2)
        logger.info(f"QM9 index created with {len(index['molecules'])} molecules")
        return True
    except Exception as e:
        logger.error(f"Failed to save QM9 index: {e}")
        return False


def load_qm9_index(cache_dir: Path) -> Optional[Dict]:
    """
    Load QM9 index for fast lookups.
    
    Args:
        cache_dir: Directory containing QM9 index
        
    Returns:
        Index dictionary, or None if not found
    """
    cache_dir = Path(cache_dir)
    index_path = cache_dir / "qm9_index.json"
    
    if not index_path.exists():
        logger.warning(f"QM9 index not found at {index_path}")
        return None
    
    try:
        with open(index_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        logger.error(f"Failed to load QM9 index: {e}")
        return None


def ensure_qm9_ready(cache_dir: Path, auto_download: bool = True) -> bool:
    """
    Ensure QM9 dataset is downloaded and extracted.
    
    This is a convenience function that handles the full workflow.
    
    Args:
        cache_dir: Directory for QM9 data
        auto_download: Automatically download if not present (default: True)
        
    Returns:
        True if QM9 is ready to use, False otherwise
    """
    cache_dir = Path(cache_dir)
    dataset_path = cache_dir / QM9_FILENAME
    extract_dir = cache_dir / "qm9_molecules"
    
    # Check if already extracted and ready
    if extract_dir.exists() and any(extract_dir.glob("*.xyz")):
        return True
    
    # Check if archive exists
    if dataset_path.exists():
        # Extract it
        logger.info("QM9 archive found, extracting...")
        return extract_qm9_dataset(cache_dir)
    
    # Need to download
    if auto_download:
        logger.info("QM9 not found. Downloading...")
        if download_qm9_dataset(cache_dir):
            # Now extract
            logger.info("Download complete. Extracting...")
            return extract_qm9_dataset(cache_dir)
        else:
            logger.error("Failed to download QM9 dataset")
            return False
    else:
        logger.warning("QM9 not found. Set auto_download=True or call download_qm9_dataset() first.")
        return False


def get_qm9_molecule_with_extraction(cache_dir: Path, mol_id: str) -> Optional[Dict]:
    """
    Load a QM9 molecule, ensuring dataset is extracted first.
    
    Args:
        cache_dir: Directory containing QM9 data
        mol_id: Molecule ID (e.g., "dsgdb9nsd_000001" or just "1")
        
    Returns:
        Dictionary with molecule data, or None if not found
    """
    # Ensure dataset is ready
    if not ensure_qm9_ready(cache_dir):
        return None
    
    # Handle different ID formats
    if isinstance(mol_id, int) or mol_id.isdigit():
        mol_id = f"dsgdb9nsd_{int(mol_id):06d}"
    
    return load_qm9_molecule(cache_dir, mol_id)


def search_qm9_with_extraction(cache_dir: Path, 
                              property_name: str,
                              min_value: Optional[float] = None,
                              max_value: Optional[float] = None,
                              limit: Optional[int] = None) -> List[Dict]:
    """
    Search QM9 molecules by property, ensuring dataset is extracted first.
    
    Args:
        cache_dir: Directory containing QM9 data
        property_name: Name of property to search
        min_value: Minimum property value
        max_value: Maximum property value
        limit: Maximum number of results
        
    Returns:
        List of molecules matching criteria
    """
    # Ensure dataset is ready
    if not ensure_qm9_ready(cache_dir):
        return []
    
    return search_qm9_by_property(cache_dir, property_name, min_value, max_value, limit)