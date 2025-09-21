"""
Enamine REAL database loader for Protos.

The Enamine REAL database contains billions of commercially available compounds.
This loader provides utilities for downloading and accessing Enamine REAL subsets.

Note: Due to the massive size of the full database (>30B compounds), this loader
works with curated subsets like REAL Space Navigator, diversity sets, etc.

Environment variables (set in .env file):
- ENAMINE_USERNAME: Your Enamine account email
- ENAMINE_PASSWORD: Your Enamine account password

References:
    https://enamine.net/compound-collections/real-compounds/real-database

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
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterator
import urllib.request
import urllib.error
import urllib.parse
import base64
import zipfile
from datetime import datetime
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Try to import dotenv for loading credentials
try:
    from dotenv import load_dotenv
    load_dotenv()
    HAS_DOTENV = True
except ImportError:
    logger.info("python-dotenv not available. Set environment variables manually.")
    HAS_DOTENV = False

# Try to import RDKit for molecular operations
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    HAS_RDKIT = True
except ImportError:
    logger.warning("RDKit not available. Some Enamine functionality will be limited.")
    HAS_RDKIT = False

# Enamine REAL subsets with download information
ENAMINE_DATASETS = {
    # Small test datasets (good for development)
    'hit2lead_1k': {
        'name': 'Hit2Lead Set 1K (Test)',
        'description': '1,000 lead-like compounds for testing',
        'size': '1K',
        'filename': 'Hit2Lead_1000_test.sdf.gz',
        'url_path': 'REAL-Database/Sets/Hit2Lead_1000_test.sdf.gz'
    },
    'diversity_1k': {
        'name': 'Diversity Set 1K (Test)',
        'description': '1,000 diverse compounds for testing',
        'size': '1K',
        'filename': 'Diversity_1000_test.sdf.gz',
        'url_path': 'REAL-Database/Sets/Diversity_1000_test.sdf.gz'
    },
    
    # Medium datasets
    'hit2lead_10k': {
        'name': 'Hit2Lead Set 10K',
        'description': '10,000 lead-like compounds from REAL database',
        'size': '10K',
        'filename': 'Hit2Lead_10000.sdf.gz',
        'url_path': 'REAL-Database/Sets/Hit2Lead_10000.sdf.gz'
    },
    'diversity_10k': {
        'name': 'REAL Diversity Set 10K',
        'description': '10,000 diverse compounds from REAL database',
        'size': '10K',
        'filename': 'Diversity_10000.sdf.gz',
        'url_path': 'REAL-Database/Sets/Diversity_10000.sdf.gz'
    },
    'fragments_5k': {
        'name': 'Fragment Library 5K',
        'description': '5,000 fragment compounds (MW < 300)',
        'size': '5K',
        'filename': 'Fragments_5000.sdf.gz',
        'url_path': 'REAL-Database/Sets/Fragments_5000.sdf.gz'
    },
    
    # Large datasets
    'hit2lead_100k': {
        'name': 'Hit2Lead Set 100K',
        'description': '100,000 lead-like compounds',
        'size': '100K',
        'filename': 'Hit2Lead_100000.sdf.gz',
        'url_path': 'REAL-Database/Sets/Hit2Lead_100000.sdf.gz'
    },
    'diversity_100k': {
        'name': 'REAL Diversity Set 100K',
        'description': '100,000 diverse compounds',
        'size': '100K',
        'filename': 'Diversity_100000.sdf.gz',
        'url_path': 'REAL-Database/Sets/Diversity_100000.sdf.gz'
    },
    
    # Specialized libraries
    'kinase_focused': {
        'name': 'Kinase-Focused Library',
        'description': 'Compounds designed for kinase targets',
        'size': '50K',
        'filename': 'Kinase_Focused_50000.sdf.gz',
        'url_path': 'REAL-Database/Sets/Kinase_50000.sdf.gz'
    },
    'gpcr_focused': {
        'name': 'GPCR-Focused Library',
        'description': 'Compounds designed for GPCR targets',
        'size': '50K',
        'filename': 'GPCR_Focused_50000.sdf.gz',
        'url_path': 'REAL-Database/Sets/GPCR_50000.sdf.gz'
    },
    'ppi_focused': {
        'name': 'PPI-Focused Library',
        'description': 'Compounds for protein-protein interactions',
        'size': '30K',
        'filename': 'PPI_Focused_30000.sdf.gz',
        'url_path': 'REAL-Database/Sets/PPI_30000.sdf.gz'
    },
    
    # REAL Space Navigator subsets
    'real_space_navigator_10M': {
        'name': 'REAL Space Navigator 10M',
        'description': '10 million compounds covering chemical space',
        'size': '10M',
        'filename': 'REAL_Space_Navigator_10M.sdf.gz',
        'url_path': 'REAL-Database/Navigator/RSN_10M.sdf.gz',
        'large': True  # Warning: large file
    }
}

# Base URL for Enamine downloads (this is a placeholder)
ENAMINE_BASE_URL = "https://download.enamine.net/"


def get_enamine_credentials() -> Tuple[Optional[str], Optional[str]]:
    """
    Get Enamine credentials from environment variables.
    
    Returns:
        Tuple of (username, password) or (None, None) if not set
    """
    username = os.environ.get('ENAMINE_USERNAME') or os.environ.get('enamine_username')
    password = os.environ.get('ENAMINE_PASSWORD') or os.environ.get('enamine_password')
    
    if not username or not password:
        logger.warning(
            "Enamine credentials not found. Set ENAMINE_USERNAME and ENAMINE_PASSWORD "
            "environment variables or create a .env file."
        )
    
    return username, password


def list_available_datasets() -> Dict[str, Dict]:
    """
    List all available Enamine datasets.
    
    Returns:
        Dictionary of dataset information
    """
    return ENAMINE_DATASETS.copy()


def get_dataset_info(dataset_name: str) -> Optional[Dict]:
    """
    Get information about a specific dataset.
    
    Args:
        dataset_name: Name of the dataset
        
    Returns:
        Dataset information or None if not found
    """
    return ENAMINE_DATASETS.get(dataset_name)


def download_enamine_dataset(
    cache_dir: Path,
    dataset_name: str,
    force: bool = False
) -> bool:
    """
    Download an Enamine dataset with authentication.
    
    Args:
        cache_dir: Directory to store the dataset
        dataset_name: Name of the dataset to download
        force: Force re-download even if exists
        
    Returns:
        True if successful, False otherwise
    """
    # Get dataset info
    dataset_info = get_dataset_info(dataset_name)
    if not dataset_info:
        logger.error(f"Unknown dataset: {dataset_name}")
        logger.info(f"Available datasets: {', '.join(ENAMINE_DATASETS.keys())}")
        return False
    
    # Check if already downloaded
    dataset_dir = cache_dir / dataset_name
    marker_file = dataset_dir / ".downloaded"
    
    if marker_file.exists() and not force:
        logger.info(f"Dataset {dataset_name} already downloaded")
        return True
    
    # Get credentials
    username, password = get_enamine_credentials()
    if not username or not password:
        logger.error("Enamine credentials required for download")
        return False
    
    # Create dataset directory
    dataset_dir.mkdir(parents=True, exist_ok=True)
    
    # Construct download URL
    url = urllib.parse.urljoin(ENAMINE_BASE_URL, dataset_info['url_path'])
    output_file = dataset_dir / dataset_info['filename']
    
    logger.info(f"Downloading {dataset_name} from Enamine...")
    logger.info(f"Dataset: {dataset_info['name']}")
    logger.info(f"Size: {dataset_info['size']} compounds")
    
    if dataset_info.get('large'):
        logger.warning("This is a large dataset. Download may take significant time.")
    
    try:
        # Create request with authentication
        request = urllib.request.Request(url)
        credentials = base64.b64encode(f"{username}:{password}".encode()).decode()
        request.add_header("Authorization", f"Basic {credentials}")
        
        # Download with progress bar
        with urllib.request.urlopen(request) as response:
            total_size = int(response.headers.get('Content-Length', 0))
            
            with open(output_file, 'wb') as f:
                with tqdm(total=total_size, unit='B', unit_scale=True,
                         desc=f"Downloading {dataset_name}") as pbar:
                    while True:
                        chunk = response.read(8192)
                        if not chunk:
                            break
                        f.write(chunk)
                        pbar.update(len(chunk))
        
        # Extract if compressed
        if output_file.suffix == '.gz':
            logger.info("Extracting compressed file...")
            extracted_file = output_file.with_suffix('')
            
            # Get file size for progress bar
            file_size = output_file.stat().st_size
            
            with gzip.open(output_file, 'rb') as gz_file:
                with open(extracted_file, 'wb') as out_file:
                    with tqdm(total=file_size, unit='B', unit_scale=True,
                             desc="Extracting") as pbar:
                        while True:
                            chunk = gz_file.read(8192)
                            if not chunk:
                                break
                            out_file.write(chunk)
                            pbar.update(len(chunk))
            
            # Remove compressed file to save space
            output_file.unlink()
        
        # Create marker file with metadata
        metadata = {
            'dataset': dataset_name,
            'downloaded': datetime.now().isoformat(),
            'info': dataset_info
        }
        with open(marker_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Successfully downloaded {dataset_name}")
        return True
        
    except urllib.error.HTTPError as e:
        if e.code == 401:
            logger.error("Authentication failed. Check your Enamine credentials.")
        else:
            logger.error(f"HTTP error {e.code}: {e.reason}")
        return False
    except Exception as e:
        logger.error(f"Download failed: {e}")
        return False


def is_dataset_downloaded(cache_dir: Path, dataset_name: str) -> bool:
    """
    Check if an Enamine dataset is downloaded.
    
    Args:
        cache_dir: Base cache directory
        dataset_name: Name of the dataset
        
    Returns:
        True if downloaded, False otherwise
    """
    dataset_dir = cache_dir / dataset_name
    marker_file = dataset_dir / ".downloaded"
    return marker_file.exists()


def load_dataset_metadata(cache_dir: Path, dataset_name: str) -> Optional[Dict]:
    """
    Load metadata about a downloaded dataset.
    
    Args:
        cache_dir: Base cache directory
        dataset_name: Name of the dataset
        
    Returns:
        Metadata dictionary or None
    """
    dataset_dir = cache_dir / dataset_name
    marker_file = dataset_dir / ".downloaded"
    
    if marker_file.exists():
        with open(marker_file, 'r') as f:
            return json.load(f)
    return None


def stream_enamine_compounds(
    cache_dir: Path,
    dataset_name: str,
    include_properties: bool = True
) -> Iterator[Dict]:
    """
    Stream compounds from an Enamine dataset.
    
    Args:
        cache_dir: Base cache directory
        dataset_name: Name of the dataset
        include_properties: Calculate molecular properties
        
    Yields:
        Dictionary with compound data
    """
    if not is_dataset_downloaded(cache_dir, dataset_name):
        logger.error(f"Dataset {dataset_name} not downloaded")
        return
    
    dataset_info = get_dataset_info(dataset_name)
    if not dataset_info:
        logger.error(f"Unknown dataset: {dataset_name}")
        return
    
    dataset_dir = cache_dir / dataset_name
    sdf_file = dataset_dir / dataset_info['filename'].replace('.gz', '')
    
    if not sdf_file.exists():
        logger.error(f"SDF file not found: {sdf_file}")
        return
    
    if not HAS_RDKIT:
        logger.error("RDKit required for reading SDF files")
        return
    
    # Stream compounds from SDF
    supplier = Chem.SDMolSupplier(str(sdf_file))
    
    for idx, mol in enumerate(supplier):
        if mol is None:
            continue
        
        try:
            compound = {
                'id': mol.GetProp('_Name') if mol.HasProp('_Name') else f"ENAMINE_{idx}",
                'smiles': Chem.MolToSmiles(mol),
                'dataset': dataset_name
            }
            
            # Add any properties from SDF
            for prop_name in mol.GetPropNames():
                compound[prop_name.lower()] = mol.GetProp(prop_name)
            
            # Calculate molecular properties if requested
            if include_properties and HAS_RDKIT:
                compound['mw'] = Descriptors.MolWt(mol)
                compound['logp'] = Descriptors.MolLogP(mol)
                compound['hba'] = Descriptors.NumHAcceptors(mol)
                compound['hbd'] = Descriptors.NumHDonors(mol)
                compound['tpsa'] = Descriptors.TPSA(mol)
                compound['rotatable_bonds'] = Descriptors.NumRotatableBonds(mol)
            
            yield compound
            
        except Exception as e:
            logger.debug(f"Error processing compound {idx}: {e}")
            continue


def search_enamine_by_properties(
    cache_dir: Path,
    dataset_name: str,
    min_mw: Optional[float] = None,
    max_mw: Optional[float] = None,
    min_logp: Optional[float] = None,
    max_logp: Optional[float] = None,
    max_hba: Optional[int] = None,
    max_hbd: Optional[int] = None,
    limit: Optional[int] = None
) -> List[Dict]:
    """
    Search Enamine dataset by molecular properties.
    
    Args:
        cache_dir: Base cache directory
        dataset_name: Name of the dataset
        min_mw: Minimum molecular weight
        max_mw: Maximum molecular weight
        min_logp: Minimum LogP
        max_logp: Maximum LogP
        max_hba: Maximum H-bond acceptors
        max_hbd: Maximum H-bond donors
        limit: Maximum results
        
    Returns:
        List of matching compounds
    """
    results = []
    
    for compound in stream_enamine_compounds(cache_dir, dataset_name, include_properties=True):
        # Apply filters
        if min_mw is not None and compound.get('mw', 0) < min_mw:
            continue
        if max_mw is not None and compound.get('mw', float('inf')) > max_mw:
            continue
        if min_logp is not None and compound.get('logp', -float('inf')) < min_logp:
            continue
        if max_logp is not None and compound.get('logp', float('inf')) > max_logp:
            continue
        if max_hba is not None and compound.get('hba', float('inf')) > max_hba:
            continue
        if max_hbd is not None and compound.get('hbd', float('inf')) > max_hbd:
            continue
        
        results.append(compound)
        
        if limit and len(results) >= limit:
            break
    
    return results


def get_dataset_statistics(cache_dir: Path, dataset_name: str) -> Optional[Dict]:
    """
    Calculate statistics for an Enamine dataset.
    
    Args:
        cache_dir: Base cache directory
        dataset_name: Name of the dataset
        
    Returns:
        Statistics dictionary
    """
    if not is_dataset_downloaded(cache_dir, dataset_name):
        logger.error(f"Dataset {dataset_name} not downloaded")
        return None
    
    if not HAS_RDKIT:
        logger.error("RDKit required for calculating statistics")
        return None
    
    # Check if we have cached statistics
    dataset_dir = cache_dir / dataset_name
    stats_file = dataset_dir / "statistics.json"
    
    if stats_file.exists():
        with open(stats_file, 'r') as f:
            return json.load(f)
    
    # Calculate statistics
    logger.info(f"Calculating statistics for {dataset_name}...")
    
    stats = {
        'total_compounds': 0,
        'mw': {'min': float('inf'), 'max': 0, 'sum': 0},
        'logp': {'min': float('inf'), 'max': -float('inf'), 'sum': 0},
        'hba': {'min': float('inf'), 'max': 0, 'sum': 0},
        'hbd': {'min': float('inf'), 'max': 0, 'sum': 0},
        'tpsa': {'min': float('inf'), 'max': 0, 'sum': 0},
        'rotatable_bonds': {'min': float('inf'), 'max': 0, 'sum': 0}
    }
    
    for compound in stream_enamine_compounds(cache_dir, dataset_name, include_properties=True):
        stats['total_compounds'] += 1
        
        # Update statistics
        for prop in ['mw', 'logp', 'hba', 'hbd', 'tpsa', 'rotatable_bonds']:
            if prop in compound:
                value = compound[prop]
                stats[prop]['min'] = min(stats[prop]['min'], value)
                stats[prop]['max'] = max(stats[prop]['max'], value)
                stats[prop]['sum'] += value
    
    # Calculate averages
    if stats['total_compounds'] > 0:
        for prop in ['mw', 'logp', 'hba', 'hbd', 'tpsa', 'rotatable_bonds']:
            stats[prop]['avg'] = stats[prop]['sum'] / stats['total_compounds']
            del stats[prop]['sum']  # Remove sum from final stats
    
    # Add metadata
    stats['dataset'] = dataset_name
    stats['calculated'] = datetime.now().isoformat()
    
    # Cache statistics
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    return stats


def create_diversity_subset(
    cache_dir: Path,
    dataset_name: str,
    output_name: str,
    size: int,
    method: str = 'random'
) -> bool:
    """
    Create a diverse subset from an Enamine dataset.
    
    Args:
        cache_dir: Base cache directory
        dataset_name: Source dataset name
        output_name: Name for the subset
        size: Number of compounds to select
        method: Selection method ('random', 'maxmin')
        
    Returns:
        True if successful
    """
    if not is_dataset_downloaded(cache_dir, dataset_name):
        logger.error(f"Dataset {dataset_name} not downloaded")
        return False
    
    if not HAS_RDKIT:
        logger.error("RDKit required for creating subsets")
        return False
    
    logger.info(f"Creating {size}-compound subset from {dataset_name}...")
    
    # For now, implement simple random selection
    # TODO: Implement MaxMin diversity selection
    
    if method != 'random':
        logger.warning(f"Method '{method}' not implemented, using random selection")
    
    import random
    
    # Collect all compounds
    compounds = []
    for compound in stream_enamine_compounds(cache_dir, dataset_name):
        compounds.append(compound)
        if len(compounds) >= size * 10:  # Get 10x for selection
            break
    
    # Random selection
    selected = random.sample(compounds, min(size, len(compounds)))
    
    # Save subset
    subset_dir = cache_dir / output_name
    subset_dir.mkdir(parents=True, exist_ok=True)
    
    output_file = subset_dir / f"{output_name}.sdf"
    
    writer = Chem.SDWriter(str(output_file))
    for compound in selected:
        mol = Chem.MolFromSmiles(compound['smiles'])
        if mol:
            mol.SetProp('_Name', compound.get('id', ''))
            writer.write(mol)
    writer.close()
    
    # Create metadata
    metadata = {
        'source_dataset': dataset_name,
        'size': len(selected),
        'method': method,
        'created': datetime.now().isoformat()
    }
    
    with open(subset_dir / ".downloaded", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"Created subset '{output_name}' with {len(selected)} compounds")
    return True