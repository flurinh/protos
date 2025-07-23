#!/usr/bin/env python3
"""
Download real PDB structures for testing.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from protos.loaders.download_structures import download_protein_structures

# List of small, reliable PDB structures for testing
TEST_PDB_IDS = [
    "1ubq",  # Ubiquitin - small, well-studied protein
    "1crn",  # Crambin - one of the smallest proteins
    "1l2y",  # Trp-cage miniprotein - very small
    "2gb1",  # Protein G B1 domain - small, stable
    "1tqn",  # Small peptide
    "1uaz",  # Small protein
    "3nir",  # Small opsin fragment
]

def main():
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/protos_data")
    structure_dir = data_dir / "structure" / "mmcif"
    structure_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Downloading test structures to {structure_dir}")
    
    # Check which structures already exist
    existing = []
    missing = []
    
    for pdb_id in TEST_PDB_IDS:
        cif_file = structure_dir / f"{pdb_id}.cif"
        if cif_file.exists() and cif_file.stat().st_size > 1000:  # Real CIF files are > 1KB
            existing.append(pdb_id)
        else:
            missing.append(pdb_id)
    
    print(f"Existing structures: {existing}")
    print(f"Missing structures: {missing}")
    
    if missing:
        print(f"\nDownloading {len(missing)} structures...")
        try:
            download_protein_structures(
                missing,
                target_folder=str(structure_dir)
            )
            print("Download complete!")
        except Exception as e:
            print(f"Error downloading structures: {e}")
            return 1
    else:
        print("All test structures already available!")
    
    # Also ensure test-data directory has these structures
    test_data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/tests/test-data/structure/mmcif")
    test_data_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\nCopying to test-data directory: {test_data_dir}")
    import shutil
    for pdb_id in TEST_PDB_IDS:
        src_file = structure_dir / f"{pdb_id}.cif"
        dst_file = test_data_dir / f"{pdb_id}.cif"
        if src_file.exists() and not dst_file.exists():
            shutil.copy2(src_file, dst_file)
            print(f"Copied {pdb_id}.cif")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())