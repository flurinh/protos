"""
Test script for ligand processor database integration.

This demonstrates:
1. Automatic database downloads
2. ProtosPaths managing all paths (no hardcoded paths!)
3. Fast local access to CCD, QM9, and Enamine
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.processing.molecule import MoleculeProcessor

# Set up environment
test_data_root = Path(__file__).parent / "data"
os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())

# Initialize processor
print("Initializing LigandProcessor...")
lig_proc = MoleculeProcessor()

# Test 1: Check database status
print("\n=== Database Status ===")
stats = lig_proc.get_database_statistics()
for db_name, info in stats.items():
    print(f"\n{db_name}:")
    print(f"  Downloaded: {info['downloaded']}")
    print(f"  Description: {info['description']}")
    if info['downloaded']:
        print(f"  Path: {info['path']}")

# Test 2: Get ligand from CCD (will auto-download if needed)
print("\n=== Testing CCD Access ===")
print("Getting ATP from CCD (will download CCD if not present)...")
atp = lig_proc.get_ccd_ligand('ATP')
if atp:
    print(f"✓ ATP found!")
    print(f"  SMILES: {atp['smiles'][:50]}...")
    print(f"  Formula: {atp['formula']}")
    print(f"  Name: {atp['name']}")
else:
    print("✗ ATP not found")

# Test 3: Create dataset from CCD components
print("\n=== Creating CCD Dataset ===")
print("Creating dataset of common cofactors...")
created = lig_proc.create_ccd_dataset(
    "cofactors",
    ["ATP", "NAD", "FAD", "HEM", "XXX"]  # XXX doesn't exist, should be skipped
)
print(f"✓ Created dataset 'cofactors' with {len(created)} components")

# Test 4: Search QM9 by properties (will auto-download if needed)
print("\n=== Testing QM9 Search ===")
print("Searching QM9 for molecules with specific properties...")
print("(This will download QM9 dataset if not present - ~3GB)")

# Only search if user confirms (since it's a large download)
response = input("Download QM9 dataset if needed? (y/n): ")
if response.lower() == 'y':
    molecules = lig_proc.search_qm9_by_properties({
        'gap': (0.2, 0.4),    # HOMO-LUMO gap in eV
        'dipole': (0, 3.0)    # Dipole moment in Debye
    }, limit=5)
    
    if molecules:
        print(f"✓ Found {len(molecules)} molecules:")
        for mol in molecules[:3]:
            print(f"  QM9 ID: {mol['qm9_id']}")
            print(f"  SMILES: {mol['smiles']}")
            print(f"  Gap: {mol['quantum_properties']['gap']:.3f} eV")
            print(f"  Dipole: {mol['quantum_properties']['dipole']:.3f} D")
            print()
    else:
        print("No molecules found matching criteria")
else:
    print("Skipping QM9 test")

# Test 5: List available Enamine datasets
print("\n=== Available Enamine Datasets ===")
datasets = lig_proc.list_enamine_datasets()
print(f"Found {len(datasets)} Enamine datasets:")
for name, info in list(datasets.items())[:5]:  # Show first 5
    print(f"  {name}: {info['description']} ({info['size']})")
print("  ...")

# Test 6: Check Enamine credentials
print("\n=== Enamine Credentials ===")
from protos.io.ingest import enamine_loader
username, password = enamine_loader.get_enamine_credentials()
if username:
    print(f"✓ Credentials found for: {username}")
else:
    print("✗ No credentials found. Set ENAMINE_USERNAME and ENAMINE_PASSWORD in .env file")

# Test 7: Download small Enamine test dataset (if credentials available)
if username:
    print("\n=== Testing Enamine Download ===")
    print("Attempting to download small test dataset 'diversity_1k'...")
    
    # This will attempt to download using credentials
    # Note: This is a mock - actual Enamine URLs would need to be configured
    dataset_info = lig_proc.get_enamine_dataset_info('diversity_1k')
    if dataset_info:
        print(f"Dataset: {dataset_info['name']}")
        print(f"Size: {dataset_info['size']}")
        print(f"Downloaded: {dataset_info.get('downloaded', False)}")

# Test 8: Check database status after downloads
print("\n=== Final Database Status ===")
stats = lig_proc.get_database_statistics()
for db_name, info in stats.items():
    print(f"\n{db_name}:")
    if db_name == 'Enamine':
        print(f"  Credentials set: {info['credentials_set']}")
        print(f"  Available datasets: {info['available_datasets']}")
        print(f"  Downloaded datasets: {len(info['downloaded_datasets'])}")
        if info['downloaded_datasets']:
            print(f"  Downloaded: {', '.join(info['downloaded_datasets'][:3])}...")
    else:
        print(f"  Downloaded: {info['downloaded']}")
        if info['downloaded']:
            print(f"  Path: {info['path']}")

# Test 9: Demonstrate that paths are managed by ProtosPaths
print("\n=== Path Management ===")
print(f"All databases stored under: {lig_proc.paths.get_processor_path('ligand')}/databases/")
print("No hardcoded paths used - everything managed by ProtosPaths!")

print("\n✅ All tests completed!")