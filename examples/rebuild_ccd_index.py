#!/usr/bin/env python3
"""
Rebuild CCD index with fixed SMILES extraction.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.io.ingest import ccd_loader_unified as ccd_loader

# Set up paths
data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# Initialize ProtosPaths
paths = ProtosPaths(data_root=str(data_dir.absolute()))

# Get CCD directory
ccd_dir = Path(paths.get_subdir_path("molecule", "cache_dir")) / "databases" / "ccd"
print(f"CCD directory: {ccd_dir}")

# Force rebuild index
print("\nRebuilding CCD index with fixed SMILES extraction...")
success = ccd_loader.build_ccd_index(ccd_dir, force=True)

if success:
    print("\n✅ Index rebuilt successfully!")
    
    # Get statistics
    stats = ccd_loader.get_ccd_statistics(ccd_dir)
    print(f"\nCCD Statistics:")
    print(f"  Total components: {stats['total_components']:,}")
    print(f"  With SMILES: {stats['with_smiles']:,}")
    print(f"  With InChI: {stats['with_inchi']:,}")
    
    # Test a component with SMILES
    test_components = ['ATP', 'HEM', 'NAD']
    print(f"\nTesting components:")
    for comp_id in test_components:
        comp = ccd_loader.get_ccd_component(ccd_dir, comp_id)
        if comp:
            print(f"\n{comp_id}:")
            print(f"  Name: {comp.get('name', 'N/A')}")
            print(f"  Formula: {comp.get('formula', 'N/A')}")
            if comp.get('smiles'):
                print(f"  SMILES: {comp['smiles'][:50]}...")
            if comp.get('smiles_canonical'):
                print(f"  SMILES (canonical): {comp['smiles_canonical'][:50]}...")
else:
    print("\n✗ Failed to rebuild index")
    sys.exit(1)