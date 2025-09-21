#!/usr/bin/env python3
"""
Debug CCD component loading.
"""

import os
import sys
from pathlib import Path
import json

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
ccd_dir = Path(paths.get_subdir_path("ligand", "cache_dir")) / "databases" / "ccd"

print(f"CCD directory: {ccd_dir}")
print(f"CCD ready: {ccd_loader.is_ccd_ready(ccd_dir)}")

# Try to load ATP
print("\nLoading ATP component...")
atp = ccd_loader.get_ccd_component(ccd_dir, 'ATP')

if atp:
    print(f"ATP found!")
    print(f"Keys: {list(atp.keys())}")
    for key, value in atp.items():
        if isinstance(value, str) and len(value) > 50:
            print(f"{key}: {value[:50]}...")
        else:
            print(f"{key}: {value}")
else:
    print("ATP not found!")

# Check the index file
print("\n\nChecking index file...")
index_path = ccd_dir / "ccd_index.json"
if index_path.exists():
    with open(index_path, 'r') as f:
        index = json.load(f)
    
    # Check ATP in index
    if 'ATP' in index['components']:
        print(f"ATP in index: {index['components']['ATP']}")
    
    # Check stats
    print(f"\nIndex stats:")
    print(f"  Total: {len(index['components'])}")
    print(f"  With SMILES: {len(index.get('with_smiles', []))}")
    
    # Show a few components with SMILES
    if index.get('with_smiles'):
        print(f"\nFirst 5 components with SMILES: {index['with_smiles'][:5]}")