#!/usr/bin/env python3
"""
Fix CCD performance by building and using an index.
"""

import os
import sys
from pathlib import Path
import time

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.ingest import ccd_loader

# CCD directory
ccd_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd")

print("=== CCD Performance Fix ===")

# Check if index exists
index_path = ccd_dir / "ccd_index.json"
if not index_path.exists():
    print("Building CCD index for fast lookups...")
    start = time.time()
    success = ccd_loader.create_ccd_index(ccd_dir)
    if success:
        print(f"✓ Index built in {time.time() - start:.1f} seconds")
    else:
        print("✗ Failed to build index")
        sys.exit(1)
else:
    print("✓ CCD index already exists")

# Load the index
print("\nLoading CCD index...")
index = ccd_loader.load_ccd_index(ccd_dir)
if index:
    print(f"✓ Index loaded with {len(index.get('components', {}))} components")
    
    # Show some stats
    print(f"\nIndex contents:")
    print(f"  Components: {len(index.get('components', {}))}")
    print(f"  By type: {len(index.get('by_type', {}))} types")
    print(f"  By formula: {len(index.get('by_formula', {}))} unique formulas")
    print(f"  With SMILES: {len(index.get('has_smiles', []))} components")
else:
    print("✗ Failed to load index")
    
print("\n✅ Index ready for fast lookups!")