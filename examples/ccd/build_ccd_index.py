#!/usr/bin/env python3
"""
Build fast CCD index using gemmi.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.loaders.ccd_index_builder import build_ccd_index_with_gemmi

# CCD paths
ccd_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd")
ccd_file = ccd_dir / "components.cif.gz"

print("=== Building Fast CCD Index ===")
print(f"CCD file: {ccd_file}")
print(f"Output directory: {ccd_dir}")

if not ccd_file.exists():
    print(f"✗ CCD file not found at {ccd_file}")
    print("Please download CCD first.")
    sys.exit(1)

print(f"\nFile size: {ccd_file.stat().st_size / 1024 / 1024:.1f} MB")

# Build the index
print("\nBuilding index with gemmi...")
success = build_ccd_index_with_gemmi(ccd_file, ccd_dir)

if success:
    print("\n✅ CCD index built successfully!")
    
    # Check created files
    index_file = ccd_dir / "ccd_index.json"
    data_file = ccd_dir / "ccd_ligands.pkl"
    
    if index_file.exists():
        print(f"\nIndex file: {index_file} ({index_file.stat().st_size / 1024 / 1024:.1f} MB)")
    if data_file.exists():
        print(f"Data file: {data_file} ({data_file.stat().st_size / 1024 / 1024:.1f} MB)")
else:
    print("\n✗ Failed to build CCD index")
    sys.exit(1)