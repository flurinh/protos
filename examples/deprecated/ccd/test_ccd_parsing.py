#!/usr/bin/env python3
"""
Test CCD parsing to understand the format.
"""

import gzip
from pathlib import Path

# CCD file path
ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd/components.cif.gz")

if ccd_file.exists():
    print(f"CCD file found: {ccd_file}")
    print(f"Size: {ccd_file.stat().st_size / 1024 / 1024:.1f} MB")
    
    # Read first few blocks to understand format
    with gzip.open(ccd_file, 'rt') as f:
        print("\n=== First 100 lines of CCD ===")
        for i, line in enumerate(f):
            print(f"{i+1:3d}: {line.rstrip()}")
            if i >= 100:
                break
else:
    print(f"CCD file not found at {ccd_file}")