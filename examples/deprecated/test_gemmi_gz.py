#!/usr/bin/env python3
"""
Test if gemmi can read gz files directly.
"""

import gemmi
from pathlib import Path

ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd/components.cif.gz")

print(f"Testing gemmi with: {ccd_file}")
print(f"File exists: {ccd_file.exists()}")
print(f"File size: {ccd_file.stat().st_size / 1024 / 1024:.1f} MB")

# Test if gemmi can read .gz files directly
try:
    print("\nTrying to read .gz file directly with gemmi...")
    doc = gemmi.cif.read(str(ccd_file))
    print(f"✓ Success! Loaded {len(doc)} blocks")
    
    # Test accessing a few components
    print("\nTesting access to components:")
    for code in ['ATP', 'NAD', 'HEM']:
        if code in doc:
            block = doc[code]
            print(f"  {code}: Found")
        else:
            print(f"  {code}: Not found")
            
except Exception as e:
    print(f"✗ Failed: {e}")
    
# Also test the integrity of the gz file
print("\n\nTesting gzip integrity...")
import gzip
try:
    with gzip.open(ccd_file, 'rb') as f:
        # Try to read in chunks
        total = 0
        while True:
            chunk = f.read(1024*1024)  # 1MB chunks
            if not chunk:
                break
            total += len(chunk)
        print(f"✓ Successfully read {total / 1024 / 1024:.1f} MB")
except Exception as e:
    print(f"✗ Failed: {e}")