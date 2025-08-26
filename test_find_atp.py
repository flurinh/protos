#!/usr/bin/env python3
"""
Find ATP in CCD file to understand its format.
"""

import gzip
from pathlib import Path

# CCD file path
ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd/components.cif.gz")

if ccd_file.exists():
    print(f"Searching for ATP in CCD file...")
    
    with gzip.open(ccd_file, 'rt') as f:
        in_atp_block = False
        atp_lines = []
        line_count = 0
        
        for line in f:
            line_count += 1
            
            # Check if we're starting ATP block
            if line.strip() == "data_ATP":
                in_atp_block = True
                atp_lines = [line]
                print(f"Found ATP at line {line_count}")
            elif in_atp_block:
                # Check if we're starting a new block
                if line.startswith("data_") and line.strip() != "data_ATP":
                    break
                atp_lines.append(line)
        
        if atp_lines:
            print(f"\n=== ATP block ({len(atp_lines)} lines) ===")
            for i, line in enumerate(atp_lines[:200]):  # First 200 lines
                print(f"{i+1:3d}: {line.rstrip()}")
        else:
            print("ATP not found in file")
else:
    print(f"CCD file not found at {ccd_file}")