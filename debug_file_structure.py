#!/usr/bin/env python3
"""
Debug CCD file structure for ATP.
"""

import gzip
from pathlib import Path

# CCD file path
ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd/components.cif.gz")

print("=== ATP Block Structure ===")

with gzip.open(ccd_file, 'rt') as f:
    in_atp_block = False
    line_num = 0
    
    for line in f:
        if line.strip() == "data_ATP":
            in_atp_block = True
            line_num = 0
        elif in_atp_block:
            if line.startswith("data_") and line.strip() != "data_ATP":
                break
            
            line_num += 1
            # Show structure markers
            if (line.startswith("_") or 
                line.startswith("#") or 
                line.startswith("loop_") or
                line.strip().startswith("ATP") or
                not line.strip()):
                print(f"{line_num:4d}: {line.rstrip()}")