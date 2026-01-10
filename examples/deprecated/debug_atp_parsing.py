#!/usr/bin/env python3
"""
Debug ATP parsing issue.
"""

import gzip
from pathlib import Path

# CCD file path
ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd/components.cif.gz")

print("=== Debugging ATP Parsing ===")

with gzip.open(ccd_file, 'rt') as f:
    in_atp_block = False
    atp_lines = []
    
    for line in f:
        if line.strip() == "data_ATP":
            in_atp_block = True
            atp_lines = [line]
        elif in_atp_block:
            if line.startswith("data_") and line.strip() != "data_ATP":
                break
            atp_lines.append(line)
    
    # Look for _chem_comp lines
    print("\n_chem_comp lines:")
    for line in atp_lines:
        if line.strip().startswith("_chem_comp."):
            print(f"  {line.rstrip()}")
    
    # Look for SMILES lines
    print("\nSMILES lines:")
    for line in atp_lines:
        if "SMILES" in line and line.strip().startswith("ATP"):
            print(f"  {line.rstrip()}")
    
    # Show how the name is stored
    print("\nName-related lines:")
    for line in atp_lines:
        if "_chem_comp.name" in line:
            idx = atp_lines.index(line)
            print(f"  Line {idx}: {line.rstrip()}")
            # Show next few lines
            for i in range(1, 5):
                if idx + i < len(atp_lines):
                    print(f"  Line {idx+i}: {atp_lines[idx+i].rstrip()}")