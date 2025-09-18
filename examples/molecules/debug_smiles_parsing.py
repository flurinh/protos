#!/usr/bin/env python3
"""
Debug SMILES parsing for ATP.
"""

import os
import sys
from pathlib import Path

# Add src to path  
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.loaders.ccd_loader import parse_ccd_cif_block
import gzip

# CCD file path
ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd/components.cif.gz")

print("=== Debugging SMILES Parsing ===")

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
    
    # Parse with our function
    parsed = parse_ccd_cif_block(atp_lines)
    
    print(f"Parsed data keys: {list(parsed.keys())}")
    print(f"\nID: {parsed.get('id')}")
    print(f"Name: {parsed.get('name')}")
    print(f"Formula: {parsed.get('formula')}")
    
    # Check for SMILES in parsed data
    print("\nSMILES-related fields:")
    for key, value in parsed.items():
        if 'smiles' in key.lower():
            print(f"  {key}: {value[:50]}...")
    
    # Show what the SMILES lines look like
    print("\nRaw SMILES lines from file:")
    for i, line in enumerate(atp_lines):
        if "SMILES" in line and line.strip().startswith("ATP"):
            print(f"  Line {i}: {repr(line.rstrip())}")
    
    # Test if the condition matches
    print(f"\nID from parsed data: {repr(parsed.get('id'))}")
    test_line = "ATP SMILES"
    print(f"Does '{test_line}' start with ID? {test_line.startswith(parsed.get('id', ''))}")