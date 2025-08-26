#!/usr/bin/env python3
"""
Test CCD file format to understand SMILES extraction.
"""

import gzip
from pathlib import Path

ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/cache/databases/ccd/components.cif.gz")

print("=== Examining CCD File Format ===\n")

# Read first few components to understand format
with gzip.open(ccd_file, 'rt') as f:
    lines_read = 0
    current_component = None
    in_descriptor_loop = False
    descriptor_headers = []
    
    for line in f:
        lines_read += 1
        
        # Track component
        if line.startswith('data_'):
            if current_component and lines_read > 1000:
                break  # Stop after examining a few components
            current_component = line.strip()
            print(f"\n{current_component}")
            in_descriptor_loop = False
            descriptor_headers = []
        
        # Look for descriptor loop
        if '_pdbx_chem_comp_descriptor.' in line:
            if not in_descriptor_loop:
                print("  Found descriptor section:")
                in_descriptor_loop = True
            
            # Collect headers
            if line.strip().startswith('_pdbx_chem_comp_descriptor.'):
                header = line.strip()
                descriptor_headers.append(header)
                print(f"    Header: {header}")
        
        # Look for SMILES in descriptor data
        elif in_descriptor_loop and line.strip() and not line.startswith('_') and not line.startswith('#'):
            # This should be data
            if 'SMILES' in line:
                print(f"    SMILES line: {line.strip()[:100]}...")
            elif current_component and current_component.replace('data_', '') in line:
                # First line of data usually contains the component ID
                print(f"    Data line: {line.strip()[:100]}...")
                
        # Exit descriptor loop
        elif in_descriptor_loop and (line.startswith('_') and '_pdbx_chem_comp_descriptor.' not in line):
            in_descriptor_loop = False
            
        if lines_read > 5000:
            break

print(f"\nExamined {lines_read} lines")