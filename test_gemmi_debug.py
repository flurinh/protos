#!/usr/bin/env python3
"""
Debug gemmi CCD extraction.
"""

import gemmi
import gzip
import tempfile
from pathlib import Path

ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/cache/databases/ccd/components.cif.gz")

# Extract just the first component for debugging
with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
    with gzip.open(ccd_file, 'rt') as f_in:
        lines_written = 0
        components = 0
        
        for line in f_in:
            tmp.write(line)
            lines_written += 1
            
            if line.startswith('data_') and lines_written > 1:
                components += 1
                if components > 2:  # Get first 2 components
                    break
    
    tmp_path = Path(tmp.name)

print(f"Extracted {lines_written} lines")

# Parse with gemmi
doc = gemmi.cif.read(str(tmp_path))
print(f"Loaded {len(doc)} blocks\n")

# Debug the first component
if len(doc) > 0:
    comp_id = '000'
    if comp_id in doc:
        block = doc[comp_id]
        print(f"Debugging component {comp_id}:")
        
        # Try to get descriptor loop
        desc_loop = block.find_loop('_pdbx_chem_comp_descriptor.comp_id')
        if desc_loop:
            print(f"\nDescriptor loop type: {type(desc_loop)}")
            print(f"Length: {len(desc_loop)}")
            
            # Check what methods/attributes are available
            print(f"\nLoop attributes: {[attr for attr in dir(desc_loop) if not attr.startswith('_')]}")
            
            # Try to access first row
            print(f"\nFirst row:")
            if len(desc_loop) > 0:
                first_row = desc_loop[0]
                print(f"  Type: {type(first_row)}")
                print(f"  Content: {first_row}")
                
                if hasattr(first_row, '__len__'):
                    print(f"  Length: {len(first_row)}")
                    print(f"  Items: {list(first_row)}")

# Also try a different approach - look for all loops
print("\n\nAll loops in block:")
if comp_id in doc:
    block = doc[comp_id]
    # Get all items in block
    for item in block:
        if hasattr(item, 'loop'):
            print(f"Found loop: {item.loop}")

# Clean up
tmp_path.unlink()