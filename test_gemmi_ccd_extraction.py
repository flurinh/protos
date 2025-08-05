#!/usr/bin/env python3
"""
Test proper gemmi extraction of CCD SMILES data.
"""

import gemmi
import gzip
import tempfile
from pathlib import Path

ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/cache/databases/ccd/components.cif.gz")

print("=== Testing Gemmi CCD Extraction ===\n")

# Extract a portion of the file for testing
print("Extracting portion of CCD file...")
with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
    with gzip.open(ccd_file, 'rt') as f_in:
        # Read until we get a few components including ATP
        lines_written = 0
        found_atp = False
        
        for line in f_in:
            tmp.write(line)
            lines_written += 1
            
            if line.startswith('data_ATP'):
                found_atp = True
            
            # Continue for a bit after ATP to complete its block
            if found_atp and line.startswith('data_') and not line.startswith('data_ATP'):
                break
                
            if lines_written > 500000:  # Safety limit
                break
    
    tmp_path = Path(tmp.name)

print(f"Extracted {lines_written} lines to {tmp_path}")
print(f"ATP found: {found_atp}")

# Now parse with gemmi
print("\nParsing with gemmi...")
doc = gemmi.cif.read(str(tmp_path))

print(f"Loaded {len(doc)} blocks")

# Test specific components
test_components = ['000', '001', 'ATP', 'NAD', 'HEM']

for comp_id in test_components:
    if comp_id in doc:
        print(f"\n{comp_id}:")
        block = doc[comp_id]
        
        # Get basic info
        name = block.find_value('_chem_comp.name')
        formula = block.find_value('_chem_comp.formula')
        comp_type = block.find_value('_chem_comp.type')
        
        print(f"  Name: {name}")
        print(f"  Formula: {formula}")
        print(f"  Type: {comp_type}")
        
        # Try to get SMILES from descriptor loop
        try:
            desc_loop = block.find_loop('_pdbx_chem_comp_descriptor.comp_id')
            if desc_loop:
                print(f"  Descriptor loop found with {len(desc_loop)} rows")
                
                # Try different approach - iterate and use find_value
                for i in range(len(desc_loop)):
                    try:
                        desc_type = desc_loop[i]['_pdbx_chem_comp_descriptor.type']
                        if 'SMILES' in desc_type:
                            descriptor = desc_loop[i]['_pdbx_chem_comp_descriptor.descriptor']
                            print(f"  {desc_type}: {descriptor[:50]}...")
                    except:
                        # Alternative approach
                        try:
                            row_data = desc_loop[i]
                            # Access by index if dict access fails
                            if hasattr(row_data, '__getitem__'):
                                if len(row_data) >= 5:  # Ensure we have enough columns
                                    desc_type = str(row_data[1])  # type is usually second column
                                    if 'SMILES' in desc_type:
                                        descriptor = str(row_data[4])  # descriptor is usually fifth column
                                        print(f"  {desc_type}: {descriptor[:50]}...")
                        except:
                            pass
            else:
                print("  No descriptor loop found")
                
        except Exception as e:
            print(f"  Error getting descriptors: {e}")
    else:
        print(f"\n{comp_id}: Not found in extracted portion")

# Clean up
tmp_path.unlink()
print("\nTest complete!")