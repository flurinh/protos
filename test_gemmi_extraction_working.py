#!/usr/bin/env python3
"""
Test proper gemmi extraction for CCD SMILES.
"""

import gemmi
import gzip
import tempfile
from pathlib import Path

ccd_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/cache/databases/ccd/components.cif.gz")

print("=== Testing Proper Gemmi CCD Extraction ===\n")

# Extract first few components
with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as tmp:
    with gzip.open(ccd_file, 'rt') as f_in:
        lines_written = 0
        components = 0
        
        for line in f_in:
            tmp.write(line)
            lines_written += 1
            
            if line.startswith('data_') and lines_written > 1:
                components += 1
                if components > 10:  # Get first 10 components
                    break
    
    tmp_path = Path(tmp.name)

print(f"Extracted {lines_written} lines ({components} components)")

# Parse with gemmi
doc = gemmi.cif.read(str(tmp_path))
print(f"Loaded {len(doc)} blocks\n")

# Process components and show the correct way to extract SMILES
for comp_id in ['000', '001', '002', 'ATP', 'NAD']:
    if comp_id in doc:
        block = doc[comp_id]
        print(f"\n{comp_id}:")
        
        # Get basic info
        name = block.find_value('_chem_comp.name')
        print(f"  Name: {name}")
        
        # The CORRECT way to get descriptor data from gemmi
        found_smiles = False
        
        # Find the descriptor loop by iterating through block items
        for item in block:
            if hasattr(item, 'loop') and item.loop:
                loop = item.loop
                # Check if this is the descriptor loop by looking at tags
                tags = [loop.tags[i] for i in range(len(loop.tags))]
                
                if '_pdbx_chem_comp_descriptor.type' in tags:
                    type_idx = tags.index('_pdbx_chem_comp_descriptor.type')
                    desc_idx = tags.index('_pdbx_chem_comp_descriptor.descriptor')
                    
                    # Iterate through rows
                    for row_idx in range(loop.length()):
                        desc_type = loop[row_idx, type_idx]  # Note: tuple indexing!
                        descriptor = loop[row_idx, desc_idx]
                        
                        # Strip quotes if present
                        if descriptor and descriptor.startswith('"') and descriptor.endswith('"'):
                            descriptor = descriptor[1:-1]
                        
                        if 'SMILES' in desc_type:
                            print(f"  {desc_type}: {descriptor[:50]}...")
                            found_smiles = True
                    break
        
        if not found_smiles:
            print("  No SMILES found")

# Clean up
tmp_path.unlink()
print("\n\nThis is how gemmi extraction should work!")