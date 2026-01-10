#!/usr/bin/env python3
"""
Test updated CCD parsing.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.ingest import ccd_loader

# CCD directory
ccd_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/databases/ccd")

print("=== Testing Updated CCD Parser ===")

# Test loading some common ligands
ligands = ['ATP', 'NAD', 'HEM', 'FAD', 'ADP', 'GTP']

for ligand_id in ligands:
    print(f"\nTesting {ligand_id}:")
    
    # Load component
    comp_data = ccd_loader.load_ccd_component(ccd_dir, ligand_id)
    if comp_data:
        print(f"  ✓ Component loaded")
        print(f"  Name: {comp_data.get('name', 'N/A')}")
        print(f"  Formula: {comp_data.get('formula', 'N/A')}")
        
        # Check for SMILES
        smiles = ccd_loader.get_ccd_smiles(ccd_dir, ligand_id)
        if smiles:
            print(f"  ✓ SMILES found: {smiles[:50]}...")
        else:
            print(f"  ✗ No SMILES found")
            # Debug: show what fields we have
            smiles_fields = [k for k in comp_data.keys() if 'smiles' in k.lower()]
            print(f"  Available SMILES fields: {smiles_fields}")
    else:
        print(f"  ✗ Failed to load component")

# Test the safe loader
print("\n=== Testing Safe CCD Loader ===")
atp = ccd_loader.get_ccd_ligand_safe(ccd_dir, "ATP")
if atp:
    print("✓ ATP loaded with safe loader:")
    for key, value in atp.items():
        if key == 'smiles':
            print(f"  {key}: {value[:50]}...")
        else:
            print(f"  {key}: {value}")
else:
    print("✗ Failed to load ATP with safe loader")