#!/usr/bin/env python3
"""
Check CCD pickle file.
"""

import pickle
from pathlib import Path

ccd_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data/ligand/cache/databases/ccd")
data_path = ccd_dir / "ccd_ligands.pkl"

print(f"Loading pickle file: {data_path}")
print(f"File size: {data_path.stat().st_size / 1024 / 1024:.1f} MB")

with open(data_path, 'rb') as f:
    ligand_data = pickle.load(f)

print(f"\nTotal components: {len(ligand_data)}")

# Check specific components
test_comps = ['ATP', 'NAD', 'HEM', '000', '001', '002']

for comp_id in test_comps:
    if comp_id in ligand_data:
        comp = ligand_data[comp_id]
        print(f"\n{comp_id}:")
        print(f"  Name: {comp.get('name', 'N/A')}")
        print(f"  SMILES: {comp.get('smiles', 'N/A')}")
        print(f"  SMILES canonical: {comp.get('smiles_canonical', 'N/A')}")
        
# Count components with SMILES
with_smiles = sum(1 for comp in ligand_data.values() if comp.get('smiles') or comp.get('smiles_canonical'))
print(f"\nComponents with SMILES: {with_smiles}")

# Show first component with SMILES
for comp_id, comp in ligand_data.items():
    if comp.get('smiles'):
        print(f"\nFirst component with SMILES: {comp_id}")
        print(f"  SMILES: {comp['smiles'][:50]}...")
        break