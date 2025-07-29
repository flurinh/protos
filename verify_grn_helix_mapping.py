#!/usr/bin/env python
"""Verify GRN mapping by examining a specific helix in detail."""

import os
from pathlib import Path
from protos.io.paths import ProtosPaths
from protos.processing.structure.structure_processor import StructureProcessor

# Setup
data_dir = Path(__file__).parent / "data"
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

paths = ProtosPaths()
struct_proc = StructureProcessor(name="verify", paths=paths)

# Load and process
struct_proc.load_structure("3SN6")
struct_proc.update_pdb_ids()
grn_assignments = struct_proc.assign_grns(protein_family='gpcr_a', verbose=False)

# Focus on helix 3 (TM3)
print("Helix 3 (TM3) GRN mapping verification:")
print("="*50)

chain_r_mask = (struct_proc.data['pdb_id'] == '3SN6') & (struct_proc.data['auth_chain_id'] == 'R')
grn_data = struct_proc.data[chain_r_mask & struct_proc.data['grn'].notna()]

# Get all helix 3 positions
helix3_data = grn_data[grn_data['grn'].str.startswith('3.')]
helix3_ca = helix3_data[helix3_data['res_atom_name'] == 'CA'].sort_values('auth_seq_id')

print(f"\nFound {len(helix3_ca)} CA atoms with helix 3 GRNs")
print("\nHelix 3 residues (auth_seq_id order):")
print(f"{'GRN':<8} {'Residue':<8} {'auth_seq_id':<12} {'Coordinates (x,y,z)'}")
print("-"*60)

for _, row in helix3_ca.iterrows():
    print(f"{row['grn']:<8} {row['res_name3l']:<8} {row['auth_seq_id']:<12} "
          f"({row['x']:6.1f}, {row['y']:6.1f}, {row['z']:6.1f})")

# Check for gaps
auth_seq_ids = helix3_ca['auth_seq_id'].values
print(f"\nauth_seq_id range: {auth_seq_ids.min()} to {auth_seq_ids.max()}")
print(f"Expected residues: {auth_seq_ids.max() - auth_seq_ids.min() + 1}")
print(f"Found residues: {len(auth_seq_ids)}")

# Key residue check
print("\nKey residue 3.50 (DRY motif):")
dry_residue = helix3_ca[helix3_ca['grn'] == '3.50']
if not dry_residue.empty:
    res = dry_residue.iloc[0]
    print(f"  {res['res_name3l']}{res['auth_seq_id']} at GRN 3.50")
    print(f"  This should be an Arginine (R) in most GPCRs")