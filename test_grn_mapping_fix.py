#!/usr/bin/env python
"""Test the fixed GRN mapping to verify auth_seq_id mapping is correct."""

import os
from pathlib import Path
from protos.io.paths import ProtosPaths
from protos.processing.structure.structure_processor import StructureProcessor

# Setup
data_dir = Path(__file__).parent / "data"
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

paths = ProtosPaths()
struct_proc = StructureProcessor(name="test_mapping", paths=paths)

# Load structure
print("Loading structure...")
struct_proc.load_structure("3SN6")
struct_proc.update_pdb_ids()

# Get sequence to see auth_seq_id mapping
print("\nChecking sequence extraction...")
sequences = struct_proc.get_seq_dict()
print(f"Found sequences: {list(sequences.keys())}")

# Check chain R backbone
chain_r_mask = (struct_proc.data['pdb_id'] == '3SN6') & (struct_proc.data['auth_chain_id'] == 'R')
backbone_mask = chain_r_mask & (struct_proc.data['res_atom_name'] == 'CA') & (struct_proc.data['group'] == 'ATOM')
chain_r_backbone = struct_proc.data[backbone_mask].sort_values(by='gen_seq_id')

print(f"\nChain R backbone has {len(chain_r_backbone)} CA atoms")
print("\nFirst 10 residues (sequence position -> auth_seq_id):")
for i, (idx, row) in enumerate(chain_r_backbone.head(10).iterrows(), start=1):
    print(f"  Position {i}: {row['res_name1l']} -> auth_seq_id {row['auth_seq_id']}")

print("\nSequence starts with:", sequences['3SN6_R'][:10])

# Assign GRNs
print("\nAssigning GRNs...")
grn_assignments = struct_proc.assign_grns(
    protein_family='gpcr_a',
    similarity_threshold=0.2,
    verbose=False
)

# Check GRN assignments
grn_data = struct_proc.data[struct_proc.data['grn'].notna() & chain_r_mask]
print(f"\nTotal residues with GRN in chain R: {len(grn_data)}")

# Check some key positions
print("\nChecking key GPCR positions:")
key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
for pos in key_positions:
    residues = grn_data[grn_data['grn'] == pos]
    if not residues.empty:
        # Get the CA atom
        ca_residues = residues[residues['res_atom_name'] == 'CA']
        if not ca_residues.empty:
            res = ca_residues.iloc[0]
            print(f"  GRN {pos}: {res['res_name3l']} at auth_seq_id {res['auth_seq_id']}")

# Check continuity - GRN positions should map to reasonable auth_seq_ids
print("\nChecking GRN position continuity (sample from helix 3):")
helix3_grns = ['3.46', '3.47', '3.48', '3.49', '3.50', '3.51']
for grn in helix3_grns:
    residues = grn_data[(grn_data['grn'] == grn) & (grn_data['res_atom_name'] == 'CA')]
    if not residues.empty:
        res = residues.iloc[0]
        print(f"  GRN {grn}: auth_seq_id {res['auth_seq_id']}")
        
# Let's also check the actual sequence positions in the GRN assignment
if '3SN6_R' in grn_assignments:
    grn_series = grn_assignments['3SN6_R']
    print(f"\nGRN assignment has {len(grn_series)} positions")
    print("\nFirst 10 GRN assignments:")
    for i, (grn, rn) in enumerate(list(grn_series.items())[:10]):
        print(f"  {grn}: {rn}")