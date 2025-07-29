#!/usr/bin/env python
"""Simple test of GRN assignment for GPCR structures."""

import os
from pathlib import Path
from protos.io.paths import ProtosPaths
from protos.processing.structure.structure_processor import StructureProcessor

# Setup
data_dir = Path(__file__).parent / "data"
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

paths = ProtosPaths()
struct_proc = StructureProcessor(name="test", paths=paths)

# Load structure
print("Loading structure...")
struct_proc.load_structure("3SN6")
struct_proc.update_pdb_ids()

print(f"Loaded {len(struct_proc.data)} atoms")
print(f"PDB IDs: {struct_proc.pdb_ids}")

# Check sequences
sequences = struct_proc.get_seq_dict()
print(f"\nFound {len(sequences)} sequences:")
for seq_id, seq in sequences.items():
    print(f"  {seq_id}: {len(seq)} residues")

# Assign GRNs
print("\nAssigning GRNs...")
grn_assignments = struct_proc.assign_grns(
    protein_family='gpcr_a',
    similarity_threshold=0.2,
    verbose=False
)

print(f"\nGRN assignments completed for {len(grn_assignments)} chains")

# Check if GRN column exists
if 'grn' in struct_proc.data.columns:
    grn_data = struct_proc.data[struct_proc.data['grn'].notna()]
    print(f"Total residues with GRN: {len(grn_data)}")
    
    # Show unique GRNs
    unique_grns = sorted(grn_data['grn'].unique())
    print(f"Unique GRN positions: {len(unique_grns)}")
    
    # Show key positions
    key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    print("\nKey GPCR positions:")
    for pos in key_positions:
        residues = grn_data[grn_data['grn'] == pos]
        if not residues.empty:
            res = residues.iloc[0]
            print(f"  {pos}: Chain {res['auth_chain_id']} residue {res['auth_seq_id']}")
else:
    print("ERROR: GRN column not created!")