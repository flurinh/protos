"""
Example script demonstrating the use of the standardized structure module.

This script shows how to load and process structure data using the new
path management system and standardized schemas.
"""

import os
import pandas as pd
import numpy as np
from pathlib import Path

# Import path management
from protos.io.paths import ProtosPaths, get_structure_path

# Import structure utilities
from protos.processing.schema.structure import load_structure

# Set up path configuration
paths = ProtosPaths()
print(f"Data root directory: {paths.data_root}")
print(f"Structure directory: {paths.get_structure_subdir_path('structure_dir')}")

# Example 1: Load structure by PDB ID
pdb_id = "4zw9"  # Example PDB ID (substitute with a valid ID)
structure_path = get_structure_path(pdb_id)
print(f"\nStructure path for {pdb_id}: {structure_path}")

# Check if the structure file exists
if os.path.exists(structure_path):
    # Load structure
    print(f"\nLoading structure {pdb_id}...")
    structure_df = load_structure(pdb_id)
    
    # Print basic information
    print(f"Structure loaded: {len(structure_df)} atoms")
    print(f"Columns: {structure_df.columns.tolist()}")
    
    # Get chain information
    chains = structure_df['auth_chain_id'].unique()
    print(f"Chains: {chains}")
    
    # Get residue count
    residue_count = structure_df[['auth_chain_id', 'auth_seq_id']].drop_duplicates().shape[0]
    print(f"Residue count: {residue_count}")
    
    # Get alpha carbons
    ca_atoms = structure_df[structure_df['atom_name'] == 'CA']
    print(f"Alpha carbon count: {len(ca_atoms)}")
else:
    print(f"Structure file not found: {structure_path}")
    print("You may need to download the structure first or specify a valid PDB ID.")

print("\nExample complete!")