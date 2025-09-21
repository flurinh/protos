import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def fix_ligand_chain(data, max_distance=5.0):
    # Group the data by PDB ID
    grouped = data.groupby('pdb_id')

    fixed_data_list = []

    for pdb_id, pdb_data in grouped:
        # Separate protein atoms and ligands for this PDB
        protein_atoms = pdb_data[pdb_data['group'] == 'ATOM']
        ligands = pdb_data[pdb_data['group'] == 'HETATM']

        if ligands.empty:
            fixed_data_list.append(pdb_data)
            continue  # No ligands to fix for this PDB

        # Create KD-tree for protein atoms
        protein_coords = protein_atoms[['x', 'y', 'z']].values
        protein_tree = cKDTree(protein_coords)

        # Find closest protein atoms for each ligand
        ligand_coords = ligands[['x', 'y', 'z']].values
        distances, indices = protein_tree.query(ligand_coords, k=5)  # Get 5 nearest neighbors

        # Assign chains based on majority vote of nearby protein atoms
        new_chains = []
        for dist, idx in zip(distances, indices):
            nearby_chains = protein_atoms.iloc[idx[dist <= max_distance]]['auth_chain_id'].values
            if len(nearby_chains) > 0:
                new_chain = pd.Series(nearby_chains).mode().iloc[0]  # Most common chain
            else:
                new_chain = ''  # No nearby atoms within max_distance
            new_chains.append(new_chain)

        # Assign the new chains to the ligands
        ligands['auth_chain_id'] = new_chains

        # Remove ligands that couldn't be assigned to a chain
        ligands = ligands[ligands['auth_chain_id'] != '']

        # Combine the updated ligands with the original protein atoms for this PDB
        updated_pdb_data = pd.concat([protein_atoms, ligands]).sort_index()

        fixed_data_list.append(updated_pdb_data)

    # Combine all the updated PDB data
    return pd.concat(fixed_data_list)

# Usage example:
# cp.data = fix_ligand_chain(cp.data, max_distance=5.0)