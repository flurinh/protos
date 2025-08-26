import numpy as np
import pandas as pd
import gc
from tqdm import tqdm
from protos.processing.structure.struct_alignment import get_structure_alignment
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment


def find_bidirectional_closest_neighbors(df1, df2):
    # Initialize lists to store the results
    closest_neighbors_df1_to_df2 = []
    closest_neighbors_df2_to_df1 = []

    # Define a function to calculate distances and find the closest neighbor
    def calculate_closest_neighbors(source_df, target_df, source_key, target_key):
        closest_neighbors = []
        for index1, row1 in source_df.iterrows():
            min_distance = float('inf')
            closest_grn = None

            for index2, row2 in target_df.iterrows():
                distance = np.linalg.norm(np.array([row1['x'], row1['y'], row1['z']]) -
                                          np.array([row2['x'], row2['y'], row2['z']]))
                if distance < min_distance:
                    min_distance = distance
                    closest_grn = row2['grn']

            closest_neighbors.append({source_key: row1['grn'], target_key: closest_grn, 'distance': min_distance})

        return pd.DataFrame(closest_neighbors)

    # Calculate closest neighbors from df1 to df2
    closest_neighbors_df1_to_df2 = calculate_closest_neighbors(df1, df2, 'grn_df1', 'grn_df2')

    # Calculate closest neighbors from df2 to df1
    closest_neighbors_df2_to_df1 = calculate_closest_neighbors(df2, df1, 'grn_df2', 'grn_df1')

    return closest_neighbors_df1_to_df2, closest_neighbors_df2_to_df1


def find_closest_sidechain_atoms_to_retinal(structure_df, retinal_xyz):
    """
    Identifies the closest atom of each sidechain to retinal for a given protein fold,
    including the minimum distance, while preserving all original columns from structure_df.
    - structure_df: DataFrame containing the full protein structure with detailed columns.
    - retinal_df: DataFrame containing retinal atoms with ['x', 'y', 'z'].
    Returns a DataFrame with all columns from structure_df for the closest atoms, plus 'min_dist'.
    """
    closest_atoms_data = []

    # Calculate pairwise distances between atoms in the structure and retinal atoms
    for grn in structure_df['grn'].unique():
        grn_atoms = structure_df[structure_df['grn'] == grn]

        if grn_atoms.empty:
            continue

        distances = cdist(grn_atoms[['x', 'y', 'z']], retinal_xyz)
        min_distance_index = np.argmin(distances, axis=None)
        atom_index, retinal_atom_index = np.unravel_index(min_distance_index, distances.shape)
        min_distance = distances[atom_index, retinal_atom_index]

        # Select the row of the closest atom and convert to a dictionary
        closest_atom_data = grn_atoms.iloc[atom_index].to_dict()

        # Add the minimum distance to the dictionary
        closest_atom_data['min_dist'] = min_distance

        closest_atoms_data.append(closest_atom_data)

    # Convert the list of dictionaries to a DataFrame
    return pd.DataFrame(closest_atoms_data)


def mapping_pipeline(df1, df2, schiff_base_residue_1, schiff_base_residue_2, penalty=7):
    # Step 1: Directly map the Schiff base residues, assuming distance 0 for a perfect match
    mapping_rows = [{'mo': schiff_base_residue_1, 'ao': schiff_base_residue_2, 'distance': 0}]

    # Exclude the Schiff base residues from further mapping
    df1_excl = df1[df1['grn'] != schiff_base_residue_1].reset_index(drop=True)
    df2_excl = df2[df2['grn'] != schiff_base_residue_2].reset_index(drop=True)

    # Step 2: Prepare the cost matrix for the Hungarian algorithm
    distances = np.linalg.norm(
        df1_excl[['x', 'y', 'z']].values[:, None, :] - df2_excl[['x', 'y', 'z']].values[None, :, :], axis=2)
    max_shape = max(distances.shape)
    cost_matrix = np.full((max_shape, max_shape), fill_value=penalty)  # Initialize with penalties
    cost_matrix[:distances.shape[0], :distances.shape[1]] = distances  # Fill with actual distances

    # Apply the Hungarian algorithm
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Process the results, include all assignments with their distances or penalties
    for i, j in zip(row_ind, col_ind):
        mo_grn = df1_excl.iloc[i]['grn'] if i < len(df1_excl) else None
        ao_grn = df2_excl.iloc[j]['grn'] if j < len(df2_excl) else None
        distance = cost_matrix[i, j]

        mapping_rows.append({'mo': mo_grn, 'ao': ao_grn, 'distance': distance})

    # Create the mapping DataFrame from the list of dictionaries
    mapping = pd.DataFrame(mapping_rows)

    return mapping


def process_structures(cp1, cp2, primary_interaction_dist=9, limit=None):
    mappings = {}
    cp1.data[['x', 'y', 'z']] = cp1.data[['x', 'y', 'z']].astype(np.float16)
    cp2.data[['x', 'y', 'z']] = cp2.data[['x', 'y', 'z']].astype(np.float16)

    # Get unique PDB IDs for both cp1 and cp2
    pdb_ids_cp1 = cp1.data['pdb_id'].unique()
    pdb_ids_cp2 = cp2.data['pdb_id'].unique()

    # Apply limit if specified
    if limit:
        pdb_ids_cp1 = pdb_ids_cp1[:limit]
        pdb_ids_cp2 = pdb_ids_cp2[:limit]

    for pdb_id1 in tqdm(pdb_ids_cp1):
        # Get all data for this PDB in cp1
        pdb_data_cp1 = cp1.data[cp1.data['pdb_id'] == pdb_id1]
        chains_cp1 = pdb_data_cp1['auth_chain_id'].unique()

        for chain_id1 in chains_cp1:
            # Extract structure and retinal from cp1
            struct_cp1 = pdb_data_cp1[(pdb_data_cp1['auth_chain_id'] == chain_id1) &
                                      (pdb_data_cp1['res_name3l'] != 'RET') &
                                      (pdb_data_cp1['group'] == 'ATOM')]
            ret_cp1 = pdb_data_cp1[(pdb_data_cp1['auth_chain_id'] == chain_id1) &
                                   (pdb_data_cp1['res_name3l'] == 'RET')]

            if ret_cp1.empty:
                print(f"No retinal found for {pdb_id1}_{chain_id1}. Skipping.")
                continue

            for pdb_id2 in pdb_ids_cp2:
                # Get all data for this PDB in cp2
                pdb_data_cp2 = cp2.data[cp2.data['pdb_id'] == pdb_id2]
                chains_cp2 = pdb_data_cp2['auth_chain_id'].unique()

                for chain_id2 in chains_cp2:
                    print(f"Comparing {pdb_id1}_{chain_id1} with {pdb_id2}_{chain_id2}")

                    # Extract structure and retinal from cp2
                    struct_cp2 = pdb_data_cp2[(pdb_data_cp2['auth_chain_id'] == chain_id2) &
                                              (pdb_data_cp2['res_name3l'] != 'RET') &
                                              (pdb_data_cp2['group'] == 'ATOM')]
                    ret_cp2 = pdb_data_cp2[(pdb_data_cp2['auth_chain_id'] == chain_id2) &
                                           (pdb_data_cp2['res_name3l'] == 'RET')]

                    if ret_cp2.empty:
                        print(f"No retinal found for {pdb_id2}_{chain_id2}. Skipping.")
                        continue

                    if len(ret_cp1) == 20 and len(ret_cp2) == 20:
                        try:
                            # Structure alignment and transformation
                            rot, tran, best_path, best_rmsd = get_structure_alignment(
                                ret_cp1[['x', 'y', 'z']].astype(np.float16),
                                ret_cp2[['x', 'y', 'z']].astype(np.float16))
                            struct_cp2[['x', 'y', 'z']] = np.dot(struct_cp2[['x', 'y', 'z']].astype(np.float16).values,
                                                                 rot) + tran
                            ret_cp2[['x', 'y', 'z']] = np.dot(ret_cp2[['x', 'y', 'z']].astype(np.float16).values,
                                                              rot) + tran

                            # Finding closest atoms to retinal
                            closest_atoms_cp1 = find_closest_sidechain_atoms_to_retinal(struct_cp1,
                                                                                        ret_cp1[['x', 'y', 'z']].astype(
                                                                                            np.float16))
                            closest_atoms_cp2 = find_closest_sidechain_atoms_to_retinal(struct_cp2,
                                                                                        ret_cp2[['x', 'y', 'z']].astype(
                                                                                            np.float16))
                            closest_atoms_cp1 = closest_atoms_cp1[
                                closest_atoms_cp1['min_dist'] < primary_interaction_dist]
                            closest_atoms_cp2 = closest_atoms_cp2[
                                closest_atoms_cp2['min_dist'] < primary_interaction_dist]

                            schiff_base_cp1 = closest_atoms_cp1.sort_values('min_dist').iloc[0]['grn']
                            schiff_base_cp2 = closest_atoms_cp2.sort_values('min_dist').iloc[0]['grn']

                            mapping = mapping_pipeline(closest_atoms_cp1, closest_atoms_cp2, schiff_base_cp1,
                                                       schiff_base_cp2)
                            mappings[f"{pdb_id1}_{chain_id1}_vs_{pdb_id2}_{chain_id2}"] = mapping

                        except Exception as e:
                            print(f"Error processing {pdb_id1}_{chain_id1} vs {pdb_id2}_{chain_id2}: {e}")
                            mappings[f"{pdb_id1}_{chain_id1}_vs_{pdb_id2}_{chain_id2}"] = None

                        # Cleanup to free memory
                        del struct_cp2, ret_cp2, closest_atoms_cp1, closest_atoms_cp2
                        gc.collect()

                    else:
                        print(
                            f"Skipping due to mismatch in expected lengths for {pdb_id1}_{chain_id1} vs {pdb_id2}_{chain_id2}.")
                        print(len(ret_cp1), len(ret_cp2))
                        mappings[f"{pdb_id1}_{chain_id1}_vs_{pdb_id2}_{chain_id2}"] = None

            # Cleanup after processing all cp2 structures for this cp1 structure
            del struct_cp1, ret_cp1
            gc.collect()

    return mappings


def create_universal_mapping_with_stats(mappings):
    # Filter out None values in mappings
    filtered_mappings = {key: value for key, value in mappings.items() if value is not None}

    # Aggregate all non-None mappings into a combined DataFrame
    combined = pd.concat(filtered_mappings.values(), ignore_index=True)

    # Calculate frequencies and average distances for each mo-ao pair
    stats = combined.groupby(['mo', 'ao'], dropna=False).agg(
        frequency=('mo', 'size'),
        avg_distance=('distance', 'mean')
    ).reset_index().sort_values(by=['frequency', 'avg_distance'], ascending=[False, True])

    # Initialize placeholders for the universal mapping
    universal_mapping_list = []
    used_mo = set()
    used_ao = set()

    # Iterate over the sorted stats to prioritize common and close pairings
    for _, row in stats.iterrows():
        mo, ao, frequency, avg_distance = row['mo'], row['ao'], row['frequency'], row['avg_distance']

        # Skip if either mo or ao has already been used, enforcing uniqueness
        if mo in used_mo or ao in used_ao:
            continue

        # Add the pairing to the universal mapping list
        universal_mapping_list.append({'mo': mo, 'ao': ao, 'frequency': frequency, 'avg_distance': avg_distance})
        used_mo.add(mo)
        used_ao.add(ao)

    # Handle unmatched residues by adding them with None for their pair, frequency of 0, and a default distance
    unmatched_mo = set(combined['mo']) - used_mo
    for mo in unmatched_mo:
        universal_mapping_list.append({'mo': mo, 'ao': None, 'frequency': 0, 'avg_distance': None})
    unmatched_ao = set(combined['ao']) - used_ao
    for ao in unmatched_ao:
        universal_mapping_list.append({'mo': None, 'ao': ao, 'frequency': 0, 'avg_distance': None})

    # Convert the list to a DataFrame
    universal_mapping = pd.DataFrame(universal_mapping_list)

    return universal_mapping
