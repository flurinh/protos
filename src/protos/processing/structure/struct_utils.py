import os
import sys
import numpy as np
import pandas as pd
import gemmi
from gemmi import cif
from math import degrees
from scipy.spatial.transform import Rotation as R
from sklearn.decomposition import PCA
from scipy.spatial.distance import cdist

from protos.processing.structure.struct_alignment import get_structure_alignment


CIF_COLUMNS = ['group_PDB', 'auth_asym_id', 'label_asym_id', 'label_seq_id', 'auth_seq_id',
                     'label_comp_id', 'id', 'label_atom_id', 'type_symbol', 'Cartn_x', 'Cartn_y', 'Cartn_z']


STRUCT_COLUMNS = ['pdb_id','group', 'auth_chain_id', 'gen_chain_id', 'gen_seq_id', 'auth_seq_id', 'res_name3l',
'atom_id', 'res_atom_name', 'atom_name', 'x', 'y', 'z', 'phi', 'omega', 'psi', 'res_name1l']

STRUCT_COLUMN_DTYPE = {
    'pdb_id': str,
    'group': str,
    'auth_chain_id': str,
    'gen_chain_id': str,
    'gen_seq_id': int,
    'auth_seq_id': int,
    'res_name3l': str,
    'atom_id': int,
    'res_atom_name': str,
    'atom_name': str,
    'x': float,
    'y': float,
    'z': float,
    'phi': float,
    'omega': float,
    'psi': float,
    'res_name1l': str,
    'grn': str
}

ALPHA_CARBON = 'CA'

SORTED_STRUCT_COLUMNS = ['pdb_id', 'group', 'auth_chain_id', 'gen_chain_id', 'auth_seq_id', 'gen_seq_id', 'res_name3l',
                         'res_name1l', 'phi', 'omega', 'psi', 'atom_id', 'res_atom_name', 'atom_name', 'x', 'y', 'z']


def get_distance_matrix(coordinates):
    """
    Calculate the pairwise distance matrix between 3D coordinates.
    
    Args:
        coordinates (np.ndarray or pd.DataFrame): Array of 3D coordinates with shape (n, 3)
                                                 If DataFrame, must have columns ['x', 'y', 'z']
    
    Returns:
        np.ndarray: Distance matrix with shape (n, n)
    """
    if isinstance(coordinates, pd.DataFrame):
        coords = coordinates[['x', 'y', 'z']].values
    else:
        coords = coordinates
    
    return cdist(coords, coords)


def load_structure(pdb_id, folder='data/mmcif/', structure_id=0):
    """
    Loads a PDB file in mmCIF format and processes its data.

    Args:
        pdb_id (str): The PDB ID of the file to load.

    Returns:
        pd.DataFrame: A DataFrame containing the processed data from the PDB file.
    """
    # Handle scientific notation IDs like 1e12
    original_pdb_id = str(pdb_id)
    search_pdb_id = original_pdb_id.lower()
    
    # Handle special cases for scientific notation
    if search_pdb_id in ['1.00e+12', '1.0e+12', '1e+12', '1e12', '1.00e12']:
        pdb_id = '1e12'
        print(f"Normalizing scientific notation PDB ID: {original_pdb_id} -> {pdb_id}")
    
    # First try with the original PDB ID - this is the most likely to work
    filename = os.path.join(folder, original_pdb_id + '.cif')
    if not os.path.exists(filename):
        # If that doesn't exist, try other possible variations for backward compatibility
        alt_filename = os.path.join(folder, pdb_id + '.cif')
        if os.path.exists(alt_filename):
            filename = alt_filename
            print(f"Using alternative filename: {filename}")
        # For uppercase variant
        elif '1.00E+12' in original_pdb_id.upper():
            alt_filename = os.path.join(folder, '1e12.cif')
            if os.path.exists(alt_filename):
                filename = alt_filename
                print(f"Using 1e12.cif instead of {original_pdb_id}.cif")

    cols = ['PDB'] + CIF_COLUMNS
    # Atom-wise features
    try:
        # copy all the data from mmCIF file
        block = cif.read_file(filename).sole_block()
        atom_rows = []
        table = block.find('_atom_site.', CIF_COLUMNS)
        for row in table:
            atom_rows.append([pdb_id] + list(row))
    except Exception as e:
        print("Encountered error while parsing file. %s" % e)
        sys.exit(1)

    atom_df = pd.DataFrame(data=atom_rows, columns=cols)
    atom_df['label_seq_id'] = atom_df.apply(lambda x: int(x.label_seq_id) if x.label_seq_id != '.' else np.nan, axis=1)

    # Residue-wise features
    st = gemmi.read_structure(filename, merge_chain_parts=True)
    print(f"{pdb_id}: Found {len(st)} model(s), using only model {structure_id}.")
    model = st[structure_id]
    res_rows = []
    for chain in model:
        for r, res in enumerate(chain.get_polymer()):
            try:
                prev_res = chain.previous_residue(res)
                next_res = chain.next_residue(res)
                phi, psi = gemmi.calculate_phi_psi(prev_res, res, next_res)
            except:
                phi, psi = np.nan, np.nan
            try:
                next_res = chain.next_residue(res)
                omega = gemmi.calculate_omega(res, next_res)
            except:
                omega = np.nan
            res_rows.append([res.label_seq, res.subchain, degrees(phi), degrees(omega), degrees(psi)])
    cols_extended = ['label_seq_id', 'label_asym_id', 'phi', 'omega', 'psi']
    res_df = pd.DataFrame(data=res_rows, columns=cols_extended)
    res_df['label_seq_id'] = res_df['label_seq_id'].astype(int)

    structure_df = pd.merge(atom_df, res_df, how='outer', on=['label_asym_id', 'label_seq_id'])

    structure_df['label_comp_sid'] = structure_df.apply(lambda x:
                                                        gemmi.find_tabulated_residue(
                                                            x.label_comp_id).one_letter_code,
                                                        axis=1)
    structure_df.columns = STRUCT_COLUMNS
    return structure_df[SORTED_STRUCT_COLUMNS]


def get_protein_structure_df(model, chain_id, only_CA = True, use_atom_type = False):
    x, y, z, atom_type = [], [], [], []
    for chain in model:
        if chain.name == chain_id:
            for residue in chain.get_polymer():
                for atom in residue:
                    if only_CA and atom.name == 'CA':  # Only consider Calpha atoms for the protein
                        x.append(atom.pos.x)
                        y.append(atom.pos.y)
                        z.append(atom.pos.z)
                    elif not only_CA:
                        x.append(atom.pos.x)
                        y.append(atom.pos.y)
                        z.append(atom.pos.z)
                    else:
                        pass
                    if use_atom_type:
                        print(atom)
                        atom_type.append(atom.name)
    if atom_type:
        return pd.DataFrame(np.array([x, y, z, atom_type]).T, columns=['x', 'y', 'z', 'atom_type'])
    else:
        return pd.DataFrame(np.array([x, y, z]).T, columns=['x', 'y', 'z'])


def calc_rotation_matrix(A, B):
    """
    Calculate the rotation matrix that rotates vector A to vector B.
    :param A:
    :param B:
    :return:
    """
    # Normalize vectors A and B
    A_normalized = A / np.linalg.norm(A)
    B_normalized = B / np.linalg.norm(B)

    # Find the rotation axis (cross product of A_normalized and B_normalized)
    rotation_axis = np.cross(A_normalized, B_normalized)

    # Find the rotation angle (use arccos of the dot product of A_normalized and B_normalized)
    cosine_angle = np.clip(np.dot(A_normalized, B_normalized), -1.0, 1.0)
    angle = np.arccos(cosine_angle)

    # Handle the case where A and B are in the same or opposite direction
    if np.isclose(angle, 0) or np.isclose(angle, np.pi):
        # Special handling might be needed, e.g., pick a perpendicular vector as the rotation axis
        # For now, we'll handle the zero angle case (no rotation needed)
        if np.isclose(angle, 0):
            return np.eye(3)
        # For the opposite direction, a 180-degree rotation around any perpendicular axis works
        if np.isclose(angle, np.pi):
            # Find a vector perpendicular to A (or B)
            perp_vector = np.array([A_normalized[1], -A_normalized[0], 0])
            perp_vector /= np.linalg.norm(perp_vector)
            rotation_axis = perp_vector
            angle = np.pi

    # Create the rotation matrix using Rodrigues' rotation formula
    rotation_vector = rotation_axis * angle
    rotation_matrix = R.from_rotvec(rotation_vector).as_matrix()

    return rotation_matrix


def calculate_membrane_normal(df_ca):
    """
    This function assumes that the membrane lies in the XY plane. The input is a structure df of a single structure
    (do not pass multiple models of a structure at once.)
    :param transmembrane_df: structure pandas dataframe of a membrane protein
    :return: the normal vector of the transmembrane protein bundle
    """
    # Convert the DataFrame to a NumPy array
    transmembrane_coords = df_ca[['x', 'y', 'z']].to_numpy()

    # Perform PCA
    pca = PCA(n_components=3)
    pca.fit(transmembrane_coords)

    # The normal vector to the membrane plane is the first principal component
    normal_vector = pca.components_[0]

    return normal_vector


def apply_rotation_matrix(df, rotation_matrix):
    """
    Helper function to apply a rotation matrix to structure dataframe
    :param df: structure dataframe
    :param rotation_matrix: 3x3 rotation matrix
    :return: rotated structure dataframe
    """
    coords = df[['x', 'y', 'z']].to_numpy()
    rotated_structure = np.dot(rotation_matrix,
                               coords.T).T  # Note the transposition to ensure correct matrix multiplication
    transformed_df = pd.DataFrame(rotated_structure, columns=['x', 'y', 'z'])
    return transformed_df


def flip_protein(df):
    """
    Rotates the protein structure around the x-axis by 180 degrees to correct its orientation,
    preserving the geometry.

    :param df: DataFrame containing the protein structure coordinates.
    :return: DataFrame with the rotated protein structure.
    """
    # Rotation matrix for 180 degrees around the x-axis
    rotation_matrix = np.array([[1, 0, 0],
                                [0, -1, 0],
                                [0, 0, -1]])

    coords = df[['x', 'y', 'z']].to_numpy()
    rotated_coords = np.dot(coords, rotation_matrix.T)  # Apply rotation
    df[['x', 'y', 'z']] = rotated_coords
    return df


def orient_structure_in_membrane(df, df_ret, normal_vector, max_steps=12, eps=0.01):
    """
    Iteratively reorients the structure such that its membrane normal aligns with the z-axis.
    :param df: DataFrame of a membrane protein structure with 'x', 'y', 'z' columns.
    :param normal_vector: The initial membrane normal vector.
    :param max_steps: Maximum number of iterations for refining the orientation.
    :param eps: Convergence criterion based on the angle difference to the z-axis.
    :return: The DataFrame with reoriented coordinates.
    """

    # Z-axis unit vector for target orientation
    z_axis = np.array([0, 0, 1.], dtype=np.float32)

    for _ in range(max_steps):
        coords = df[['x', 'y', 'z']].astype(np.float32).to_numpy()
        coords_ret = df_ret[['x', 'y', 'z']].astype(np.float32).to_numpy()
        # Normalize the normal vector to ensure it's a unit vector
        normal_vector = normal_vector / np.linalg.norm(normal_vector)

        # Calculate the rotation matrix for current normal to z-axis
        rotation_matrix = calc_rotation_matrix(normal_vector, z_axis)

        # Apply the rotation to the structure
        rotated_structure = np.dot(rotation_matrix, coords.T).T
        rotate_ret = np.dot(rotation_matrix, coords_ret.T).T
        df[['x', 'y', 'z']] = rotated_structure
        df_ret[['x', 'y', 'z']] = rotate_ret

        # Recalculate the normal vector to check for convergence
        new_normal = calculate_membrane_normal(df)
        new_normal /= np.linalg.norm(new_normal)  # Normalize the new normal vector

        # Check if the rotation has sufficiently aligned the normal with the z-axis
        angle_diff = np.arccos(np.clip(np.dot(new_normal, z_axis), -1.0, 1.0))
        if angle_diff < eps:
            break  # Stop if the new normal is sufficiently aligned with the z-axis

        normal_vector = new_normal  # Update the normal vector for the next iteration

    return df


def orient_structures_n_terminus_up(processed_structures):
    """
    Ensures all structures and their ligands have the correct orientation, with the N-terminus
    pointing towards the extracellular domain (positive z-axis) by rotating around the x-axis.

    :param processed_structures: Dictionary of protein structures with their data.
    :return: Updated dictionary with correctly oriented structures and ligands.
    """
    for pdb_id, structure in processed_structures.items():
        df_norm = structure['df_norm']
        df_ret = structure['ret_norm']

        # Determine orientation based on the N-terminus z-coordinate
        n_terminus_z = df_norm[~df_norm['auth_seq_id'].isna()].sort_values('auth_seq_id')['z'].iloc[0]

        # If N-terminus is on the negative z-axis, rotate the protein and the ligand
        if n_terminus_z < 0:
            df_norm = flip_protein(df_norm)
            df_ret = flip_protein(df_ret)  # Assuming df_ret has x, y, z columns

        # Update the processed structures with the potentially rotated data
        processed_structures[pdb_id].update({
            'df_norm': df_norm,
            'ret_norm': df_ret
        })

    return processed_structures


def annotate_helix_numbers(df, k, w):
    """This function assumes that the input structure is already oriented in the membrane
    and that the z-coordinate is the normal to the membrane. The input should be a
    membrane protein with transmembrane helices.
    :param df: structure dataframe
    :param k:
    :param w:
    """
    # Initialize variables
    current_helix = 0
    last_direction = None
    direction_window = []  # to store recent direction history

    # Add a new column for helix numbers, initializing with None
    df['helix'] = None

    # Iterate through the DataFrame
    for i in range(len(df) - k):
        current_z = float(df.at[i, 'z'])
        next_z = float(df.at[i + k, 'z'])

        # Determine the direction
        if abs(next_z) > abs(current_z):
            direction = 'away'
        else:
            direction = 'towards'

        # Update direction window
        direction_window.append(direction)
        if len(direction_window) > w:
            direction_window.pop(0)

        # Check for helix transition
        if len(set(direction_window)) == 1 and direction_window[-1] != last_direction and direction_window[
            -1] == 'towards':
            if current_helix <= 8:
                current_helix += 1
        last_direction = direction_window[0]

        # Annotate the residue with the current helix number
        df.at[i, 'helix'] = current_helix

    # Assign the helix number to the last few residues
    for i in range(len(df) - k, len(df)):
        df.at[i, 'helix'] = current_helix

    return df


def get_ca_ret_coords(processed_structures, pdb_id):
    pdb_ids = list(processed_structures.keys())
    chain_id = processed_structures[pdb_id]['chain_id']
    df = processed_structures[pdb_id]['df']
    df_ca = df[(df['res_atom_name'] == 'CA') &
               (df['auth_chain_id'] == chain_id)][['x', 'y', 'z']].reset_index(drop=True).astype(np.double)
    df_ret = pd.DataFrame(processed_structures[pdb_id]['ret_coords'], columns=['x', 'y', 'z'])
    return df_ca, df_ret


def get_chain(processed_structures, pdb_id):
    chain_id = processed_structures[pdb_id]['chain_id']
    df = processed_structures[pdb_id]['df']
    df_all = df[df['auth_chain_id'] == chain_id].reset_index(drop=True)
    return df_all


# Load and preprocess structures (placeholder for actual data loading)
def normalize_structures(processed_structures):
    for pdb_id, data in processed_structures.items():
        df_ca, df_ret = get_ca_ret_coords(processed_structures, pdb_id)

        # Get membrane normal (we only use the backbone)
        mean_backbone = df_ca.mean(axis=0)
        df_ca = df_ca[['x', 'y', 'z']].astype(np.float16) - mean_backbone
        df_ret = df_ret[['x', 'y', 'z']].astype(np.float16) - mean_backbone

        membrane_normal = calculate_membrane_normal(df_ca)

        # Normalize full structure
        df_norm = get_chain(processed_structures, pdb_id)
        df_norm[['x', 'y', 'z']] = df_norm[['x', 'y', 'z']].astype(np.float16) - mean_backbone
        orient_structure_in_membrane(df_norm, df_ret, membrane_normal)

        # lets check
        df_norm_ca = df_norm[df_norm['res_atom_name'] == 'CA']
        calculate_membrane_normal(df_norm_ca)

        # Update processed structures with oriented data
        processed_structures[pdb_id].update({
            'df_norm': df_norm,
            'ret_norm': df_ret,
            'membrane_normal': membrane_normal
        })


def annotate_processed_structures(processed_structures, k=7, w=5):
    """

    :param processed_structures:
    :param k:
    :param w:
    :return:
    """
    pdb_ids = list(processed_structures.keys())
    for pdb_id in pdb_ids:
        df = processed_structures[pdb_id]['df_norm']
        df_ca = df[df['res_atom_name'] == 'CA']
        gen_seq_id = df_ca['gen_seq_id'].tolist()
        xyz_df_ca = df_ca[['x', 'y', 'z']].astype(np.float16).reset_index()
        df_ca_annotated = annotate_helix_numbers(xyz_df_ca, k=k, w=w)
        df_ca_annotated['gen_seq_id'] = gen_seq_id
        # now merge on gen_seq_id with our df
        merged_df = df.merge(df_ca_annotated[['gen_seq_id', 'helix']], on='gen_seq_id', how='left')

        # Replace processed_structures entry with the updated DataFrame
        processed_structures[pdb_id]['df_norm'] = merged_df

    return processed_structures


def apply_rotation_translation_to_ret(ret_df, rotation_matrix, translation_vector):
    # Convert DataFrame to a NumPy array
    coords = ret_df[['x', 'y', 'z']].to_numpy()

    # Apply rotation
    rotated_coords = np.dot(coords, rotation_matrix)

    # Apply translation
    translated_coords = rotated_coords + translation_vector

    # Update the DataFrame with the transformed coordinates
    transformed_df = pd.DataFrame(translated_coords, columns=['x', 'y', 'z'])
    return transformed_df


def align_proteins_on_retinal(processed_structures, reference_id='4PXK', ret_id='RET'):
    # Ensure the reference structure is in the processed_structures
    if reference_id not in processed_structures:
        raise ValueError(f"Reference structure {reference_id} not found in processed structures.")

    # Extract retinal coordinates of the reference structure
    reference_coords = pd.DataFrame(processed_structures[reference_id]['ret_norm'], columns=['x', 'y', 'z'])

    aligned_results = {}

    # Align each structure based on the reference retinal coordinates
    for pdb_id, structure in processed_structures.items():
        current_coords = pd.DataFrame(structure['ret_norm'], columns=['x', 'y', 'z'])

        if pdb_id == reference_id:
            aligned_results[pdb_id] = {'rotation_matrix': np.eye(3), 'translation_vector': np.zeros(3), 'rmsd': 0}
            continue

        rot, tran, _, best_rmsd = get_structure_alignment(reference_coords, current_coords)
        aligned_results[pdb_id] = {
            'rotation_matrix': rot,
            'translation_vector': tran,
            'rmsd': best_rmsd
        }

    # Apply transformations and update processed_structures
    for pdb_id, structure in processed_structures.items():
        rot = aligned_results[pdb_id]['rotation_matrix']
        tran = aligned_results[pdb_id]['translation_vector']

        # Apply to retinal coordinates
        ret_df = pd.DataFrame(structure['ret_norm'], columns=['x', 'y', 'z'])
        transformed_ret = apply_rotation_translation_to_ret(ret_df, rot, tran)
        processed_structures[pdb_id]['ret_aligned'] = transformed_ret

        # Apply to full structure (assuming df_norm represents the normalized structure DataFrame)
        transformed_structure = apply_rotation_translation_to_ret(structure['df_norm'], rot, tran)
        processed_structures[pdb_id]['df_ret_aligned'] = transformed_structure

    return processed_structures