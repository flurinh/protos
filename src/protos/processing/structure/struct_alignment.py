import numpy as np
import pandas as pd
from Bio.PDB.ccealign import run_cealign
from Bio.PDB.qcprot import QCPSuperimposer


def align_structures(coords1, coords2, window_size=8, max_gap=30):
    """
    Aligns two structures using CEalign and returns the rotation, translation, and RMSD.
    Returns the rotated coordinates of the second structure.
    :param coords1: This is the reference structure (Pandas DataFrame or numpy array).
    :param coords2: This is the structure to be aligned (Pandas DataFrame or numpy array).
    :param window_size: The window size for CEalign.
    :param max_gap: The maximum gap size for CEalign.
    :return: coords2_rot_df, rot, tran, best_path, best_rmsd
    """
    # Convert to list format that CEalign expects
    if isinstance(coords1, pd.DataFrame):
        coords1 = coords1.values.tolist()
    elif isinstance(coords1, np.ndarray):
        coords1 = coords1.tolist()
        
    if isinstance(coords2, pd.DataFrame):
        coords2 = coords2.values.tolist()
    elif isinstance(coords2, np.ndarray):
        coords2 = coords2.tolist()
    # Get CEAlignment object from run_cealign
    alignment = run_cealign(coords1, coords2, window_size, max_gap)[0]
    
    # Extract path information from the alignment object
    if hasattr(alignment, 'path'):
        paths = [(alignment.path[0], alignment.path[1])]
    else:
        raise RuntimeError("Failed to find alignment path in CEAlignment object.")
    
    # Create unique paths set
    unique_paths = {(tuple(pA), tuple(pB)) for pA, pB in paths}

    best_rmsd, best_u = 1e6, None
    for u_path in unique_paths:
        idxA, idxB = u_path

        coordsA = np.array([coords1[i] for i in idxA])
        coordsB = np.array([coords2[i] for i in idxB])

        aln = QCPSuperimposer()
        aln.set(coordsA, coordsB)
        aln.run()
        if aln.rms < best_rmsd:
            best_rmsd = aln.rms
            best_u = (aln.rot, aln.tran)
            best_path = u_path

    if best_u is None:
        raise RuntimeError("Failed to find a suitable alignment.")

    rot, tran = best_u

    # Apply the rotation and translation to all atoms in the second structure
    coords2_rot = np.dot(coords2, rot) + tran
    coords2_rot_df = pd.DataFrame(coords2_rot, columns=['x', 'y', 'z'])

    return coords2_rot_df, rot, tran, best_path, best_rmsd


def get_structure_alignment(coords1, coords2, window_size=8, max_gap=30):
    """
    Aligns two structures using CEalign and returns the rotation, translation, and RMSD.
    :param coords1: Pandas DataFrame or numpy array with coordinates
    :param coords2: Pandas DataFrame or numpy array with coordinates
    :param window_size: Window size for CEalign
    :param max_gap: Maximum gap size for CEalign
    :return: Rotation matrix, translation vector, best alignment path, RMSD
    """
    # Convert to list format that CEalign expects
    if isinstance(coords1, pd.DataFrame):
        coords1 = coords1.values.tolist()
    elif isinstance(coords1, np.ndarray):
        coords1 = coords1.tolist()
        
    if isinstance(coords2, pd.DataFrame):
        coords2 = coords2.values.tolist()
    elif isinstance(coords2, np.ndarray):
        coords2 = coords2.tolist()

    # Get CEAlignment object from run_cealign
    alignment = run_cealign(coords1, coords2, window_size, max_gap)[0]
    
    # Extract path information from the alignment object
    if hasattr(alignment, 'path'):
        # Make sure path has the correct format (two lists)
        if len(alignment.path) == 2:
            path_a, path_b = alignment.path
            paths = [(path_a, path_b)]
        else:
            raise RuntimeError("Unexpected path format in CEAlignment object.")
    else:
        raise RuntimeError("Failed to find alignment path in CEAlignment object.")
    
    # Create unique paths set
    unique_paths = {(tuple(pA), tuple(pB)) for pA, pB in paths}

    best_rmsd, best_u = 1e6, None
    for u_path in unique_paths:
        idxA, idxB = u_path

        coordsA = np.array([coords1[i] for i in idxA])
        coordsB = np.array([coords2[i] for i in idxB])

        aln = QCPSuperimposer()
        aln.set(coordsA, coordsB)
        aln.run()
        if aln.rms < best_rmsd:
            best_rmsd = aln.rms
            best_u = (aln.rot, aln.tran)
            best_path = u_path

    if best_u is None:
        raise RuntimeError("Failed to find a suitable alignment.")

    rot, tran = best_u
    return rot, tran, best_path, best_rmsd


def structure_comparison_ava(processed_structures):
    """
    Perform all-vs-all structure comparison using 'get_structure_alignment'.
    :param processed_structures: Dictionary of structure IDs mapped to their respective 'df_norm' DataFrames.
    :return: Numpy array containing RMSD values for all structure pairs, and a list of structure IDs for indexing.
    """
    structures_ids = list(processed_structures.keys())
    n = len(structures_ids)
    rmsd_matrix = np.zeros((n, n))

    for i in range(n):
        ref_df = processed_structures[structures_ids[i]]['df_norm']
        ref_df_ca = ref_df[ref_df['res_atom_name']=='CA']
        reference_structure = ref_df_ca[['x', 'y', 'z']].dropna()
        for j in range(i + 1, n):  # To avoid redundant calculations
            df = processed_structures[structures_ids[j]]['df_norm']
            df_ca = df[df['res_atom_name'] == 'CA']
            current_structure = df_ca[['x', 'y', 'z']].dropna()
            _, _, _, best_rmsd = get_structure_alignment(reference_structure, current_structure)
            rmsd_matrix[i, j] = best_rmsd
            rmsd_matrix[j, i] = best_rmsd  # Symmetric matrix

    return rmsd_matrix, structures_ids


def structure_comparison_1va(processed_structures):
    """
    Perform one-vs-all structure comparison using 'get_structure_alignment'.
    The first entry in the processed_structures dictionary is considered as the reference.
    :param processed_structures: Dictionary of structure IDs mapped to their respective 'df_norm' DataFrames.
    :return: List containing RMSD values between the reference structure and all other structures, and their IDs.
    """
    structures_ids = list(processed_structures.keys())
    reference_structure = processed_structures[structures_ids[0]]['df_norm'][['x', 'y', 'z']].dropna()
    rmsd_list = []

    for i in range(1, len(structures_ids)):
        current_structure = processed_structures[structures_ids[i]]['df_norm'][['x', 'y', 'z']].dropna()
        _, _, _, best_rmsd = get_structure_alignment(reference_structure, current_structure)
        rmsd_list.append(best_rmsd)

    return rmsd_list, structures_ids[1:]