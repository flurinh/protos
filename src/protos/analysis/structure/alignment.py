"""
Structure alignment and superposition functions.

This module provides functions for aligning protein structures,
including specialized alignment based on ligands like retinal.
"""

import numpy as np
import pandas as pd
from typing import Dict, Tuple, List, Optional, Union
try:
    from Bio.PDB.ccealign import run_cealign
    from Bio.PDB.qcprot import QCPSuperimposer
    BIOPYTHON_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    BIOPYTHON_AVAILABLE = False


def _to_coordinate_array(coords: Union[pd.DataFrame, np.ndarray, List[List[float]]]) -> List[List[float]]:
    """Convert various coordinate containers to a plain list of XYZ triplets."""

    if isinstance(coords, pd.DataFrame):
        if not {'x', 'y', 'z'}.issubset(coords.columns):
            raise ValueError("DataFrame must contain 'x', 'y', and 'z' columns for alignment")
        return coords[['x', 'y', 'z']].to_numpy().tolist()

    if isinstance(coords, np.ndarray):
        if coords.shape[-1] != 3:
            raise ValueError("NumPy coordinate array must have shape (N, 3)")
        return coords.tolist()

    return coords  # assume already list-like


def align_structures(coords1: Union[pd.DataFrame, np.ndarray, List[List[float]]],
                    coords2: Union[pd.DataFrame, np.ndarray, List[List[float]]],
                    window_size: int = 8,
                    max_gap: int = 30) -> Tuple[pd.DataFrame, np.ndarray, np.ndarray, Tuple, float]:
    """
    Align two structures using CEalign algorithm.
    
    Args:
        coords1: Reference structure coordinates
        coords2: Structure to be aligned
        window_size: Window size for CEalign
        max_gap: Maximum gap size for CEalign
        
    Returns:
        Tuple of (rotated coords2, rotation matrix, translation, alignment path, RMSD)
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "CEalign-based alignment requires Biopython. Install 'biopython' to use this function."
        )

    # Convert to list format for CEalign
    coords1 = _to_coordinate_array(coords1)
    coords2 = _to_coordinate_array(coords2)
    
    # Get CEAlignment object
    alignment = run_cealign(coords1, coords2, window_size, max_gap)[0]
    
    # Extract alignment path
    if hasattr(alignment, 'path'):
        paths = [(alignment.path[0], alignment.path[1])]
    else:
        raise RuntimeError("Failed to find alignment path in CEAlignment object.")
    
    # Create unique paths set
    unique_paths = {(tuple(pA), tuple(pB)) for pA, pB in paths}
    
    # Find best alignment
    best_rmsd, best_u = 1e6, None
    best_path = None
    
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
    
    # Apply transformation to all atoms
    coords2_array = np.array(coords2)
    coords2_rot = np.dot(coords2_array, rot) + tran
    coords2_rot_df = pd.DataFrame(coords2_rot, columns=['x', 'y', 'z'])
    
    return coords2_rot_df, rot, tran, best_path, best_rmsd


def get_structure_alignment(coords1: Union[pd.DataFrame, np.ndarray, List[List[float]]],
                           coords2: Union[pd.DataFrame, np.ndarray, List[List[float]]],
                           window_size: int = 8, 
                           max_gap: int = 30) -> Tuple[np.ndarray, np.ndarray, Tuple, float]:
    """
    Get alignment transformation without applying it.
    
    Args:
        coords1: Reference structure coordinates
        coords2: Structure to be aligned
        window_size: Window size for CEalign
        max_gap: Maximum gap size for CEalign
        
    Returns:
        Tuple of (rotation matrix, translation vector, alignment path, RMSD)
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError(
            "CEalign-based alignment requires Biopython. Install 'biopython' to use this function."
        )

    # Convert to list format for CEalign
    coords1 = _to_coordinate_array(coords1)
    coords2 = _to_coordinate_array(coords2)
    
    # Get CEAlignment object
    alignment = run_cealign(coords1, coords2, window_size, max_gap)[0]
    
    # Extract alignment path
    if hasattr(alignment, 'path'):
        if len(alignment.path) == 2:
            path_a, path_b = alignment.path
            paths = [(path_a, path_b)]
        else:
            raise RuntimeError("Unexpected path format in CEAlignment object.")
    else:
        raise RuntimeError("Failed to find alignment path in CEAlignment object.")
    
    # Create unique paths set
    unique_paths = {(tuple(pA), tuple(pB)) for pA, pB in paths}
    
    # Find best alignment
    best_rmsd, best_u = 1e6, None
    best_path = None
    
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


def extract_ca_and_ligand_coords(df: pd.DataFrame, 
                                chain_id: str,
                                ligand_id: str = 'RET') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract CA atoms and ligand coordinates from structure.
    
    Args:
        df: Structure DataFrame
        chain_id: Chain identifier
        ligand_id: Ligand residue name (default: 'RET' for retinal)
        
    Returns:
        Tuple of (CA coordinates DataFrame, ligand coordinates DataFrame)
    """
    # Extract CA atoms from specified chain
    df_ca = df[(df['chain_id'] == chain_id) & (df['res_atom_name'] == 'CA')].copy()
    
    # Extract ligand coordinates
    df_ligand = df[df['res_name'] == ligand_id].copy()
    
    # Return coordinate columns only
    if not df_ca.empty:
        df_ca = df_ca[['x', 'y', 'z']]
    
    if not df_ligand.empty:
        df_ligand = df_ligand[['x', 'y', 'z']]
    
    return df_ca, df_ligand


def align_on_retinal(structures: Dict[str, pd.DataFrame], 
                    reference_id: str = '4PXK',
                    ligand_id: str = 'RET',
                    chain_id: str = 'A') -> Dict[str, pd.DataFrame]:
    """
    Align structures based on retinal (or other ligand) position.
    
    Args:
        structures: Dictionary of {structure_id: DataFrame}
        reference_id: ID of reference structure
        ligand_id: Ligand residue name for alignment
        chain_id: Chain to use for alignment
        
    Returns:
        Dictionary with aligned structures
    """
    if reference_id not in structures:
        raise ValueError(f"Reference structure '{reference_id}' not found")
    
    # Get reference structure
    ref_df = structures[reference_id]
    ref_ca, ref_ligand = extract_ca_and_ligand_coords(ref_df, chain_id, ligand_id)
    
    if ref_ca.empty:
        raise ValueError(f"No CA atoms found in reference structure chain {chain_id}")
    
    aligned_structures = {reference_id: ref_df}  # Reference stays unchanged
    
    # Align all other structures
    for struct_id, df in structures.items():
        if struct_id == reference_id:
            continue
        
        # Extract CA and ligand coordinates
        struct_ca, struct_ligand = extract_ca_and_ligand_coords(df, chain_id, ligand_id)
        
        if struct_ca.empty:
            print(f"Warning: No CA atoms found for {struct_id}, skipping")
            continue
        
        try:
            # Get alignment transformation
            rot, tran, _, rmsd = get_structure_alignment(ref_ca, struct_ca)
            
            # Apply transformation to entire structure
            coords = df[['x', 'y', 'z']].values
            aligned_coords = np.dot(coords, rot) + tran
            
            # Create aligned structure
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = aligned_coords
            
            aligned_structures[struct_id] = df_aligned
            
            # If ligand present, apply same transformation
            if not struct_ligand.empty:
                ligand_coords = struct_ligand.values
                aligned_ligand_coords = np.dot(ligand_coords, rot) + tran
                # Update ligand coordinates in the aligned structure
                ligand_mask = df_aligned['res_name'] == ligand_id
                df_aligned.loc[ligand_mask, ['x', 'y', 'z']] = aligned_ligand_coords
                
        except Exception as e:
            print(f"Warning: Failed to align {struct_id}: {e}")
            aligned_structures[struct_id] = df  # Keep original if alignment fails
    
    return aligned_structures


def calculate_rmsd(coords1: Union[np.ndarray, pd.DataFrame], 
                   coords2: Union[np.ndarray, pd.DataFrame],
                   weights: Optional[np.ndarray] = None) -> float:
    """
    Calculate RMSD between two sets of coordinates.
    
    Args:
        coords1: First set of coordinates (N x 3)
        coords2: Second set of coordinates (N x 3)
        weights: Optional weights for each atom
        
    Returns:
        RMSD value in Angstroms
    """
    # Convert to numpy arrays if needed
    if isinstance(coords1, pd.DataFrame):
        coords1 = coords1[['x', 'y', 'z']].values
    if isinstance(coords2, pd.DataFrame):
        coords2 = coords2[['x', 'y', 'z']].values
    
    # Check shapes match
    if coords1.shape != coords2.shape:
        raise ValueError(f"Coordinate shapes don't match: {coords1.shape} vs {coords2.shape}")
    
    # Calculate squared differences
    diff = coords1 - coords2
    squared_diff = np.sum(diff ** 2, axis=1)
    
    # Apply weights if provided
    if weights is not None:
        if weights.shape[0] != squared_diff.shape[0]:
            raise ValueError("Weights shape doesn't match number of atoms")
        rmsd = np.sqrt(np.sum(weights * squared_diff) / np.sum(weights))
    else:
        rmsd = np.sqrt(np.mean(squared_diff))
    
    return rmsd


def kabsch_alignment(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Calculate optimal rotation matrix using Kabsch algorithm.
    
    Args:
        coords1: Reference coordinates (N x 3)
        coords2: Coordinates to align (N x 3)
        
    Returns:
        Tuple of (rotation matrix, translation vector, RMSD after alignment)
    """
    # Center both sets of coordinates
    center1 = np.mean(coords1, axis=0)
    center2 = np.mean(coords2, axis=0)
    
    coords1_centered = coords1 - center1
    coords2_centered = coords2 - center2
    
    # Calculate covariance matrix
    H = coords2_centered.T @ coords1_centered
    
    # Singular Value Decomposition
    U, S, Vt = np.linalg.svd(H)
    
    # Calculate rotation matrix
    R = Vt.T @ U.T
    
    # Ensure proper rotation (det(R) = 1)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Calculate translation
    t = center1 - R @ center2
    
    # Apply transformation and calculate RMSD
    coords2_aligned = (R @ coords2.T).T + t
    rmsd = calculate_rmsd(coords1, coords2_aligned)
    
    return R, t, rmsd


def simple_align_structures(df1: pd.DataFrame, df2: pd.DataFrame, 
                           atom_selection: str = 'CA',
                           chain_id: Optional[str] = None) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Simple structure alignment based on atom selection.
    
    Args:
        df1: Reference structure DataFrame
        df2: Structure to align DataFrame  
        atom_selection: Atom type to use for alignment ('CA', 'backbone', 'all')
        chain_id: Specific chain to use (if None, uses all chains)
        
    Returns:
        Tuple of (rotation matrix, translation vector, RMSD)
    """
    # Reset index to access all columns
    df1_reset = df1.reset_index() if isinstance(df1.index, pd.MultiIndex) else df1
    df2_reset = df2.reset_index() if isinstance(df2.index, pd.MultiIndex) else df2
    
    # Filter by chain if specified
    if chain_id:
        df1_reset = df1_reset[df1_reset['auth_chain_id'] == chain_id]
        df2_reset = df2_reset[df2_reset['auth_chain_id'] == chain_id]
    
    # Filter by atom selection
    if atom_selection == 'CA':
        mask1 = df1_reset['atom_name'] == 'CA'
        mask2 = df2_reset['atom_name'] == 'CA'
    elif atom_selection == 'backbone':
        backbone_atoms = ['N', 'CA', 'C', 'O']
        mask1 = df1_reset['atom_name'].isin(backbone_atoms)
        mask2 = df2_reset['atom_name'].isin(backbone_atoms)
    elif atom_selection == 'all':
        mask1 = df1_reset['group'] == 'ATOM'
        mask2 = df2_reset['group'] == 'ATOM'
    else:
        raise ValueError(f"Unknown atom selection: {atom_selection}")
    
    df1_filtered = df1_reset[mask1]
    df2_filtered = df2_reset[mask2]
    
    # Match residues between structures
    # Create residue identifiers
    df1_filtered['res_id'] = (df1_filtered['auth_chain_id'] + '_' + 
                              df1_filtered['auth_seq_id'].astype(str) + '_' +
                              df1_filtered['atom_name'])
    df2_filtered['res_id'] = (df2_filtered['auth_chain_id'] + '_' + 
                              df2_filtered['auth_seq_id'].astype(str) + '_' +
                              df2_filtered['atom_name'])
    
    # Find common residues
    common_res = set(df1_filtered['res_id']) & set(df2_filtered['res_id'])
    
    if not common_res:
        raise ValueError("No common residues found between structures")
    
    # Sort by res_id to ensure matching order
    df1_common = df1_filtered[df1_filtered['res_id'].isin(common_res)].sort_values('res_id')
    df2_common = df2_filtered[df2_filtered['res_id'].isin(common_res)].sort_values('res_id')
    
    # Extract coordinates
    coords1 = df1_common[['x', 'y', 'z']].values
    coords2 = df2_common[['x', 'y', 'z']].values
    
    # Perform Kabsch alignment
    rotation, translation, rmsd = kabsch_alignment(coords1, coords2)
    
    return rotation, translation, rmsd
