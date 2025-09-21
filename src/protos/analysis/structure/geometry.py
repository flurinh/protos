"""
Geometric analysis functions for protein structures.

This module provides functions for calculating distances, rotations,
and other geometric properties of protein structures.
"""

import numpy as np
import pandas as pd
from typing import Union, Optional


def calculate_distance_matrix(coordinates: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
    """
    Calculate pairwise distance matrix for a set of coordinates.
    
    Args:
        coordinates: N x 3 array or DataFrame with x, y, z columns
        
    Returns:
        N x N symmetric distance matrix
    """
    if isinstance(coordinates, pd.DataFrame):
        coords = coordinates[['x', 'y', 'z']].values
    else:
        coords = coordinates
    
    # Efficient vectorized distance calculation
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=-1))
    
    return distances


def calculate_rotation_matrix(vec_from: np.ndarray, vec_to: np.ndarray) -> np.ndarray:
    """
    Calculate rotation matrix to rotate vec_from to vec_to using Rodrigues' formula.
    
    Args:
        vec_from: Source vector (3D)
        vec_to: Target vector (3D)
        
    Returns:
        3x3 rotation matrix
    """
    # Normalize vectors
    vec_from = vec_from / np.linalg.norm(vec_from)
    vec_to = vec_to / np.linalg.norm(vec_to)
    
    # Check if vectors are already aligned
    if np.allclose(vec_from, vec_to):
        return np.eye(3)
    
    # Check if vectors are opposite
    if np.allclose(vec_from, -vec_to):
        # Find an orthogonal vector
        orthogonal = np.array([1, 0, 0]) if abs(vec_from[0]) < 0.9 else np.array([0, 1, 0])
        v = np.cross(vec_from, orthogonal)
        v = v / np.linalg.norm(v)
        return 2 * np.outer(v, v) - np.eye(3)
    
    # Rodrigues' rotation formula
    v = np.cross(vec_from, vec_to)
    s = np.linalg.norm(v)
    c = np.dot(vec_from, vec_to)
    
    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])
    
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
    
    return R


def apply_rotation(df: pd.DataFrame, rotation_matrix: np.ndarray, 
                  center: Optional[np.ndarray] = None) -> pd.DataFrame:
    """
    Apply rotation matrix to coordinates in DataFrame.
    
    Args:
        df: DataFrame with x, y, z columns
        rotation_matrix: 3x3 rotation matrix
        center: Optional center of rotation. If None, rotates around origin
        
    Returns:
        DataFrame with rotated coordinates
    """
    df_copy = df.copy()
    coords = df_copy[['x', 'y', 'z']].values
    
    if center is not None:
        coords = coords - center
    
    rotated_coords = np.dot(coords, rotation_matrix.T)
    
    if center is not None:
        rotated_coords = rotated_coords + center
    
    df_copy[['x', 'y', 'z']] = rotated_coords
    return df_copy


def apply_translation(df: pd.DataFrame, translation: np.ndarray) -> pd.DataFrame:
    """
    Apply translation to coordinates in DataFrame.
    
    Args:
        df: DataFrame with x, y, z columns
        translation: 3D translation vector
        
    Returns:
        DataFrame with translated coordinates
    """
    df_copy = df.copy()
    df_copy[['x', 'y', 'z']] += translation
    return df_copy


def flip_structure(df: pd.DataFrame, axis: str = 'x') -> pd.DataFrame:
    """
    Flip structure by rotating 180° around specified axis.
    
    Args:
        df: DataFrame with x, y, z columns
        axis: Axis to rotate around ('x', 'y', or 'z')
        
    Returns:
        DataFrame with flipped coordinates
    """
    # Create 180° rotation matrix for specified axis
    if axis == 'x':
        R = np.array([[1, 0, 0], [0, -1, 0], [0, 0, -1]])
    elif axis == 'y':
        R = np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]])
    elif axis == 'z':
        R = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])
    else:
        raise ValueError(f"Invalid axis '{axis}'. Must be 'x', 'y', or 'z'")
    
    # Calculate center of mass
    center = df[['x', 'y', 'z']].mean().values
    
    return apply_rotation(df, R, center=center)


def apply_rotation_translation_to_coords(coords: np.ndarray, 
                                       rotation: np.ndarray, 
                                       translation: np.ndarray) -> np.ndarray:
    """
    Apply rotation and translation to coordinates.
    
    Args:
        coords: N x 3 array of coordinates
        rotation: 3x3 rotation matrix
        translation: 3D translation vector
        
    Returns:
        Transformed coordinates
    """
    return np.dot(coords, rotation.T) + translation


def calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> Union[float, np.ndarray]:
    """
    Calculate Euclidean distance between coordinates.
    
    Args:
        coord1: First coordinate(s) - can be single point (3,) or array of points (N, 3)
        coord2: Second coordinate(s) - can be single point (3,) or array of points (N, 3)
        
    Returns:
        Distance(s) between coordinates
    """
    return np.linalg.norm(coord1 - coord2, axis=-1)


def find_contacts(df: pd.DataFrame, cutoff: float = 5.0, 
                  exclude_same_residue: bool = True,
                  exclude_neighbors: int = 0) -> pd.DataFrame:
    """
    Find all atom pairs within cutoff distance.
    
    Args:
        df: Structure DataFrame with coordinates
        cutoff: Distance cutoff in Angstroms
        exclude_same_residue: If True, exclude contacts within same residue
        exclude_neighbors: Exclude contacts between residues within this sequence distance
        
    Returns:
        DataFrame with columns: atom1_idx, atom2_idx, chain1, chain2, 
        resid1, resid2, resname1, resname2, atom1, atom2, distance
    """
    # Reset index to have access to all columns
    df_reset = df.reset_index()
    
    # Calculate distance matrix
    coords = df_reset[['x', 'y', 'z']].values
    dist_matrix = calculate_distance_matrix(coords)
    
    # Find pairs within cutoff
    i_indices, j_indices = np.where((dist_matrix < cutoff) & (dist_matrix > 0))
    
    # Only keep upper triangle to avoid duplicates
    mask = i_indices < j_indices
    i_indices = i_indices[mask]
    j_indices = j_indices[mask]
    
    contacts = []
    for i, j in zip(i_indices, j_indices):
        row_i = df_reset.iloc[i]
        row_j = df_reset.iloc[j]
        
        # Check exclusion criteria
        if exclude_same_residue:
            if (row_i['auth_chain_id'] == row_j['auth_chain_id'] and 
                row_i['auth_seq_id'] == row_j['auth_seq_id']):
                continue
        
        if exclude_neighbors > 0:
            if (row_i['auth_chain_id'] == row_j['auth_chain_id'] and 
                abs(row_i['auth_seq_id'] - row_j['auth_seq_id']) <= exclude_neighbors):
                continue
        
        contacts.append({
            'atom1_idx': i,
            'atom2_idx': j,
            'structure_id': row_i['structure_id'],
            'chain1': row_i['auth_chain_id'],
            'chain2': row_j['auth_chain_id'],
            'resid1': row_i['auth_seq_id'],
            'resid2': row_j['auth_seq_id'],
            'resname1': row_i.get('res_name', row_i.get('res_name3l', 'UNK')),
            'resname2': row_j.get('res_name', row_j.get('res_name3l', 'UNK')),
            'atom1': row_i['atom_name'],
            'atom2': row_j['atom_name'],
            'distance': dist_matrix[i, j]
        })
    
    return pd.DataFrame(contacts)


def calculate_center_of_mass(df: pd.DataFrame, mass_column: Optional[str] = None) -> np.ndarray:
    """
    Calculate center of mass for structure.
    
    Args:
        df: Structure DataFrame with coordinates
        mass_column: Column name containing atomic masses. 
                    If None, uses uniform weights
        
    Returns:
        3D center of mass coordinates
    """
    coords = df[['x', 'y', 'z']].values
    
    if mass_column and mass_column in df.columns:
        masses = df[mass_column].values
    else:
        # Default atomic masses for common atoms
        element_masses = {
            'C': 12.011, 'N': 14.007, 'O': 15.999, 'S': 32.065,
            'P': 30.974, 'H': 1.008, 'CA': 40.078, 'FE': 55.845,
            'ZN': 65.38, 'MG': 24.305, 'CL': 35.453, 'NA': 22.990
        }
        
        if 'element' in df.columns:
            masses = df['element'].map(element_masses).fillna(12.0).values
        else:
            # Extract element from atom name (first 1-2 characters)
            df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
            elements = df_reset['atom_name'].str.extract(r'^([A-Z][a-z]?)')[0]
            masses = elements.map(element_masses).fillna(12.0).values
    
    # Calculate center of mass
    total_mass = masses.sum()
    com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
    
    return com


def calculate_radius_of_gyration(df: pd.DataFrame, mass_column: Optional[str] = None) -> float:
    """
    Calculate radius of gyration for structure.
    
    The radius of gyration is a measure of the compactness of a structure,
    defined as the root-mean-square distance of atoms from the center of mass.
    
    Args:
        df: Structure DataFrame with coordinates
        mass_column: Column name containing atomic masses.
                    If None, uses uniform weights
        
    Returns:
        Radius of gyration in Angstroms
    """
    # Calculate center of mass
    com = calculate_center_of_mass(df, mass_column)
    
    # Get coordinates
    coords = df[['x', 'y', 'z']].values
    
    # Get masses (same logic as center of mass)
    if mass_column and mass_column in df.columns:
        masses = df[mass_column].values
    else:
        # Use uniform masses for simplicity
        masses = np.ones(len(df))
    
    # Calculate squared distances from center of mass
    distances_squared = np.sum((coords - com) ** 2, axis=1)
    
    # Calculate radius of gyration
    total_mass = masses.sum()
    rg_squared = np.sum(masses * distances_squared) / total_mass
    
    return np.sqrt(rg_squared)