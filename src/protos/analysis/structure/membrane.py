"""
Membrane protein analysis functions.

This module provides specialized functions for analyzing and orienting
membrane proteins, including helix annotation and membrane normal calculation.
"""

import numpy as np
import pandas as pd
from typing import Dict, Optional, Tuple
from sklearn.decomposition import PCA

from .geometry import calculate_rotation_matrix, apply_rotation


def calculate_membrane_normal(df_ca: pd.DataFrame) -> np.ndarray:
    """
    Calculate membrane normal vector using PCA on CA atoms.
    
    Args:
        df_ca: DataFrame with CA atoms (must have x, y, z columns)
        
    Returns:
        3D unit vector representing the membrane normal
    """
    # Extract CA coordinates
    ca_coords = df_ca[['x', 'y', 'z']].values
    
    # Apply PCA
    pca = PCA(n_components=3)
    pca.fit(ca_coords)
    
    # The third principal component (smallest variance) is the membrane normal
    membrane_normal = pca.components_[2]
    
    # Normalize to unit vector
    membrane_normal = membrane_normal / np.linalg.norm(membrane_normal)
    
    return membrane_normal


def orient_in_membrane(df: pd.DataFrame, max_iterations: int = 12) -> pd.DataFrame:
    """
    Orient structure with membrane normal along z-axis.
    
    Uses iterative refinement to find optimal membrane orientation.
    
    Args:
        df: DataFrame with protein structure
        max_iterations: Maximum refinement iterations
        
    Returns:
        DataFrame with oriented structure
    """
    df_oriented = df.copy()
    
    for iteration in range(max_iterations):
        # Extract CA atoms
        df_ca = df_oriented[df_oriented['res_atom_name'] == 'CA'].copy()
        
        if df_ca.empty:
            print("Warning: No CA atoms found for membrane orientation")
            return df_oriented
        
        # Calculate current membrane normal
        membrane_normal = calculate_membrane_normal(df_ca)
        
        # Target is z-axis
        target_normal = np.array([0, 0, 1])
        
        # Check if already aligned
        alignment = np.dot(membrane_normal, target_normal)
        if abs(alignment) > 0.999:  # Close enough
            break
        
        # Calculate rotation to align with z-axis
        if alignment < 0:  # Pointing opposite direction
            membrane_normal = -membrane_normal
        
        R = calculate_rotation_matrix(membrane_normal, target_normal)
        
        # Apply rotation
        df_oriented = apply_rotation(df_oriented, R)
    
    return df_oriented


def orient_n_terminus_up(structures: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Orient structures with N-terminus pointing up (positive z).
    
    Args:
        structures: Dictionary of {structure_id: DataFrame}
        
    Returns:
        Dictionary with oriented structures
    """
    oriented_structures = {}
    
    for struct_id, df in structures.items():
        df_copy = df.copy()
        
        # Find N-terminus (minimum residue number)
        n_term_residue = df_copy['res_seq_id'].min()
        n_term_df = df_copy[df_copy['res_seq_id'] == n_term_residue]
        
        # Find C-terminus (maximum residue number)
        c_term_residue = df_copy['res_seq_id'].max()
        c_term_df = df_copy[df_copy['res_seq_id'] == c_term_residue]
        
        if not n_term_df.empty and not c_term_df.empty:
            # Get average z-coordinate for N and C terminus
            n_term_z = n_term_df['z'].mean()
            c_term_z = c_term_df['z'].mean()
            
            # If N-terminus is below C-terminus, flip the structure
            if n_term_z < c_term_z:
                from .geometry import flip_structure
                df_copy = flip_structure(df_copy, axis='x')
        
        oriented_structures[struct_id] = df_copy
    
    return oriented_structures


def annotate_transmembrane_helices(df: pd.DataFrame, k: int = 7, w: int = 5) -> pd.DataFrame:
    """
    Annotate transmembrane helix numbers using z-coordinate clustering.
    
    Args:
        df: DataFrame with protein structure
        k: Number of helices to identify
        w: Window size for smoothing
        
    Returns:
        DataFrame with 'helix_annotation' column added
    """
    df_copy = df.copy()
    
    # Extract CA atoms
    df_ca = df_copy[df_copy['res_atom_name'] == 'CA'].copy()
    
    if df_ca.empty:
        print("Warning: No CA atoms found for helix annotation")
        df_copy['helix_annotation'] = -1
        return df_copy
    
    # Sort by residue sequence
    df_ca = df_ca.sort_values('res_seq_id')
    
    # Get z-coordinates
    z_coords = df_ca['z'].values
    residue_ids = df_ca['res_seq_id'].values
    
    # Simple peak detection based on z-coordinate
    # Find regions with high absolute z-values (membrane spanning)
    z_abs_smooth = pd.Series(np.abs(z_coords)).rolling(window=w, center=True).mean()
    
    # Find peaks (local maxima)
    peaks = []
    for i in range(1, len(z_abs_smooth) - 1):
        if pd.notna(z_abs_smooth.iloc[i]) and \
           z_abs_smooth.iloc[i] > z_abs_smooth.iloc[i-1] and \
           z_abs_smooth.iloc[i] > z_abs_smooth.iloc[i+1]:
            peaks.append(i)
    
    # Select top k peaks
    if len(peaks) > k:
        peak_values = [z_abs_smooth.iloc[p] for p in peaks]
        sorted_peaks = sorted(zip(peak_values, peaks), reverse=True)
        peaks = [p[1] for p in sorted_peaks[:k]]
        peaks.sort()  # Sort by position
    
    # Initialize all residues with no helix assignment
    df_copy['helix_annotation'] = -1
    
    # Assign helix numbers based on proximity to peaks
    if peaks:
        # Create helix regions around peaks
        for helix_num, peak_idx in enumerate(peaks, 1):
            peak_residue = residue_ids[peak_idx]
            
            # Find extent of helix (simplified approach)
            # Assign residues within a window around the peak
            window = 10  # Approximate helix length
            
            mask = (df_copy['res_seq_id'] >= peak_residue - window) & \
                   (df_copy['res_seq_id'] <= peak_residue + window)
            
            df_copy.loc[mask, 'helix_annotation'] = helix_num
    
    return df_copy


def annotate_processed_structures(processed_structures: Dict[str, Dict]) -> Dict[str, Dict]:
    """
    Apply helix annotation to multiple processed structures.
    
    Args:
        processed_structures: Dictionary of {structure_id: {'df_norm': DataFrame}}
        
    Returns:
        Updated dictionary with helix annotations added
    """
    for struct_id, struct_data in processed_structures.items():
        if 'df_norm' in struct_data:
            df = struct_data['df_norm']
            df_annotated = annotate_transmembrane_helices(df)
            struct_data['df_norm'] = df_annotated
    
    return processed_structures