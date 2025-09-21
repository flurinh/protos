"""
Structure comparison and analysis functions.

This module provides functions for comparing multiple protein structures,
including RMSD calculations and structure normalization.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional

from .alignment import get_structure_alignment
from .membrane import orient_in_membrane, annotate_transmembrane_helices


def compare_all_vs_all(structures: Dict[str, pd.DataFrame]) -> Tuple[np.ndarray, List[str]]:
    """
    Perform all-vs-all structure comparison.
    
    Args:
        structures: Dictionary of {structure_id: DataFrame with CA coordinates}
        
    Returns:
        Tuple of (RMSD matrix, list of structure IDs)
    """
    structure_ids = list(structures.keys())
    n = len(structure_ids)
    rmsd_matrix = np.zeros((n, n))
    
    for i in range(n):
        # Get reference CA atoms
        ref_df = structures[structure_ids[i]]
        if 'res_atom_name' in ref_df.columns:
            ref_ca = ref_df[ref_df['res_atom_name'] == 'CA'][['x', 'y', 'z']].dropna()
        else:
            ref_ca = ref_df[['x', 'y', 'z']].dropna()
        
        for j in range(i + 1, n):  # Upper triangle only
            # Get target CA atoms
            target_df = structures[structure_ids[j]]
            if 'res_atom_name' in target_df.columns:
                target_ca = target_df[target_df['res_atom_name'] == 'CA'][['x', 'y', 'z']].dropna()
            else:
                target_ca = target_df[['x', 'y', 'z']].dropna()
            
            try:
                # Calculate RMSD
                _, _, _, rmsd = get_structure_alignment(ref_ca, target_ca)
                rmsd_matrix[i, j] = rmsd
                rmsd_matrix[j, i] = rmsd  # Symmetric matrix
            except Exception as e:
                print(f"Warning: Failed to align {structure_ids[i]} vs {structure_ids[j]}: {e}")
                rmsd_matrix[i, j] = np.nan
                rmsd_matrix[j, i] = np.nan
    
    return rmsd_matrix, structure_ids


def compare_one_vs_all(structures: Dict[str, pd.DataFrame], 
                      reference_id: Optional[str] = None) -> Tuple[List[float], List[str]]:
    """
    Compare one structure against all others.
    
    Args:
        structures: Dictionary of {structure_id: DataFrame}
        reference_id: ID of reference structure (if None, uses first)
        
    Returns:
        Tuple of (RMSD values, structure IDs excluding reference)
    """
    structure_ids = list(structures.keys())
    
    if reference_id is None:
        reference_id = structure_ids[0]
    elif reference_id not in structure_ids:
        raise ValueError(f"Reference structure '{reference_id}' not found")
    
    # Get reference CA atoms
    ref_df = structures[reference_id]
    if 'res_atom_name' in ref_df.columns:
        ref_ca = ref_df[ref_df['res_atom_name'] == 'CA'][['x', 'y', 'z']].dropna()
    else:
        ref_ca = ref_df[['x', 'y', 'z']].dropna()
    
    rmsd_values = []
    compared_ids = []
    
    for struct_id in structure_ids:
        if struct_id == reference_id:
            continue
        
        # Get target CA atoms
        target_df = structures[struct_id]
        if 'res_atom_name' in target_df.columns:
            target_ca = target_df[target_df['res_atom_name'] == 'CA'][['x', 'y', 'z']].dropna()
        else:
            target_ca = target_df[['x', 'y', 'z']].dropna()
        
        try:
            # Calculate RMSD
            _, _, _, rmsd = get_structure_alignment(ref_ca, target_ca)
            rmsd_values.append(rmsd)
            compared_ids.append(struct_id)
        except Exception as e:
            print(f"Warning: Failed to align {reference_id} vs {struct_id}: {e}")
            rmsd_values.append(np.nan)
            compared_ids.append(struct_id)
    
    return rmsd_values, compared_ids


def normalize_structures(structures: Dict[str, pd.DataFrame], 
                        orient_membrane: bool = True,
                        annotate_helices: bool = True,
                        center: bool = True) -> Dict[str, Dict]:
    """
    Normalize structure coordinates and optionally orient/annotate.
    
    Args:
        structures: Dictionary of {structure_id: DataFrame}
        orient_membrane: Whether to orient membrane proteins
        annotate_helices: Whether to annotate transmembrane helices
        center: Whether to center structures at origin
        
    Returns:
        Dictionary with normalized structures and metadata
    """
    normalized_structures = {}
    
    for struct_id, df in structures.items():
        df_norm = df.copy()
        metadata = {'structure_id': struct_id}
        
        # Center structure at origin
        if center:
            center_of_mass = df_norm[['x', 'y', 'z']].mean()
            df_norm[['x', 'y', 'z']] -= center_of_mass.values
            metadata['center_of_mass'] = center_of_mass.to_dict()
        
        # Orient membrane proteins
        if orient_membrane:
            try:
                df_norm = orient_in_membrane(df_norm)
                metadata['membrane_oriented'] = True
            except Exception as e:
                print(f"Warning: Failed to orient {struct_id} in membrane: {e}")
                metadata['membrane_oriented'] = False
        
        # Annotate helices
        if annotate_helices:
            try:
                df_norm = annotate_transmembrane_helices(df_norm)
                metadata['helices_annotated'] = True
            except Exception as e:
                print(f"Warning: Failed to annotate helices for {struct_id}: {e}")
                metadata['helices_annotated'] = False
        
        # Normalize B-factors to [0, 1] range
        if 'b_factor' in df_norm.columns:
            b_min = df_norm['b_factor'].min()
            b_max = df_norm['b_factor'].max()
            if b_max > b_min:
                df_norm['b_factor_normalized'] = (df_norm['b_factor'] - b_min) / (b_max - b_min)
            else:
                df_norm['b_factor_normalized'] = 0.5
        
        normalized_structures[struct_id] = {
            'df_norm': df_norm,
            'metadata': metadata
        }
    
    return normalized_structures