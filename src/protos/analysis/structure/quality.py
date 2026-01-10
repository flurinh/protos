"""
Structure quality assessment and validation functions.

This module provides functions for validating protein structure quality,
including bond lengths, clashes, missing atoms, and chirality checks.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple


# Standard bond lengths (in Angstroms) with tolerances
STANDARD_BONDS = {
    ('N', 'CA'): (1.458, 0.019),    # peptide N-CA
    ('CA', 'C'): (1.525, 0.021),    # CA-C
    ('C', 'O'): (1.231, 0.020),     # carbonyl
    ('C', 'N'): (1.329, 0.014),     # peptide bond
    ('CA', 'CB'): (1.530, 0.020),   # CA-CB
    ('C', 'S'): (1.810, 0.030),     # C-S (cysteine)
    ('S', 'S'): (2.030, 0.030),     # disulfide
}

# Van der Waals radii (in Angstroms)
VDW_RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'P': 1.80, 'H': 1.20, 'FE': 1.40, 'ZN': 1.39,
    'CA': 1.76, 'MG': 1.73, 'CL': 1.75, 'NA': 2.27
}

# Standard amino acid atom sets
AA_ATOMS = {
    'ALA': {'N', 'CA', 'C', 'O', 'CB'},
    'ARG': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'},
    'ASN': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2'},
    'ASP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2'},
    'CYS': {'N', 'CA', 'C', 'O', 'CB', 'SG'},
    'GLN': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2'},
    'GLU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'},
    'GLY': {'N', 'CA', 'C', 'O'},
    'HIS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'},
    'ILE': {'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1'},
    'LEU': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'},
    'LYS': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'},
    'MET': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE'},
    'PHE': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
    'PRO': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD'},
    'SER': {'N', 'CA', 'C', 'O', 'CB', 'OG'},
    'THR': {'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2'},
    'TRP': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'},
    'TYR': {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'},
    'VAL': {'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2'},
}


def check_missing_atoms(df: pd.DataFrame, 
                       residue_templates: Optional[Dict[str, set]] = None) -> pd.DataFrame:
    """
    Check for missing atoms in each residue.
    
    Args:
        df: Structure DataFrame
        residue_templates: Optional custom residue templates.
                          If None, uses standard amino acid templates
        
    Returns:
        DataFrame with columns: chain_id, residue_id, residue_name, missing_atoms
    """
    templates = residue_templates or AA_ATOMS
    
    # Reset index to access all columns
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Group by residue
    missing_info = []
    
    for (chain_id, res_id, res_name), group in df_reset.groupby(
        ['auth_chain_id', 'auth_seq_id', 'res_name3l']):
        
        # Get expected atoms for this residue type
        if res_name in templates:
            expected_atoms = templates[res_name]
            present_atoms = set(group['atom_name'])
            missing_atoms = expected_atoms - present_atoms
            
            if missing_atoms:
                missing_info.append({
                    'structure_id': group['structure_id'].iloc[0],
                    'chain_id': chain_id,
                    'residue_id': res_id,
                    'residue_name': res_name,
                    'missing_atoms': list(missing_atoms),
                    'n_missing': len(missing_atoms),
                    'completeness': len(present_atoms) / len(expected_atoms)
                })
    
    return pd.DataFrame(missing_info)


def validate_bond_lengths(df: pd.DataFrame, 
                         tolerance: float = 3.0,
                         custom_bonds: Optional[Dict[Tuple[str, str], Tuple[float, float]]] = None) -> pd.DataFrame:
    """
    Check if bond lengths are within expected ranges.
    
    Args:
        df: Structure DataFrame
        tolerance: Number of standard deviations to allow
        custom_bonds: Optional custom bond definitions
        
    Returns:
        DataFrame with columns: atom1_idx, atom2_idx, atom1, atom2, 
                               bond_type, length, expected, deviation
    """
    bonds = STANDARD_BONDS.copy()
    if custom_bonds:
        bonds.update(custom_bonds)
    
    # Reset index for easier access
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    bond_issues = []
    
    # Group by residue to check intraresidue bonds
    for (chain_id, res_id), group in df_reset.groupby(['auth_chain_id', 'auth_seq_id']):
        atoms = group.set_index('atom_name')
        
        # Check each bond type
        for (atom1_name, atom2_name), (expected, std) in bonds.items():
            if atom1_name in atoms.index and atom2_name in atoms.index:
                # Calculate distance
                coord1 = atoms.loc[atom1_name, ['x', 'y', 'z']].values
                coord2 = atoms.loc[atom2_name, ['x', 'y', 'z']].values
                
                # Handle multiple atoms with same name (alternative conformations)
                if coord1.ndim > 1:
                    coord1 = coord1[0]
                if coord2.ndim > 1:
                    coord2 = coord2[0]
                    
                distance = np.linalg.norm(coord1 - coord2)
                deviation = abs(distance - expected) / std
                
                if deviation > tolerance:
                    bond_issues.append({
                        'structure_id': group['structure_id'].iloc[0],
                        'chain_id': chain_id,
                        'residue_id': res_id,
                        'atom1': atom1_name,
                        'atom2': atom2_name,
                        'bond_type': f"{atom1_name}-{atom2_name}",
                        'length': distance,
                        'expected': expected,
                        'std_dev': std,
                        'deviation_sigma': deviation
                    })
    
    # Also check peptide bonds between consecutive residues
    for chain_id, chain_df in df_reset.groupby('auth_chain_id'):
        # Sort by residue number
        chain_sorted = chain_df.sort_values('auth_seq_id')
        residues = chain_sorted.groupby('auth_seq_id')
        
        res_list = list(residues)
        for i in range(len(res_list) - 1):
            res1_id, res1_atoms = res_list[i]
            res2_id, res2_atoms = res_list[i + 1]
            
            # Check if consecutive in sequence
            if res2_id - res1_id == 1:
                # Find C in res1 and N in res2
                c_atoms = res1_atoms[res1_atoms['atom_name'] == 'C']
                n_atoms = res2_atoms[res2_atoms['atom_name'] == 'N']
                
                if not c_atoms.empty and not n_atoms.empty:
                    c_coord = c_atoms.iloc[0][['x', 'y', 'z']].values
                    n_coord = n_atoms.iloc[0][['x', 'y', 'z']].values
                    
                    distance = np.linalg.norm(c_coord - n_coord)
                    expected, std = bonds[('C', 'N')]
                    deviation = abs(distance - expected) / std
                    
                    if deviation > tolerance:
                        bond_issues.append({
                            'structure_id': df_reset['structure_id'].iloc[0],
                            'chain_id': chain_id,
                            'residue_id': f"{res1_id}-{res2_id}",
                            'atom1': f"C({res1_id})",
                            'atom2': f"N({res2_id})",
                            'bond_type': "peptide",
                            'length': distance,
                            'expected': expected,
                            'std_dev': std,
                            'deviation_sigma': deviation
                        })
    
    return pd.DataFrame(bond_issues)


def check_clashes(df: pd.DataFrame, 
                 clash_cutoff: float = 2.0,
                 vdw_scale: float = 0.8) -> pd.DataFrame:
    """
    Find steric clashes (atoms too close).
    
    Args:
        df: Structure DataFrame
        clash_cutoff: Absolute minimum distance cutoff
        vdw_scale: Scale factor for VdW radii (0.8 = 80% of sum)
        
    Returns:
        DataFrame with clashing atom pairs
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Extract element from atom name
    df_reset['element'] = df_reset['atom_name'].str.extract(r'^([A-Z][a-z]?)')[0]
    
    # Get coordinates
    coords = df_reset[['x', 'y', 'z']].values
    n_atoms = len(coords)
    
    clashes = []
    
    # Check all atom pairs
    for i in range(n_atoms - 1):
        for j in range(i + 1, n_atoms):
            # Skip if same residue
            if (df_reset.iloc[i]['auth_chain_id'] == df_reset.iloc[j]['auth_chain_id'] and
                df_reset.iloc[i]['auth_seq_id'] == df_reset.iloc[j]['auth_seq_id']):
                continue
            
            # Calculate distance
            dist = np.linalg.norm(coords[i] - coords[j])
            
            # Get VdW radii
            elem1 = df_reset.iloc[i]['element']
            elem2 = df_reset.iloc[j]['element']
            r1 = VDW_RADII.get(elem1, 1.5)
            r2 = VDW_RADII.get(elem2, 1.5)
            
            # Check for clash
            min_allowed = max(clash_cutoff, vdw_scale * (r1 + r2))
            
            if dist < min_allowed:
                clashes.append({
                    'structure_id': df_reset.iloc[i]['structure_id'],
                    'atom1_idx': i,
                    'atom2_idx': j,
                    'chain1': df_reset.iloc[i]['auth_chain_id'],
                    'chain2': df_reset.iloc[j]['auth_chain_id'],
                    'resid1': df_reset.iloc[i]['auth_seq_id'],
                    'resid2': df_reset.iloc[j]['auth_seq_id'],
                    'resname1': df_reset.iloc[i].get('res_name3l', 'UNK'),
                    'resname2': df_reset.iloc[j].get('res_name3l', 'UNK'),
                    'atom1': df_reset.iloc[i]['atom_name'],
                    'atom2': df_reset.iloc[j]['atom_name'],
                    'distance': dist,
                    'min_allowed': min_allowed,
                    'severity': (min_allowed - dist) / min_allowed
                })
    
    return pd.DataFrame(clashes)


def calculate_b_factor_statistics(df: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate B-factor statistics.
    
    Args:
        df: Structure DataFrame with b_factor column
        
    Returns:
        Dictionary with statistics: mean, std, percentiles, high_b_factor info
    """
    # Check if b_factor column exists
    if 'b_factor' not in df.columns:
        return {
            'error': 'No b_factor column found',
            'mean': np.nan,
            'std': np.nan
        }
    
    # Reset index for analysis
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Overall statistics
    b_factors = df_reset['b_factor'].dropna()
    
    if b_factors.empty:
        return {
            'error': 'No valid b_factor values',
            'mean': np.nan,
            'std': np.nan
        }
    
    stats = {
        'mean': b_factors.mean(),
        'std': b_factors.std(),
        'median': b_factors.median(),
        'min': b_factors.min(),
        'max': b_factors.max(),
        'q25': b_factors.quantile(0.25),
        'q75': b_factors.quantile(0.75),
        'q95': b_factors.quantile(0.95),
    }
    
    # Statistics by chain
    chain_stats = {}
    for chain_id, chain_df in df_reset.groupby('auth_chain_id'):
        chain_b = chain_df['b_factor'].dropna()
        if not chain_b.empty:
            chain_stats[chain_id] = {
                'mean': chain_b.mean(),
                'std': chain_b.std(),
                'n_atoms': len(chain_b)
            }
    stats['by_chain'] = chain_stats
    
    # Statistics by residue type
    restype_stats = {}
    for res_name, res_df in df_reset.groupby('res_name3l'):
        res_b = res_df['b_factor'].dropna()
        if len(res_b) >= 5:  # Only include if enough atoms
            restype_stats[res_name] = {
                'mean': res_b.mean(),
                'std': res_b.std(),
                'n_atoms': len(res_b)
            }
    stats['by_residue_type'] = restype_stats
    
    # Identify high B-factor regions (>95th percentile)
    high_b_threshold = stats['q95']
    high_b_residues = df_reset[df_reset['b_factor'] > high_b_threshold].groupby(
        ['auth_chain_id', 'auth_seq_id', 'res_name3l']
    )['b_factor'].mean().sort_values(ascending=False)
    
    stats['high_b_residues'] = [
        {
            'chain': chain,
            'residue_id': res_id,
            'residue_name': res_name,
            'mean_b_factor': mean_b
        }
        for (chain, res_id, res_name), mean_b in high_b_residues.head(10).items()
    ]
    
    return stats


def validate_chirality(df: pd.DataFrame) -> pd.DataFrame:
    """
    Check amino acid chirality (L vs D).
    
    Uses the improper dihedral angle defined by N-CA-C-CB to determine chirality.
    L-amino acids should have a positive improper dihedral.
    
    Args:
        df: Structure DataFrame
        
    Returns:
        DataFrame with residues that have unusual chirality
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    chirality_issues = []
    
    # Group by residue
    for (chain_id, res_id, res_name), group in df_reset.groupby(
        ['auth_chain_id', 'auth_seq_id', 'res_name3l']):
        
        # Skip glycine (no chirality)
        if res_name == 'GLY':
            continue
        
        # Get required atoms
        atoms = group.set_index('atom_name')
        required = ['N', 'CA', 'C', 'CB']
        
        if all(atom in atoms.index for atom in required):
            # Get coordinates - extract and ensure float64 dtype
            def get_coords(atoms_df, atom_name):
                data = atoms_df.loc[atom_name, ['x', 'y', 'z']]
                if isinstance(data, pd.DataFrame):
                    # Multiple conformations - take first
                    return data.iloc[0].values.astype(np.float64)
                else:
                    # Single atom
                    return data.values.astype(np.float64)

            n_coord = get_coords(atoms, 'N')
            ca_coord = get_coords(atoms, 'CA')
            c_coord = get_coords(atoms, 'C')
            cb_coord = get_coords(atoms, 'CB')
            
            # Calculate improper dihedral angle
            # Vectors from CA
            v1 = n_coord - ca_coord
            v2 = c_coord - ca_coord
            v3 = cb_coord - ca_coord
            
            # Cross product gives normal to N-CA-C plane
            normal = np.cross(v1, v2)
            
            # Dot product with CB vector
            chirality_sign = np.dot(normal, v3)
            
            # L-amino acids should have positive sign
            if chirality_sign < 0:
                chirality_issues.append({
                    'structure_id': group['structure_id'].iloc[0],
                    'chain_id': chain_id,
                    'residue_id': res_id,
                    'residue_name': res_name,
                    'chirality': 'D',
                    'improper_dihedral': np.degrees(np.arcsin(chirality_sign / 
                                                              (np.linalg.norm(normal) * np.linalg.norm(v3))))
                })
    
    return pd.DataFrame(chirality_issues)


def check_chain_breaks(df: pd.DataFrame, max_ca_distance: float = 4.2) -> pd.DataFrame:
    """
    Check for chain breaks based on CA-CA distances.
    
    Args:
        df: Structure DataFrame  
        max_ca_distance: Maximum allowed CA-CA distance between consecutive residues
        
    Returns:
        DataFrame with detected chain breaks
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    chain_breaks = []
    
    # Check each chain
    for chain_id, chain_df in df_reset.groupby('auth_chain_id'):
        # Get CA atoms only
        ca_atoms = chain_df[chain_df['atom_name'] == 'CA'].sort_values('auth_seq_id')
        
        if len(ca_atoms) < 2:
            continue
        
        # Check consecutive CA distances
        for i in range(len(ca_atoms) - 1):
            res1 = ca_atoms.iloc[i]
            res2 = ca_atoms.iloc[i + 1]
            
            # Calculate distance
            coord1 = res1[['x', 'y', 'z']].values
            coord2 = res2[['x', 'y', 'z']].values
            distance = np.linalg.norm(coord1 - coord2)
            
            # Check if too far (indicating missing residues)
            if distance > max_ca_distance:
                chain_breaks.append({
                    'structure_id': res1['structure_id'],
                    'chain_id': chain_id,
                    'residue1_id': res1['auth_seq_id'],
                    'residue1_name': res1.get('res_name3l', 'UNK'),
                    'residue2_id': res2['auth_seq_id'],
                    'residue2_name': res2.get('res_name3l', 'UNK'),
                    'ca_distance': distance,
                    'gap_size': res2['auth_seq_id'] - res1['auth_seq_id'] - 1
                })
    
    return pd.DataFrame(chain_breaks)


def validate_structure_integrity(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """
    Run comprehensive structure validation.
    
    Performs all quality checks and returns results.
    
    Args:
        df: Structure DataFrame
        
    Returns:
        Dictionary with validation results:
        - missing_atoms: DataFrame of missing atoms
        - bond_issues: DataFrame of problematic bonds
        - clashes: DataFrame of steric clashes
        - chirality_issues: DataFrame of D-amino acids
        - chain_breaks: DataFrame of chain discontinuities
        - b_factor_stats: B-factor statistics dict
    """
    results = {}
    
    # Check missing atoms
    results['missing_atoms'] = check_missing_atoms(df)
    
    # Validate bond lengths
    results['bond_issues'] = validate_bond_lengths(df)
    
    # Check for clashes
    results['clashes'] = check_clashes(df)
    
    # Validate chirality
    results['chirality_issues'] = validate_chirality(df)
    
    # Check chain breaks
    results['chain_breaks'] = check_chain_breaks(df)
    
    # Calculate B-factor statistics
    results['b_factor_stats'] = calculate_b_factor_statistics(df)
    
    # Summary statistics
    results['summary'] = {
        'n_missing_atoms': len(results['missing_atoms']),
        'n_bond_issues': len(results['bond_issues']),
        'n_clashes': len(results['clashes']),
        'n_chirality_issues': len(results['chirality_issues']),
        'n_chain_breaks': len(results['chain_breaks']),
        'mean_b_factor': results['b_factor_stats'].get('mean', np.nan)
    }
    
    return results