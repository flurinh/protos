"""
Structural property calculation functions.

This module provides functions for calculating various structural properties
including secondary structure, solvent accessibility, hydrophobic moments,
and binding site identification.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union


# Hydrophobicity scales (Kyte-Doolittle)
HYDROPHOBICITY = {
    'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5,
    'CYS': 2.5, 'GLN': -3.5, 'GLU': -3.5, 'GLY': -0.4,
    'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9,
    'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8,
    'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2
}


def calculate_secondary_structure(df: pd.DataFrame) -> pd.Series:
    """
    Assign secondary structure elements using DSSP-like algorithm.
    
    This is a simplified version that uses phi/psi angles and hydrogen bonding
    patterns to assign helix (H), sheet (E), and coil (C) states.
    
    Args:
        df: Structure DataFrame with backbone atoms
        
    Returns:
        Series with (chain_id, residue_id) -> ss_type mapping
    """
    # Reset index to access all data
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    ss_assignments = {}
    
    # Process each chain
    for chain_id, chain_df in df_reset.groupby('auth_chain_id'):
        # Get backbone atoms
        backbone = chain_df[chain_df['atom_name'].isin(['N', 'CA', 'C', 'O'])]
        
        # Group by residue and sort
        residues = backbone.groupby('auth_seq_id').filter(
            lambda x: set(x['atom_name']) == {'N', 'CA', 'C', 'O'}
        ).groupby('auth_seq_id')
        
        res_list = sorted(residues.groups.keys())
        
        # Calculate phi/psi angles if possible
        for i, res_id in enumerate(res_list):
            # Default to coil
            ss_type = 'C'
            
            # Need at least 3 consecutive residues for phi/psi
            if i > 0 and i < len(res_list) - 1:
                prev_res = res_list[i-1]
                next_res = res_list[i+1]
                
                # Check if consecutive
                if res_id - prev_res == 1 and next_res - res_id == 1:
                    # Get atoms for dihedral calculation
                    try:
                        # Phi: C(-1) - N - CA - C
                        c_prev = backbone[(backbone['auth_seq_id'] == prev_res) & 
                                        (backbone['atom_name'] == 'C')].iloc[0]
                        n_curr = backbone[(backbone['auth_seq_id'] == res_id) & 
                                        (backbone['atom_name'] == 'N')].iloc[0]
                        ca_curr = backbone[(backbone['auth_seq_id'] == res_id) & 
                                         (backbone['atom_name'] == 'CA')].iloc[0]
                        c_curr = backbone[(backbone['auth_seq_id'] == res_id) & 
                                        (backbone['atom_name'] == 'C')].iloc[0]
                        
                        # Psi: N - CA - C - N(+1)
                        n_next = backbone[(backbone['auth_seq_id'] == next_res) & 
                                        (backbone['atom_name'] == 'N')].iloc[0]
                        
                        # Calculate dihedral angles
                        phi = calculate_dihedral(
                            c_prev[['x', 'y', 'z']].values,
                            n_curr[['x', 'y', 'z']].values,
                            ca_curr[['x', 'y', 'z']].values,
                            c_curr[['x', 'y', 'z']].values
                        )
                        
                        psi = calculate_dihedral(
                            n_curr[['x', 'y', 'z']].values,
                            ca_curr[['x', 'y', 'z']].values,
                            c_curr[['x', 'y', 'z']].values,
                            n_next[['x', 'y', 'z']].values
                        )
                        
                        # Assign secondary structure based on Ramachandran regions
                        if -180 <= phi <= -45 and -60 <= psi <= 60:
                            ss_type = 'E'  # Beta sheet
                        elif -80 <= phi <= -35 and -50 <= psi <= -20:
                            ss_type = 'H'  # Alpha helix
                        
                    except (IndexError, KeyError):
                        pass
            
            ss_assignments[(chain_id, res_id)] = ss_type
    
    return pd.Series(ss_assignments)


def calculate_dihedral(p1: np.ndarray, p2: np.ndarray, 
                      p3: np.ndarray, p4: np.ndarray) -> float:
    """
    Calculate dihedral angle between four points.
    
    Args:
        p1, p2, p3, p4: 3D coordinates of four atoms
        
    Returns:
        Dihedral angle in degrees (-180 to 180)
    """
    # Vectors
    v1 = p2 - p1
    v2 = p3 - p2
    v3 = p4 - p3
    
    # Normal vectors to planes
    n1 = np.cross(v1, v2)
    n2 = np.cross(v2, v3)
    
    # Normalize
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    
    # Angle between planes
    cos_angle = np.dot(n1, n2)
    angle = np.arccos(np.clip(cos_angle, -1, 1))
    
    # Determine sign
    if np.dot(n1, v3) < 0:
        angle = -angle
    
    return np.degrees(angle)


def calculate_solvent_accessibility(df: pd.DataFrame, 
                                   probe_radius: float = 1.4,
                                   n_sphere_points: int = 100) -> pd.Series:
    """
    Calculate solvent accessible surface area per residue.
    
    This is a simplified SASA calculation using a rolling ball algorithm.
    
    Args:
        df: Structure DataFrame
        probe_radius: Radius of water probe in Angstroms
        n_sphere_points: Number of points to sample on sphere surface
        
    Returns:
        Series with (chain_id, residue_id) -> SASA value mapping
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Van der Waals radii
    vdw_radii = {
        'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
        'P': 1.80, 'H': 1.20, 'FE': 1.40, 'ZN': 1.39
    }
    
    # Generate uniform sphere points
    sphere_points = generate_sphere_points(n_sphere_points)
    
    sasa_by_residue = {}
    
    # Get all atom coordinates and radii
    coords = df_reset[['x', 'y', 'z']].values
    elements = df_reset['atom_name'].str.extract(r'^([A-Z][a-z]?)')[0]
    radii = elements.map(vdw_radii).fillna(1.5).values + probe_radius
    
    # Build neighbor list for efficiency
    # (simplified - in practice would use spatial data structure)
    
    # Calculate SASA for each residue
    for (chain_id, res_id), res_atoms in df_reset.groupby(['auth_chain_id', 'auth_seq_id']):
        res_sasa = 0
        
        for idx, atom in res_atoms.iterrows():
            atom_coord = atom[['x', 'y', 'z']].values
            atom_radius = radii[idx]
            
            # Test sphere points
            exposed_points = 0
            for sphere_point in sphere_points:
                test_point = atom_coord + atom_radius * sphere_point
                
                # Check if blocked by other atoms
                blocked = False
                for j, other_coord in enumerate(coords):
                    if j != idx:
                        dist = np.linalg.norm(test_point - other_coord)
                        if dist < radii[j]:
                            blocked = True
                            break
                
                if not blocked:
                    exposed_points += 1
            
            # Calculate exposed surface area
            atom_sasa = 4 * np.pi * atom_radius**2 * (exposed_points / n_sphere_points)
            res_sasa += atom_sasa
        
        sasa_by_residue[(chain_id, res_id)] = res_sasa
    
    return pd.Series(sasa_by_residue)


def generate_sphere_points(n_points: int) -> np.ndarray:
    """
    Generate uniformly distributed points on unit sphere.
    
    Uses golden angle spiral method.
    
    Args:
        n_points: Number of points to generate
        
    Returns:
        Array of shape (n_points, 3) with unit vectors
    """
    indices = np.arange(0, n_points, dtype=float) + 0.5
    
    # Golden angle
    theta = np.arccos(1 - 2 * indices / n_points)
    phi = np.pi * (1 + 5**0.5) * indices
    
    # Convert to Cartesian
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    
    return np.column_stack([x, y, z])


def calculate_hydrophobic_moment(df: pd.DataFrame, 
                                window_size: int = 11,
                                angle: float = 100.0) -> pd.Series:
    """
    Calculate hydrophobic moment for helical regions.
    
    The hydrophobic moment is a measure of amphipathicity in helices,
    useful for identifying membrane-spanning regions.
    
    Args:
        df: Structure DataFrame
        window_size: Size of sliding window (default 11 for ~3 helix turns)
        angle: Angle between consecutive residues in helix (default 100°)
        
    Returns:
        Series with (chain_id, residue_id) -> hydrophobic moment
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    moments = {}
    
    # Process each chain
    for chain_id, chain_df in df_reset.groupby('auth_chain_id'):
        # Get CA atoms with residue names
        ca_atoms = chain_df[chain_df['atom_name'] == 'CA'].sort_values('auth_seq_id')
        
        if len(ca_atoms) < window_size:
            continue
        
        res_ids = ca_atoms['auth_seq_id'].values
        res_names = ca_atoms['res_name3l'].values
        
        # Convert angle to radians
        angle_rad = np.radians(angle)
        
        # Sliding window
        for i in range(len(ca_atoms) - window_size + 1):
            window_res_ids = res_ids[i:i + window_size]
            window_res_names = res_names[i:i + window_size]
            
            # Calculate hydrophobic moment vector
            hx = 0
            hy = 0
            
            for j, res_name in enumerate(window_res_names):
                h_value = HYDROPHOBICITY.get(res_name, 0)
                theta = j * angle_rad
                hx += h_value * np.cos(theta)
                hy += h_value * np.sin(theta)
            
            # Magnitude of hydrophobic moment
            h_moment = np.sqrt(hx**2 + hy**2) / window_size
            
            # Assign to central residue
            central_idx = i + window_size // 2
            central_res_id = res_ids[central_idx]
            moments[(chain_id, central_res_id)] = h_moment
    
    return pd.Series(moments)


def identify_binding_sites(df: pd.DataFrame, 
                          ligand_df: Optional[pd.DataFrame] = None,
                          cavity_min_volume: float = 100.0) -> List[Dict]:
    """
    Identify potential binding pockets using cavity detection.
    
    This is a simplified version that identifies pockets based on:
    1. Known ligand positions (if provided)
    2. Cavities in the protein surface
    
    Args:
        df: Structure DataFrame
        ligand_df: Optional DataFrame with ligand atoms
        cavity_min_volume: Minimum cavity volume in Å³
        
    Returns:
        List of binding site descriptions with:
        - site_id: Unique identifier
        - center: Center coordinates
        - residues: List of pocket residues
        - volume: Estimated volume
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    binding_sites = []
    
    # 1. Identify sites from known ligands
    if ligand_df is not None and not ligand_df.empty:
        ligand_reset = ligand_df.reset_index() if isinstance(ligand_df.index, pd.MultiIndex) else ligand_df
        
        # Group ligands by residue (each ligand is a potential site)
        for (res_name, res_id), lig_atoms in ligand_reset.groupby(['res_name', 'auth_seq_id']):
            # Skip water and ions
            if res_name in ['HOH', 'WAT', 'NA', 'CL', 'CA', 'ZN', 'MG']:
                continue
            
            # Calculate ligand center
            lig_center = lig_atoms[['x', 'y', 'z']].mean().values
            
            # Find protein residues within 5Å of ligand
            protein_coords = df_reset[['x', 'y', 'z']].values
            distances = np.linalg.norm(protein_coords - lig_center, axis=1)
            
            nearby_indices = np.where(distances <= 5.0)[0]
            nearby_atoms = df_reset.iloc[nearby_indices]
            
            # Get unique residues
            pocket_residues = nearby_atoms.groupby(['auth_chain_id', 'auth_seq_id']).size()
            
            if len(pocket_residues) >= 3:  # At least 3 residues
                binding_sites.append({
                    'site_id': f'ligand_{res_name}_{res_id}',
                    'type': 'ligand_binding',
                    'ligand': res_name,
                    'center': lig_center.tolist(),
                    'residues': [
                        {'chain': chain, 'residue_id': res_id}
                        for (chain, res_id) in pocket_residues.index
                    ],
                    'n_residues': len(pocket_residues),
                    'volume': estimate_pocket_volume(nearby_atoms)
                })
    
    # 2. Identify cavities using alpha shapes (simplified)
    # Group atoms by spatial proximity
    coords = df_reset[['x', 'y', 'z']].values
    
    # Find potential cavity centers using grid-based approach
    grid_spacing = 3.0
    x_range = np.arange(coords[:, 0].min(), coords[:, 0].max(), grid_spacing)
    y_range = np.arange(coords[:, 1].min(), coords[:, 1].max(), grid_spacing)
    z_range = np.arange(coords[:, 2].min(), coords[:, 2].max(), grid_spacing)
    
    cavity_id = 1
    for x in x_range:
        for y in y_range:
            for z in z_range:
                grid_point = np.array([x, y, z])
                
                # Check if point is inside protein but not too close to atoms
                distances = np.linalg.norm(coords - grid_point, axis=1)
                
                # Inside protein: some atoms within 8Å but none within 3Å
                if distances.min() > 3.0 and (distances < 8.0).sum() > 20:
                    # Potential cavity center
                    nearby_indices = np.where(distances <= 8.0)[0]
                    cavity_atoms = df_reset.iloc[nearby_indices]
                    
                    # Check if forms enclosed pocket
                    if is_enclosed_pocket(cavity_atoms, grid_point):
                        pocket_residues = cavity_atoms.groupby(['auth_chain_id', 'auth_seq_id']).size()
                        
                        volume = estimate_pocket_volume(cavity_atoms)
                        if volume >= cavity_min_volume:
                            binding_sites.append({
                                'site_id': f'cavity_{cavity_id}',
                                'type': 'predicted_cavity',
                                'center': grid_point.tolist(),
                                'residues': [
                                    {'chain': chain, 'residue_id': res_id}
                                    for (chain, res_id) in pocket_residues.index
                                ],
                                'n_residues': len(pocket_residues),
                                'volume': volume
                            })
                            cavity_id += 1
    
    return binding_sites


def estimate_pocket_volume(pocket_atoms: pd.DataFrame) -> float:
    """
    Estimate volume of a binding pocket.
    
    Uses convex hull approximation.
    
    Args:
        pocket_atoms: DataFrame of atoms forming the pocket
        
    Returns:
        Estimated volume in Å³
    """
    coords = pocket_atoms[['x', 'y', 'z']].values
    
    if len(coords) < 4:
        return 0.0
    
    # Simple volume estimation using bounding box
    # (In practice, would use convex hull or alpha shapes)
    min_coords = coords.min(axis=0)
    max_coords = coords.max(axis=0)
    
    # Bounding box volume reduced by packing factor
    bbox_volume = np.prod(max_coords - min_coords)
    estimated_volume = bbox_volume * 0.5  # Rough estimate
    
    return estimated_volume


def is_enclosed_pocket(atoms: pd.DataFrame, center: np.ndarray, 
                      n_directions: int = 26) -> bool:
    """
    Check if a point is in an enclosed pocket.
    
    Tests if rays from the center hit atoms in multiple directions.
    
    Args:
        atoms: DataFrame of nearby atoms
        center: Center point to test
        n_directions: Number of directions to test
        
    Returns:
        True if pocket appears enclosed
    """
    coords = atoms[['x', 'y', 'z']].values
    
    # Test directions (26 = faces, edges, and corners of cube)
    directions = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                directions.append(np.array([dx, dy, dz], dtype=float))
    
    # Normalize directions
    directions = [d / np.linalg.norm(d) for d in directions]
    
    # Count how many directions are blocked
    blocked_directions = 0
    
    for direction in directions:
        # Check if ray hits any atom
        for coord in coords:
            # Vector from center to atom
            to_atom = coord - center
            
            # Project onto direction
            projection = np.dot(to_atom, direction)
            
            # Check if in front and close to ray
            if projection > 0:
                # Distance from ray
                ray_dist = np.linalg.norm(to_atom - projection * direction)
                if ray_dist < 3.0:  # Within 3Å of ray
                    blocked_directions += 1
                    break
    
    # Pocket is enclosed if most directions are blocked
    return blocked_directions >= len(directions) * 0.7


def calculate_electrostatic_potential(df: pd.DataFrame, 
                                     grid_spacing: float = 1.0,
                                     distance_cutoff: float = 20.0) -> np.ndarray:
    """
    Calculate simple electrostatic potential on a grid.
    
    Uses Coulomb's law with partial charges.
    
    Args:
        df: Structure DataFrame
        grid_spacing: Grid resolution in Angstroms
        distance_cutoff: Maximum distance for calculations
        
    Returns:
        3D grid of potential values with grid info as attributes
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Partial charges by atom type (simplified)
    charges = {
        'N': -0.3, 'O': -0.5, 'C': 0.0, 'S': 0.0,
        'NZ': 1.0,  # Lysine
        'NH1': 0.5, 'NH2': 0.5,  # Arginine
        'OD1': -0.5, 'OD2': -0.5,  # Aspartate
        'OE1': -0.5, 'OE2': -0.5,  # Glutamate
    }
    
    # Get atom charges
    atom_charges = df_reset['atom_name'].map(charges).fillna(0.0).values
    atom_coords = df_reset[['x', 'y', 'z']].values
    
    # Define grid
    min_coords = atom_coords.min(axis=0) - 5
    max_coords = atom_coords.max(axis=0) + 5
    
    x_grid = np.arange(min_coords[0], max_coords[0], grid_spacing)
    y_grid = np.arange(min_coords[1], max_coords[1], grid_spacing)
    z_grid = np.arange(min_coords[2], max_coords[2], grid_spacing)
    
    # Calculate potential on grid
    potential = np.zeros((len(x_grid), len(y_grid), len(z_grid)))
    
    # Coulomb constant (simplified units)
    k_e = 332.0  # kcal·Å/mol·e²
    
    for i, x in enumerate(x_grid):
        for j, y in enumerate(y_grid):
            for k, z in enumerate(z_grid):
                grid_point = np.array([x, y, z])
                
                # Sum contributions from all charged atoms
                for atom_idx, charge in enumerate(atom_charges):
                    if charge != 0:
                        distance = np.linalg.norm(grid_point - atom_coords[atom_idx])
                        
                        if distance < distance_cutoff and distance > 0.1:
                            # Coulomb potential with distance-dependent dielectric
                            dielectric = 4.0 * distance  # Simple distance-dependent
                            potential[i, j, k] += k_e * charge / (dielectric * distance)
    
    # Add grid information as attributes (store separately)
    grid_info = {
        'x_grid': x_grid,
        'y_grid': y_grid, 
        'z_grid': z_grid,
        'spacing': grid_spacing
    }
    
    return potential, grid_info


def identify_surface_residues(df: pd.DataFrame, 
                            sasa_threshold: float = 20.0) -> pd.Series:
    """
    Identify surface-exposed residues based on SASA.
    
    Args:
        df: Structure DataFrame
        sasa_threshold: Minimum SASA to be considered surface-exposed
        
    Returns:
        Series with (chain_id, residue_id) -> is_surface mapping
    """
    # Calculate SASA
    sasa_values = calculate_solvent_accessibility(df, n_sphere_points=50)
    
    # Identify surface residues
    surface_residues = sasa_values > sasa_threshold
    
    return surface_residues