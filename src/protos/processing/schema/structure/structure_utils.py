"""
Structure utilities with standardized schemas.

This module provides refactored versions of structure processing functions
that use the standardized schemas defined in schema_definitions.py.
"""

import os
import sys
import logging
import numpy as np
import pandas as pd
import gemmi
from gemmi import cif
from math import degrees
from pathlib import Path
from typing import Dict, List, Tuple, Union, Optional, Any

# Import path management
from protos.io.paths import get_structure_path

# Import schema definitions and interfaces
from protos.processing.schema.schema_definitions import (
    STRUCTURE_CORE_COLUMNS,
    STRUCTURE_EXTENDED_COLUMNS,
    STRUCTURE_ALL_COLUMNS,
    STRUCTURE_COLUMN_ORDER,
    create_empty_structure_df,
    validate_structure_df
)
from protos.processing.schema.interface_definitions import StructureInterface
from protos.processing.schema.conversion_utilities import (
    cif_to_structure_df,
    three_to_one_letter_code
)

# Configure logger
logger = logging.getLogger(__name__)


def load_structure(file_path: Union[str, Path], structure_id: Optional[str] = None) -> pd.DataFrame:
    """
    Load a structure from a file.
    
    This is an implementation of the StructureInterface.load_structure method
    that uses standardized schemas and the path management system.
    
    Args:
        file_path: Path to the structure file or PDB ID
        structure_id: Identifier for the structure (defaults to filename without extension)
        
    Returns:
        DataFrame containing the structure data according to the standard schema
        
    Raises:
        FileNotFoundError: If the structure file is not found
        ValueError: If the structure cannot be loaded
    """
    # Convert to Path object for consistent handling
    file_path = Path(str(file_path))
    
    # If the file doesn't exist, it might be a PDB ID
    if not file_path.exists() and not os.path.isabs(file_path):
        # Try to interpret as a PDB ID and get the path
        pdb_id = str(file_path).split('.')[0]  # Strip extension if any
        file_path = get_structure_path(pdb_id, create_if_missing=False)
        
        # Use PDB ID as structure_id if not provided
        if structure_id is None:
            structure_id = pdb_id
    
    # Check if file exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Structure file not found: {file_path}")
    
    # If no structure_id provided, use filename without extension
    if structure_id is None:
        structure_id = file_path.stem
    
    # Load atom data from mmCIF file
    try:
        # Get CIF block
        block = cif.read_file(str(file_path)).sole_block()
        
        # Convert to DataFrame
        atom_df = cif_block_to_df(block, structure_id)
        
        # Load residue-wise features (torsion angles)
        torsion_df = extract_torsion_angles(file_path, structure_id=0)
        
        # Merge atom and residue data
        structure_df = pd.merge(
            atom_df, 
            torsion_df, 
            how='outer', 
            on=['gen_chain_id', 'gen_seq_id']
        )
        
        # Add 1-letter residue codes
        structure_df['res_name1l'] = structure_df['res_name3l'].apply(three_to_one_letter_code)
        
        # Create res_atom_name column if not present
        if 'res_atom_name' not in structure_df.columns:
            structure_df['res_atom_name'] = structure_df['res_name3l'] + '.' + structure_df['atom_name']
        
        # Validate the resulting DataFrame
        if not validate_structure_df(structure_df):
            logger.warning(f"Structure validation failed for {structure_id}")
        
        # Ensure columns are in the standard order
        for col in STRUCTURE_COLUMN_ORDER:
            if col not in structure_df.columns:
                structure_df[col] = np.nan
        
        # Return a DataFrame with standard column order
        return structure_df[STRUCTURE_COLUMN_ORDER]
        
    except Exception as e:
        logger.error(f"Error loading structure {structure_id}: {e}")
        raise ValueError(f"Failed to load structure: {e}")


def cif_block_to_df(block: gemmi.cif.Block, structure_id: str) -> pd.DataFrame:
    """
    Convert a CIF block to a DataFrame with standard structure schema.
    
    Args:
        block: CIF block to convert
        structure_id: Identifier for the structure
        
    Returns:
        DataFrame with structure data
    """
    # CIF columns to extract
    cif_columns = [
        'group_PDB', 'auth_asym_id', 'label_asym_id', 'label_seq_id', 'auth_seq_id',
        'label_comp_id', 'id', 'label_atom_id', 'type_symbol', 'Cartn_x', 'Cartn_y', 'Cartn_z'
    ]
    
    # Extract data from CIF block
    atom_rows = []
    table = block.find('_atom_site.', cif_columns)
    for row in table:
        atom_rows.append([structure_id] + list(row))
    
    # Create DataFrame
    cols = ['pdb_id'] + cif_columns
    atom_df = pd.DataFrame(data=atom_rows, columns=cols)
    
    # Perform column mapping to standard schema
    column_mapping = {
        'group_PDB': 'group',
        'label_asym_id': 'gen_chain_id',
        'auth_asym_id': 'auth_chain_id',
        'label_seq_id': 'gen_seq_id',
        'auth_seq_id': 'auth_seq_id',
        'label_comp_id': 'res_name3l',
        'label_atom_id': 'atom_name',
        'Cartn_x': 'x',
        'Cartn_y': 'y',
        'Cartn_z': 'z',
        'id': 'atom_id',
    }
    
    # Rename columns
    atom_df = atom_df.rename(columns=column_mapping)
    
    # Convert label_seq_id to integer
    atom_df['gen_seq_id'] = atom_df['gen_seq_id'].apply(
        lambda x: int(x) if x != '.' else np.nan
    )
    
    # Convert coordinates to float
    for col in ['x', 'y', 'z']:
        atom_df[col] = atom_df[col].astype(float)
    
    # Convert atom_id to integer
    atom_df['atom_id'] = atom_df['atom_id'].astype(int)
    
    return atom_df


def extract_torsion_angles(file_path: Union[str, Path], structure_id: int = 0) -> pd.DataFrame:
    """
    Extract torsion angles (phi, psi, omega) from a structure file.
    
    Args:
        file_path: Path to the structure file
        structure_id: Model ID to use (default: 0)
        
    Returns:
        DataFrame with torsion angles
    """
    # Load structure with gemmi
    st = gemmi.read_structure(str(file_path), merge_chain_parts=True)
    logger.debug(f"Found {len(st)} model(s), using model {structure_id}.")
    
    # Use specified model
    model = st[structure_id]
    
    # Extract residue-wise features
    res_rows = []
    for chain in model:
        for r, res in enumerate(chain.get_polymer()):
            # Calculate phi and psi angles
            try:
                prev_res = chain.previous_residue(res)
                next_res = chain.next_residue(res)
                phi, psi = gemmi.calculate_phi_psi(prev_res, res, next_res)
            except Exception:
                phi, psi = np.nan, np.nan
            
            # Calculate omega angle
            try:
                next_res = chain.next_residue(res)
                omega = gemmi.calculate_omega(res, next_res)
            except Exception:
                omega = np.nan
            
            # Add row with residue data
            res_rows.append([
                res.label_seq, 
                res.subchain, 
                degrees(phi), 
                degrees(omega), 
                degrees(psi)
            ])
    
    # Create DataFrame
    cols = ['gen_seq_id', 'gen_chain_id', 'phi', 'omega', 'psi']
    res_df = pd.DataFrame(data=res_rows, columns=cols)
    
    # Convert seq_id to integer
    res_df['gen_seq_id'] = res_df['gen_seq_id'].astype(int)
    
    return res_df


def normalize_structures(structure_dfs: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """
    Normalize multiple structure DataFrames for consistent comparison.
    
    This function centers and aligns structures to ensure they are in a
    consistent coordinate frame. It is a standardized version of the original
    function that uses the standard schema for structure DataFrames.
    
    Args:
        structure_dfs: Dictionary mapping structure IDs to structure DataFrames
        
    Returns:
        Dictionary mapping structure IDs to normalized structure DataFrames
    """
    # Create a copy of the input dictionary to avoid modifying the original
    normalized_structures = {}
    
    # Process each structure
    for pdb_id, structure_df in structure_dfs.items():
        # Get CA and retinal coordinates
        df_ca, df_ret = get_ca_ret_coords(structure_df)
        
        # Center structure on backbone
        mean_backbone = df_ca[['x', 'y', 'z']].mean(axis=0)
        
        # Create a copy of the structure DataFrame
        df_norm = structure_df.copy()
        
        # Center the structure
        for col in ['x', 'y', 'z']:
            df_norm[col] = df_norm[col] - mean_backbone[col]
        
        # Center retinal coordinates
        df_ret_norm = None
        if not df_ret.empty:
            df_ret_norm = df_ret.copy()
            for col in ['x', 'y', 'z']:
                df_ret_norm[col] = df_ret_norm[col] - mean_backbone[col]
        
        # Store normalized structure and retinal coordinates
        normalized_structures[pdb_id] = {
            'structure': df_norm,
            'ca_coords': df_ca,
            'retinal': df_ret_norm,
            'center': mean_backbone
        }
    
    return normalized_structures


def get_ca_ret_coords(structure_df: pd.DataFrame, chain_id: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Extract alpha carbon and retinal coordinates from a structure.
    
    This is a standardized version of the original function that uses the
    standard schema for structure DataFrames.
    
    Args:
        structure_df: Structure DataFrame in standard format
        chain_id: Chain ID to use (default: first chain in the structure)
        
    Returns:
        Tuple of (alpha_carbon_coords, retinal_coords) as DataFrames
    """
    # If chain_id is not provided, use the first chain in the structure
    if chain_id is None:
        chain_id = structure_df['auth_chain_id'].iloc[0]
    
    # Get alpha carbon atoms for the specified chain
    df_ca = structure_df[
        (structure_df['atom_name'] == 'CA') & 
        (structure_df['auth_chain_id'] == chain_id)
    ].reset_index(drop=True)
    
    # Get coordinates
    ca_coords = df_ca[['x', 'y', 'z']].astype(np.float64)
    
    # Get retinal atoms (RET residue)
    df_ret = structure_df[
        structure_df['res_name3l'] == 'RET'
    ].reset_index(drop=True)
    
    # If no retinal is found, return empty DataFrame
    if df_ret.empty:
        ret_coords = pd.DataFrame(columns=['x', 'y', 'z'])
    else:
        ret_coords = df_ret[['x', 'y', 'z']].astype(np.float64)
    
    return ca_coords, ret_coords