"""
Structure data schemas and constants.

This module defines the standard column names, data types, and validation
functions for protein structure DataFrames in Protos.
"""

import pandas as pd
from typing import Dict, List, Optional

# CIF atom site columns - from mmCIF/PDBx format
CIF_COLUMNS = ['group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id', 'label_comp_id', 'label_asym_id',
               'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code', 'Cartn_x', 'Cartn_y', 'Cartn_z', 
               'occupancy', 'B_iso_or_equiv', 'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id', 
               'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num']

# Standard internal column names for structure DataFrames
STRUCT_COLUMNS = ['structure_id', 'group', 'auth_chain_id', 'label_chain_id', 'gen_seq_id', 'auth_seq_id', 'res_name',
                  'atom_id', 'atom_name', 'x', 'y', 'z', 'phi', 'omega', 'psi', 'res_name3l', 'res_name1l', 'grn']

# Data type definitions for structure columns
STRUCT_COLUMN_DTYPE = {
    'structure_id': str,  # Unique identifier for the structure
    'group': str,         # ATOM or HETATM
    'atom_id': int,       # Atom serial number
    'atom_name': str,     # Atom name (CA, CB, etc.)
    'element': str,       # Element symbol
    'alt_id': str,        # Alternate location indicator
    
    # Residue information
    'res_name': str,      # Residue name (3-letter code)
    'res_name3l': str,    # Same as res_name (for compatibility)
    'res_name1l': str,    # 1-letter amino acid code
    'auth_seq_id': int,   # Author residue sequence number
    'label_seq_id': int,  # Label residue sequence number
    'gen_seq_id': int,    # Generic sequence ID (sequential)
    'insertion': str,     # Insertion code
    'grn': str,           # Generic residue number (e.g., '3.50')
    
    # Chain information
    'auth_chain_id': str,   # Author chain ID
    'label_chain_id': str,  # Label chain ID
    'entity_id': int,       # Entity identifier
    
    # Coordinates and properties
    'x': float,           # X coordinate
    'y': float,           # Y coordinate
    'z': float,           # Z coordinate
    'occupancy': float,   # Occupancy
    'b_factor': float,    # B-factor/temperature factor
    'charge': str,        # Formal charge
    
    # Torsion angles (calculated)
    'phi': float,         # Phi angle (degrees)
    'psi': float,         # Psi angle (degrees)
    'omega': float,       # Omega angle (degrees)
    
    # Other
    'model_num': int,     # Model number
    'res_id': str,        # Composite residue ID
    'auth_comp_id': str   # Author compound ID
}

# Atom type constants
ALPHA_CARBON = 'CA'

# Column ordering for standardized output
SORTED_STRUCT_COLUMNS = ['structure_id', 'group', 'atom_id', 'atom_name', 'element', 
                         'auth_chain_id', 'auth_seq_id', 'res_name', 'res_name1l',
                         'x', 'y', 'z', 'occupancy', 'b_factor',
                         'phi', 'psi', 'omega', 'grn']


def validate_structure_df(df: pd.DataFrame, required_columns: Optional[List[str]] = None) -> bool:
    """
    Validate structure DataFrame against schema.
    
    Args:
        df: DataFrame to validate
        required_columns: Optional list of required columns. If None, uses STRUCT_COLUMNS
        
    Returns:
        True if valid, raises ValueError otherwise
    """
    if required_columns is None:
        required_columns = ['x', 'y', 'z', 'res_atom_name']  # Minimal requirements
    
    missing_cols = [col for col in required_columns if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Check data types for numeric columns
    numeric_cols = ['x', 'y', 'z']
    for col in numeric_cols:
        if col in df.columns and not pd.api.types.is_numeric_dtype(df[col]):
            raise ValueError(f"Column '{col}' must be numeric")
    
    return True


def apply_structure_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply standard data types to structure DataFrame columns.
    
    Args:
        df: DataFrame to process
        
    Returns:
        DataFrame with corrected data types
    """
    df_copy = df.copy()
    
    for col, dtype in STRUCT_COLUMN_DTYPE.items():
        if col in df_copy.columns:
            try:
                if dtype == float:
                    df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce')
                elif dtype == int:
                    # For integer columns, convert to nullable int to handle NaN
                    df_copy[col] = pd.to_numeric(df_copy[col], errors='coerce').astype('Int64')
                elif dtype == str:
                    df_copy[col] = df_copy[col].astype(str)
            except Exception as e:
                # Log but don't fail - some columns might have mixed types
                pass
    
    return df_copy