"""
Schema definitions for standard DataFrame formats used across protos.

This module defines the standard column names, data types, and schemas for
DataFrames used in different components of the protos package. These schemas
ensure consistent data exchange between processors and utilities.

The module includes schemas for:
- Structure data (atom-level structure information)
- GRN data (GRN assignments and tables)
- Sequence data (sequence alignments and annotations)
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Union, Optional, Any

# -----------------------------------------------------------------------------
# Structure Schema Definitions
# -----------------------------------------------------------------------------

# Core structure columns that all structure DataFrames must contain
STRUCTURE_CORE_COLUMNS = {
    # Identifiers
    'pdb_id': str,        # PDB/mmCIF identifier (e.g., "1abc")
    'auth_chain_id': str, # Author-assigned chain identifier (e.g., "A")
    'gen_chain_id': str,  # Software-generated chain identifier
    
    # Residue information
    'auth_seq_id': int,   # Author-assigned residue number
    'gen_seq_id': int,    # Software-generated residue number (sequential)
    'res_name3l': str,    # 3-letter residue name (e.g., "ALA")
    'res_name1l': str,    # 1-letter residue name (e.g., "A")
    
    # Atom information
    'atom_id': int,       # Atom identifier within structure
    'atom_name': str,     # Atom name (e.g., "CA")
    'res_atom_name': str, # Combined residue and atom name (e.g., "ALA.CA")
    
    # Coordinates
    'x': float,           # X coordinate in Angstroms
    'y': float,           # Y coordinate in Angstroms
    'z': float,           # Z coordinate in Angstroms
}

# Extended structure columns for additional information
STRUCTURE_EXTENDED_COLUMNS = {
    # Backbone torsion angles
    'phi': float,         # Phi torsion angle in radians
    'psi': float,         # Psi torsion angle in radians
    'omega': float,       # Omega torsion angle in radians
    
    # Secondary structure
    'ss_type': str,       # Secondary structure type (e.g., "H", "E", "C")
    
    # B-factor and occupancy
    'b_factor': float,    # B-factor (temperature factor)
    'occupancy': float,   # Occupancy value
    
    # Additional classification
    'group': str,         # Group classification (e.g., "ATOM", "HETATM")
    'entity_id': str,     # Entity identifier
}

# All structure columns
STRUCTURE_ALL_COLUMNS = {
    **STRUCTURE_CORE_COLUMNS,
    **STRUCTURE_EXTENDED_COLUMNS
}

# Helper function to create a new empty structure DataFrame
def create_empty_structure_df() -> pd.DataFrame:
    """
    Create an empty DataFrame with the standard structure schema.
    
    Returns:
        Empty DataFrame with correct columns and data types
    """
    df = pd.DataFrame(columns=STRUCTURE_COLUMN_ORDER)
    for col, dtype in STRUCTURE_ALL_COLUMNS.items():
        df[col] = df[col].astype(dtype)
    return df

# Preferred order for structure columns
STRUCTURE_COLUMN_ORDER = [
    'pdb_id', 'group', 'auth_chain_id', 'gen_chain_id',
    'gen_seq_id', 'auth_seq_id', 'res_name3l', 'atom_id',
    'res_atom_name', 'atom_name', 'x', 'y', 'z',
    'phi', 'psi', 'omega', 'ss_type', 'b_factor', 'occupancy',
    'res_name1l', 'entity_id'
]

# Standard multi-index for structure DataFrames
STRUCTURE_INDEX_NAMES = ['structure_idx', 'atom_idx']

# Common atom selectors
ALPHA_CARBON = 'CA'
BACKBONE_ATOMS = ['N', 'CA', 'C', 'O']
SIDECHAIN_ATOMS = 'sidechain'  # Special value for sidechain selection

# -----------------------------------------------------------------------------
# GRN Schema Definitions
# -----------------------------------------------------------------------------

# Core GRN table columns
# Note: GRN tables are typically organized with GRNs as columns
# and protein/gene identifiers as rows
GRN_CORE_COLUMNS = {
    'id': str,            # Protein/gene identifier
    'name': str,          # Protein/gene name
    'species': str,       # Species name
    'family': str,        # Protein family
    'grn_system': str,    # GRN numbering system used
}

# Special values for GRN tables
GRN_GAP_SYMBOL = '-'      # Symbol for gaps in GRN tables
GRN_UNKNOWN_SYMBOL = 'X'  # Symbol for unknown residues

# GRN format specifications
GRN_PATTERNS = {
    'standard': r'^(\d+)x(\d+)$',  # e.g., 1x50
    'n_term': r'^n\.(\d+)$',       # e.g., n.10
    'c_term': r'^c\.(\d+)$',       # e.g., c.5
    'loop': r'^([1-8])([1-8])\.(\d{3})$'  # e.g., 12.003, 65.011
}

# Documentation for GRN formats
GRN_FORMAT_DOCS = {
    'standard': "Standard GRN format: <helix>x<position> (e.g., 1x50)",
    'n_term': "N-terminal format: n.<position> (e.g., n.10)",
    'c_term': "C-terminal format: c.<position> (e.g., c.5)",
    'loop': """Loop region format: <closer helix><further helix>.<distance> where:
            - First digit: Closer helix (1-8)
            - Second digit: Further helix (1-8)
            - Three-digit decimal: Distance from closer helix (001-999)
            Examples: 12.003 (between helix 1-2, closer to 1, distance 3)
                     65.011 (between helix 5-6, closer to 6, distance 11)"""
}

# -----------------------------------------------------------------------------
# Sequence Schema Definitions
# -----------------------------------------------------------------------------

# Standard columns for sequence alignment results
SEQUENCE_ALIGNMENT_COLUMNS = {
    'query_id': str,           # Query sequence identifier
    'target_id': str,          # Target sequence identifier
    'sequence_identity': float, # Sequence identity percentage
    'alignment_length': int,   # Length of the alignment
    'mismatches': int,         # Number of mismatches
    'gap_openings': int,       # Number of gap openings
    'query_start': int,        # Start position in query
    'query_end': int,          # End position in query
    'target_start': int,       # Start position in target
    'target_end': int,         # End position in target
    'e_value': float,          # E-value of alignment
    'bit_score': float,        # Bit score of alignment
}

# Standard columns for sequence feature tables
SEQUENCE_FEATURE_COLUMNS = {
    'seq_id': str,             # Sequence identifier
    'position': int,           # Position in sequence (1-indexed)
    'residue': str,            # Residue at position
    'feature_type': str,       # Type of feature
    'feature_value': Any,      # Value of feature
    'confidence': float,       # Confidence score (0-1)
}

# -----------------------------------------------------------------------------
# Function Definitions for Schema Validation
# -----------------------------------------------------------------------------

def validate_structure_df(df: pd.DataFrame, required_columns: List[str] = None) -> bool:
    """
    Validate a structure DataFrame against the standard schema.
    
    Args:
        df: DataFrame to validate
        required_columns: List of required columns (defaults to STRUCTURE_CORE_COLUMNS)
        
    Returns:
        bool: True if valid, False otherwise
    
    Raises:
        ValueError: If the DataFrame is missing required columns or has incorrect types
    """
    if required_columns is None:
        required_columns = list(STRUCTURE_CORE_COLUMNS.keys())
    
    # Check for required columns
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Structure DataFrame missing required columns: {missing_columns}")
    
    # Check column types for core columns
    for col, expected_type in STRUCTURE_CORE_COLUMNS.items():
        if col in df.columns:
            # Skip validation for empty DataFrames
            if df.empty:
                continue
                
            if col in ['x', 'y', 'z'] and df[col].dtype != np.float64:
                # Convert coordinate columns to float if they're not already
                df[col] = df[col].astype(np.float64)
            elif col in ['gen_seq_id', 'auth_seq_id', 'atom_id'] and not pd.api.types.is_integer_dtype(df[col].dtype):
                # Convert integer columns to int if they're not already
                df[col] = df[col].astype(int)
            elif col in ['pdb_id', 'auth_chain_id', 'gen_chain_id', 'res_name3l', 'res_name1l', 'atom_name', 'res_atom_name']:
                # Convert string columns to string if they're not already
                if not pd.api.types.is_string_dtype(df[col].dtype):
                    df[col] = df[col].astype(str)
    
    return True

def validate_grn_table(df: pd.DataFrame) -> bool:
    """
    Validate a GRN table DataFrame against the standard schema.
    
    Args:
        df: DataFrame to validate
        
    Returns:
        bool: True if valid, False otherwise
        
    Raises:
        ValueError: If the GRN table has invalid format
    """
    # Check for index (should be protein/gene identifiers)
    if df.index.name is None:
        df.index.name = 'id'
    
    # Check that columns are GRNs (should match GRN patterns)
    import re
    
    valid_grn_patterns = [re.compile(pattern) for pattern in GRN_PATTERNS.values()]
    invalid_grns = []
    
    for col in df.columns:
        # Skip metadata columns
        if col in GRN_CORE_COLUMNS:
            continue
            
        # Check if the column is a valid GRN
        if not any(pattern.match(str(col)) for pattern in valid_grn_patterns):
            invalid_grns.append(col)
    
    if invalid_grns:
        raise ValueError(f"GRN table contains invalid GRN columns: {invalid_grns}")
    
    return True

def validate_sequence_alignment_df(df: pd.DataFrame) -> bool:
    """
    Validate a sequence alignment DataFrame against the standard schema.
    
    Args:
        df: DataFrame to validate
        
    Returns:
        bool: True if valid, False otherwise
        
    Raises:
        ValueError: If the alignment DataFrame is missing required columns
    """
    required_cols = ['query_id', 'target_id', 'sequence_identity', 
                     'alignment_length', 'query_start', 'query_end',
                     'target_start', 'target_end']
    
    missing_columns = [col for col in required_cols if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Alignment DataFrame missing required columns: {missing_columns}")
    
    return True

def create_empty_structure_df() -> pd.DataFrame:
    """
    Create an empty structure DataFrame with the standard schema.
    
    Returns:
        pd.DataFrame: Empty DataFrame with the standard structure schema
    """
    # Create DataFrame with all columns
    df = pd.DataFrame(columns=STRUCTURE_COLUMN_ORDER)
    
    # Set column types
    for col, dtype in STRUCTURE_ALL_COLUMNS.items():
        if col in df.columns:
            df[col] = df[col].astype(dtype)
    
    # Create multi-index
    idx = pd.MultiIndex.from_arrays([[], []], names=STRUCTURE_INDEX_NAMES)
    df.index = idx
    
    return df

def create_empty_grn_table() -> pd.DataFrame:
    """
    Create an empty GRN table with the standard schema.
    
    Returns:
        pd.DataFrame: Empty DataFrame with the standard GRN table schema
    """
    # Create DataFrame with core columns
    df = pd.DataFrame(columns=list(GRN_CORE_COLUMNS.keys()))
    
    # Set index name
    df.index.name = 'id'
    
    return df

def create_empty_sequence_alignment_df() -> pd.DataFrame:
    """
    Create an empty sequence alignment DataFrame with the standard schema.
    
    Returns:
        pd.DataFrame: Empty DataFrame with the standard sequence alignment schema
    """
    # Create DataFrame with all columns
    df = pd.DataFrame(columns=list(SEQUENCE_ALIGNMENT_COLUMNS.keys()))
    
    # Set column types
    for col, dtype in SEQUENCE_ALIGNMENT_COLUMNS.items():
        df[col] = df[col].astype(dtype)
    
    return df