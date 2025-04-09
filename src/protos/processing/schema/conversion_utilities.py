"""
Conversion utilities for transforming data between different formats.

This module provides utility functions for converting between different data formats
used in the protos package. These utilities ensure consistent data exchange and
compatibility between components.

The module includes functions for converting:
- Between DataFrames and dictionaries
- Between different DataFrame schemas
- Between different file formats
- Between different data representations
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Union, Optional, Any
from pathlib import Path
import os
import re

# Import standard schemas
from protos.processing.schema.schema_definitions import (
    STRUCTURE_CORE_COLUMNS,
    STRUCTURE_EXTENDED_COLUMNS,
    STRUCTURE_ALL_COLUMNS,
    STRUCTURE_COLUMN_ORDER,
    GRN_PATTERNS,
    SEQUENCE_ALIGNMENT_COLUMNS
)

# -----------------------------------------------------------------------------
# Structure Data Conversions
# -----------------------------------------------------------------------------

def cif_to_structure_df(cif_df: pd.DataFrame, structure_id: Optional[str] = None) -> pd.DataFrame:
    """
    Convert a DataFrame from CIF format to the standard structure format.
    
    Args:
        cif_df: DataFrame in CIF format
        structure_id: Identifier for the structure (defaults to first PDB ID found)
        
    Returns:
        DataFrame in standard structure format
    """
    # Create an empty DataFrame with the correct columns
    structure_df = pd.DataFrame(columns=list(STRUCTURE_ALL_COLUMNS.keys()))
    
    # Set PDB ID if not provided
    if structure_id is None and 'label_entity_id' in cif_df.columns:
        structure_id = cif_df['label_entity_id'].iloc[0]
    elif structure_id is None:
        structure_id = 'UNK'
    
    # Map CIF columns to structure columns
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
    
    # Copy columns from CIF to structure DataFrame
    for cif_col, struct_col in column_mapping.items():
        if cif_col in cif_df.columns:
            structure_df[struct_col] = cif_df[cif_col]
    
    # Add PDB ID column
    structure_df['pdb_id'] = structure_id
    
    # Create res_atom_name column
    if 'res_name3l' in structure_df.columns and 'atom_name' in structure_df.columns:
        structure_df['res_atom_name'] = structure_df['res_name3l'] + '.' + structure_df['atom_name']
    
    # Convert 3-letter residue codes to 1-letter codes
    if 'res_name3l' in structure_df.columns:
        structure_df['res_name1l'] = structure_df['res_name3l'].apply(three_to_one_letter_code)
    
    # Set column data types
    for col, dtype in STRUCTURE_ALL_COLUMNS.items():
        if col in structure_df.columns:
            structure_df[col] = structure_df[col].astype(dtype)
    
    return structure_df

def three_to_one_letter_code(three_letter: str) -> str:
    """
    Convert a three-letter amino acid code to one-letter code.
    
    Args:
        three_letter: Three-letter amino acid code (e.g., 'ALA')
        
    Returns:
        One-letter amino acid code (e.g., 'A')
    """
    # Standard amino acid mapping
    aa_map = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        # Non-standard amino acids
        'MSE': 'M',  # Selenomethionine
        'UNK': 'X',  # Unknown
        'HYP': 'P',  # Hydroxyproline
        'PCA': 'E',  # Pyroglutamic acid
        'SEP': 'S',  # Phosphoserine
        'TPO': 'T',  # Phosphothreonine
        'PTR': 'Y',  # Phosphotyrosine
        'ACE': 'X',  # Acetyl group
        'NH2': 'X',  # Amide group
        # Add more mappings as needed
    }
    
    # Handle None, NaN, or empty strings
    if three_letter is None or pd.isna(three_letter) or three_letter == '':
        return 'X'
    
    # Convert to uppercase for consistency
    three_letter = three_letter.upper()
    
    # Return mapped value or 'X' for unknown
    return aa_map.get(three_letter, 'X')


def structure_df_to_cif(structure_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert a DataFrame from standard structure format to CIF format.
    
    Args:
        structure_df: DataFrame in standard structure format
        
    Returns:
        DataFrame in CIF format
    """
    # Create an empty DataFrame for CIF format
    cif_df = pd.DataFrame()
    
    # Map structure columns to CIF columns
    column_mapping = {
        'group': 'group_PDB',
        'gen_chain_id': 'label_asym_id',
        'auth_chain_id': 'auth_asym_id',
        'gen_seq_id': 'label_seq_id',
        'auth_seq_id': 'auth_seq_id',
        'res_name3l': 'label_comp_id',
        'atom_name': 'label_atom_id',
        'x': 'Cartn_x',
        'y': 'Cartn_y',
        'z': 'Cartn_z',
        'atom_id': 'id',
    }
    
    # Copy columns from structure to CIF DataFrame
    for struct_col, cif_col in column_mapping.items():
        if struct_col in structure_df.columns:
            cif_df[cif_col] = structure_df[struct_col]
    
    # Add required CIF columns that might not be present
    if 'type_symbol' not in cif_df.columns and 'atom_name' in structure_df.columns:
        cif_df['type_symbol'] = structure_df['atom_name'].apply(lambda x: x[0] if len(x) > 0 else 'X')
    
    return cif_df

def dict_to_structure_df(structure_dict: Dict[str, Any]) -> pd.DataFrame:
    """
    Convert a dictionary to a standard structure DataFrame.
    
    Args:
        structure_dict: Dictionary with structure data
        
    Returns:
        DataFrame in standard structure format
    """
    # Create an empty DataFrame with the correct columns
    structure_df = pd.DataFrame(columns=list(STRUCTURE_ALL_COLUMNS.keys()))
    
    # Copy data from dictionary to DataFrame
    for key, value in structure_dict.items():
        if key in STRUCTURE_ALL_COLUMNS and len(value) > 0:
            structure_df[key] = value
    
    # Set column data types
    for col, dtype in STRUCTURE_ALL_COLUMNS.items():
        if col in structure_df.columns:
            structure_df[col] = structure_df[col].astype(dtype)
    
    return structure_df

def structure_df_to_dict(structure_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Convert a standard structure DataFrame to a dictionary.
    
    Args:
        structure_df: DataFrame in standard structure format
        
    Returns:
        Dictionary with structure data
    """
    # Create a dictionary from the DataFrame
    structure_dict = {}
    
    # Copy data from DataFrame to dictionary
    for col in structure_df.columns:
        structure_dict[col] = structure_df[col].tolist()
    
    return structure_dict

def three_to_one_letter_code(three_letter_code: str) -> str:
    """
    Convert a three-letter amino acid code to a one-letter code.
    
    Args:
        three_letter_code: Three-letter amino acid code
        
    Returns:
        One-letter amino acid code
    """
    # Define the mapping from three-letter to one-letter codes
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
        'ASX': 'B', 'GLX': 'Z', 'XAA': 'X', 'UNK': 'X', '': '-'
    }
    
    # Convert the three-letter code to uppercase
    three_letter_code = str(three_letter_code).upper()
    
    # Return the one-letter code or 'X' if not found
    return three_to_one.get(three_letter_code, 'X')

def one_to_three_letter_code(one_letter_code: str) -> str:
    """
    Convert a one-letter amino acid code to a three-letter code.
    
    Args:
        one_letter_code: One-letter amino acid code
        
    Returns:
        Three-letter amino acid code
    """
    # Define the mapping from one-letter to three-letter codes
    one_to_three = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
        'B': 'ASX', 'Z': 'GLX', 'X': 'UNK', '-': ''
    }
    
    # Convert the one-letter code to uppercase
    one_letter_code = str(one_letter_code).upper()
    
    # Return the three-letter code or 'UNK' if not found
    return one_to_three.get(one_letter_code, 'UNK')

# -----------------------------------------------------------------------------
# GRN Data Conversions
# -----------------------------------------------------------------------------

def grn_str_to_float(grn_str: str) -> float:
    """
    Convert a GRN string to a float for sorting.
    
    Args:
        grn_str: GRN string (e.g., '1x50')
        
    Returns:
        Float representation (e.g., 1.5)
    """
    # Check for the different GRN patterns
    for pattern_name, pattern in GRN_PATTERNS.items():
        match = re.match(pattern, grn_str)
        if match:
            if pattern_name == 'standard':
                # Format: 1x50 -> 1.5
                helix = int(match.group(1))
                pos = int(match.group(2))
                return helix + pos / 100.0
            elif pattern_name == 'n_term':
                # Format: n.10 -> -0.1
                pos = int(match.group(1))
                return -pos / 100.0
            elif pattern_name == 'c_term':
                # Format: c.5 -> 8.05 (assuming 8 helices)
                pos = int(match.group(1))
                return 8.0 + pos / 100.0
            elif pattern_name == 'loop':
                # Format: 2.45 -> 2.045
                helix = int(match.group(1))
                pos = int(match.group(2))
                return helix + pos / 1000.0
    
    # Default value for unrecognized formats
    return 0.0

def grn_float_to_str(grn_float: float) -> str:
    """
    Convert a GRN float to a string.
    
    Args:
        grn_float: Float representation (e.g., 1.5)
        
    Returns:
        GRN string (e.g., '1x50')
    """
    if grn_float < 0:
        # N-terminal region: -0.1 -> n.10
        pos = int(abs(grn_float) * 100)
        return f"n.{pos}"
    elif grn_float > 8:
        # C-terminal region: 8.05 -> c.5
        pos = int((grn_float - 8.0) * 100)
        return f"c.{pos}"
    else:
        # Check if it's a loop or standard notation
        helix = int(grn_float)
        decimal_part = grn_float - helix
        
        if decimal_part < 0.1:
            # Loop notation: 2.045 -> 2.45
            pos = int(decimal_part * 1000)
            return f"{helix}.{pos}"
        else:
            # Standard notation: 1.5 -> 1x50
            pos = int(decimal_part * 100)
            return f"{helix}x{pos}"

def sort_grns(grn_list: List[str]) -> List[str]:
    """
    Sort a list of GRNs in standard order.
    
    Args:
        grn_list: List of GRN strings
        
    Returns:
        Sorted list of GRNs
    """
    # Convert GRNs to floats for sorting
    grn_floats = [(grn, grn_str_to_float(grn)) for grn in grn_list]
    
    # Sort by float values
    sorted_grns = sorted(grn_floats, key=lambda x: x[1])
    
    # Return the sorted GRN strings
    return [grn for grn, _ in sorted_grns]

def grn_mapping_to_df(grn_mapping: Dict[int, str], sequence: Optional[str] = None) -> pd.DataFrame:
    """
    Convert a GRN mapping dictionary to a DataFrame.
    
    Args:
        grn_mapping: Dictionary mapping sequence positions to GRNs
        sequence: The amino acid sequence (optional)
        
    Returns:
        DataFrame with GRN mapping
    """
    # Create DataFrame from the mapping
    df = pd.DataFrame(list(grn_mapping.items()), columns=['position', 'grn'])
    
    # Sort by position
    df = df.sort_values('position')
    
    # Add amino acid information if sequence is provided
    if sequence is not None:
        df['residue'] = df['position'].apply(lambda pos: sequence[pos-1] if 1 <= pos <= len(sequence) else 'X')
    
    return df

def df_to_grn_mapping(df: pd.DataFrame) -> Dict[int, str]:
    """
    Convert a DataFrame with GRN mapping to a dictionary.
    
    Args:
        df: DataFrame with columns 'position' and 'grn'
        
    Returns:
        Dictionary mapping sequence positions to GRNs
    """
    # Check that required columns exist
    if 'position' not in df.columns or 'grn' not in df.columns:
        raise ValueError("DataFrame must have 'position' and 'grn' columns")
    
    # Create mapping dictionary
    grn_mapping = dict(zip(df['position'], df['grn']))
    
    return grn_mapping

# -----------------------------------------------------------------------------
# Sequence Data Conversions
# -----------------------------------------------------------------------------

def fasta_to_dict(fasta_content: str) -> Dict[str, str]:
    """
    Convert FASTA content to a dictionary of sequences.
    
    Args:
        fasta_content: Content of a FASTA file
        
    Returns:
        Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    # Parse FASTA content
    for line in fasta_content.strip().split('\n'):
        if line.startswith('>'):
            # Save previous sequence if any
            if current_id is not None:
                sequences[current_id] = ''.join(current_seq)
            
            # Start new sequence
            current_id = line[1:].strip()
            current_seq = []
        elif line.strip():
            # Add to current sequence
            current_seq.append(line.strip())
    
    # Add the last sequence
    if current_id is not None:
        sequences[current_id] = ''.join(current_seq)
    
    return sequences

def dict_to_fasta(sequences: Dict[str, str], width: int = 80) -> str:
    """
    Convert a dictionary of sequences to FASTA content.
    
    Args:
        sequences: Dictionary mapping sequence IDs to sequences
        width: Line width for the FASTA file
        
    Returns:
        FASTA content as a string
    """
    fasta_lines = []
    
    # Generate FASTA content
    for seq_id, sequence in sequences.items():
        # Add header
        fasta_lines.append(f">{seq_id}")
        
        # Add sequence with line wrapping
        for i in range(0, len(sequence), width):
            fasta_lines.append(sequence[i:i+width])
    
    return '\n'.join(fasta_lines)

def alignment_result_to_df(query_id: str,
                         target_id: str,
                         aligned_query: str,
                         aligned_target: str,
                         score: float) -> pd.DataFrame:
    """
    Convert alignment results to a DataFrame.
    
    Args:
        query_id: Query sequence ID
        target_id: Target sequence ID
        aligned_query: Aligned query sequence
        aligned_target: Aligned target sequence
        score: Alignment score
        
    Returns:
        DataFrame with alignment results
    """
    # Calculate alignment statistics
    alignment_length = len(aligned_query)
    identity = sum(1 for q, t in zip(aligned_query, aligned_target) if q == t and q != '-' and t != '-')
    gaps = aligned_query.count('-') + aligned_target.count('-')
    
    if alignment_length > 0:
        sequence_identity = (identity / alignment_length) * 100
    else:
        sequence_identity = 0
    
    # Find start and end positions
    query_start = aligned_query.find(next(c for c in aligned_query if c != '-')) + 1
    query_end = len(aligned_query) - aligned_query[::-1].find(next(c for c in aligned_query[::-1] if c != '-'))
    
    target_start = aligned_target.find(next(c for c in aligned_target if c != '-')) + 1
    target_end = len(aligned_target) - aligned_target[::-1].find(next(c for c in aligned_target[::-1] if c != '-'))
    
    # Create DataFrame
    df = pd.DataFrame({
        'query_id': [query_id],
        'target_id': [target_id],
        'sequence_identity': [sequence_identity],
        'alignment_length': [alignment_length],
        'mismatches': [alignment_length - identity - gaps],
        'gap_openings': [sum(1 for i in range(1, alignment_length) if (aligned_query[i] == '-' and aligned_query[i-1] != '-') or (aligned_target[i] == '-' and aligned_target[i-1] != '-'))],
        'query_start': [query_start],
        'query_end': [query_end],
        'target_start': [target_start],
        'target_end': [target_end],
        'e_value': [0.0],  # Not applicable for simple alignments
        'bit_score': [score]
    })
    
    return df

def df_to_alignment_result(df: pd.DataFrame, include_aligned_sequences: bool = False) -> Dict[str, Any]:
    """
    Convert an alignment DataFrame to a dictionary with alignment results.
    
    Args:
        df: DataFrame with alignment results
        include_aligned_sequences: Whether to include aligned sequences (not available in standard format)
        
    Returns:
        Dictionary with alignment results
    """
    # Check that the DataFrame has at least one row
    if df.empty:
        return {}
    
    # Extract the first row
    row = df.iloc[0]
    
    # Create result dictionary
    result = {
        'query_id': row['query_id'],
        'target_id': row['target_id'],
        'sequence_identity': row['sequence_identity'],
        'alignment_length': row['alignment_length'],
        'query_start': row['query_start'],
        'query_end': row['query_end'],
        'target_start': row['target_start'],
        'target_end': row['target_end'],
        'score': row['bit_score']
    }
    
    # Add aligned sequences if requested (not available in standard format)
    if include_aligned_sequences:
        result['aligned_query'] = ''
        result['aligned_target'] = ''
    
    return result