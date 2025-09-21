"""
Sequence extraction and analysis functions for protein structures.

This module provides functions for extracting sequences from structures,
mapping between structure and sequence numbering, and identifying missing residues.
"""

import pandas as pd
from typing import Dict, List, Optional, Tuple, Set


# Three-letter to one-letter amino acid code mapping
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X', 'MSE': 'M',  # Selenomethionine
    'SEC': 'U',  # Selenocysteine
    'PYL': 'O',  # Pyrrolysine
}

# One-letter to three-letter mapping
AA_1TO3 = {v: k for k, v in AA_3TO1.items() if k not in ['UNK', 'MSE']}


def extract_sequence(df: pd.DataFrame, 
                    chain_id: Optional[str] = None,
                    one_letter: bool = True,
                    include_gaps: bool = True) -> str:
    """
    Extract amino acid sequence from structure.
    
    Args:
        df: Structure DataFrame
        chain_id: Chain identifier (if None, extracts first chain)
        one_letter: If True, returns one-letter code, else three-letter
        include_gaps: If True, includes '-' for missing residues
        
    Returns:
        Amino acid sequence string
    """
    # Reset index to access all columns
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Filter by chain if specified
    if chain_id is not None:
        df_chain = df_reset[df_reset['auth_chain_id'] == chain_id]
    else:
        # Use first chain
        chain_id = df_reset['auth_chain_id'].iloc[0]
        df_chain = df_reset[df_reset['auth_chain_id'] == chain_id]
    
    if df_chain.empty:
        return ""
    
    # Get CA atoms only (one per residue)
    ca_atoms = df_chain[df_chain['atom_name'] == 'CA'].copy()
    
    if ca_atoms.empty:
        # Fallback to any backbone atom
        ca_atoms = df_chain[df_chain['atom_name'].isin(['N', 'C', 'O'])].copy()
        if not ca_atoms.empty:
            # Remove duplicates, keep first atom per residue
            ca_atoms = ca_atoms.drop_duplicates(subset=['auth_seq_id'], keep='first')
    
    if ca_atoms.empty:
        return ""
    
    # Sort by residue number
    ca_atoms = ca_atoms.sort_values('auth_seq_id')
    
    # Get residue names and numbers
    residues = ca_atoms[['auth_seq_id', 'res_name3l']].values
    
    if include_gaps:
        # Build sequence with gaps
        sequence_parts = []
        prev_resid = None
        
        for resid, resname in residues:
            # Add gaps if residues are missing
            if prev_resid is not None and resid - prev_resid > 1:
                n_gaps = resid - prev_resid - 1
                if one_letter:
                    sequence_parts.append('-' * n_gaps)
                else:
                    sequence_parts.extend(['---'] * n_gaps)
            
            # Add residue
            if one_letter:
                aa_code = AA_3TO1.get(resname.upper(), 'X')
                sequence_parts.append(aa_code)
            else:
                sequence_parts.append(resname)
            
            prev_resid = resid
        
        if one_letter:
            return ''.join(sequence_parts)
        else:
            return ' '.join(sequence_parts)
    else:
        # No gaps
        if one_letter:
            return ''.join(AA_3TO1.get(resname.upper(), 'X') for _, resname in residues)
        else:
            return ' '.join(resname for _, resname in residues)


def extract_all_sequences(df: pd.DataFrame, 
                         one_letter: bool = True,
                         min_length: int = 10) -> Dict[str, str]:
    """
    Extract sequences for all chains.
    
    Args:
        df: Structure DataFrame
        one_letter: If True, returns one-letter codes
        min_length: Minimum sequence length to include
        
    Returns:
        Dictionary mapping chain_id to sequence
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    sequences = {}
    
    # Get unique chains
    chains = df_reset['auth_chain_id'].unique()
    
    for chain_id in chains:
        seq = extract_sequence(df, chain_id, one_letter=one_letter, include_gaps=False)
        
        if len(seq) >= min_length:
            sequences[chain_id] = seq
    
    return sequences


def map_structure_to_sequence(df: pd.DataFrame, 
                            sequence: str,
                            chain_id: Optional[str] = None) -> Dict[int, int]:
    """
    Map structure residue numbers to sequence positions.
    
    Handles cases where structure numbering doesn't match sequence positions
    due to missing residues, insertions, or non-standard numbering.
    
    Args:
        df: Structure DataFrame
        sequence: Target sequence (one-letter code)
        chain_id: Chain to map (if None, uses first chain)
        
    Returns:
        Dictionary mapping structure residue number to sequence position (0-based)
    """
    # Extract structure sequence
    struct_seq = extract_sequence(df, chain_id, one_letter=True, include_gaps=False)
    
    if not struct_seq:
        return {}
    
    # Simple sequence alignment to find best mapping
    # This is a basic implementation - could use more sophisticated alignment
    best_offset = 0
    best_matches = 0
    
    # Try different offsets to find best alignment
    for offset in range(max(0, -len(struct_seq)), len(sequence)):
        matches = 0
        for i, aa in enumerate(struct_seq):
            seq_pos = i + offset
            if 0 <= seq_pos < len(sequence) and sequence[seq_pos] == aa:
                matches += 1
        
        if matches > best_matches:
            best_matches = matches
            best_offset = offset
    
    # Build mapping using best offset
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    if chain_id is not None:
        df_chain = df_reset[df_reset['auth_chain_id'] == chain_id]
    else:
        chain_id = df_reset['auth_chain_id'].iloc[0]
        df_chain = df_reset[df_reset['auth_chain_id'] == chain_id]
    
    # Get CA atoms
    ca_atoms = df_chain[df_chain['atom_name'] == 'CA'].sort_values('auth_seq_id')
    
    mapping = {}
    for i, (_, atom) in enumerate(ca_atoms.iterrows()):
        struct_resid = atom['auth_seq_id']
        seq_pos = i + best_offset
        
        if 0 <= seq_pos < len(sequence):
            resname = atom['res_name3l']
            expected_aa = AA_3TO1.get(resname.upper(), 'X')
            
            # Only map if amino acids match
            if sequence[seq_pos] == expected_aa:
                mapping[struct_resid] = seq_pos
    
    return mapping


def identify_missing_residues(df: pd.DataFrame, 
                            expected_sequence: str,
                            chain_id: Optional[str] = None) -> List[Tuple[int, str]]:
    """
    Compare structure to expected sequence and identify missing residues.
    
    Args:
        df: Structure DataFrame
        expected_sequence: Expected full sequence (one-letter code)
        chain_id: Chain to analyze
        
    Returns:
        List of (position, residue) tuples for missing residues
    """
    # Get structure to sequence mapping
    mapping = map_structure_to_sequence(df, expected_sequence, chain_id)
    
    if not mapping:
        # No mapping found, all residues are missing
        return [(i, aa) for i, aa in enumerate(expected_sequence)]
    
    # Find missing positions
    mapped_positions = set(mapping.values())
    missing = []
    
    for pos, aa in enumerate(expected_sequence):
        if pos not in mapped_positions:
            missing.append((pos, aa))
    
    return missing


def get_sequence_segments(df: pd.DataFrame, 
                         chain_id: Optional[str] = None) -> List[Tuple[int, int, str]]:
    """
    Identify continuous sequence segments in structure.
    
    Useful for finding structured regions vs loops/missing regions.
    
    Args:
        df: Structure DataFrame
        chain_id: Chain to analyze
        
    Returns:
        List of (start_resid, end_resid, sequence) tuples
    """
    # Reset index
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    
    # Filter by chain
    if chain_id is not None:
        df_chain = df_reset[df_reset['auth_chain_id'] == chain_id]
    else:
        chain_id = df_reset['auth_chain_id'].iloc[0]
        df_chain = df_reset[df_reset['auth_chain_id'] == chain_id]
    
    # Get CA atoms
    ca_atoms = df_chain[df_chain['atom_name'] == 'CA'].sort_values('auth_seq_id')
    
    if ca_atoms.empty:
        return []
    
    segments = []
    current_start = None
    current_seq = []
    prev_resid = None
    
    for _, atom in ca_atoms.iterrows():
        resid = atom['auth_seq_id']
        resname = atom['res_name3l']
        aa = AA_3TO1.get(resname.upper(), 'X')
        
        if current_start is None:
            # Start new segment
            current_start = resid
            current_seq = [aa]
            prev_resid = resid
        elif resid - prev_resid == 1:
            # Continue segment
            current_seq.append(aa)
            prev_resid = resid
        else:
            # Gap found, save current segment and start new
            segments.append((current_start, prev_resid, ''.join(current_seq)))
            current_start = resid
            current_seq = [aa]
            prev_resid = resid
    
    # Add final segment
    if current_start is not None:
        segments.append((current_start, prev_resid, ''.join(current_seq)))
    
    return segments


def compare_chain_sequences(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compare sequences between chains to identify similar/identical chains.
    
    Args:
        df: Structure DataFrame
        
    Returns:
        DataFrame with pairwise sequence comparisons
    """
    sequences = extract_all_sequences(df, one_letter=True, min_length=5)
    
    if len(sequences) < 2:
        return pd.DataFrame()
    
    comparisons = []
    chains = sorted(sequences.keys())
    
    for i, chain1 in enumerate(chains):
        for j in range(i + 1, len(chains)):
            chain2 = chains[j]
            
            seq1 = sequences[chain1]
            seq2 = sequences[chain2]
            
            # Calculate sequence identity
            if len(seq1) == len(seq2):
                # Same length - direct comparison
                matches = sum(a == b for a, b in zip(seq1, seq2))
                identity = matches / len(seq1)
                aligned_length = len(seq1)
            else:
                # Different lengths - find best alignment
                identity, aligned_length = calculate_sequence_identity(seq1, seq2)
            
            comparisons.append({
                'chain1': chain1,
                'chain2': chain2,
                'len1': len(seq1),
                'len2': len(seq2),
                'identity': identity,
                'aligned_length': aligned_length,
                'is_identical': identity > 0.99
            })
    
    return pd.DataFrame(comparisons)


def calculate_sequence_identity(seq1: str, seq2: str) -> Tuple[float, int]:
    """
    Calculate sequence identity between two sequences.
    
    Simple implementation - for production use, would use proper alignment algorithm.
    
    Args:
        seq1: First sequence
        seq2: Second sequence
        
    Returns:
        Tuple of (identity fraction, aligned length)
    """
    if not seq1 or not seq2:
        return 0.0, 0
    
    # Try different alignments and find best
    best_identity = 0.0
    best_length = 0
    
    # Slide shorter sequence along longer
    if len(seq1) <= len(seq2):
        short, long = seq1, seq2
    else:
        short, long = seq2, seq1
    
    for offset in range(len(long) - len(short) + 1):
        matches = 0
        for i, aa in enumerate(short):
            if long[offset + i] == aa:
                matches += 1
        
        identity = matches / len(short)
        if identity > best_identity:
            best_identity = identity
            best_length = len(short)
    
    return best_identity, best_length


def annotate_sequence_conservation(df: pd.DataFrame, 
                                 alignment: Dict[str, str]) -> pd.Series:
    """
    Annotate residues with sequence conservation scores.
    
    Args:
        df: Structure DataFrame
        alignment: Dictionary mapping chain_id to aligned sequences
                  (must have gaps '-' to maintain alignment)
        
    Returns:
        Series with (chain_id, residue_id) -> conservation score
    """
    if len(alignment) < 2:
        # Need at least 2 sequences for conservation
        return pd.Series(dtype=float)
    
    # Calculate conservation at each position
    sequences = list(alignment.values())
    seq_length = len(sequences[0])
    
    # Ensure all sequences have same length
    if not all(len(s) == seq_length for s in sequences):
        raise ValueError("All sequences must be aligned (same length)")
    
    # Calculate position-wise conservation
    conservation_scores = []
    for pos in range(seq_length):
        # Get amino acids at this position (excluding gaps)
        aas = [s[pos] for s in sequences if s[pos] != '-']
        
        if not aas:
            conservation_scores.append(0.0)
        else:
            # Calculate frequency of most common amino acid
            aa_counts = pd.Series(aas).value_counts()
            max_count = aa_counts.iloc[0]
            conservation = max_count / len(aas)
            conservation_scores.append(conservation)
    
    # Map conservation to structure residues
    df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df
    residue_conservation = {}
    
    for chain_id, aligned_seq in alignment.items():
        # Get structure sequence for mapping
        chain_df = df_reset[df_reset['auth_chain_id'] == chain_id]
        ca_atoms = chain_df[chain_df['atom_name'] == 'CA'].sort_values('auth_seq_id')
        
        # Map structure positions to alignment positions
        struct_pos = 0
        for align_pos, aa in enumerate(aligned_seq):
            if aa != '-':
                if struct_pos < len(ca_atoms):
                    resid = ca_atoms.iloc[struct_pos]['auth_seq_id']
                    residue_conservation[(chain_id, resid)] = conservation_scores[align_pos]
                    struct_pos += 1
    
    return pd.Series(residue_conservation)