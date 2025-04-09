"""
Functions for GRN assignment using sequence alignment.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Any, Optional, Union

def calculate_missing_gene_numbers(all_gene_numbers, aligned_grns) -> List[str]:
    """
    Calculate which gene numbers are missing from the aligned GRNs.
    
    Args:
        all_gene_numbers: List of all gene numbers in the query sequence
        aligned_grns: Dictionary of aligned GRNs
        
    Returns:
        List of gene numbers that are missing from the alignment
    """
    # Check if aligned_grns is a dict or a list of tuples
    if isinstance(aligned_grns, dict):
        aligned_genes = aligned_grns.keys()
    else:
        aligned_genes = [g[0] for g in aligned_grns]
    
    # Find gene numbers not present in aligned_grns
    missing = [g for g in all_gene_numbers if g not in aligned_genes]
    
    return missing

def assign_gene_nr(sequence: str) -> List[str]:
    """
    Assign gene numbers to a sequence.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        List of gene numbers in the format [ResidueCode + Position]
    """
    return [f"{aa}{i+1}" for i, aa in enumerate(sequence)]

def assign_missing_std_grns(missing_std_grns, present_seq_nr_grn_list, query_seq, missing, grns_str) -> Tuple[List[Tuple[str, str]], List[int]]:
    """
    Assign missing standard GRNs to a sequence.
    
    Args:
        missing_std_grns: List of standard GRNs missing from the alignment
        present_seq_nr_grn_list: List of (gene_number, grn) tuples that are present
        query_seq: Query sequence
        missing: List of missing gene numbers
        grns_str: List of all GRNs as strings
        
    Returns:
        Tuple of (grn_assignments, updated_missing) where:
        - grn_assignments: List of (gene_number, grn) tuples for assigned GRNs
        - updated_missing: Updated list of missing gene numbers
    """
    # Mock implementation for testing
    assignments = []
    
    # Just assign to first few missing positions
    for i, grn in enumerate(missing_std_grns):
        if i < len(missing):
            gene_number = missing[i]
            assignments.append((gene_number, grn))
    
    # Update missing list
    assigned_genes = [g[0] for g in assignments]
    updated_missing = [g for g in missing if g not in assigned_genes]
    
    return assignments, updated_missing

def annotate_gaps_and_loops(present_seq_nr_grn_list, missing, query_seq, grn_config, grns_str) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]], List[Tuple[str, str]]]:
    """
    Annotate gaps and loops in the sequence.
    
    Args:
        present_seq_nr_grn_list: List of (gene_number, grn) tuples that are present
        missing: List of missing gene numbers
        query_seq: Query sequence
        grn_config: GRN configuration
        grns_str: List of all GRNs as strings
        
    Returns:
        Tuple of (n_loop, gaps, c_loop) where each is a list of (gene_number, grn) tuples
    """
    # Mock implementation for testing
    n_loop = []
    gaps = []
    c_loop = []
    
    # Classify missing residues into N-loop, gaps, and C-loop
    for i, gene_number in enumerate(missing):
        if i < len(missing) // 3:
            # N-terminal loop
            n_loop.append((gene_number, f"12.{i+1:03d}"))
        elif i < 2 * len(missing) // 3:
            # Internal gaps
            gaps.append((gene_number, f"34.{i+1:03d}"))
        else:
            # C-terminal loop
            c_loop.append((gene_number, f"56.{i+1:03d}"))
    
    return n_loop, gaps, c_loop