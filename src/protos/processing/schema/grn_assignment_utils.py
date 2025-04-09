"""
Utility functions for GRN assignment with improved loop region handling.

This module provides functions for assigning GRNs to protein sequences,
with a specific focus on correctly handling loop regions according to the
standardized format: <closer helix><further helix>.<distance>
"""

import re
import logging
from typing import Dict, List, Tuple, Union, Optional, Any, Set

from protos.processing.schema.grn_utils_updated import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string
)


def assign_gene_nr(seq: str) -> List[str]:
    """
    Assign gene numbers to the full sequence.
    
    Args:
        seq: Amino acid sequence
        
    Returns:
        List of amino acids with position numbers (e.g., ['M1', 'E2', 'L3', ...])
    """
    return [seq[x] + str(x + 1) for x in range(len(seq))]


def get_closest_present_seqnr(missing_seqnr: int, 
                            present_seq_nr_grn_list: List[Tuple[str, str]], 
                            loop_side: str = 'n') -> Tuple[Tuple[str, str], int]:
    """
    Given a list of sequence numbers and GRNs pairs, return the closest pair for a missing seqnr.
    
    Args:
        missing_seqnr: Missing sequence number
        present_seq_nr_grn_list: List of (seq_nr, grn) pairs
        loop_side: Which side of the loop to search ('n' for N-terminal, 'c' for C-terminal)
        
    Returns:
        Tuple of (closest_pair, min_distance)
    """
    min_dist = 1000
    closest = None
    
    for g in present_seq_nr_grn_list:
        seq_nr = int(g[0][1:])  # Extract number from seq_id (e.g., 'A15' -> 15)
        dist = seq_nr - missing_seqnr

        # Skip if wrong direction
        if loop_side == 'n' and dist >= 0:
            continue
        if loop_side == 'c' and dist <= 0:
            continue

        abs_dist = abs(dist)
        if abs_dist < abs(min_dist):
            min_dist = dist
            closest = g
            
    return closest, min_dist


def annotate_loop_region(interval: List[int], 
                       present_seq_nr_grn_list: List[Tuple[str, str]], 
                       query_seq: str) -> List[Tuple[str, str]]:
    """
    Annotate residues in a loop region between two helices.
    
    Args:
        interval: List of sequence numbers in the loop
        present_seq_nr_grn_list: List of (seq_nr, grn) pairs for already annotated residues
        query_seq: Amino acid sequence
        
    Returns:
        List of (seq_id, grn) pairs for the loop region
    """
    # Get the helices bordering this loop
    known_grns = [g[1] for g in present_seq_nr_grn_list 
                if 'x' in g[1] and len(g[1].split('x')[0]) == 1]
    
    if not known_grns:
        return []  # No helices found to anchor the loop
        
    # Find helix numbers by looking at standard GRNs (e.g., '1x50')
    helix_numbers = sorted(set([int(g.split('x')[0]) for g in known_grns]))
    
    # If there are not at least 2 helices, we can't annotate a loop
    if len(helix_numbers) < 2:
        return []
        
    # Find the two helices that this interval is between
    # Get sequence numbers for known helices
    helix_seqnrs = {}
    for helix in helix_numbers:
        helix_grns = [g for g in present_seq_nr_grn_list 
                    if 'x' in g[1] and g[1].split('x')[0] == str(helix)]
        if helix_grns:
            helix_seqnrs[helix] = [int(g[0][1:]) for g in helix_grns]
    
    # Find the helices this interval is between
    interval_min = min(interval)
    interval_max = max(interval)
    
    helix_before = None
    helix_after = None
    
    for helix, seqnrs in helix_seqnrs.items():
        if max(seqnrs) < interval_min:
            if helix_before is None or max(helix_seqnrs[helix_before]) < max(seqnrs):
                helix_before = helix
        if min(seqnrs) > interval_max:
            if helix_after is None or min(helix_seqnrs[helix_after]) > min(seqnrs):
                helix_after = helix
    
    # If we can't determine both bordering helices, return empty list
    if helix_before is None or helix_after is None:
        return []
        
    # Now annotate the loop residues
    loop_grns = []
    
    # For each residue in the interval
    for seqnr in interval:
        # Find distance to the closest helix
        dist_before = min(abs(seqnr - nr) for nr in helix_seqnrs[helix_before])
        dist_after = min(abs(seqnr - nr) for nr in helix_seqnrs[helix_after])
        
        # Determine which helix is closer
        if dist_before <= dist_after:
            closer_helix = helix_before
            further_helix = helix_after
            distance = dist_before
        else:
            closer_helix = helix_after
            further_helix = helix_before
            distance = dist_after
        
        # Format the loop GRN with the correct format: <closer helix><further helix>.<distance>
        grn = f"{closer_helix}{further_helix}.{distance:03d}"
        
        # Create the seq_id
        seq_id = query_seq[seqnr - 1] + str(seqnr)
        
        # Add to the list
        loop_grns.append((seq_id, grn))
    
    return loop_grns


def calculate_missing_ntail_grns(aligned_grns: Dict[str, str], 
                               missing_gene_numbers: List[str], 
                               grns_float: List[float]) -> Tuple[List[Tuple[str, str]], int]:
    """
    Calculate GRNs for missing N-terminal residues.
    
    Args:
        aligned_grns: Dictionary mapping sequence positions to GRNs
        missing_gene_numbers: List of missing gene numbers
        grns_float: List of float representations of standard GRNs
        
    Returns:
        Tuple of (n_tail_list, first_gene_number_int)
    """
    beginning_tm1_float = grns_float[0]
    first_grn = list(aligned_grns.values())[0]
    first_grn_float = parse_grn_str2float(first_grn)
    first_gene_number = list(aligned_grns.keys())[0]

    first_gene_number_int = int(first_gene_number[1:]) - 1

    missing_grns_tm1_ = int(100 * (first_grn_float - beginning_tm1_float))
    missing_grns_tm1 = min(first_gene_number_int, missing_grns_tm1_)
    missing_tm1 = min(missing_grns_tm1, missing_grns_tm1_)
    missing_ntail = max(0, first_gene_number_int - missing_tm1)

    n_tail_float = grns_float[:grns_float.index(first_grn_float)]
    n_tail_float += [-(i + 1) for i in range(missing_ntail)]
    n_tail_float = sorted(n_tail_float)

    n_tail_str = [parse_grn_float2str(x) for x in n_tail_float]
    n_tail_list = list(zip(missing_gene_numbers[:first_gene_number_int][::-1], n_tail_str[::-1]))[::-1]
    return n_tail_list, first_gene_number_int


def calculate_missing_ctail_grns(aligned_grns: Dict[str, str], 
                               missing_gene_numbers: List[str], 
                               query_gene_len: int, 
                               grns_float: List[float]) -> Tuple[List[Tuple[str, str]], Optional[int]]:
    """
    Calculate GRNs for missing C-terminal residues.
    
    Args:
        aligned_grns: Dictionary mapping sequence positions to GRNs
        missing_gene_numbers: List of missing gene numbers
        query_gene_len: Length of the query gene
        grns_float: List of float representations of standard GRNs
        
    Returns:
        Tuple of (c_tail_list, last_gene_number_int)
    """
    # Check if we have any aligned GRNs
    if not aligned_grns:
        return [], None

    # Extract the float value for the end of the last section from the standardized GRNs
    ending_tmx_float = grns_float[-1]

    # Get information about the last known GRN in the alignment
    last_grn = list(aligned_grns.values())[-1]
    last_grn_float = parse_grn_str2float(last_grn)
    last_gene_number = list(aligned_grns.keys())[-1]
    last_gene_number_int = int(last_gene_number[1:])
    missing_grns_last_section_ = int(100 * (ending_tmx_float - last_grn_float))
    missing_grns_last_section = min(query_gene_len - last_gene_number_int, missing_grns_last_section_)
    missing_last_section = max(0, query_gene_len - last_gene_number_int - missing_grns_last_section)

    c_tail_float = grns_float[grns_float.index(last_grn_float) + 1:]
    c_tail_float += [100 + i + 1 for i in range(missing_last_section)]
    c_tail_float = sorted(c_tail_float)

    c_tail_str = [parse_grn_float2str(x) for x in c_tail_float]
    c_tail_list = list(zip(missing_gene_numbers[-(query_gene_len - last_gene_number_int):], c_tail_str))

    return c_tail_list, last_gene_number_int


def valid_jump(prev_ref_grn: Optional[str], 
              curr_ref_grn: str, 
              prev_query_key: Optional[str], 
              curr_query_key: str, 
              max_alignment_gap: int = 1) -> bool:
    """
    Check if a jump from one aligned position to another is valid.
    
    Args:
        prev_ref_grn: Previous reference GRN
        curr_ref_grn: Current reference GRN
        prev_query_key: Previous query key
        curr_query_key: Current query key
        max_alignment_gap: Maximum allowed alignment gap
        
    Returns:
        True if the jump is valid, False otherwise
    """
    # If there is no previous reference GRN, return True
    if prev_ref_grn is None:
        return True
    
    # For standard GRNs, check TM region
    if 'x' in prev_ref_grn and 'x' in curr_ref_grn:
        prev_grn_tm = int(prev_ref_grn.split('x')[0])
        curr_grn_tm = int(curr_ref_grn.split('x')[0])

        # If the TM numbers are not equal, return True
        if prev_grn_tm != curr_grn_tm:
            return True

        # Parse the GRNs to float values
        prev_grn_float = parse_grn_str2float(prev_ref_grn)
        curr_grn_float = parse_grn_str2float(curr_ref_grn)

        # Extract the residue numbers from the query keys
        prev_query_num = int(prev_query_key[1:])
        curr_query_num = int(curr_query_key[1:])

        # Calculate the differences between GRNs and residue numbers
        grn_diff = abs(curr_grn_float - prev_grn_float)
        rn_diff = abs(curr_query_num - prev_query_num)

        # Check if the differences meet the criteria, return True if they do
        if (grn_diff <= .1) and (rn_diff <= max_alignment_gap):
            return True
    
    # For other GRN types (N-terminal, C-terminal, loop)
    # We could add more specific rules here
    
    return False


def get_correctly_aligned_grns(all_query_gene_numbers: List[str], 
                              reference_grn_dict: Dict[str, str], 
                              alignment: Tuple[str, str, str], 
                              max_alignment_gap: int = 1) -> Dict[str, str]:
    """
    Get correctly aligned GRNs from an alignment.
    
    Args:
        all_query_gene_numbers: List of all query gene numbers
        reference_grn_dict: Dictionary mapping reference residues to GRNs
        alignment: Tuple of (query_seq, match_line, ref_seq)
        max_alignment_gap: Maximum allowed alignment gap
        
    Returns:
        Dictionary mapping query keys to GRNs
    """
    query_seq, match_line, ref_seq = alignment
    pointer_B = 0
    pointer_A = 0
    ref_keys = list(reference_grn_dict.keys())
    result_dict = {}

    # Initialize previous values for reference and query keys and pairs
    prev_pair = None
    prev_query_key = None

    # Iterate through the match_line
    for i in range(len(match_line)):
        if (pointer_B < len(all_query_gene_numbers)) and (pointer_A < len(ref_keys)):
            query_key = all_query_gene_numbers[pointer_B]
            ref_key = ref_keys[pointer_A]

            # Check if the current position in the match_line is an alignment match or mismatch
            if match_line[i] == '|' or match_line[i] == '.':
                if ref_seq[i] != '-' and query_seq[i] != '-':
                    curr_pair = reference_grn_dict[ref_key]
                    
                    # Normalize the GRN to ensure consistent format
                    curr_pair = normalize_grn_format(curr_pair)

                    if valid_jump(prev_pair, curr_pair, prev_query_key, query_key, max_alignment_gap):
                        result_dict[query_key] = curr_pair

                        # Update the previous values for the next iteration
                        prev_pair = curr_pair
                        prev_query_key = query_key

            # Update the pointers based on the gaps in the sequences
            if query_seq[i] != '-':
                pointer_B += 1
            if ref_seq[i] != '-':
                pointer_A += 1
                
    return result_dict