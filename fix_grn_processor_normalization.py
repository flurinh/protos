"""
Patch for GRN processor to fix normalization issues in expand_annotation.

This module provides fixed versions of key functions that ensure:
1. Dense assignment in standard regions (no gaps)
2. All residues are assigned (no missing)
3. Proper ordering (both GRN and residue numbers)
4. No duplicate assignments
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from protos.processing.grn.grn_utils import (
    parse_grn_str2float, parse_grn_float2str, sort_grns_str, get_grn_interval
)
from protos.processing.grn.grn_assignment import (
    assign_gene_nr, calculate_missing_gene_numbers, sort_grn_rn_pairs
)


def expand_annotation_fixed(new_row, query_seq, alignment, max_alignment_gap=1, 
                           protein_family: str = 'gpcr_a', verbose=0):
    """
    Fixed version of expand_annotation that ensures proper normalization.
    
    Key fixes:
    1. Ensures all residues 1 to len(query_seq) are assigned
    2. Compresses gaps in standard regions  
    3. Maintains proper ordering
    4. Prevents duplicates
    """
    from protos.processing.grn.grn_utils import GRNConfigManager
    
    # Create reference_grn_dict from the input row
    aligned_grns = {v: k for (k, v) in new_row.to_dict().items() if (v != '-') & ('.' in k)}
    
    if verbose >= 1:
        print(f"\n=== Fixed Expand Annotation Started ===")
        print(f"Query length: {len(query_seq)}, Initial GRNs: {len(aligned_grns)}")
    
    # Calculate query length and gene numbers
    query_gene_len = len(query_seq)
    
    if isinstance(query_seq, str):
        all_query_gene_numbers = assign_gene_nr(query_seq)
    elif isinstance(query_seq, list):
        all_query_gene_numbers = query_seq
    else:
        raise ValueError("Query sequence must be string or list")

    # Initialize GRN configuration
    config = GRNConfigManager()
    grn_config_std = config.get_config(protein_family=protein_family, strict=False)
    
    # Generate all standard GRNs
    grns_str_std = []
    if grn_config_std:
        for region_name, (start_grn, end_grn) in grn_config_std.items():
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_std.extend(region_grns)
    
    grns_str_std = list(set(grns_str_std))
    grns_str_std = sort_grns_str(grns_str_std)
    
    # Filter to single TM GRNs
    grns_str_std_filtered = []
    for grn in grns_str_std:
        if '.' in grn:
            tm_part = grn.split('.')[0]
            if len(tm_part) == 1 and tm_part.isdigit():
                grns_str_std_filtered.append(grn)

    grns_str_std = grns_str_std_filtered
    grns_float_std = [parse_grn_str2float(x) for x in grns_str_std]
    
    # Build complete assignment mapping
    grn_assignments = {}  # residue_number -> GRN
    
    # Step 1: Add initial aligned positions
    for gene_num, grn in aligned_grns.items():
        rn_num = int(gene_num[1:])
        grn_assignments[rn_num] = grn
    
    if verbose >= 1:
        print(f"\nStep 1: Initial assignments: {len(grn_assignments)}")
    
    # Step 2: Process N-terminal
    if grn_assignments:
        first_assigned = min(grn_assignments.keys())
        first_grn = grn_assignments[first_assigned]
        first_grn_float = parse_grn_str2float(first_grn)
        
        # Find the start of TM1 in our standard GRNs
        tm1_start = None
        for grn in grns_str_std:
            if grn.startswith('1.'):
                tm1_start = parse_grn_str2float(grn)
                break
        
        if tm1_start and first_grn_float >= tm1_start:
            # Assign N-terminal positions
            n_term_positions = first_assigned - 1
            for i in range(n_term_positions):
                rn_num = i + 1
                grn_assignments[rn_num] = f"n.{n_term_positions - i}"
                
            if verbose >= 1:
                print(f"Step 2: Added {n_term_positions} N-terminal positions")
    
    # Step 3: Process C-terminal
    if grn_assignments:
        last_assigned = max(grn_assignments.keys())
        last_grn = grn_assignments[last_assigned]
        
        # Assign C-terminal positions
        c_term_positions = query_gene_len - last_assigned
        for i in range(c_term_positions):
            rn_num = last_assigned + i + 1
            grn_assignments[rn_num] = f"c.{i + 1}"
            
        if verbose >= 1:
            print(f"Step 3: Added {c_term_positions} C-terminal positions")
    
    # Step 4: Fill internal gaps
    all_positions = set(range(1, query_gene_len + 1))
    assigned_positions = set(grn_assignments.keys())
    missing_positions = sorted(all_positions - assigned_positions)
    
    if verbose >= 1:
        print(f"\nStep 4: Filling {len(missing_positions)} internal gaps")
    
    for missing_pos in missing_positions:
        # Find nearest neighbors
        prev_pos = None
        next_pos = None
        
        for pos in sorted(grn_assignments.keys()):
            if pos < missing_pos:
                prev_pos = pos
            elif pos > missing_pos and next_pos is None:
                next_pos = pos
                break
        
        if prev_pos and next_pos:
            prev_grn = grn_assignments[prev_pos]
            next_grn = grn_assignments[next_pos]
            
            # Determine appropriate assignment
            new_grn = _interpolate_grn(prev_grn, next_grn, prev_pos, next_pos, missing_pos, grn_config_std)
            grn_assignments[missing_pos] = new_grn
    
    # Step 5: Compress gaps in standard regions
    if verbose >= 1:
        print(f"\nStep 5: Compressing gaps in standard regions")
    
    grn_assignments = _compress_standard_region_gaps(grn_assignments, grn_config_std, grns_str_std, verbose)
    
    # Step 6: Convert to output format
    grn_list = []
    rn_list = []
    
    for rn_num in sorted(grn_assignments.keys()):
        grn = grn_assignments[rn_num]
        rn = query_seq[rn_num - 1] + str(rn_num)
        grn_list.append(grn)
        rn_list.append(rn)
    
    # Verify no missing positions
    missing = [i for i in range(1, query_gene_len + 1) if i not in grn_assignments]
    
    if verbose >= 1:
        if missing:
            print(f"\n   WARNING: {len(missing)} positions still missing!")
        else:
            print(f"\n   SUCCESS: All {query_gene_len} positions annotated")
        print(f"\n=== Fixed Expand Annotation Completed ===\n")
    
    return grn_list, rn_list, missing


def _interpolate_grn(prev_grn: str, next_grn: str, 
                    prev_pos: int, next_pos: int, missing_pos: int,
                    grn_config: dict) -> str:
    """Interpolate GRN for missing position between two assigned positions."""
    
    prev_float = parse_grn_str2float(prev_grn)
    next_float = parse_grn_str2float(next_grn)
    
    # Check if in same TM region
    prev_tm = _get_tm_number(prev_grn)
    next_tm = _get_tm_number(next_grn)
    
    if prev_tm and next_tm and prev_tm == next_tm:
        # Within same TM - use standard position
        # Linear interpolation
        position_ratio = (missing_pos - prev_pos) / (next_pos - prev_pos)
        interpolated = prev_float + position_ratio * (next_float - prev_float)
        return parse_grn_float2str(round(interpolated, 2))
    
    elif prev_tm and next_tm and prev_tm != next_tm:
        # Between different TMs - assign to loop
        distance_from_prev = missing_pos - prev_pos
        if prev_tm < next_tm:
            loop_id = f"{prev_tm}{next_tm}"
        else:
            loop_id = f"{next_tm}{prev_tm}"
        return f"{loop_id}.{distance_from_prev:03d}"
    
    else:
        # One or both are not standard TMs
        # Simple interpolation
        position_ratio = (missing_pos - prev_pos) / (next_pos - prev_pos)
        interpolated = prev_float + position_ratio * (next_float - prev_float)
        return parse_grn_float2str(interpolated)


def _get_tm_number(grn: str) -> Optional[int]:
    """Extract TM number from a GRN string."""
    if not grn or grn == '-':
        return None
    
    if grn.startswith('n.') or grn.startswith('c.'):
        return None
    
    if '.' in grn:
        parts = grn.split('.')
        if len(parts[0]) == 1 and parts[0].isdigit():
            return int(parts[0])
    elif 'x' in grn:
        parts = grn.split('x')
        if parts[0].isdigit():
            return int(parts[0])
    
    return None


def _compress_standard_region_gaps(grn_assignments: Dict[int, str], 
                                  grn_config: dict,
                                  grns_str_std: List[str],
                                  verbose: int) -> Dict[int, str]:
    """Compress gaps within standard TM regions to ensure dense assignment."""
    
    # Process each TM region
    for region_name, (start_grn, end_grn) in grn_config.items():
        # Handle both 'TM1' and 'tm1' formats
        if not (region_name.upper().startswith('TM') and region_name[2:].isdigit()):
            continue
            
        start_float = parse_grn_str2float(start_grn)
        end_float = parse_grn_str2float(end_grn)
        tm_num = _get_tm_number(start_grn)
        
        if not tm_num:
            continue
        
        # Find all positions assigned to this TM
        tm_positions = []
        for pos, grn in grn_assignments.items():
            grn_float = parse_grn_str2float(grn)
            if start_float <= grn_float <= end_float and _get_tm_number(grn) == tm_num:
                tm_positions.append(pos)
        
        if not tm_positions:
            continue
            
        tm_positions.sort()
        
        # Check if positions are consecutive
        is_consecutive = all(tm_positions[i+1] - tm_positions[i] == 1 
                           for i in range(len(tm_positions) - 1))
        
        if not is_consecutive:
            if verbose >= 2:
                print(f"  Compressing gaps in {region_name}")
            
            # Get available GRN positions in this region
            available_grns = [g for g in grns_str_std 
                            if start_float <= parse_grn_str2float(g) <= end_float]
            
            # Reassign consecutively
            for i, pos in enumerate(tm_positions):
                if i < len(available_grns):
                    old_grn = grn_assignments[pos]
                    new_grn = available_grns[i]
                    if old_grn != new_grn and verbose >= 3:
                        print(f"    Position {pos}: {old_grn} -> {new_grn}")
                    grn_assignments[pos] = new_grn
    
    return grn_assignments


# Monkey patch the original function if needed
def patch_grn_processor():
    """Apply the fixed expand_annotation to GRN processor."""
    import protos.processing.grn.grn_table_utils
    protos.processing.grn.grn_table_utils.expand_annotation = expand_annotation_fixed
    print("GRN processor patched with normalization fixes")


if __name__ == "__main__":
    # Test the fixed function
    print("Testing fixed expand_annotation...")
    
    # Create test data
    test_row = pd.Series({
        '1.28': 'T10',
        '1.35': 'V15', 
        '3.28': 'L25',
        '7.47': 'Y90'
    })
    
    test_seq = 'A' * 100  # 100 residue sequence
    test_alignment = (test_seq, '|' * 100, 'B' * 100)  # Dummy alignment
    
    grn_list, rn_list, missing = expand_annotation_fixed(
        test_row, test_seq, test_alignment, 
        protein_family='gpcr_a', verbose=2
    )
    
    print(f"\nResults:")
    print(f"Assigned: {len(grn_list)} positions")
    print(f"Missing: {len(missing)} positions")
    
    # Check for issues
    rn_nums = [int(rn[1:]) for rn in rn_list]
    print(f"\nValidation:")
    print(f"Sequential: {rn_nums == sorted(rn_nums)}")
    print(f"No gaps: {max(rn_nums) - min(rn_nums) + 1 == len(rn_nums)}")
    print(f"Complete: {len(rn_nums) == len(test_seq)}")
    
    # Show first few assignments
    print(f"\nFirst 10 assignments:")
    for i in range(min(10, len(grn_list))):
        print(f"  {rn_list[i]} -> {grn_list[i]}")