"""
Fixed versions of GRN assignment functions that ensure:
1. No gaps in standard regions (dense assignment)
2. No missing residues  
3. Proper ordering of residues
4. No duplicate residues
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
from protos.processing.grn.grn_utils import (
    parse_grn_str2float, parse_grn_float2str, sort_grns_str
)

def normalize_grn_assignment(grn_rn_pairs: List[Tuple[str, str]], 
                           query_length: int,
                           grn_config: dict,
                           verbose: int = 0) -> List[Tuple[str, str]]:
    """
    Normalize GRN assignments to ensure:
    - All residues are assigned (no missing)
    - Residues are in sequential order
    - No duplicates
    - Gaps are compressed in standard regions
    
    Args:
        grn_rn_pairs: List of (GRN, residue_number) tuples
        query_length: Total length of the query sequence
        grn_config: GRN configuration with TM regions
        verbose: Verbosity level
        
    Returns:
        Normalized list of (GRN, residue_number) tuples
    """
    if verbose >= 1:
        print(f"\n=== Normalizing GRN Assignment ===")
        print(f"Input: {len(grn_rn_pairs)} assignments for {query_length} residues")
    
    # Step 1: Convert to dictionaries for easier manipulation
    rn_to_grn = {}  # residue number -> GRN
    grn_to_rn = {}  # GRN -> residue number
    
    for grn, rn in grn_rn_pairs:
        rn_num = int(rn[1:])
        
        # Check for duplicates
        if rn_num in rn_to_grn:
            if verbose >= 2:
                print(f"  WARNING: Duplicate assignment for {rn}: {rn_to_grn[rn_num]} vs {grn}")
            # Keep the assignment with better GRN (standard over loop)
            old_grn = rn_to_grn[rn_num]
            if _is_standard_grn(grn) and not _is_standard_grn(old_grn):
                rn_to_grn[rn_num] = grn
                grn_to_rn[grn] = rn_num
        else:
            rn_to_grn[rn_num] = grn
            grn_to_rn[grn] = rn_num
    
    # Step 2: Check for missing residues
    all_rn_nums = set(range(1, query_length + 1))
    assigned_rn_nums = set(rn_to_grn.keys())
    missing_rn_nums = all_rn_nums - assigned_rn_nums
    
    if verbose >= 1 and missing_rn_nums:
        print(f"  Found {len(missing_rn_nums)} missing residues")
        if verbose >= 2 and len(missing_rn_nums) < 20:
            print(f"    Missing: {sorted(missing_rn_nums)}")
    
    # Step 3: Fill missing residues
    if missing_rn_nums:
        rn_to_grn = _fill_missing_residues(rn_to_grn, missing_rn_nums, grn_config, verbose)
    
    # Step 4: Compress gaps in standard regions
    rn_to_grn = _compress_standard_gaps(rn_to_grn, grn_config, verbose)
    
    # Step 5: Convert back to list format and ensure ordering
    normalized_pairs = []
    for rn_num in sorted(rn_to_grn.keys()):
        rn = f"{chr(65 + (rn_num - 1) % 26)}{rn_num}"  # A1, B2, etc.
        grn = rn_to_grn[rn_num]
        normalized_pairs.append((grn, rn))
    
    if verbose >= 1:
        print(f"  Output: {len(normalized_pairs)} assignments")
        print(f"  All residues assigned: {len(normalized_pairs) == query_length}")
        
    return normalized_pairs


def _is_standard_grn(grn: str) -> bool:
    """Check if a GRN is a standard TM position (not loop or terminal)."""
    if not grn or grn == '-':
        return False
    
    # N/C terminals
    if grn.startswith('n.') or grn.startswith('c.'):
        return False
    
    # Loops (e.g., 12.003)
    if '.' in grn:
        parts = grn.split('.')
        if len(parts[0]) == 2:  # Loop format
            return False
        elif len(parts[0]) == 1:  # Standard format (e.g., 1.50)
            return True
    
    # Legacy x format
    if 'x' in grn:
        return True
        
    return False


def _get_tm_region(grn: str) -> Optional[int]:
    """Get the TM region number from a GRN (1-8)."""
    if not _is_standard_grn(grn):
        return None
        
    if '.' in grn:
        return int(grn.split('.')[0])
    elif 'x' in grn:
        return int(grn.split('x')[0])
    
    return None


def _fill_missing_residues(rn_to_grn: Dict[int, str], 
                          missing_rn_nums: set,
                          grn_config: dict,
                          verbose: int) -> Dict[int, str]:
    """Fill missing residues with appropriate GRN assignments."""
    
    if verbose >= 2:
        print(f"\n  Filling {len(missing_rn_nums)} missing residues...")
    
    # Sort all residue numbers
    all_rn_nums = sorted(list(rn_to_grn.keys()) + list(missing_rn_nums))
    
    for missing_rn in sorted(missing_rn_nums):
        # Find nearest assigned neighbors
        prev_rn = None
        next_rn = None
        
        for rn in sorted(rn_to_grn.keys()):
            if rn < missing_rn:
                prev_rn = rn
            elif rn > missing_rn and next_rn is None:
                next_rn = rn
                break
        
        # Determine appropriate GRN
        if prev_rn is None and next_rn is not None:
            # Missing at N-terminus
            next_grn = rn_to_grn[next_rn]
            next_float = parse_grn_str2float(next_grn)
            
            if next_float < 0:  # Already in N-term
                new_grn = f"n.{int(abs(next_float)) + (next_rn - missing_rn)}"
            else:
                # Assign to N-terminal
                new_grn = f"n.{next_rn - missing_rn}"
                
        elif prev_rn is not None and next_rn is None:
            # Missing at C-terminus
            prev_grn = rn_to_grn[prev_rn]
            prev_float = parse_grn_str2float(prev_grn)
            
            if prev_float >= 100:  # Already in C-term
                offset = int(prev_float - 100)
                new_grn = f"c.{offset + (missing_rn - prev_rn)}"
            else:
                # Assign to C-terminal
                new_grn = f"c.{missing_rn - prev_rn}"
                
        elif prev_rn is not None and next_rn is not None:
            # Missing in the middle
            prev_grn = rn_to_grn[prev_rn]
            next_grn = rn_to_grn[next_rn]
            
            # Check if between different TM regions
            prev_tm = _get_tm_region(prev_grn)
            next_tm = _get_tm_region(next_grn)
            
            if prev_tm and next_tm and prev_tm != next_tm:
                # In a loop region
                distance_from_prev = missing_rn - prev_rn
                if prev_tm < next_tm:
                    # Normal loop (e.g., between TM2 and TM3)
                    loop_id = f"{prev_tm}{next_tm}"
                else:
                    # Reverse notation
                    loop_id = f"{next_tm}{prev_tm}"
                    
                new_grn = f"{loop_id}.{distance_from_prev:03d}"
            else:
                # Within same region or between non-standard regions
                # Interpolate based on position
                position_ratio = (missing_rn - prev_rn) / (next_rn - prev_rn)
                prev_float = parse_grn_str2float(prev_grn)
                next_float = parse_grn_str2float(next_grn)
                
                # For standard regions, use discrete positions
                if _is_standard_grn(prev_grn) and _is_standard_grn(next_grn):
                    # Round to nearest position
                    interpolated = prev_float + position_ratio * (next_float - prev_float)
                    new_grn = parse_grn_float2str(round(interpolated, 2))
                else:
                    # For loops/terminals, assign to loop
                    if prev_tm:
                        distance_from_prev = missing_rn - prev_rn
                        loop_id = f"{prev_tm}{prev_tm}"
                        new_grn = f"{loop_id}.{distance_from_prev:03d}"
                    else:
                        # Fallback to simple interpolation
                        interpolated = prev_float + position_ratio * (next_float - prev_float)
                        new_grn = parse_grn_float2str(interpolated)
        else:
            # No neighbors (shouldn't happen)
            new_grn = f"n.{missing_rn}"
            
        rn_to_grn[missing_rn] = new_grn
        
        if verbose >= 3:
            print(f"    Filled residue {missing_rn} with {new_grn}")
    
    return rn_to_grn


def _compress_standard_gaps(rn_to_grn: Dict[int, str],
                           grn_config: dict,
                           verbose: int) -> Dict[int, str]:
    """Compress gaps in standard TM regions to ensure dense assignment."""
    
    if verbose >= 2:
        print(f"\n  Compressing gaps in standard regions...")
    
    # Process each TM region
    for region_name, (start_grn, end_grn) in grn_config.items():
        if not region_name.startswith('TM'):
            continue
            
        start_float = parse_grn_str2float(start_grn)
        end_float = parse_grn_str2float(end_grn)
        tm_num = _get_tm_region(start_grn)
        
        if not tm_num:
            continue
        
        # Find all residues assigned to this TM region
        tm_residues = []
        for rn_num, grn in rn_to_grn.items():
            grn_float = parse_grn_str2float(grn)
            if start_float <= grn_float <= end_float and _is_standard_grn(grn):
                tm_residues.append(rn_num)
        
        if not tm_residues:
            continue
            
        # Sort residues
        tm_residues.sort()
        
        # Check for gaps
        has_gap = False
        for i in range(len(tm_residues) - 1):
            if tm_residues[i+1] - tm_residues[i] > 1:
                has_gap = True
                break
        
        if has_gap and verbose >= 2:
            print(f"    Found gap in {region_name}, compressing...")
        
        # Reassign GRNs sequentially within the region
        if has_gap:
            # Get all standard positions in this region
            positions = []
            pos = start_float
            while pos <= end_float:
                positions.append(pos)
                pos = round(pos + 0.01, 2)
            
            # Assign positions sequentially
            for i, rn_num in enumerate(tm_residues):
                if i < len(positions):
                    new_grn = parse_grn_float2str(positions[i])
                    if verbose >= 3 and rn_to_grn[rn_num] != new_grn:
                        print(f"      {rn_num}: {rn_to_grn[rn_num]} -> {new_grn}")
                    rn_to_grn[rn_num] = new_grn
    
    return rn_to_grn


def validate_grn_assignment(grn_rn_pairs: List[Tuple[str, str]], 
                           query_length: int,
                           verbose: int = 0) -> Tuple[bool, List[str]]:
    """
    Validate GRN assignments for common issues.
    
    Returns:
        Tuple of (is_valid, list_of_issues)
    """
    issues = []
    
    # Check 1: All residues assigned
    rn_nums = [int(rn[1:]) for _, rn in grn_rn_pairs]
    missing = set(range(1, query_length + 1)) - set(rn_nums)
    if missing:
        issues.append(f"Missing {len(missing)} residues: {sorted(missing)[:10]}...")
    
    # Check 2: Duplicates
    duplicates = [x for x in rn_nums if rn_nums.count(x) > 1]
    if duplicates:
        issues.append(f"Duplicate residues: {set(duplicates)}")
    
    # Check 3: Sequential ordering
    if rn_nums != sorted(rn_nums):
        issues.append("Residues not in sequential order")
    
    # Check 4: GRN ordering
    grn_floats = [parse_grn_str2float(grn) for grn, _ in grn_rn_pairs]
    if grn_floats != sorted(grn_floats):
        issues.append("GRNs not in sequential order")
    
    # Check 5: Gaps in standard regions
    for i in range(len(grn_rn_pairs) - 1):
        grn1, rn1 = grn_rn_pairs[i]
        grn2, rn2 = grn_rn_pairs[i+1]
        
        if _is_standard_grn(grn1) and _is_standard_grn(grn2):
            tm1 = _get_tm_region(grn1)
            tm2 = _get_tm_region(grn2)
            
            if tm1 == tm2:  # Same TM region
                rn_diff = int(rn2[1:]) - int(rn1[1:])
                grn_diff = parse_grn_str2float(grn2) - parse_grn_str2float(grn1)
                
                if rn_diff == 1 and grn_diff > 0.02:
                    issues.append(f"Gap in TM{tm1}: {grn1}->{grn2} (positions {grn_diff:.2f} apart)")
    
    is_valid = len(issues) == 0
    
    if verbose >= 1:
        print(f"\n=== Validation Results ===")
        print(f"Valid: {is_valid}")
        if issues:
            print("Issues found:")
            for issue in issues:
                print(f"  - {issue}")
    
    return is_valid, issues


if __name__ == "__main__":
    # Test the normalization
    print("Testing GRN normalization...")
    
    # Example with gaps and missing residues
    test_pairs = [
        ('n.5', 'A1'),
        ('n.3', 'B2'),
        ('1.28', 'E5'),   # Gap in N-term
        ('1.30', 'F6'),   # Gap in standard region
        ('1.35', 'G7'),   # Large gap
        ('23.005', 'K11'), # Loop
        ('3.28', 'O15'),  # Missing residues 12-14
        ('3.29', 'P16'),
        ('7.55', 'T20'),  # Missing residues 17-19
    ]
    
    query_length = 25
    grn_config = {
        'TM1': ['1.28', '1.60'],
        'TM2': ['2.38', '2.68'], 
        'TM3': ['3.22', '3.56'],
        'TM7': ['7.30', '7.60']
    }
    
    print(f"\nInput: {len(test_pairs)} assignments for {query_length} residues")
    
    # Validate before
    is_valid, issues = validate_grn_assignment(test_pairs, query_length, verbose=1)
    
    # Normalize
    normalized = normalize_grn_assignment(test_pairs, query_length, grn_config, verbose=2)
    
    # Validate after
    is_valid, issues = validate_grn_assignment(normalized, query_length, verbose=1)
    
    print(f"\nNormalized to {len(normalized)} assignments")
    print("First 10 assignments:")
    for grn, rn in normalized[:10]:
        print(f"  {rn} -> {grn}")