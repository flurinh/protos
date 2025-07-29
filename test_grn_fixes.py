"""Test script to analyze GRN assignment issues in detail."""

import pandas as pd
import numpy as np
from pathlib import Path
import os

# Set up paths
os.environ["PROTOS_DATA_ROOT"] = "/mnt/c/Users/hidbe/PycharmProjects/protos/data"
os.environ["PROTOS_REF_DATA_ROOT"] = "/mnt/c/Users/hidbe/PycharmProjects/protos/data"

from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_table_utils import annotate_grnp
from protos.processing.grn.grn_utils import parse_grn_str2float, parse_grn_float2str, sort_grns_str
from protos.processing.grn.grn_assignment import (
    calculate_missing_ntail_grns, 
    calculate_missing_ctail_grns,
    annotate_gaps_and_loops,
    assign_missing_std_grns
)
from protos.processing.grn.grn_table_utils import is_sequential

def analyze_grn_issues():
    """Analyze potential issues in GRN assignment."""
    
    print("=== GRN Assignment Analysis ===\n")
    
    # Load GRN processor
    print("1. Loading GRN processor...")
    grnp = GRNProcessor(
        name="test_grn"
    )
    
    # Load reference table
    print("2. Loading reference table...")
    ref_data = grnp.load_reference_table("mo_ref")
    grnp.data = ref_data
    print(f"   Loaded {len(grnp.data)} sequences")
    
    # Get a test sequence
    test_seq_name = grnp.data.index[0]
    test_row = grnp.data.iloc[0]
    print(f"   Using test sequence: {test_seq_name}")
    
    # Extract sequence
    seq_parts = []
    grn_positions = []
    for grn, res in test_row.items():
        if res != '-':
            seq_parts.append(res[0])
            grn_positions.append(grn)
    
    test_seq = ''.join(seq_parts)
    print(f"   Sequence length: {len(test_seq)}")
    print(f"   GRN positions: {len(grn_positions)}")
    
    # Test specific problematic scenarios
    print("\n3. Testing GRN assignment edge cases...")
    
    # Test 1: N-terminal assignment
    print("\n   Test 1: N-terminal assignment")
    aligned_grns = {'T10': '1.28', 'V15': '1.33', 'L20': '1.38'}  # Sparse initial alignment
    missing_gene_numbers = [f'X{i}' for i in range(1, 10)]  # Missing N-terminal
    grns_float = [-5.0, -4.0, -3.0, -2.0, -1.0, 1.28, 1.29, 1.30]  # Available GRNs
    
    n_tail, first_pos = calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers, grns_float)
    print(f"      N-tail assignments: {len(n_tail)}")
    if n_tail:
        print(f"      First few: {n_tail[:3]}")
    
    # Test 2: C-terminal assignment
    print("\n   Test 2: C-terminal assignment")
    aligned_grns = {'T280': '7.47', 'V285': '7.52', 'L290': '7.57'}  # End positions
    missing_gene_numbers = [f'X{i}' for i in range(291, 301)]  # Missing C-terminal
    grns_float = [7.47, 7.52, 7.57, 7.58, 7.59, 100.01, 100.02]  # Available GRNs
    
    c_tail, last_pos = calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, 300, grns_float)
    print(f"      C-tail assignments: {len(c_tail)}")
    if c_tail:
        print(f"      First few: {c_tail[:3]}")
    
    # Test 3: Gap/Loop assignment
    print("\n   Test 3: Gap and loop assignment")
    present_list = [('V50', '2.50'), ('L60', '3.28')]  # Gap between TM2 and TM3
    missing = list(range(51, 60))  # Missing positions
    test_seq_gap = 'X' * 100  # Dummy sequence
    
    # Create a simple config for testing
    grn_config = {
        'TM2': ['2.38', '2.68'],
        'TM3': ['3.22', '3.56']
    }
    grns_str = ['2.50', '2.51', '2.52', '3.27', '3.28', '3.29']
    
    nloop, gaps, cloop = annotate_gaps_and_loops(present_list, missing, test_seq_gap, grn_config, grns_str)
    print(f"      N-loops: {len(nloop)}, Gaps: {len(gaps)}, C-loops: {len(cloop)}")
    
    # Test 4: Check for sequential issues
    print("\n   Test 4: Sequential residue checking")
    test_rn_lists = [
        ['A1', 'B2', 'C3', 'D4'],  # Sequential
        ['A1', 'B2', 'D4', 'C3'],  # Out of order
        ['A1', 'B2', 'C4', 'D5'],  # Gap
        ['A1', 'B2', 'C2', 'D3']   # Duplicate
    ]
    
    for i, rn_list in enumerate(test_rn_lists):
        is_seq = is_sequential(rn_list)
        print(f"      List {i+1}: {rn_list} -> Sequential: {is_seq}")
    
    # Test 5: Full annotation on a real sequence
    print("\n4. Testing full annotation pipeline...")
    
    # Create a test sequence with known issues
    test_query = "TESTPROTEIN"
    test_sequence = "MGSHHHHHHGSGLVPRGSHMASMTGGQQMGRGSMGSHHHHHH" + test_seq[:50]  # Add extra N-term
    
    print(f"   Query: {test_query}")
    print(f"   Length: {len(test_sequence)}")
    
    try:
        result = annotate_grnp(
            grnp=grnp,
            query_name=test_query,
            query_seq=test_sequence,
            add_to_GRNP=False,
            verbose=2,
            protein_family='microbial_opsins'
        )
        
        print("\n   Result summary:")
        print(f"   Annotated positions: {len(result)}")
        
        # Check for issues
        grn_list = list(result.index)
        rn_list = list(result.values)
        
        # Check if residues are sequential
        seq_nums = [int(rn[1:]) for rn in rn_list if rn != '-']
        is_sorted = all(seq_nums[i] <= seq_nums[i+1] for i in range(len(seq_nums)-1))
        print(f"   Residues ordered: {is_sorted}")
        
        # Check for gaps
        expected_nums = set(range(1, len(test_sequence) + 1))
        actual_nums = set(seq_nums)
        missing_nums = expected_nums - actual_nums
        print(f"   Missing residues: {len(missing_nums)}")
        if missing_nums and len(missing_nums) < 10:
            print(f"   Missing: {sorted(missing_nums)}")
        
        # Check for duplicates
        duplicates = [x for x in seq_nums if seq_nums.count(x) > 1]
        print(f"   Duplicate residues: {len(set(duplicates))}")
        if duplicates:
            print(f"   Duplicates: {set(duplicates)}")
        
        # Check GRN ordering
        grn_floats = [parse_grn_str2float(g) for g in grn_list]
        grn_sorted = all(grn_floats[i] <= grn_floats[i+1] for i in range(len(grn_floats)-1))
        print(f"   GRNs ordered: {grn_sorted}")
        
        # Check for gap compression
        print("\n   Checking gap regions...")
        for i in range(len(grn_list) - 1):
            grn1 = grn_list[i]
            grn2 = grn_list[i+1]
            
            # Check if both are in standard regions (not loops)
            if '.' in grn1 and '.' in grn2:
                tm1 = grn1.split('.')[0]
                tm2 = grn2.split('.')[0]
                
                if len(tm1) == 1 and len(tm2) == 1 and tm1.isdigit() and tm2.isdigit():
                    # Both are standard positions
                    float1 = parse_grn_str2float(grn1)
                    float2 = parse_grn_str2float(grn2)
                    
                    diff = float2 - float1
                    if diff > 0.02:  # More than 2 positions apart
                        rn1 = int(rn_list[i][1:])
                        rn2 = int(rn_list[i+1][1:])
                        rn_diff = rn2 - rn1
                        
                        if rn_diff == 1:  # Sequential residues with GRN gap
                            print(f"   Gap detected: {grn1} -> {grn2} (GRN gap: {diff:.3f}, RN gap: {rn_diff})")
        
    except Exception as e:
        print(f"\n   ERROR during annotation: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n=== Analysis Complete ===")

if __name__ == "__main__":
    analyze_grn_issues()