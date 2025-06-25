#!/usr/bin/env python3
"""
Test expand_annotation to identify the list index error
"""

import sys
from pathlib import Path
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.grn.grn_table_utils import expand_annotation
from protos.processing.grn.grn_assignment import assign_gene_nr
from protos.processing.sequence.seq_alignment import init_aligner, align_blosum62, format_alignment
from protos.io.fasta_utils import read_fasta

def test_expand_annotation():
    """Debug expand_annotation issues."""
    
    print("Testing expand_annotation")
    print("=" * 60)
    
    # Test data
    query_seq = "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVK"
    
    # Create a simple new_row with GRN assignments
    new_row_data = {
        '1x44': 'G11',
        '1x45': 'V12', 
        '1x46': 'S13',
        '1x47': 'Q14',
        '1x48': 'A15',
        '1x49': 'Q16',
        '1x50': 'I17'
    }
    new_row = pd.Series(new_row_data)
    
    print("Input row:")
    print(new_row)
    
    # Create alignment
    aligner = init_aligner()
    ref_seq = "GVSQAQITGRPEWIWLALGTALMGLGTLYFLVK"
    alignment = align_blosum62(query_seq, ref_seq, aligner)
    formatted = format_alignment(alignment)
    
    print("\nAlignment:")
    print(f"Query:  {formatted[0]}")
    print(f"Match:  {formatted[1]}")
    print(f"Ref:    {formatted[2]}")
    
    # Test expand_annotation
    print("\nTesting expand_annotation...")
    try:
        grn_list, rn_list, missing = expand_annotation(
            new_row,
            query_seq,
            formatted,
            max_alignment_gap=1,
            protein_family='microbial_opsins',
            verbose=1  # Enable verbose output
        )
        
        print(f"\nSuccess!")
        print(f"GRN list: {len(grn_list)} positions")
        print(f"RN list: {len(rn_list)} positions")
        print(f"Missing: {len(missing)} positions")
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        
        # Debug the specific error
        print("\n\nDebugging the error...")
        
        # Check what happens in expand_annotation

        aligned_grns = {v: k for (k, v) in new_row.to_dict().items() if (v != '-') & ('x' in k)}
        print(f"\naligned_grns: {aligned_grns}")
        
        # Check assign_gene_nr
        all_query_gene_numbers = assign_gene_nr(query_seq)
        print(f"\nall_query_gene_numbers (first 10): {all_query_gene_numbers[:10]}")
        
        # Check the reference dict format
        reference_grn_dict = {v: k for k, v in aligned_grns.items()}
        print(f"\nreference_grn_dict: {reference_grn_dict}")
        
        # Try get_correctly_aligned_grns
        from protos.processing.grn.grn_assignment import get_correctly_aligned_grns
        
        try:
            result = get_correctly_aligned_grns(
                all_query_gene_numbers,
                reference_grn_dict,
                formatted,
                max_alignment_gap=1
            )
            print(f"\nget_correctly_aligned_grns result: {result}")
        except Exception as e2:
            print(f"\nError in get_correctly_aligned_grns: {e2}")
            traceback.print_exc()

if __name__ == "__main__":
    test_expand_annotation()