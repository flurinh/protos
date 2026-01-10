#!/usr/bin/env python3
"""
Final test of GRN annotation with all fixes
"""

import sys
from pathlib import Path
import pandas as pd

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.grn.grn_table_utils import (
    init_row_from_alignment,
    expand_annotation
)
from protos.processing.grn.grn_utils import get_seq
from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    format_alignment
)
from protos.processing.grn.grn_assignment import parse_grn_float2str

def test_final():
    """Test final GRN annotation."""
    
    print("FINAL GRN ANNOTATION TEST")
    print("=" * 60)
    
    # Test parse_grn_float2str
    print("\n1. Testing parse_grn_float2str:")
    test_values = [1.50, 2.50, 3.50, 7.50, 1.44, 1.01]
    for val in test_values:
        result = parse_grn_float2str(val)
        print(f"  {val} -> {result}")
    
    project_root = Path(__file__).parent.parent
    
    # Load reference
    ref_file = project_root / "src/protos/reference_data/grn/ref/mo_ref.csv"
    ref_table = pd.read_csv(ref_file, index_col=0)
    ref_table = ref_table.fillna('-')
    
    # Test sequence
    test_seq = "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG"
    
    # Get reference
    ref_id = 'BR' if 'BR' in ref_table.index else ref_table.index[0]
    ref_seq = get_seq(ref_id, ref_table)
    
    # Align
    aligner = init_aligner()
    alignment = align_blosum62(test_seq, ref_seq, aligner)
    formatted = format_alignment(alignment)
    
    # Create initial annotation
    ref_row = ref_table.loc[ref_id]
    ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
    seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
    
    new_row = init_row_from_alignment(formatted, seq_pos2grn)
    
    print(f"\n2. Initial annotation: {len(new_row)} positions")
    
    # Check for .5 positions
    dot_5_positions = [k for k in new_row.index if '.5' in k]
    print(f"\n3. Positions with '.5': {dot_5_positions}")
    for pos in dot_5_positions:
        print(f"   {pos}: {new_row[pos]}")
    
    # Expand
    try:
        grn_list, rn_list, missing = expand_annotation(
            new_row,
            test_seq,
            formatted,
            max_alignment_gap=1,
            protein_family='microbial_opsins',
            verbose=0
        )
        
        print(f"\n4. After expansion: {len(grn_list)} positions")
        
        # Check .5 positions in result
        dot_5_in_result = [(grn, rn) for grn, rn in zip(grn_list, rn_list) if '.5' in grn]
        print(f"\n5. Positions with '.5' after expansion: {len(dot_5_in_result)}")
        for grn, rn in dot_5_in_result[:10]:
            print(f"   {grn}: {rn}")
        
        # Create final row
        final_row = pd.Series(dict(zip(grn_list, rn_list)))
        
        # Add missing columns
        for col in ref_table.columns:
            if col not in final_row.index:
                final_row[col] = '-'
        
        # Check specific .50 columns
        print("\n6. Checking specific columns in final output:")
        test_cols = ['1.50', '2.50', '3.50', '7.50']
        for col in test_cols:
            if col in final_row.index:
                print(f"   {col}: {final_row[col]}")
            else:
                print(f"   {col}: NOT IN INDEX")
                
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_final()