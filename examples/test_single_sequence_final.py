#!/usr/bin/env python3
"""
Test single sequence annotation to verify .50 positions
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
from protos.io.fasta_utils import read_fasta

def test_single():
    """Test single sequence."""
    
    print("TEST SINGLE SEQUENCE ANNOTATION")
    print("=" * 60)
    
    project_root = Path(__file__).parent.parent
    
    # Load reference
    ref_file = project_root / "src/protos/reference_data/grn/ref/mo_ref.csv"
    ref_table = pd.read_csv(ref_file, index_col=0)
    ref_table = ref_table.fillna('-')
    
    # Load AR1_A sequence
    sequences = read_fasta(str(project_root / "src/protos/reference_data/sequence/fasta/opsin_sequences_from_yaml.fasta"))
    test_seq = sequences['AR1_A']
    
    print(f"Test sequence: AR1_A")
    print(f"Length: {len(test_seq)}")
    
    # Use BR as reference
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
    
    # Expand
    grn_list, rn_list, missing = expand_annotation(
        new_row,
        test_seq,
        formatted,
        max_alignment_gap=1,
        protein_family='microbial_opsins',
        verbose=0
    )
    
    # Create final row
    final_row = pd.Series(dict(zip(grn_list, rn_list)))
    
    # Add missing columns
    for col in ref_table.columns:
        if col not in final_row.index:
            final_row[col] = '-'
    
    # Reorder
    final_row = final_row[ref_table.columns]
    
    # Check .50 positions
    print("\n.50 positions in AR1_A:")
    dot_50_cols = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    for col in dot_50_cols:
        print(f"  {col}: {final_row[col]}")
    
    # Save to proper GRN tables location
    output_file = project_root / "src/protos/reference_data/grn/tables" / "ar1_test.csv"
    output_file.parent.mkdir(exist_ok=True, parents=True)
    
    # Create DataFrame with just this sequence
    output_df = pd.DataFrame([final_row], index=['AR1_A'])
    output_df.to_csv(output_file)
    
    print(f"\nSaved to: {output_file}")
    
    # Also check the first 10 values
    print("\nFirst 10 positions:")
    for i, col in enumerate(ref_table.columns[:10]):
        print(f"  {col}: {final_row[col]}")

if __name__ == "__main__":
    test_single()