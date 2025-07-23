#!/usr/bin/env python
"""
Test script for GRN assignment extracted from proc.ipynb
This script demonstrates the GRN assignment process for GPCR structures.
"""

import pandas as pd
from pathlib import Path
import os

# Protos imports
from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.grn import GRNProcessor
from protos.processing.grn.grn_table_utils import (
    init_row_from_alignment,
    expand_annotation,
)
from protos.processing.grn.grn_utils import (
    get_seq, sort_grns_str, GRNConfigManager, parse_grn_str2float
)
from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    mmseqs2_align2,
    mmseqs2_align,
    format_alignment
)
from protos.processing.sequence.mmseqs_helper import *
from protos.processing.sequence.mmseqs_utils import *
from protos.processing.grn.grn_assignment import parse_grn_float2str

from protos.processing.grn.grn_utils import *

def main():
    # Setup data directory
    datadir = 'C:\\Users\\hidbe\\PycharmProjects\\protos'
    test_data_root = Path(datadir) / "data"
    test_data_root.mkdir(exist_ok=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    # Initialize processors
    print("Initializing processors...")
    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    
    # Load dataset
    print("\nLoading GPCR dataset...")
    struct_proc.load_dataset('gpcr_agonist_inverse_agonist')
    struct_proc.update_pdb_ids()
    struct_proc.update_pdb_chain_ids()
    
    # Extract sequences
    print("\nExtracting sequences...")
    sequence = struct_proc.get_seq_dict()
    
    # Load reference GRN table
    print("\nLoading reference GRN table...")
    print(f"ProtosPaths data_root: {paths.data_root}")
    print(f"GRN processor ref_dir: {grn_proc.path_ref_dir}")
    
    # Check if file exists
    ref_file = grn_proc.path_ref_dir / "gpcrdb_ref.csv"
    print(f"Looking for reference file at: {ref_file}")
    print(f"File exists: {ref_file.exists()}")
    
    grn_proc.data = grn_proc.load_reference_table('gpcrdb_ref')
    grn_proc.ids = grn_proc.data.index.tolist()
    
    print(f"Loaded {len(grn_proc.data)} sequences from reference table")
    print(f"First 5 IDs: {grn_proc.ids[:5]}")
    print(f"First 10 columns: {list(grn_proc.data.columns[:10])}")
    print(f"Sample GRN data for first sequence:")
    print(grn_proc.data.iloc[0, :10])
    
    # Perform alignment
    print("\nPerforming MMseqs2 alignment...")
    alignment_df = mmseqs2_align2(struct_proc.chain_dict, grn_proc.get_seq_dict())
    
    # Filter by sequence identity
    alignment_df = alignment_df[alignment_df['sequence_identity'] > .25]
    print(f"Found {len(alignment_df)} alignments with > 25% sequence identity")
    
    # Get unique queries
    queries = alignment_df['query_id'].unique().tolist()
    print(f"\nProcessing {len(queries)} unique queries")
    
    # Initialize aligner and config
    aligner = init_aligner()
    config = GRNConfigManager(paths=paths)
    grn_config_std = config.get_config(protein_family='gpcr_a', strict=False)

    grns_str_std = []
    if grn_config_std:
        for region_name, (start_grn, end_grn) in grn_config_std.items():
            # Generate GRNs for this interval (auto-generate since we don't have a predefined list)
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_std.extend(region_grns)

    # Remove duplicates and sort
    grns_str_std = list(set(grns_str_std))
    grns_str_std = sort_grns_str(grns_str_std)
    
    # Process each query
    rows = {}
    print("\nAnnotating sequences...")
    
    for i, q in enumerate(queries):
        print(f"\n{'='*60}")
        print(f"Processing {i+1}/{len(queries)}: {q}")
        print(f"{'='*60}")
        
        # Get best sequence match
        query_alignments = alignment_df[alignment_df['query_id'] == q]
        best_alignment = query_alignments.loc[query_alignments['sequence_identity'].idxmax()]
        ref_id = best_alignment['target_id']
        
        print(f"Best match: {ref_id} (identity: {best_alignment['sequence_identity']:.1%})")
        
        # Get sequences
        test_seq = struct_proc.chain_dict[q]
        grn_proc.apply_interval(grns_str_std)
        ref_seq = grn_proc.get_seq_dict()[ref_id]

        
        # Align sequences
        alignment = align_blosum62(test_seq, ref_seq, aligner)
        formatted = format_alignment(alignment)
        print(f"\nAlignment:\n{formatted}")
        
        # Create initial annotation
        ref_row = grn_proc.data.loc[ref_id]
        ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
        seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
        
        new_row = init_row_from_alignment(formatted, seq_pos2grn)
        
        # Expand annotation
        print(f"\nExpanding annotation...")
        grn_list, rn_list, missing = expand_annotation(
            new_row,
            test_seq.replace('-', ''),
            formatted,
            max_alignment_gap=1,
            protein_family='gpcr_a',
            verbose=2
        )
        
        # Create final row
        final_row = pd.Series(dict(zip(grn_list, rn_list)))
        rows[q] = final_row
        
        print(f"\nAnnotated {len(grn_list)} positions, {len(missing)} missing")
    
    # Create final dataframe
    print("\n\nCreating final GRN table...")
    df = pd.DataFrame.from_dict(rows, orient='index')
    cols = df.columns.tolist()
    df = df[sort_grns_str(cols)].fillna('-')
    
    # Save results
    output_file = "grn_assignments_gpcr.csv"
    df.to_csv(output_file)
    print(f"\nSaved GRN assignments to {output_file}")
    
    # Display summary
    print(f"\nSummary:")
    print(f"- Processed {len(df)} structures")
    print(f"- Annotated {len(df.columns)} GRN positions")
    print(f"- Coverage: {(df != '-').sum().sum() / (len(df) * len(df.columns)):.1%}")
    
    return df


if __name__ == "__main__":
    df = main()
    print("\nDone!")