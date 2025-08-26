#!/usr/bin/env python
"""
Annotate opsin sequences from opsin_sequences_from_yaml.fasta.
This script uses the SequenceProcessor to load sequences and annotate them with GRN.
"""

import pandas as pd
from pathlib import Path
import os

# Protos imports
from protos.io.paths import ProtosPaths
from protos.processing.sequence.sequence_processor import SequenceProcessor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_table_utils import (
    init_row_from_alignment,
    expand_annotation,
)
from protos.processing.grn.grn_utils import (
    get_seq, sort_grns_str, GRNConfigManager, parse_grn_str2float,
    parse_grn_float2str, get_grn_interval
)
from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    mmseqs2_align2,
    format_alignment
)
from protos.processing.sequence.mmseqs_helper import *
from protos.processing.sequence.mmseqs_utils import *

def main():
    # Setup data directory
    datadir = Path(__file__).parent.absolute()
    test_data_root = datadir / "data"
    test_data_root.mkdir(exist_ok=True, parents=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    # Initialize processors
    print("Initializing processors...")
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    
    # Load sequences using SequenceProcessor
    print("\nLoading opsin sequences from YAML...")
    # load_entity will automatically handle multi-sequence FASTA files
    loaded_sequences = seq_proc.load_entity('opsin_sequences_from_yaml')
    
    # Handle both single and multi-sequence files
    if isinstance(loaded_sequences, str):
        # Single sequence
        test_sequences = {'opsin_sequences_from_yaml': loaded_sequences}
        print(f"Loaded 1 sequence of {len(loaded_sequences)} residues")
    else:
        # Multi-sequence file
        test_sequences = loaded_sequences
        print(f"Loaded {len(test_sequences)} sequences")
    
    # Display first 10 sequences
    print("\nFirst 10 sequences:")
    for i, (seq_id, sequence) in enumerate(test_sequences.items()):
        if i >= 10:
            print(f"  ... and {len(test_sequences) - 10} more sequences")
            break
        print(f"  {seq_id}: {len(sequence)} residues")
    
    # Load microbial opsin reference GRN table
    print("\nLoading microbial opsin reference GRN table...")
    ref_file = grn_proc.path_ref_dir / "mo_ref.csv"
    print(f"Looking for reference file at: {ref_file}")
    print(f"File exists: {ref_file.exists()}")
    
    # Load and clean the reference table
    grn_proc.data = pd.read_csv(ref_file, index_col=0)
    # Replace NaN values with '-' to avoid float subscript errors
    grn_proc.data = grn_proc.data.fillna('-')
    grn_proc.ids = grn_proc.data.index.tolist()
    
    print(f"Loaded {len(grn_proc.data)} reference sequences")
    print(f"First 5 Reference IDs: {grn_proc.ids[:5]}")
    print(f"Number of GRN positions: {len(grn_proc.data.columns)}")
    
    # Build reference sequences from cleaned data
    ref_sequences = {}
    for idx in grn_proc.ids:
        row = grn_proc.data.loc[idx]
        # Build sequence from non-gap positions
        seq_parts = []
        for val in row.values:
            if val not in ['-', 'X', '.', ''] and not pd.isna(val):
                # Handle both single letters and longer annotations
                if isinstance(val, str) and len(val) > 0:
                    seq_parts.append(val[0])
        ref_sequences[idx] = ''.join(seq_parts)
    
    print(f"Built {len(ref_sequences)} reference sequences")
    
    # Clean up temp directory before alignment
    import shutil
    import platform
    mmseqs_tmp_dir = Path("temp/mmseqs_tmp")
    if mmseqs_tmp_dir.exists():
        print("\nCleaning up existing MMseqs2 temp directory...")
        if platform.system() == "Windows":
            try:
                for item in mmseqs_tmp_dir.iterdir():
                    if item.is_symlink():
                        item.unlink()
                shutil.rmtree(mmseqs_tmp_dir)
            except Exception as e:
                print(f"Warning: Could not fully clean temp directory: {e}")
                import subprocess
                try:
                    subprocess.run(['wsl', 'rm', '-rf', 'temp/mmseqs_tmp'], check=True)
                    print("Cleaned up using WSL")
                except:
                    print("Could not clean up temp directory, continuing anyway...")
        else:
            shutil.rmtree(mmseqs_tmp_dir)
    
    # Perform alignment
    print("\nPerforming MMseqs2 alignment...")
    alignment_df = mmseqs2_align2(test_sequences, ref_sequences)
    
    if alignment_df is None:
        print("ERROR: MMseqs2 alignment failed, returning None")
        return None
    
    # Filter by sequence identity - use lower threshold for diverse opsins
    alignment_df = alignment_df[alignment_df['sequence_identity'] > .15]
    print(f"Found {len(alignment_df)} alignments with > 15% sequence identity")
    
    # Get unique queries
    queries = alignment_df['query_id'].unique().tolist()
    print(f"\nWill process {len(queries)} unique queries")
    
    # Initialize aligner and config
    aligner = init_aligner()
    config = GRNConfigManager(paths=paths)
    
    # Get microbial opsin config
    grn_config_str = config.get_config(protein_family='microbial_opsins', strict=True)
    
    if not grn_config_str:
        print("No specific microbial opsin config found, using all GRN positions from reference")
        grns_str_str = list(grn_proc.data.columns)
    else:
        grns_str_str = []
        for region_name, (start_grn, end_grn) in grn_config_str.items():
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_str.extend(region_grns)
        grns_str_str = list(set(grns_str_str))
    
    grns_str_str = sort_grns_str(grns_str_str)
    print(f"Using {len(grns_str_str)} GRN positions for annotation")
    
    # Process each query
    rows = {}
    print("\nAnnotating sequences...")
    
    # Process in batches to show progress
    batch_size = 10
    for batch_start in range(0, len(queries), batch_size):
        batch_end = min(batch_start + batch_size, len(queries))
        print(f"\nProcessing batch {batch_start//batch_size + 1} (sequences {batch_start+1}-{batch_end})...")
        
        for i in range(batch_start, batch_end):
            q = queries[i]
            
            # Get best sequence match
            query_alignments = alignment_df[alignment_df['query_id'] == q]
            if query_alignments.empty:
                print(f"  Skipping {q}: No alignments found")
                continue
                
            best_alignment = query_alignments.loc[query_alignments['e_value'].idxmin()]
            ref_id = best_alignment['target_id']
            
            print(f"  {q}: Best match {ref_id} (identity: {best_alignment['sequence_identity']:.1%})")
            
            # Get sequences
            test_seq = test_sequences[q]
            ref_seq = ref_sequences[ref_id]
            
            # Align sequences
            alignment = align_blosum62(test_seq, ref_seq, aligner)
            formatted = format_alignment(alignment)
            
            # Create initial annotation
            ref_row = grn_proc.data.loc[ref_id]
            ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
            seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
            
            new_row = init_row_from_alignment(formatted, seq_pos2grn)
            grns = [grn for grn in grns_str_str if grn in new_row.index]
            new_row = new_row[grns]
            
            # Expand annotation with error handling
            try:
                grn_list, rn_list, missing = expand_annotation(
                    new_row,
                    test_seq.replace('-', ''),
                    formatted,
                    max_alignment_gap=1,
                    protein_family='microbial_opsins',
                    verbose=0  # Quiet mode for batch processing
                )
                
                # Create final row
                final_row = pd.Series(dict(zip(grn_list, rn_list)))
                rows[q] = final_row
                
            except Exception as e:
                print(f"    Warning: Failed to expand annotation for {q}: {e}")
                # Use initial annotation if expansion fails
                rows[q] = new_row
    
    # Create final dataframe
    print("\n\nCreating final GRN table...")
    df = pd.DataFrame.from_dict(rows, orient='index')
    if not df.empty:
        cols = df.columns.tolist()
        df = df[sort_grns_str(cols)].fillna('-')
    
    # Save results following Protos conventions
    output_name = "opsin_sequences_from_yaml_grn_annotations"
    
    # Set the processor's data to our results
    grn_proc.data = df
    grn_proc.dataset = output_name
    grn_proc.ids = df.index.tolist()
    grn_proc.grns = df.columns.tolist()
    
    # Save the GRN table using the processor
    grn_proc.save_grn_table(dataset_id=output_name, normalize_formats=False)
    print(f"\nSaved GRN table using GRNProcessor.save_grn_table()")
    
    # Also save a local copy for reference
    output_file = f"grn_assignments_{output_name}.csv"
    df.to_csv(output_file)
    print(f"Also saved local copy to {output_file}")
    
    # Display summary
    print(f"\nSummary:")
    print(f"- Processed {len(df)} sequences out of {len(test_sequences)} total")
    print(f"- Annotated {len(df.columns)} GRN positions")
    if not df.empty:
        coverage = (df != '-').sum().sum() / (len(df) * len(df.columns))
        print(f"- Average coverage: {coverage:.1%}")
    
    # Show sample of the results
    if not df.empty:
        print(f"\nSample of GRN assignments (first 5 sequences, first 10 positions):")
        print(df.iloc[:5, :10])
        
        # Show key positions for first few sequences
        print(f"\nKey positions (x.50) for first 5 sequences:")
        key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
        for seq_id in df.index[:5]:
            print(f"\n{seq_id}:")
            for pos in key_positions:
                if pos in df.columns and df.loc[seq_id, pos] != '-':
                    print(f"  {pos}: {df.loc[seq_id, pos]}")
    
    return df


if __name__ == "__main__":
    df = main()
    print("\nDone!")