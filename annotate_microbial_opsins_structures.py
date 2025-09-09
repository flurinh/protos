#!/usr/bin/env python
"""
Annotate microbial opsin sequences using GRN assignment.
This script follows Protos data management principles by:
1. Loading structures from a dataset
2. Extracting sequences from structures
3. Aligning against microbial opsin reference GRN table
"""

import pandas as pd
from pathlib import Path
import os

# Protos imports
from protos.io.paths import ProtosPaths
from protos.processing.structure.structure_processor import StructureProcessor
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
    mmseqs2_align,
    format_alignment
)
from protos.processing.sequence.mmseqs_helper import *
from protos.processing.sequence.mmseqs_utils import *

def main():
    # Setup data directory - use current directory
    datadir = Path(__file__).parent.absolute()
    test_data_root = datadir / "data"
    test_data_root.mkdir(exist_ok=True, parents=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    # Initialize processors
    print("Initializing processors...")
    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    
    # Option 1: Load from structure dataset if available
    try:
        print("\nAttempting to load microbial opsin structures...")
        struct_proc.load_dataset('test_mo')
        struct_proc.update_pdb_ids()
        struct_proc.update_pdb_chain_ids()
        
        # Extract sequences from structures
        print("\nExtracting sequences from structures...")
        test_sequences = struct_proc.get_seq_dict()
        print(f"Loaded {len(test_sequences)} sequences from structures")
        
    except Exception as e:
        print(f"Could not load structure dataset: {e}")
        print("\nFalling back to loading sequences directly...")
        
        # Option 2: Load sequences directly from FASTA
        # This is less ideal but works if no structure dataset exists
        from protos.io.fasta_utils import read_fasta
        
        # Check for test_mo.fasta
        fasta_path = seq_proc.path_fasta_dir / "test_mo.fasta"
        if fasta_path.exists():
            test_sequences = read_fasta(str(fasta_path))
            print(f"Loaded {len(test_sequences)} sequences from {fasta_path}")
        else:
            # Try loading individual sequence files for known MO structures
            test_sequences = {}
            mo_ids = ['1UAZ', '3DDL', '4PXK', '1C3W', '1VGO', '1XIO', '1JGJ']
            for mo_id in mo_ids:
                try:
                    seq = seq_proc.load_entity(mo_id)
                    if seq:
                        if isinstance(seq, str):
                            test_sequences[mo_id] = seq
                        elif isinstance(seq, dict):
                            test_sequences.update(seq)
                except:
                    pass
            print(f"Loaded {len(test_sequences)} sequences from individual files")
    
    if not test_sequences:
        print("ERROR: No sequences found to process!")
        return None
    
    # Display loaded sequences
    for seq_id, sequence in test_sequences.items():
        print(f"  {seq_id}: {len(sequence)} residues")
    
    # Load microbial opsin reference GRN table
    print("\nLoading microbial opsin reference GRN table...")
    print(f"GRN processor ref_dir: {grn_proc.path_ref_dir}")
    
    # Load reference table with proper handling of missing values
    ref_file = grn_proc.path_ref_dir / "mo_ref.csv"
    print(f"Looking for reference file at: {ref_file}")
    print(f"File exists: {ref_file.exists()}")
    
    # Load and clean the reference table
    grn_proc.data = pd.read_csv(ref_file, index_col=0)
    # Replace NaN values with '-' to avoid the float subscript error
    grn_proc.data = grn_proc.data.fillna('-')
    grn_proc.ids = grn_proc.data.index.tolist()
    
    print(f"Loaded {len(grn_proc.data)} reference sequences")
    print(f"First 5 Reference IDs: {grn_proc.ids[:5]}")
    print(f"Number of GRN positions: {len(grn_proc.data.columns)}")
    
    # Get reference sequences - need to handle the cleaned data
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
        print("Cleaning up existing MMseqs2 temp directory...")
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
    
    # Check if alignment failed
    if alignment_df is None:
        print("ERROR: MMseqs2 alignment failed, returning None")
        return None
    
    # Filter by sequence identity - use lower threshold for MO
    alignment_df = alignment_df[alignment_df['sequence_identity'] > .20]
    print(f"Found {len(alignment_df)} alignments with > 20% sequence identity")
    
    # Get unique queries
    queries = alignment_df['query_id'].unique().tolist()
    print(f"\nProcessing {len(queries)} unique queries")
    
    # Initialize aligner and config
    aligner = init_aligner()
    config = GRNConfigManager(paths=paths)


    # Try to get microbial opsin config
    grn_config_str = config.get_config(protein_family='microbial_opsins', strict=True)

    # If no specific config, use all GRN positions from reference
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
    print(f"strict {len(grns_str_str)} GRN positions for annotation:", grns_str_str)

    # Try to get microbial opsin config
    grn_config_str = config.get_config(protein_family='microbial_opsins', strict=False)
    
    # If no specific config, use all GRN positions from reference
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
    print(f"standard {len(grns_str_str)} GRN positions for annotation:", grns_str_str)

    
    # Process each query
    rows = {}
    print("\nAnnotating sequences...")
    
    for i, q in enumerate(queries):
        print(f"\n{'='*60}")
        print(f"Processing {i+1}/{len(queries)}: {q}")
        print(f"{'='*60}")
        
        # Get best sequence match
        query_alignments = alignment_df[alignment_df['query_id'] == q]
        print(f"query {q} alignments", query_alignments.head(5))
        best_alignment = query_alignments.loc[query_alignments['e_value'].idxmin()]
        ref_id = best_alignment['target_id']
        
        print(f"Best match: {ref_id} (identity: {best_alignment['sequence_identity']:.1%})")
        
        # Get sequences
        test_seq = test_sequences[q]
        ref_seq = ref_sequences[ref_id]
        
        # Align sequences
        alignment = align_blosum62(test_seq, ref_seq, aligner)
        formatted = format_alignment(alignment)
        print(f"\nAlignment preview (first 100 chars):")
        print(f"Query:  {formatted[0][:100]}...")
        print(f"Match:  {formatted[1][:100]}...")
        print(f"Ref:    {formatted[2][:100]}...")
        
        # Create initial annotation
        ref_row = grn_proc.data.loc[ref_id]
        ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
        seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
        
        new_row = init_row_from_alignment(formatted, seq_pos2grn)
        grns = [grn for grn in grns_str_str if grn in new_row.index]
        new_row = new_row[grns]

        print("new_row", new_row.head(30))
        
        # Check key positions
        print(f"\nInitial annotation has {len(new_row)} positions")
        key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
        print(f"Key positions (x.50) in initial annotation:")
        for pos in key_positions:
            if pos in new_row.index:
                print(f"  {pos} -> {new_row[pos]}")
        
        # Expand annotation
        print(f"\nExpanding annotation...")
        grn_list, rn_list, missing = expand_annotation(
            new_row,
            test_seq.replace('-', ''),
            formatted,
            max_alignment_gap=1,
            protein_family='microbial_opsins',
            verbose=1
        )
        
        # Check key positions in final annotation
        key_in_final = [(g, r) for g, r in zip(grn_list, rn_list) if g.endswith('.50')]
        print(f"\nKey positions (x.50) in final annotation: {len(key_in_final)}")
        for grn, rn in key_in_final:
            print(f"  {grn} -> {rn}")
        
        # Create final row
        final_row = pd.Series(dict(zip(grn_list, rn_list)))
        rows[q] = final_row
        
        print(f"\nSummary: Annotated {len(grn_list)} positions, {len(missing)} missing")
    
    # Create final dataframe
    print("\n\nCreating final GRN table...")
    df = pd.DataFrame.from_dict(rows, orient='index')
    cols = df.columns.tolist()
    df = df[sort_grns_str(cols)].fillna('-')
    
    # Save results following Protos conventions
    output_name = "microbial_opsins_grn_annotations"
    
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
    print(f"- Processed {len(df)} sequences")
    print(f"- Annotated {len(df.columns)} GRN positions")
    print(f"- Coverage: {(df != '-').sum().sum() / (len(df) * len(df.columns)):.1%}")
    
    # Show sample of the results
    print(f"\nSample of GRN assignments:")
    print(df.iloc[:, :10])  # First 10 columns
    
    return df


if __name__ == "__main__":
    df = main()
    print("\nDone!")