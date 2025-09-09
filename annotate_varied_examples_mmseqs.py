#!/usr/bin/env python
"""
Annotate opsin sequences using MMseqs2 alignments directly (without re-aligning with Biopython).
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
    format_alignment,
    load_alignment_file,
    detect_mmseqs2,
    windows_to_wsl_path
)
from protos.processing.sequence.mmseqs_helper import *
from protos.processing.sequence.mmseqs_utils import *
import subprocess


def mmseqs2_align_sensitive(query_seqs, ref_seqs, temp_folder='temp'):
    """
    Modified version of mmseqs2_align2 with HIGHER SENSITIVITY for divergent sequences.
    """
    def write_fasta_file(seqs, filename):
        with open(filename, 'w') as fasta_file:
            for key, value in seqs.items():
                fasta_file.write(f'>{key}\n{value}\n')

    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    if not os.path.exists(os.path.join(temp_folder, "mmseqs_tmp")):
        os.makedirs(os.path.join(temp_folder, "mmseqs_tmp"))

    # Detect MMseqs2
    path_mmseqs, use_wsl = detect_mmseqs2()

    if not path_mmseqs:
        print("MMseqs2 not found. Please install it or set MMSEQS_PATH environment variable.")
        return None

    write_fasta_file(ref_seqs, os.path.join(temp_folder, 'ref_seqs.fasta'))
    write_fasta_file(query_seqs, os.path.join(temp_folder, 'query_seqs.fasta'))

    # Set command prefix based on WSL usage
    if use_wsl:
        cmd_prefix = ['wsl']
    else:
        cmd_prefix = []

    try:
        # Convert paths to WSL format if needed
        if use_wsl:
            ref_fasta = windows_to_wsl_path(os.path.join(temp_folder, 'ref_seqs.fasta'))
            query_fasta = windows_to_wsl_path(os.path.join(temp_folder, 'query_seqs.fasta'))
            sequences_db = windows_to_wsl_path(os.path.join(temp_folder, 'mmseqs_tmp', 'sequences_db'))
            query_db = windows_to_wsl_path(os.path.join(temp_folder, 'mmseqs_tmp', 'query_db'))
            results_db = windows_to_wsl_path(os.path.join(temp_folder, 'mmseqs_tmp', 'results'))
            tmp_dir = windows_to_wsl_path(os.path.join(temp_folder, 'mmseqs_tmp'))
            alignment_tsv = windows_to_wsl_path(os.path.join(temp_folder, 'alignment_results.tsv'))
        else:
            ref_fasta = f"{temp_folder}/ref_seqs.fasta"
            query_fasta = f"{temp_folder}/query_seqs.fasta"
            sequences_db = f"{temp_folder}/mmseqs_tmp/sequences_db"
            query_db = f"{temp_folder}/mmseqs_tmp/query_db"
            results_db = f"{temp_folder}/mmseqs_tmp/results"
            tmp_dir = f"{temp_folder}/mmseqs_tmp"
            alignment_tsv = f"{temp_folder}/alignment_results.tsv"

        # Create databases
        subprocess.run(cmd_prefix + [path_mmseqs, 'createdb', ref_fasta,
                                     sequences_db], check=True)
        subprocess.run(cmd_prefix + [path_mmseqs, 'createdb', query_fasta,
                                     query_db], check=True)
        
        # Run search with INCREASED SENSITIVITY
        print("Running MMseqs2 with high sensitivity settings...")
        subprocess.run(cmd_prefix + [
            path_mmseqs, 'search', 
            query_db,
            sequences_db, 
            results_db,
            tmp_dir,
            '-s', '7.5',           # High sensitivity
            '--max-seqs', '10000',  # Get more results
            '-e', '10',            # Higher E-value threshold
            '--min-seq-id', '0.1', # Lower identity threshold (10%)
            '-a'                   # Include all hits
        ], check=True)
        
        # Convert results to readable format with alignment sequences
        subprocess.run(cmd_prefix + [
            path_mmseqs, 'convertalis', query_db,
            sequences_db, results_db,
            alignment_tsv,
            '--format-output', 'query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln'
        ], check=True)

        # Clean up
        import shutil
        if use_wsl:
            subprocess.run(['wsl', 'rm', '-rf', tmp_dir])
        else:
            shutil.rmtree(os.path.join(temp_folder, "mmseqs_tmp"))

        # Load alignment results
        alignment_df = load_alignment_file(os.path.join(temp_folder, 'alignment_results.tsv'))
        
        # Check if we have alignment sequences
        if 'qaln' not in alignment_df.columns:
            # Try to load as custom format
            column_names = ['query_id', 'target_id', 'sequence_identity', 'alignment_length', 
                          'mismatches', 'gap_opens', 'query_start', 'query_end', 
                          'target_start', 'target_end', 'e_value', 'bit_score', 'qaln', 'taln']
            alignment_df = pd.read_csv(os.path.join(temp_folder, 'alignment_results.tsv'), 
                                     sep='\t', header=None, names=column_names)
        
        print(f"\nAlignment results: {len(alignment_df)} total alignments")
        print(f"Unique queries aligned: {alignment_df['query_id'].nunique()}")
        if 'qaln' in alignment_df.columns:
            print("Alignment sequences available")
        
        return alignment_df

    except subprocess.CalledProcessError as e:
        print(f"MMseqs2 error: {e}")
        return None
    except Exception as e:
        print(f"Error running MMseqs2: {e}")
        return None


def transfer_grn_from_mmseqs_alignment(query_seq, ref_seq, ref_grn_row, alignment_info):
    """
    Transfer GRN annotations using MMseqs2 alignment with actual alignment strings.
    
    Args:
        query_seq: Full query sequence string
        ref_seq: Full reference sequence string  
        ref_grn_row: Pandas Series with GRN annotations for reference
        alignment_info: Dict with alignment details including 'qaln' and 'taln'
    
    Returns:
        Pandas Series with transferred GRN annotations
    """
    # Get aligned sequences from MMseqs2
    query_aln = alignment_info.get('qaln', '')
    target_aln = alignment_info.get('taln', '')
    
    if not query_aln or not target_aln:
        # Fallback to coordinate-based transfer if no alignment strings
        return transfer_grn_from_coordinates(query_seq, ref_seq, ref_grn_row, alignment_info)
    
    # Get start positions (1-based from MMseqs2)
    t_start = int(alignment_info['target_start']) - 1  # Convert to 0-based
    
    # Create mapping from reference sequence positions to GRN
    ref_dict = {grn: res for grn, res in ref_grn_row.to_dict().items() if res != '-'}
    ref_grns = list(ref_dict.keys())
    
    # Map reference sequence positions to GRN positions
    ref_pos_to_grn = {}
    ref_pos = 0
    for i, grn in enumerate(ref_grns):
        ref_pos_to_grn[ref_pos] = grn
        ref_pos += 1
    
    # Transfer annotations using the alignment
    transferred = {}
    ref_seq_pos = t_start  # Current position in reference sequence
    
    for i, (q_char, t_char) in enumerate(zip(query_aln, target_aln)):
        if t_char != '-':  # Not a gap in target
            if ref_seq_pos in ref_pos_to_grn:
                grn = ref_pos_to_grn[ref_seq_pos]
                if q_char != '-':  # Not a gap in query
                    transferred[grn] = q_char
            ref_seq_pos += 1
    
    return pd.Series(transferred)


def transfer_grn_from_coordinates(query_seq, ref_seq, ref_grn_row, alignment_info):
    """
    Fallback method using only coordinates (original implementation).
    """
    # Convert to 0-based indexing
    q_start = int(alignment_info['query_start']) - 1
    q_end = int(alignment_info['query_end'])
    t_start = int(alignment_info['target_start']) - 1
    t_end = int(alignment_info['target_end'])
    
    # Get the aligned portions
    query_aligned = query_seq[q_start:q_end]
    ref_aligned = ref_seq[t_start:t_end]
    
    # Create a mapping from reference positions to GRN
    ref_dict = {grn: res for grn, res in ref_grn_row.to_dict().items() if res != '-'}
    ref_grns = list(ref_dict.keys())
    
    # Map reference sequence positions to GRN positions
    ref_pos_to_grn = {}
    ref_pos = 0
    for i, grn in enumerate(ref_grns):
        ref_pos_to_grn[ref_pos] = grn
        ref_pos += 1
    
    # Transfer annotations for aligned region
    transferred = {}
    
    # Simple transfer assuming no gaps
    for i, aa in enumerate(query_aligned):
        ref_pos = t_start + i
        if ref_pos in ref_pos_to_grn:
            grn = ref_pos_to_grn[ref_pos]
            transferred[grn] = aa
    
    return pd.Series(transferred)


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
    print("\nLoading opsin sequences...")
    loaded_sequences = seq_proc.load_entity('mo_small')
    
    # Handle both single and multi-sequence files
    if isinstance(loaded_sequences, str):
        test_sequences = {'opsin_sequences_from_yaml': loaded_sequences}
        print(f"Loaded 1 sequence of {len(loaded_sequences)} residues")
    else:
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
    
    # Load and clean the reference table
    grn_proc.data = pd.read_csv(ref_file, index_col=0)
    grn_proc.data = grn_proc.data.fillna('-')
    grn_proc.ids = grn_proc.data.index.tolist()
    
    print(f"Loaded {len(grn_proc.data)} reference sequences")
    print(f"Number of GRN positions: {len(grn_proc.data.columns)}")
    
    # Build reference sequences from cleaned data
    ref_sequences = {}
    for idx in grn_proc.ids:
        row = grn_proc.data.loc[idx]
        seq_parts = []
        for val in row.values:
            if val not in ['-', 'X', '.', ''] and not pd.isna(val):
                if isinstance(val, str) and len(val) > 0:
                    seq_parts.append(val[0])
        seq = ''.join(seq_parts)
        ref_sequences[idx] = seq
    
    print(f"Built {len(ref_sequences)} reference sequences")
    
    # Clean up temp directory before alignment
    import shutil
    import platform
    mmseqs_tmp_dir = Path("temp/mmseqs_tmp")
    if mmseqs_tmp_dir.exists():
        print("\nCleaning up existing MMseqs2 temp directory...")
        if platform.system() == "Windows":
            try:
                subprocess.run(['wsl', 'rm', '-rf', 'temp/mmseqs_tmp'], check=True)
            except:
                pass
        else:
            shutil.rmtree(mmseqs_tmp_dir)
    
    # Perform alignment with SENSITIVE settings
    print("\nPerforming MMseqs2 alignment with HIGH SENSITIVITY...")
    alignment_df = mmseqs2_align_sensitive(test_sequences, ref_sequences)
    
    if alignment_df is None:
        print("ERROR: MMseqs2 alignment failed")
        return None
    
    # Filter by sequence identity
    alignment_df = alignment_df[alignment_df['sequence_identity'] > .15]
    print(f"Found {len(alignment_df)} alignments with > 15% sequence identity")
    
    # Get unique queries
    queries = alignment_df['query_id'].unique().tolist()
    print(f"\nWill process {len(queries)} unique queries")
    
    # Show sequences that failed to align
    all_sequence_ids = list(test_sequences.keys())
    aligned_ids = set(queries)
    failed_ids = [seq_id for seq_id in all_sequence_ids if seq_id not in aligned_ids]
    
    if failed_ids:
        print(f"\n{len(failed_ids)} sequences failed to align with > 15% identity:")
        for seq_id in failed_ids:
            print(f"  - {seq_id}")
    
    # Get GRN config
    config = GRNConfigManager(paths=paths)
    grn_config_str = config.get_config(protein_family='microbial_opsins', strict=True)
    
    if not grn_config_str:
        grns_str_str = list(grn_proc.data.columns)
    else:
        grns_str_str = []
        for region_name, (start_grn, end_grn) in grn_config_str.items():
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_str.extend(region_grns)
        grns_str_str = list(set(grns_str_str))
    
    grns_str_str = sort_grns_str(grns_str_str)
    print(f"Using {len(grns_str_str)} GRN positions for annotation")
    
    # Process each query using MMseqs2 alignment coordinates
    rows = {}
    print("\nAnnotating sequences using MMseqs2 alignments directly...")
    
    for i, q in enumerate(queries):
        if i % 10 == 0:
            print(f"Processing sequence {i+1}/{len(queries)}...")
        
        # Get best alignment
        query_alignments = alignment_df[alignment_df['query_id'] == q]
        if query_alignments.empty:
            continue
            
        best_alignment = query_alignments.loc[query_alignments['e_value'].idxmin()]
        ref_id = best_alignment['target_id']
        
        # Get sequences
        test_seq = test_sequences[q]
        ref_seq = ref_sequences[ref_id]
        
        # Get reference GRN row
        ref_row = grn_proc.data.loc[ref_id]
        
        # Transfer GRN using alignment coordinates
        alignment_info = {
            'query_start': int(best_alignment['query_start']),
            'query_end': int(best_alignment['query_end']),
            'target_start': int(best_alignment['target_start']),
            'target_end': int(best_alignment['target_end'])
        }
        
        transferred_grn = transfer_grn_from_mmseqs_alignment(
            test_seq, ref_seq, ref_row, alignment_info
        )
        
        # Filter to requested GRN positions
        final_row = pd.Series(index=grns_str_str, dtype=str)
        for grn in grns_str_str:
            if grn in transferred_grn.index:
                final_row[grn] = transferred_grn[grn]
            else:
                final_row[grn] = '-'
        
        rows[q] = final_row
        
        if (i + 1) % 10 == 0:
            print(f"  Processed: {q} -> {ref_id} ({best_alignment['sequence_identity']:.1%} identity)")
    
    # Create final dataframe
    print("\n\nCreating final GRN table...")
    df = pd.DataFrame.from_dict(rows, orient='index')
    if not df.empty:
        cols = df.columns.tolist()
        df = df[sort_grns_str(cols)].fillna('-')
    
    # Save results
    output_name = "mo_small_mmseqs_direct"
    
    grn_proc.data = df
    grn_proc.dataset = output_name
    grn_proc.ids = df.index.tolist()
    grn_proc.grns = df.columns.tolist()
    
    grn_proc.save_grn_table(dataset_id=output_name, normalize_formats=False)
    print(f"\nSaved GRN table using GRNProcessor.save_grn_table()")
    
    # Also save a local copy
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
    
    # Show sample of results
    if not df.empty:
        print(f"\nSample of GRN assignments (first 5 sequences, first 10 positions):")
        print(df.iloc[:5, :10])
    
    return df


if __name__ == "__main__":
    df = main()
    print("\nDone!")