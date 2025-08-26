#!/usr/bin/env python
"""
Annotate mo_small dataset sequences using GRN assignment with MORE SENSITIVE alignment.
This script is modified to use higher sensitivity for MMseqs2 alignment to catch divergent sequences.
"""

import pandas as pd
from pathlib import Path
import os
import sys
import subprocess
import shutil

# Setup paths
project_root = Path(__file__).parent.absolute()
sys.path.insert(0, str(project_root))

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
    format_alignment,
    load_alignment_file,
    detect_mmseqs2,
    windows_to_wsl_path
)
from protos.processing.sequence.mmseqs_helper import *
from protos.processing.sequence.mmseqs_utils import *
from protos.io.fasta_utils import read_fasta


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
        # -s 7.5 increases sensitivity (default is 5.7)
        # --max-seqs 1000 to ensure we get all hits
        # -e 10 to allow higher E-values
        # --min-seq-id 0.1 to allow low identity matches (10%)
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
        
        # Convert results to readable format
        subprocess.run(cmd_prefix + [path_mmseqs, 'convertalis', query_db,
                                     sequences_db, results_db,
                                     alignment_tsv], check=True)

        # Clean up
        if use_wsl:
            subprocess.run(['wsl', 'rm', '-rf', tmp_dir])
        else:
            shutil.rmtree(os.path.join(temp_folder, "mmseqs_tmp"))

        # Load alignment results
        alignment_df = load_alignment_file(os.path.join(temp_folder, 'alignment_results.tsv'))
        
        print(f"\nAlignment results: {len(alignment_df)} total alignments")
        print(f"Unique queries aligned: {alignment_df['query_id'].nunique()}")
        print(f"Query IDs found: {sorted(alignment_df['query_id'].unique())}")
        
        return alignment_df

    except subprocess.CalledProcessError as e:
        print(f"MMseqs2 error: {e}")
        return None
    except Exception as e:
        print(f"Error running MMseqs2: {e}")
        return None


def main():
    # Setup data directory
    data_dir = project_root / "data"
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize ProtosPaths
    paths = ProtosPaths()
    
    # Initialize processors
    print("Initializing processors...")
    seq_proc = SequenceProcessor(paths=paths)
    grn_proc = GRNProcessor(paths=paths)
    
    # Step 1: Read sequences from mo_small.fasta (with simplified names)
    print("\nReading sequences from mo_small.fasta...")
    fasta_path = seq_proc.path_fasta_dir / "mo_small.fasta"
    
    if not fasta_path.exists():
        print(f"ERROR: {fasta_path} not found!")
        return None
        
    # Read sequences from FASTA file
    sequences = read_fasta(str(fasta_path))
    print(f"Loaded {len(sequences)} sequences from mo_small.fasta")
    
    # Display loaded sequences
    print("\nSequences loaded:")
    for seq_id, sequence in sequences.items():
        print(f"  {seq_id}: {len(sequence)} residues")
    
    # Step 2: Register sequences and create dataset
    print("\nRegistering sequences and creating mo_small dataset...")
    
    # Get list of sequence IDs for dataset
    sequence_ids = list(sequences.keys())
    
    # Create the mo_small dataset
    seq_proc.create_dataset(
        dataset_name="mo_small",
        entity_names=sequence_ids,
        metadata={
            "description": "Small set of microbial opsin sequences for GRN annotation (simplified names)",
            "source": "mo_small.fasta",
            "organism": "various microbial species",
            "naming_map": "mo_small_name_mapping.txt",
            "note": "Names simplified to MO_001-MO_028 for compatibility"
        }
    )
    print(f"Created dataset 'mo_small' with {len(sequence_ids)} sequences")
    
    # Step 3: Load microbial opsin reference GRN table
    print("\nLoading microbial opsin reference GRN table...")
    
    # Load reference table with proper handling of missing values
    ref_file = grn_proc.path_ref_dir / "mo_ref.csv"
    
    if not ref_file.exists():
        print(f"ERROR: Reference file {ref_file} not found!")
        return None
        
    print(f"Loading reference from: {ref_file}")
    
    # Load and clean the reference table
    grn_proc.data = pd.read_csv(ref_file, index_col=0)
    # Replace NaN values with '-' to avoid the float subscript error
    grn_proc.data = grn_proc.data.fillna('-')
    grn_proc.ids = grn_proc.data.index.tolist()
    
    print(f"Loaded {len(grn_proc.data)} reference sequences")
    print(f"Number of GRN positions: {len(grn_proc.data.columns)}")
    
    # Get reference sequences
    ref_sequences = grn_proc.get_seq_dict()
    print(f"Built {len(ref_sequences)} reference sequences")
    
    # Clean up temp directory before alignment
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
                try:
                    subprocess.run(['wsl', 'rm', '-rf', 'temp/mmseqs_tmp'], check=True)
                    print("Cleaned up using WSL")
                except:
                    print("Could not clean up temp directory, continuing anyway...")
        else:
            shutil.rmtree(mmseqs_tmp_dir)
    
    # Step 4: Perform alignment with HIGHER SENSITIVITY
    print("\nPerforming MMseqs2 alignment with HIGH SENSITIVITY...")
    print(f"Query sequences: {len(sequences)}")
    print(f"Reference sequences: {len(ref_sequences)}")
    
    alignment_df = mmseqs2_align_sensitive(sequences, ref_sequences)
    
    # Check if alignment failed
    if alignment_df is None:
        print("ERROR: MMseqs2 alignment failed, returning None")
        return None
    
    print(f"\nTotal alignments found: {len(alignment_df)}")
    
    # Show identity distribution
    print("\nSequence identity distribution:")
    print(alignment_df['sequence_identity'].describe())
    
    # Filter by sequence identity - use VERY LOW threshold for MO
    min_identity = 0.15  # 15% identity
    filtered_df = alignment_df[alignment_df['sequence_identity'] > min_identity]
    print(f"\nAlignments with > {min_identity*100}% identity: {len(filtered_df)}")
    
    # Get unique queries
    queries = filtered_df['query_id'].unique().tolist()
    print(f"\nProcessing {len(queries)} unique queries (out of {len(sequences)} total)")
    
    # Show which sequences didn't align
    missing_seqs = set(sequences.keys()) - set(queries)
    if missing_seqs:
        print(f"\nSequences with no alignment above {min_identity*100}% identity:")
        for seq in sorted(missing_seqs):
            print(f"  - {seq}")
    
    # Initialize aligner and config
    aligner = init_aligner(match_score=2, mismatch_score=-1, extend_gap_score=-.05, open_gap_score=-15)
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
    print(f"Using {len(grns_str_str)} GRN positions for annotation")
    
    # Step 5: Process each query and annotate
    rows = {}
    print(f"\nAnnotating {len(queries)} sequences...")
    
    for i, q in enumerate(queries):
        print(f"\n{'='*60}")
        print(f"Processing {i+1}/{len(queries)}: {q}")
        print(f"{'='*60}")
        
        # Get best sequence match
        query_alignments = filtered_df[filtered_df['query_id'] == q]
        best_alignment = query_alignments.loc[query_alignments['e_value'].idxmin()]
        ref_id = best_alignment['target_id']
        
        print(f"Best match: {ref_id} (identity: {best_alignment['sequence_identity']:.1%}, E-value: {best_alignment['e_value']:.2e})")
        
        # Get sequences
        test_seq = sequences[q]
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
        print(new_row)
        new_row = new_row[grns]
        
        # Check key positions
        print(f"\nInitial annotation has {len(new_row)} positions")
        key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
        print(f"Key positions (x.50) in initial annotation:")
        for pos in key_positions:
            if pos in new_row.index:
                print(f"  {pos} -> {new_row[pos]}")
        
        # Expand annotation
        try:
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
        except Exception as e:
            print(f"ERROR expanding annotation: {e}")
            # Use initial annotation if expansion fails
            rows[q] = new_row
    
    # Step 6: Create and save final GRN table
    print("\n\nCreating final GRN table...")
    df = pd.DataFrame.from_dict(rows, orient='index')
    
    if df.empty:
        print("ERROR: No sequences could be annotated!")
        return None
        
    cols = df.columns.tolist()
    df = df[sort_grns_str(cols)].fillna('-')
    
    # Save results following Protos conventions
    output_name = "mo_small_grn_annotations"
    
    # Set the processor's data to our results
    grn_proc.data = df
    grn_proc.dataset = output_name
    grn_proc.ids = df.index.tolist()
    grn_proc.grns = df.columns.tolist()
    
    # Save the GRN table using the processor
    grn_proc.save_grn_table(dataset_id=output_name, normalize_formats=False)
    print(f"\nSaved GRN table to: data/grn/tables/{output_name}.csv")
    print(f"Created dataset metadata at: data/grn/datasets/{output_name}.json")
    
    # Display summary
    print(f"\nSummary:")
    print(f"- Total sequences in dataset: {len(sequences)}")
    print(f"- Sequences with alignments: {len(queries)}")
    print(f"- Sequences annotated: {len(df)}")
    print(f"- GRN positions annotated: {len(df.columns)}")
    print(f"- Coverage: {(df != '-').sum().sum() / (len(df) * len(df.columns)):.1%}")
    
    # Show sample of the results
    print(f"\nSample of GRN assignments (first 10 columns):")
    print(df.iloc[:, :10])
    
    # Show key positions across all sequences
    print(f"\nKey positions (x.50) across all sequences:")
    key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    available_keys = [pos for pos in key_positions if pos in df.columns]
    if available_keys:
        print(df[available_keys])
    
    return df


if __name__ == "__main__":
    df = main()
    print("\nDone!")