#!/usr/bin/env python
"""
Annotate mo_small dataset sequences using GRN assignment with section-wise alignment.

This script uses a new approach where:
1. MMseqs2 finds the best reference sequence with high sensitivity
2. For each GRN section (TM1-TM7), we align with high gap penalties to detect and align each section
3. We update the target sequence pivot after each section alignment
"""

import pandas as pd
from pathlib import Path
import os
import sys
import subprocess
import shutil
from typing import Dict, List, Tuple
import json

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


def section_wise_grn_transfer(query_seq: str, ref_seq: str, ref_grn_row: pd.Series, 
                            grn_config: Dict, aligner, verbose: bool = False) -> pd.Series:
    """
    Transfer GRN annotations by aligning each section separately with high gap penalties.
    
    This approach:
    1. For each section (TM1-TM7), extract the reference section sequence
    2. Align it to the remaining query sequence with very high gap opening penalty
    3. Transfer the GRN annotations for that section
    4. Update the query position to the end of the aligned section
    
    Args:
        query_seq: Full query sequence
        ref_seq: Full reference sequence 
        ref_grn_row: GRN annotations for the reference sequence
        grn_config: Dictionary with section definitions (e.g., {'tm1': ['1.44', '1.53'], ...})
        aligner: Biopython aligner (will create modified versions with high gap penalty)
        verbose: Print detailed alignment info
        
    Returns:
        pd.Series with GRN annotations for the query sequence
    """
    # Create a high gap penalty aligner for section detection
    # Very high open gap penalty to prevent gaps within sections
    section_aligner = init_aligner(
        match_score=2, 
        mismatch_score=-1, 
        extend_gap_score=-0.5,  # Less penalty for extending gaps
        open_gap_score=-50      # Very high penalty for opening gaps
    )
    
    # Initialize result dictionary
    query_grn_annotations = {}
    
    # Build reference GRN to position mapping
    ref_grn_to_pos = {}
    pos = 0
    for grn in ref_grn_row.index:
        res = ref_grn_row[grn]
        if res != '-':
            ref_grn_to_pos[grn] = pos
            pos += 1
    
    # Track position in query sequence
    query_offset = 0
    
    # Process each section in order
    section_order = ['tm1', 'tm2', 'tm3', 'tm4', 'tm5', 'tm6', 'tm7']
    
    for section_name in section_order:
        if section_name not in grn_config:
            if verbose:
                print(f"Warning: Section {section_name} not in config, skipping")
            continue
            
        start_grn, end_grn = grn_config[section_name]
        
        # Get GRNs for this section
        section_grns = get_grn_interval(start_grn, end_grn)
        
        # Extract reference section sequence
        ref_section_seq = ""
        ref_section_grns = []
        
        for grn in section_grns:
            if grn in ref_grn_row.index:
                res = ref_grn_row[grn]
                if res != '-' and res != 'X' and res != '.':
                    ref_section_seq += res
                    ref_section_grns.append(grn)
        
        if not ref_section_seq:
            if verbose:
                print(f"Section {section_name}: No sequence found in reference")
            continue
            
        # Get remaining query sequence
        remaining_query = query_seq[query_offset:]
        
        if not remaining_query:
            if verbose:
                print(f"Section {section_name}: No more query sequence to align")
            break
        
        if verbose:
            print(f"\n{'='*60}")
            print(f"Processing section: {section_name} ({start_grn}-{end_grn})")
            print(f"Reference section: {ref_section_seq} ({len(ref_section_seq)} aa)")
            print(f"Query remaining: {len(remaining_query)} aa from position {query_offset}")
        
        # Align section to remaining query with high gap penalty
        try:
            section_alignment = align_blosum62(ref_section_seq, remaining_query, section_aligner)
            section_formatted = format_alignment(section_alignment)
            
            if verbose:
                print(f"\nSection alignment (first 100 chars):")
                print(f"Ref:   {section_formatted[0][:100]}")
                print(f"Match: {section_formatted[1][:100]}")  
                print(f"Query: {section_formatted[2][:100]}")
            
            # Transfer annotations for this section
            ref_aln = section_formatted[0]
            query_aln = section_formatted[2]
            
            ref_pos = 0
            query_pos = query_offset
            
            for i in range(len(ref_aln)):
                ref_char = ref_aln[i]
                query_char = query_aln[i]
                
                # Track position in reference section
                if ref_char != '-':
                    if ref_pos < len(ref_section_grns):
                        current_grn = ref_section_grns[ref_pos]
                        
                        # If query has a residue at this position, transfer the annotation
                        if query_char != '-':
                            query_grn_annotations[current_grn] = query_char
                            
                            if verbose and current_grn.endswith('.50'):
                                print(f"  Annotated key position {current_grn} -> {query_char}")
                    
                    ref_pos += 1
                
                # Track position in query
                if query_char != '-':
                    query_pos += 1
            
            # Update query offset to end of aligned region
            # Find the last aligned position in the query
            last_query_pos = query_offset
            for i in range(len(query_aln)):
                if query_aln[i] != '-':
                    last_query_pos = query_offset + sum(1 for c in query_aln[:i+1] if c != '-')
            
            query_offset = last_query_pos
            
            if verbose:
                print(f"Section {section_name}: Annotated {len([g for g in ref_section_grns if g in query_grn_annotations])} positions")
                print(f"Query position updated to: {query_offset}")
                
        except Exception as e:
            if verbose:
                print(f"Error aligning section {section_name}: {e}")
            continue
    
    # Convert to pandas Series with all GRN positions
    all_grns = sort_grns_str(list(ref_grn_row.index))
    result = pd.Series([query_grn_annotations.get(grn, '-') for grn in all_grns], index=all_grns)
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"Total positions annotated: {(result != '-').sum()}")
        key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
        print("Key positions (x.50) annotated:")
        for pos in key_positions:
            if pos in result.index:
                print(f"  {pos} -> {result[pos]}")
    
    return result


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
    
    # Step 1: Read sequences from mo_small.fasta
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
            "description": "Small set of microbial opsin sequences for GRN annotation (section-wise)",
            "source": "mo_small.fasta",
            "organism": "various microbial species",
            "naming_map": "mo_small_name_mapping.txt",
            "method": "section-wise alignment with high gap penalties"
        }
    )
    print(f"Created dataset 'mo_small' with {len(sequence_ids)} sequences")
    
    # Step 3: Load microbial opsin reference GRN table
    print("\nLoading microbial opsin reference GRN table...")
    
    # Load reference table
    ref_file = grn_proc.path_ref_dir / "mo_ref.csv"
    
    if not ref_file.exists():
        print(f"ERROR: Reference file {ref_file} not found!")
        return None
        
    print(f"Loading reference from: {ref_file}")
    
    # Load and clean the reference table
    grn_proc.data = pd.read_csv(ref_file, index_col=0)
    grn_proc.data = grn_proc.data.fillna('-')
    grn_proc.ids = grn_proc.data.index.tolist()
    
    print(f"Loaded {len(grn_proc.data)} reference sequences")
    print(f"Number of GRN positions: {len(grn_proc.data.columns)}")
    
    # Get reference sequences
    ref_sequences = grn_proc.get_seq_dict()
    print(f"Built {len(ref_sequences)} reference sequences")
    
    # Load GRN config for microbial opsins
    config_manager = GRNConfigManager(paths=paths)
    grn_config = config_manager.get_config(protein_family='microbial_opsins', strict=True)
    
    if not grn_config:
        print("ERROR: No microbial opsin GRN config found!")
        return None
    
    print(f"\nLoaded GRN config with sections: {list(grn_config.keys())}")
    
    # Clean up temp directory
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
    
    # Step 4: Perform MMseqs2 alignment with high sensitivity
    print("\nPerforming MMseqs2 alignment with HIGH SENSITIVITY...")
    print(f"Query sequences: {len(sequences)}")
    print(f"Reference sequences: {len(ref_sequences)}")
    
    alignment_df = mmseqs2_align_sensitive(sequences, ref_sequences)
    
    if alignment_df is None:
        print("ERROR: MMseqs2 alignment failed")
        return None
    
    print(f"\nTotal alignments found: {len(alignment_df)}")
    
    # Filter by sequence identity
    min_identity = 0.15  # 15% identity threshold
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
    
    # Initialize aligner
    aligner = init_aligner(match_score=2, mismatch_score=-1, extend_gap_score=-.05, open_gap_score=-15)
    
    # Step 5: Process each query with section-wise alignment
    rows = {}
    print(f"\nAnnotating {len(queries)} sequences with section-wise alignment...")
    
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
        
        # Get reference GRN row
        ref_row = grn_proc.data.loc[ref_id]
        
        # Perform section-wise GRN transfer
        verbose = (i < 3)  # Verbose for first few sequences
        grn_annotations = section_wise_grn_transfer(
            test_seq, ref_seq, ref_row, grn_config, aligner, verbose=verbose
        )
        
        rows[q] = grn_annotations
        
        # Summary for this sequence
        n_annotated = (grn_annotations != '-').sum()
        print(f"\nSummary for {q}: Annotated {n_annotated} positions")
    
    # Step 6: Create and save final GRN table
    print("\n\nCreating final GRN table...")
    df = pd.DataFrame.from_dict(rows, orient='index')
    
    if df.empty:
        print("ERROR: No sequences could be annotated!")
        return None
        
    cols = df.columns.tolist()
    df = df[sort_grns_str(cols)].fillna('-')
    
    # Save results
    output_name = "mo_small_grn_section_wise"
    
    # Set the processor's data to our results
    grn_proc.data = df
    grn_proc.dataset = output_name
    grn_proc.ids = df.index.tolist()
    grn_proc.grns = df.columns.tolist()
    
    # Save the GRN table
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
    
    # Show key positions across all sequences
    print(f"\nKey positions (x.50) across all sequences:")
    key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    available_keys = [pos for pos in key_positions if pos in df.columns]
    if available_keys:
        print(df[available_keys])
    
    # Compare with original approach by loading the previous results
    original_csv = project_root / "grn_assignments_mo_small.csv"
    if original_csv.exists():
        print("\nComparing with original approach...")
        original_df = pd.read_csv(original_csv, index_col=0)
        
        # Compare coverage
        original_coverage = (original_df != '-').sum().sum() / (len(original_df) * len(original_df.columns))
        new_coverage = (df != '-').sum().sum() / (len(df) * len(df.columns))
        
        print(f"Original coverage: {original_coverage:.1%}")
        print(f"Section-wise coverage: {new_coverage:.1%}")
        
        # Compare key positions
        print("\nKey position comparison:")
        for pos in available_keys:
            if pos in original_df.columns:
                orig_annotated = (original_df[pos] != '-').sum()
                new_annotated = (df[pos] != '-').sum()
                print(f"  {pos}: Original={orig_annotated}, Section-wise={new_annotated}")
    
    return df


if __name__ == "__main__":
    df = main()
    print("\nDone!")