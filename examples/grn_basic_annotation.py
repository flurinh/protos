#!/usr/bin/env python3
"""
Basic GRN Annotation Example
============================

This script demonstrates basic GRN (Generic Residue Numbering) annotation 
for microbial opsin sequences using the Protos framework.

Key Features:
- Zero configuration setup using ProtosPaths
- Loads real microbial opsin sequences 
- Uses mo_ref as reference GRN table
- Assigns GRN positions to query sequences
- Saves annotated GRN table

Usage:
    python examples/grn_basic_annotation.py
    
Output:
    - GRN table saved to: data/grn/tables/opsin_grn_annotations.csv
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from protos.processing.grn import GRNProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.sequence import (
    init_aligner, align_blosum62, format_alignment, mmseqs2_align2
)
from protos.processing.grn.grn_utils import (
    get_seq, sort_grns_str, GRNConfigManager, parse_grn_str2float
)
from protos.processing.grn.grn_table_utils import (
    init_row_from_alignment, expand_annotation
)
import pandas as pd
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Run GRN annotation with real microbial opsin data."""
    
    print("\n" + "="*80)
    print("GRN ANNOTATION WITH REAL MICROBIAL OPSIN DATA")
    print("="*80)
    
    # Enable debug logging
    logging.basicConfig(level=logging.DEBUG, 
                       format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # 1. Initialize processors with zero configuration
    print("\n1. Initializing Processors (Zero Configuration)")
    print("-" * 60)
    
    # Set up ProtosPaths to use the Windows data directory from WSL
    from protos.io.paths.path_config import ProtosPaths
    
    # Windows path converted to WSL format
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir)
    paths = ProtosPaths(data_root=str(data_dir))
    
    grn_proc = GRNProcessor(name="opsin_grn_demo", paths=paths)
    seq_proc = SequenceProcessor(name="opsin_seq_demo", paths=paths)
    
    print(f"✅ GRN Processor initialized")
    print(f"✅ Sequence Processor initialized")
    
    # Initialize GRN configuration for microbial opsins
    config_manager = GRNConfigManager(paths=paths)
    
    # 2. Load real opsin sequences
    print("\n2. Loading Real Opsin Sequences")
    print("-" * 60)
    
    # The processor will look for sequences in its managed directories
    # No hardcoded paths needed!
    sequence_name = "opsin_sequences_from_yaml"
    
    # Try to load the sequences using the processor
    try:
        # First check if it exists as an entity
        if seq_proc.entity_exists(sequence_name):
            sequences = seq_proc.load_entity(sequence_name)
            if isinstance(sequences, str):
                # Single sequence, wrap in dict
                sequences = {sequence_name: sequences}
        else:
            # Try loading as a dataset file
            sequences = seq_proc.load_sequences(sequence_name)
    except FileNotFoundError:
        print(f"❌ Sequence file '{sequence_name}' not found")
        print("Please ensure the opsin sequences are in the sequence/fasta/ directory")
        print(f"Expected location: {seq_proc.path_fasta_dir}/{sequence_name}.fasta")
        return
    print(f"✅ Loaded {len(sequences)} opsin sequences")
    
    # Show first few sequences
    print("\nOpsin sequences loaded:")
    for i, (seq_id, seq) in enumerate(sequences.items()):
        if i < 5:
            print(f"  - {seq_id}: {len(seq)} AA")
    
    # 3. Load reference GRN table for microbial opsins
    print("\n3. Loading Reference GRN Table")
    print("-" * 60)
    
    # Load mo_ref as the reference table
    try:
        # Load mo_ref.csv directly as it has a different format
        mo_ref_path = data_dir / "grn" / "ref" / "mo_ref.csv"
        if not mo_ref_path.exists():
            raise FileNotFoundError(f"mo_ref.csv not found at {mo_ref_path}")
            
        # Read mo_ref with first column as index
        mo_ref_df = pd.read_csv(mo_ref_path, index_col=0)
        
        # Set it in the GRN processor
        grn_proc.data = mo_ref_df
        grn_proc.ids = list(mo_ref_df.index)
        grn_proc.grns = list(mo_ref_df.columns)
        
        print(f"✅ Loaded reference table 'mo_ref' with {len(mo_ref_df)} proteins")
        print(f"   Reference proteins: {list(mo_ref_df.index[:5])}...")
        print(f"   Number of GRN positions: {len(mo_ref_df.columns)}")
        
        # Show TM coverage in reference
        tm_coverage = {}
        for col in mo_ref_df.columns:
            if '.' in col and col[0].isdigit():
                try:
                    tm = int(col.split('.')[0])
                    if 1 <= tm <= 7:
                        tm_coverage[tm] = tm_coverage.get(tm, 0) + 1
                except:
                    pass
        print(f"   TM coverage: {sorted(tm_coverage.items())}")
        
    except Exception as e:
        print(f"❌ Could not load mo_ref reference table: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # 4. Process each opsin sequence
    print("\n4. Processing Opsin Sequences")
    print("-" * 60)
    
    # Get GRN intervals for microbial opsins
    grn_intervals = config_manager.get_intervals("microbial_opsins")
    if not grn_intervals:
        print("⚠️  No GRN intervals found for microbial_opsins, using all positions")
    else:
        print(f"✅ Loaded GRN intervals for microbial opsins")
    
    # Get reference sequences from mo_ref table
    # The mo_ref table contains residue identifiers like "T144" at each GRN position
    # We need to extract just the amino acid codes to build sequences
    ref_sequences = {}
    
    for ref_id in grn_proc.data.index:
        # Extract amino acid sequence from the GRN positions
        grn_row = grn_proc.data.loc[ref_id]
        # Get only the amino acid codes (first letter of residue identifier)
        aa_seq = ""
        for grn_pos in sorted(grn_row.index, key=lambda x: parse_grn_str2float(x)):
            residue = grn_row[grn_pos]
            if residue != '-' and isinstance(residue, str) and len(residue) > 0:
                # Check if it's a valid residue identifier (e.g., "T144")
                if residue[0] in 'ACDEFGHIKLMNPQRSTVWY':
                    aa_seq += residue[0]
        
        if aa_seq:
            ref_sequences[ref_id] = aa_seq
    
    print(f"✅ Extracted {len(ref_sequences)} reference sequences from mo_ref")
    if ref_sequences:
        print(f"   Reference sequence lengths: {[len(seq) for seq in list(ref_sequences.values())[:3]]}...")
    
    # Initialize aligner once for all sequences
    aligner = init_aligner()
    
    # Process all opsin sequences
    results = {}
    num_to_process = len(sequences)  # Process ALL sequences
    
    print(f"\nProcessing all {num_to_process} sequences")
    
    # Process sequences in batches for better performance
    for i, (seq_id, query_seq) in enumerate(sequences.items()):
        if i >= num_to_process:
            break
        
        # Show progress every 10 sequences
        if i % 10 == 0:
            print(f"\nProgress: {i}/{num_to_process} sequences processed...")
            
        print(f"[{i+1}/{num_to_process}] {seq_id} ({len(query_seq)} AA)", end=" ")
        
        # Find best matching reference
        query_dict = {seq_id: query_seq}
        
        # Skip MMseqs2 for faster analysis - just use direct alignment
        best_ref, best_score = None, -float('inf')
        
        for ref_name, ref_seq in ref_sequences.items():
            alignment = align_blosum62(query_seq, ref_seq, aligner)
            if alignment.score > best_score:
                best_score = alignment.score
                best_ref = ref_name
        
        print(f"-> {best_ref}", end=" ")
        
        # Perform pairwise alignment
        ref_seq = ref_sequences[best_ref]
        
        # Transfer GRN annotations using the reference from mo_ref
        if best_ref not in grn_proc.data.index:
            print(f"WARNING: Reference {best_ref} not found in GRN table")
            continue
            
        ref_row = grn_proc.data.loc[best_ref]
        
        # Create a mapping from sequence position to GRN
        ref_grn_dict = {grn: res for grn, res in ref_row.to_dict().items() if res not in ['-', 'X']}
        
        # Only show detailed info for first few sequences
        if i < 5:
            print(f"  Reference {best_ref} has {len(ref_grn_dict)} annotated positions")
        
        # Build position mapping from alignment
        # Extract just the amino acids from the reference GRN positions
        seq_pos2grn = {}
        
        # Build a sequence from just the annotated GRN positions
        ref_seq_from_grn = ""
        grn_positions = []
        
        for grn_pos in sorted(grn_proc.data.columns, key=lambda x: parse_grn_str2float(x)):
            residue = ref_row[grn_pos]
            if residue != '-' and isinstance(residue, str) and len(residue) > 0:
                # Check if it's a valid residue identifier
                if residue[0] in 'ACDEFGHIKLMNPQRSTVWY':
                    ref_seq_from_grn += residue[0]
                    grn_positions.append(grn_pos)
                    # Map position in this sequence to GRN position
                    seq_pos2grn[len(ref_seq_from_grn)] = grn_pos
        
        if not ref_seq_from_grn:
            print("WARNING: Could not extract reference sequence from GRN table")
            continue
            
        # Re-align with the GRN reference sequence (only annotated positions)
        alignment_grn = align_blosum62(query_seq, ref_seq_from_grn, aligner)
        formatted_grn = format_alignment(alignment_grn)
        
        # Initialize GRN row from GRN alignment
        new_row = init_row_from_alignment(formatted_grn, seq_pos2grn)
        
        print(f"  Transferred {len(new_row)} GRN positions")
        
        # For now, use the initial annotation without expansion
        # The expand_annotation function expects a different config format
        # TODO: Update to use proper expansion once config format is aligned
        final_grn = new_row if isinstance(new_row, pd.Series) else pd.Series(new_row)
        
        # Optional: Add neighboring positions based on microbial opsin intervals
        if len(final_grn) > 10:  # Only expand if we have reasonable initial coverage
            try:
                # Get microbial opsin intervals from config
                mo_config = config_manager.configs.get('config', {}).get('microbial_opsins', {})
                standard_intervals = mo_config.get('standard', {})
                
                # Check if we can add positions at the edges of TM regions
                expanded = final_grn.copy()
                for tm_name, (start_grn, end_grn) in standard_intervals.items():
                    if tm_name.startswith('tm'):
                        # Check if we have any positions in this TM
                        tm_num = tm_name[2]  # Extract TM number
                        tm_positions = [grn for grn in final_grn.index if grn.startswith(f"{tm_num}.")]
                        
                        if tm_positions:
                            # We have some positions in this TM, potentially expand edges
                            min_pos = min([parse_grn_str2float(p) for p in tm_positions])
                            max_pos = max([parse_grn_str2float(p) for p in tm_positions])
                            
                            # Add a few positions at the edges if they're within the standard range
                            start_float = parse_grn_str2float(start_grn)
                            end_float = parse_grn_str2float(end_grn)
                            
                            # This is a simplified expansion - in practice you'd check the alignment
                            print(f"    TM{tm_num}: positions {len(tm_positions)}, range {start_grn}-{end_grn}")
                
                final_grn = expanded
                
            except Exception as e:
                print(f"  Note: Simplified expansion not applied: {str(e)}")
        
        coverage = len(final_grn) / len(query_seq) * 100
        print(f"  Coverage: {coverage:.1f}% ({len(final_grn)}/{len(query_seq)} positions)")
        
        # Store results
        results[seq_id] = final_grn
        
        # Only show detailed output for first few sequences
        if i < 5:
            # Define conserved positions for microbial opsins
            conserved_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
            print("  Key GRN positions:")
            for pos in conserved_positions:
                if hasattr(final_grn, 'index') and pos in final_grn.index and final_grn[pos] != '-':
                    print(f"    {pos}: {final_grn[pos]}")
            
            # Show TM coverage
            tm_coverage = {}
            if hasattr(final_grn, 'index'):
                for grn in final_grn.index:
                    if '.' in grn and not grn.startswith('n.') and not grn.startswith('c.'):
                        try:
                            tm = int(grn.split('.')[0])
                            if 1 <= tm <= 7:
                                tm_coverage[tm] = tm_coverage.get(tm, 0) + 1
                        except:
                            pass
            
            tm_summary = ", ".join([f"TM{tm}:{count}" for tm, count in sorted(tm_coverage.items())])
            print(f"  TM coverage: {tm_summary}")
    
    # 5. Create output GRN table
    print("\n" + "="*60)
    print("5. Creating Output GRN Table")
    print("="*60)
    
    if results:
        # Combine all results into a DataFrame
        all_grns = set()
        for grn_series in results.values():
            all_grns.update(grn_series.index)
        
        # Create DataFrame with all GRN columns
        output_data = []
        for seq_id, grn_series in results.items():
            row_data = {grn: grn_series.get(grn, '-') for grn in all_grns}
            row_data['entity_id'] = seq_id
            output_data.append(row_data)
        
        output_df = pd.DataFrame(output_data)
        output_df = output_df.set_index('entity_id')
        
        # Sort columns - don't use sort_grns_str since it converts to x notation
        grn_cols = [col for col in all_grns if '.' in col]
        def sort_key(grn):
            try:
                parts = grn.split('.')
                if len(parts) == 2:
                    return (int(parts[0]), int(parts[1]))
                return (999, 999)
            except:
                return (999, 999)
        sorted_cols = sorted(grn_cols, key=sort_key)
        output_df = output_df[sorted_cols]
        
        # Save results with timestamp to avoid permission issues
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_name = f"opsin_grn_annotations_{timestamp}"
        grn_proc.data = output_df
        grn_proc.ids = output_df.index.tolist()
        grn_proc.grns = output_df.columns.tolist()
        
        saved_path = grn_proc.save_grn_table(output_name)
        print(f"✅ Saved GRN table to: {saved_path}")
        print(f"✅ Annotated {len(output_df)} sequences")
        
        # Summary statistics
        non_gap_counts = output_df.apply(lambda x: (x != '-').sum(), axis=1)
        
        print("\nAnnotation Summary:")
        print(f"  - Input sequences: {len(sequences)}")
        print(f"  - Processed: {len(results)}")
        print(f"  - Average positions per sequence: {non_gap_counts.mean():.1f}")
        print(f"  - Min positions: {non_gap_counts.min()}")
        print(f"  - Max positions: {non_gap_counts.max()}")
        
        # TM coverage summary
        all_tm_coverage = {tm: 0 for tm in range(1, 8)}
        for _, row in output_df.iterrows():
            for grn in output_df.columns:
                if '.' in grn and row[grn] != '-':
                    try:
                        tm = int(grn.split('.')[0])
                        if 1 <= tm <= 7:
                            all_tm_coverage[tm] += 1
                    except:
                        pass
        
        print("\nOverall TM coverage:")
        for tm in range(1, 8):
            avg_coverage = all_tm_coverage[tm] / len(output_df)
            print(f"  - TM{tm}: {avg_coverage:.1f} positions per sequence")
    
    print("\n" + "="*80)
    print("GRN ANNOTATION COMPLETE")
    print("="*80)




if __name__ == "__main__":
    main()