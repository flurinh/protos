"""
Direct GRN Annotation Example
=============================

This script annotates a test sequence using the GRN system,
working directly with the data to avoid notation conversion issues.

Usage:
    python examples/grn.py
"""

import os
import sys
import pandas as pd
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from protos.processing.grn.grn_table_utils import (
    annotate_gpcr,
    init_row_from_alignment,
    expand_annotation,
    GRNConfigManager,
    init_grn_intervals
)
from protos.processing.grn.grn_utils import get_seq, sort_grns_str
from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    format_alignment,
    mmseqs2_align
)
from protos.io.fasta_utils import read_fasta


class DirectGRNAnnotator:
    """Direct GRN annotation without GRNProcessor to avoid notation issues."""
    
    def __init__(self):
        self.project_root = Path(__file__).parent.parent
        self.ref_data_dir = self.project_root / "src" / "protos" / "reference_data"
        self.test_data_dir = self.project_root / "tests" / "test-data"
        
    def load_reference_table(self):
        """Load reference GRN table directly."""
        grn_file = self.ref_data_dir / "grn" / "ref" / "mo_ref.csv"
        print(f"Loading reference from: {grn_file}")
        
        df = pd.read_csv(grn_file, index_col=0)
        df = df.fillna('-')
        
        print(f"Loaded {len(df)} reference proteins")
        print(f"Columns: {len(df.columns)} GRN positions")
        
        return df
    
    def get_sequences_from_table(self, grn_table):
        """Extract sequences from GRN table."""
        sequences = {}
        
        for protein_id in grn_table.index:
            sequence = get_seq(protein_id, grn_table)
            sequences[protein_id] = sequence
            
        return sequences
    
    def simple_annotation_demo(self, query_seq, ref_table):
        """Demonstrate simple GRN annotation process."""
        print("\n3. Simple Annotation Demonstration")
        print("-" * 40)
        
        # Find a good reference - use Bacteriorhodopsin if available
        if 'BR' in ref_table.index:
            ref_id = 'BR'
        else:
            ref_id = ref_table.index[0]
        
        print(f"Using reference: {ref_id}")
        ref_seq = get_seq(ref_id, ref_table)
        print(f"Reference length: {len(ref_seq)}")
        
        # Perform alignment
        print("\nPerforming pairwise alignment...")
        aligner = init_aligner()
        alignment = align_blosum62(query_seq, ref_seq, aligner)
        formatted = format_alignment(alignment)
        
        print("Alignment preview:")
        print(f"Query:  {formatted[0][:60]}...")
        print(f"        {formatted[1][:60]}...")
        print(f"Ref:    {formatted[2][:60]}...")
        
        # Create initial GRN mapping from alignment
        ref_row = ref_table.loc[ref_id]
        ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
        seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
        
        # Initialize row from alignment
        new_row = init_row_from_alignment(formatted, seq_pos2grn)
        
        print(f"\nInitial GRN assignments: {len(new_row)} positions")
        
        # Show some assignments
        print("\nExample GRN assignments:")
        count = 0
        for grn, res in new_row.items():
            if count < 10 and res != '-':
                print(f"  {grn}: {res}")
                count += 1
        
        return new_row
    
    def create_full_annotation(self, query_name, query_seq, ref_table):
        """Create full GRN annotation for a sequence."""
        print("\n4. Full Annotation Process")
        print("-" * 40)
        
        # Get reference sequences
        ref_sequences = self.get_sequences_from_table(ref_table)
        
        # Try to use MMseqs2 for best match
        try:
            print("Searching for best reference match...")
            hits = mmseqs2_align(query_seq, ref_sequences)
            best_match = hits['target_id'].iloc[0]
            print(f"Best match: {best_match}")
        except Exception as e:
            print(f"MMseqs2 not available: {e}")
            print("Using BR as reference...")
            best_match = 'BR' if 'BR' in ref_table.index else ref_table.index[0]
        
        # Get reference sequence
        ref_seq = ref_sequences[best_match]
        
        # Align
        aligner = init_aligner()
        alignment = align_blosum62(query_seq, ref_seq, aligner)
        formatted = format_alignment(alignment)
        
        # Create initial annotation
        ref_row = ref_table.loc[best_match]
        ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
        seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
        
        new_row = init_row_from_alignment(formatted, seq_pos2grn)
        
        # Get strict GRN positions from config
        config = GRNConfigManager(protein_family='microbial_opsins')
        grn_config_strict = config.get_config(strict=True)
        
        # Filter to keep only strict GRNs that are present in new_row
        # Note: Don't use init_grn_intervals as it may generate positions not in our data
        strict_grns_in_row = []
        for grn in new_row.index:
            # Handle both 'x' and '.' notation
            if 'x' in grn:
                parts = grn.split('x')
            elif '.' in grn:
                parts = grn.split('.')
            else:
                continue
                
            if len(parts[0]) == 1:
                # Check if this GRN falls within strict boundaries
                tm_num = parts[0]
                grn_num = int(parts[1])
                
                # Check each TM region
                for tm_key, bounds in grn_config_strict.items():
                    if bounds and len(bounds) == 2:
                        # Extract TM number from config key (e.g., 'tm1' -> '1')
                        if tm_key.startswith('tm'):
                            config_tm = tm_key[2:]
                            if config_tm == tm_num:
                                # Check if GRN is within bounds
                                # Handle both notations for bounds
                                if 'x' in bounds[0]:
                                    start_num = int(bounds[0].split('x')[1])
                                elif '.' in bounds[0]:
                                    start_num = int(bounds[0].split('.')[1])
                                else:
                                    start_num = 0
                                    
                                if 'x' in bounds[1]:
                                    end_num = int(bounds[1].split('x')[1])
                                elif '.' in bounds[1]:
                                    end_num = int(bounds[1].split('.')[1])
                                else:
                                    end_num = 100
                                if start_num <= grn_num <= end_num:
                                    strict_grns_in_row.append(grn)
                                    break
        
        if strict_grns_in_row:
            new_row_strict = new_row[strict_grns_in_row]
        else:
            # If no strict GRNs found, keep all for now
            new_row_strict = new_row
            
        print(f"Strict GRN positions: {len(new_row_strict)}")
        
        # Expand annotation
        try:
            print("\nExpanding annotation...")
            # Get sequence from strict row
            new_row_seq = ''.join([x[0] for x in new_row_strict.tolist() if x != '-']).replace('-', '')
            
            if len(new_row_seq) > 0:
                alignment2 = align_blosum62(query_seq, new_row_seq, aligner)
                formatted2 = format_alignment(alignment2)
                
                grn_list, rn_list, missing = expand_annotation(
                    new_row_strict,  # Use strict row
                    query_seq,
                    formatted2,
                    max_alignment_gap=1,
                    protein_family='microbial_opsins',
                    verbose=0
                )
            else:
                # If no sequence from strict, use full annotation
                print("No strict sequence found, using full annotation...")
                grn_list, rn_list, missing = expand_annotation(
                    new_row,  # Use full row
                    query_seq,
                    formatted,  # Use original alignment
                    max_alignment_gap=1,
                    protein_family='microbial_opsins',
                    verbose=0
                )
            
            print(f"Expanded to {len(grn_list)} positions")
            if missing:
                print(f"Missing residues: {len(missing)}")
            
            # Create final row
            final_row = pd.Series(dict(zip(grn_list, rn_list)))
            
            # Ensure all columns from reference table are present
            for col in ref_table.columns:
                if col not in final_row.index:
                    final_row[col] = '-'
            
            # Reorder to match reference table
            final_row = final_row[ref_table.columns]
            
            return final_row
            
        except Exception as e:
            print(f"Error in expansion: {e}")
            # Return initial annotation
            for col in ref_table.columns:
                if col not in new_row.index:
                    new_row[col] = '-'
            return new_row[ref_table.columns]


def main():
    """Run direct GRN annotation."""
    print("\n" + "="*60)
    print("DIRECT GRN ANNOTATION")
    print("="*60)
    
    annotator = DirectGRNAnnotator()
    
    # Load sequences from the specified file
    input_fasta = annotator.project_root / "input" / "opsin_sequences_from_yaml.fasta"
    print(f"\n1. Loading sequences from: {input_fasta}")
    
    if not input_fasta.exists():
        # Fallback to test sequence
        input_fasta = annotator.test_data_dir / "sequence" / "fasta" / "test_mo.fasta"
        print(f"File not found, using test sequence from: {input_fasta}")
    
    sequences = read_fasta(str(input_fasta))
    print(f"Loaded {len(sequences)} sequences")
    
    # Show all sequences
    for seq_id, seq in sequences.items():
        print(f"\n{seq_id}: {len(seq)} AA")
    
    # Load reference table
    print("\n2. Loading reference GRN table...")
    ref_table = annotator.load_reference_table()
    
    # Process each sequence
    all_annotations = {}
    key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.48', '6.50', '7.50']
    
    for seq_id, query_seq in sequences.items():
        print(f"\n{'='*60}")
        print(f"Processing: {seq_id}")
        print(f"{'='*60}")
        
        # Skip if sequence is too short
        if len(query_seq) < 50:
            print(f"Skipping {seq_id} - sequence too short ({len(query_seq)} AA)")
            continue
            
        # Simple annotation demo for first sequence only
        if seq_id == list(sequences.keys())[0]:
            initial_annotation = annotator.simple_annotation_demo(query_seq, ref_table)
        
        # Full annotation
        try:
            full_annotation = annotator.create_full_annotation(seq_id, query_seq, ref_table)
            all_annotations[seq_id] = full_annotation
            
            # Display key positions for this sequence
            print(f"\nKey positions for {seq_id}:")
            for pos in key_positions:
                if pos in full_annotation.index:
                    value = full_annotation[pos]
                    if value != '-':
                        print(f"  {pos}: {value}")
                        
        except Exception as e:
            print(f"Error annotating {seq_id}: {e}")
            continue
    
    # Create output table with only new annotations
    print("\n5. Creating output table with newly annotated sequences...")
    
    # Create DataFrame with only the new annotations
    if all_annotations:
        output_df = pd.DataFrame(all_annotations).T
        
        # Ensure columns are in the same order as reference table
        output_df = output_df.reindex(columns=ref_table.columns, fill_value='-')
    else:
        # Empty DataFrame with same columns as reference
        output_df = pd.DataFrame(columns=ref_table.columns)
    
    # Save only the new annotations to proper GRN tables location
    output_dir = annotator.project_root / "src/protos/reference_data/grn/tables"
    output_dir.mkdir(exist_ok=True, parents=True)  # Create directory if it doesn't exist
    output_file = output_dir / "grn_annotated_opsins_new_only_v2.csv"
    output_df.to_csv(output_file)
    print(f"Saved to: {output_file}")
    
    # Summary of annotations
    print("\n6. Summary of all annotations:")
    print("-" * 80)
    print(f"{'Sequence ID':<20} {'1x50':<6} {'2x50':<6} {'3x50':<6} {'7x50 (Schiff)':<15}")
    print("-" * 80)
    
    for seq_id, annotation in all_annotations.items():
        values = []
        for pos in ['1.50', '2.50', '3.50', '7.50']:
            if pos in annotation.index and annotation[pos] != '-':
                values.append(annotation[pos][:6])
            else:
                values.append('-')
        print(f"{seq_id:<20} {values[0]:<6} {values[1]:<6} {values[2]:<6} {values[3]:<15}")
    
    print("\n" + "="*60)
    print("BATCH ANNOTATION COMPLETE")
    print("="*60)
    print(f"\nSummary:")
    print(f"  - Input sequences: {len(sequences)}")
    print(f"  - Successfully annotated: {len(all_annotations)}")
    print(f"  - Reference proteins used: {len(ref_table)}")
    print(f"  - New sequences in output: {len(output_df)}")
    print(f"  - Output file: {output_file}")
    print(f"\nNote: Output contains ONLY newly annotated sequences, not reference proteins")


if __name__ == "__main__":
    main()