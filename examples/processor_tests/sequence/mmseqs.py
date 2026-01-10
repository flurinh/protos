#!/usr/bin/env python
"""
Convert MMseqs2 alignment output to dense alignment format like Biopython.

This script takes MMseqs2 alignment coordinates and sequences to produce:
- Full sequences with gaps inserted at appropriate positions
- Alignment midline with '|' for matches, '.' for mismatches
"""

import pandas as pd
from pathlib import Path
import os

# Protos imports
from protos.io.paths import ProtosPaths
from protos.processing.sequence.sequence_processor import SequenceProcessor
from protos.processing.grn.grn_processor import GRNProcessor


def mmseqs_to_dense_alignment(query_seq, target_seq, query_start, query_end, 
                             target_start, target_end, query_aln, target_aln):
    """
    Convert MMseqs2 local alignment to dense (global-like) alignment format.
    
    Args:
        query_seq: Full query sequence
        target_seq: Full target sequence
        query_start: 1-based start position in query
        query_end: 1-based end position in query
        target_start: 1-based start position in target
        target_end: 1-based end position in target
        query_aln: Aligned query sequence from MMseqs2
        target_aln: Aligned target sequence from MMseqs2
        
    Returns:
        dict with:
            - 'query_aligned': Full query sequence with gaps
            - 'target_aligned': Full target sequence with gaps
            - 'midline': Alignment midline with | and .
            - 'query_start': 1-based start in dense alignment
            - 'target_start': 1-based start in dense alignment
    """
    # Convert to 0-based indexing
    q_start = query_start - 1
    q_end = query_end
    t_start = target_start - 1
    t_end = target_end
    
    # Initialize dense alignment sequences
    query_dense = []
    target_dense = []
    midline = []
    
    # Add leading gaps/residues before alignment starts
    if q_start > 0 or t_start > 0:
        # Determine which sequence starts first
        if q_start < t_start:
            # Query starts first
            # Add query residues with gaps in target
            for i in range(q_start):
                query_dense.append(query_seq[i])
                target_dense.append('-')
                midline.append(' ')
            
            # Add remaining target residues with gaps in query
            for i in range(t_start - q_start):
                query_dense.append('-')
                target_dense.append(target_seq[i])
                midline.append(' ')
        else:
            # Target starts first
            # Add target residues with gaps in query
            for i in range(t_start):
                query_dense.append('-')
                target_dense.append(target_seq[i])
                midline.append(' ')
            
            # Add remaining query residues with gaps in target
            for i in range(q_start - t_start):
                query_dense.append(query_seq[i])
                target_dense.append('-')
                midline.append(' ')
    
    # Add the aligned region
    for q_char, t_char in zip(query_aln, target_aln):
        query_dense.append(q_char)
        target_dense.append(t_char)
        
        # Create midline
        if q_char == '-' or t_char == '-':
            midline.append(' ')
        elif q_char == t_char:
            midline.append('|')
        else:
            midline.append('.')
    
    # Add trailing gaps/residues after alignment ends
    q_remaining = len(query_seq) - q_end
    t_remaining = len(target_seq) - t_end
    
    if q_remaining > 0 or t_remaining > 0:
        # Add remaining residues
        max_remaining = max(q_remaining, t_remaining)
        
        for i in range(max_remaining):
            if i < q_remaining:
                query_dense.append(query_seq[q_end + i])
            else:
                query_dense.append('-')
            
            if i < t_remaining:
                target_dense.append(target_seq[t_end + i])
            else:
                target_dense.append('-')
            
            midline.append(' ')
    
    # Convert lists to strings
    query_aligned = ''.join(query_dense)
    target_aligned = ''.join(target_dense)
    midline_str = ''.join(midline)
    
    # Determine actual start positions in dense alignment (1-based)
    dense_query_start = 1
    dense_target_start = 1
    
    # Find first non-gap position for query
    for i, char in enumerate(query_aligned):
        if char != '-':
            dense_query_start = i + 1
            break
    
    # Find first non-gap position for target
    for i, char in enumerate(target_aligned):
        if char != '-':
            dense_target_start = i + 1
            break
    
    return {
        'query_aligned': query_aligned,
        'target_aligned': target_aligned,
        'midline': midline_str,
        'query_start': dense_query_start,
        'target_start': dense_target_start,
        'alignment_length': len(query_aligned)
    }


def format_alignment_display(alignment_dict, query_id='Query', target_id='Target', 
                           line_length=60):
    """
    Format alignment for display similar to Biopython's format.
    
    Args:
        alignment_dict: Output from mmseqs_to_dense_alignment
        query_id: Name of query sequence
        target_id: Name of target sequence
        line_length: Number of characters per line
        
    Returns:
        Formatted string for display
    """
    query_aln = alignment_dict['query_aligned']
    target_aln = alignment_dict['target_aligned']
    midline = alignment_dict['midline']
    
    lines = []
    
    # Process alignment in chunks
    query_pos = alignment_dict['query_start'] - 1
    target_pos = alignment_dict['target_start'] - 1
    
    for start in range(0, len(query_aln), line_length):
        end = min(start + line_length, len(query_aln))
        
        # Count non-gap residues in this chunk
        query_chunk = query_aln[start:end]
        target_chunk = target_aln[start:end]
        midline_chunk = midline[start:end]
        
        # Update positions
        query_residues = sum(1 for c in query_chunk if c != '-')
        target_residues = sum(1 for c in target_chunk if c != '-')
        
        # Format query line
        query_end = query_pos + query_residues
        lines.append(f"{query_id:<8} {query_pos+1:>6} {query_chunk} {query_end}")
        
        # Format midline
        lines.append(f"{'':8} {'':>6} {midline_chunk}")
        
        # Format target line
        target_end = target_pos + target_residues
        lines.append(f"{target_id:<8} {target_pos+1:>6} {target_chunk} {target_end}")
        
        lines.append("")  # Empty line between blocks
        
        # Update positions
        query_pos = query_end
        target_pos = target_end
    
    return '\n'.join(lines)


def test_alignment():
    """Test the alignment conversion with example data."""
    # Example from actual MMseqs2 output
    query_seq = "FVYAIVSMLFLMFNSFGVVQLVQMSKYVWCSRCCRCTDDNEKFNSAIELAYTVLSLVAKTLLGWMMYANVLA"
    target_seq = "LRRYNVVAAVVHAAQAVAVLAIATFLALSALFHVIVRWVEYSLSSSLMIVIIGLFGVNASMILFGWLQWLPFVFGCLAGAVPWVGIVFVYGIIVSLFLLFNVFALVQERLYITLSLVAKSSALA"
    
    # MMseqs2 alignment info
    query_start = 14
    query_end = 73
    target_start = 84
    target_end = 116
    query_aln = "FVYAIVSMLFLMFNSFGVVQLVQMSKYVWCSRCCRCTDDNEKFNSAIELAYTVLSLVAKT"
    target_aln = "FVYGIIISLFLLFNSFALVQ---------------------------ERAYIVLSLVAKS"
    
    # Convert to dense alignment
    result = mmseqs_to_dense_alignment(
        query_seq, target_seq,
        query_start, query_end,
        target_start, target_end,
        query_aln, target_aln
    )
    
    # Display formatted alignment
    formatted = format_alignment_display(result, 'MO_024', '6SU3')
    print("Dense Alignment Format:")
    print("=" * 80)
    print(formatted)
    
    # Show summary
    print("\nAlignment Summary:")
    print(f"Query length: {len(query_seq)}")
    print(f"Target length: {len(target_seq)}")
    print(f"Alignment length: {result['alignment_length']}")
    print(f"Query coverage: {(query_end - query_start + 1) / len(query_seq):.1%}")
    print(f"Target coverage: {(target_end - target_start + 1) / len(target_seq):.1%}")
    
    # Calculate identity
    matches = sum(1 for c in result['midline'] if c == '|')
    aligned_positions = sum(1 for a, b in zip(query_aln, target_aln) if a != '-' and b != '-')
    print(f"Sequence identity: {matches}/{aligned_positions} ({matches/aligned_positions:.1%})")


def process_mmseqs_results(alignment_file):
    """
    Process MMseqs2 alignment results with full sequences.
    
    Args:
        alignment_file: Path to MMseqs2 alignment TSV file with qaln,taln columns
        
    Returns:
        List of alignment dictionaries
    """
    # Load alignments
    column_names = ['query_id', 'target_id', 'sequence_identity', 'alignment_length', 
                    'mismatches', 'gap_opens', 'query_start', 'query_end', 
                    'target_start', 'target_end', 'e_value', 'bit_score', 'qaln', 'taln']
    
    df = pd.read_csv(alignment_file, sep='\t', header=None, names=column_names)
    
    print(f"Loaded {len(df)} alignments")
    
    # For each alignment, we would need the full sequences
    # This is just a placeholder - in practice you'd load these from FASTA files
    alignments = []
    
    for _, row in df.iterrows():
        # Here you would look up the full sequences
        # For now, we'll just store the alignment info
        alignments.append({
            'query_id': row['query_id'],
            'target_id': row['target_id'],
            'query_aln': row['qaln'],
            'target_aln': row['taln'],
            'query_start': row['query_start'],
            'query_end': row['query_end'],
            'target_start': row['target_start'],
            'target_end': row['target_end'],
            'e_value': row['e_value'],
            'identity': row['sequence_identity']
        })
    
    return alignments


if __name__ == "__main__":
    # Run test with example data
    test_alignment()
    
    # If alignment file exists, process it
    alignment_file = Path("temp/alignment_results.tsv")
    if alignment_file.exists():
        print(f"\n\nProcessing {alignment_file}...")
        alignments = process_mmseqs_results(alignment_file)
        print(f"Found {len(alignments)} alignments")
        
        # Show first few alignments
        for i, aln in enumerate(alignments[:3]):
            print(f"\nAlignment {i+1}: {aln['query_id']} vs {aln['target_id']}")
            print(f"Identity: {aln['identity']:.1%}, E-value: {aln['e_value']}")