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
    # We need to handle both sequences completely
    
    # First, determine the total length needed for the N-terminal region
    # This should be the maximum of the two starting positions
    n_terminal_length = max(q_start, t_start)
    
    # Fill in the N-terminal region
    for i in range(n_terminal_length):
        # Add query residue or gap
        if i < q_start:
            query_dense.append(query_seq[i])
        else:
            query_dense.append('-')
        
        # Add target residue or gap  
        if i < t_start:
            target_dense.append(target_seq[i])
        else:
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
    query_pos = 0  # Track position in original sequence
    target_pos = 0
    
    # Count initial positions
    for i in range(alignment_dict['query_start'] - 1):
        if query_aln[i] != '-':
            query_pos += 1
    for i in range(alignment_dict['target_start'] - 1):
        if target_aln[i] != '-':
            target_pos += 1
    
    for start in range(0, len(query_aln), line_length):
        end = min(start + line_length, len(query_aln))
        
        # Get chunks
        query_chunk = query_aln[start:end]
        target_chunk = target_aln[start:end]
        midline_chunk = midline[start:end]
        
        # Count non-gap residues in this chunk
        query_residues = sum(1 for c in query_chunk if c != '-')
        target_residues = sum(1 for c in target_chunk if c != '-')
        
        # Format query line
        if query_residues > 0:
            query_start_pos = query_pos + 1
            query_pos += query_residues
            lines.append(f"{query_id:<8} {query_start_pos:>6} {query_chunk} {query_pos}")
        else:
            lines.append(f"{query_id:<8} {'':>6} {query_chunk}")
        
        # Format midline
        lines.append(f"{'':8} {'':>6} {midline_chunk}")
        
        # Format target line
        if target_residues > 0:
            target_start_pos = target_pos + 1
            target_pos += target_residues
            lines.append(f"{target_id:<8} {target_start_pos:>6} {target_chunk} {target_pos}")
        else:
            lines.append(f"{target_id:<8} {'':>6} {target_chunk}")
        
        lines.append("")  # Empty line between blocks
    
    return '\n'.join(lines)


def load_sequences():
    """Load query and reference sequences using Protos."""
    # Setup paths
    datadir = Path(__file__).parent.absolute()
    test_data_root = datadir / "data"
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Initialize processors
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    
    # Load query sequences (mo_small dataset)
    print("Loading query sequences...")
    query_sequences = seq_proc.load_entity('mo_small')
    print(f"Loaded {len(query_sequences)} query sequences")
    
    # Load reference sequences from GRN table
    print("Loading reference sequences...")
    ref_file = grn_proc.path_ref_dir / "mo_ref.csv"
    grn_proc.data = pd.read_csv(ref_file, index_col=0)
    grn_proc.data = grn_proc.data.fillna('-')
    
    # Build reference sequences
    ref_sequences = {}
    for idx in grn_proc.data.index:
        row = grn_proc.data.loc[idx]
        seq_parts = []
        for val in row.values:
            if val not in ['-', 'X', '.', ''] and not pd.isna(val):
                if isinstance(val, str) and len(val) > 0:
                    seq_parts.append(val[0])
        seq = ''.join(seq_parts)
        ref_sequences[idx] = seq
    
    print(f"Built {len(ref_sequences)} reference sequences")
    return query_sequences, ref_sequences


def test_alignment_with_protos():
    """Test the alignment conversion with real sequences from Protos."""
    # Load real sequences
    try:
        query_sequences, ref_sequences = load_sequences()
        
        # Use MO_024 and 6SU3 as example
        query_seq = query_sequences.get('MO_024', '')
        target_seq = ref_sequences.get('6SU3', '')
        
        if query_seq and target_seq:
            print(f"\nUsing Protos sequences:")
            print(f"Query (MO_024): {len(query_seq)} aa")
            print(f"Target (6SU3): {len(target_seq)} aa")
        else:
            raise ValueError("Could not load sequences")
            
    except Exception as e:
        print(f"Could not load from Protos: {e}")
        # Fallback to hardcoded sequences
        query_seq = "IAIAATVLHFLNVVGMLAIYWGGFDGETTRDSCFQTTVSYLKWTQIAEGANTEGLVTFNTSSGLYGVKHGTVAGNTLRLHWLIVAFHALSFLFQGVVLIPGWYRYENRVKEGSNPMRFIEYSISASIMLVSVALISGIIDENELITISVLCGATQMCGLVSETIVSNVEKLKNEVRGDIRCIVHNLHIAATVAHLSGWVMVMVGYGVIWRYFILSTSESYSSPPEFVYAIVSMLFLMFNSFGVVQLVQMSKYVWCSRCCRCTDDNEKFNSAIELAYTVLSLVAKTLLGWMMYANVLA"
        target_seq = "LQNFNRIAGVFHLLQMAVALFLGLSALFHFIVRWVEYSLSSSVMIVLIAIFGVNASMILFGWLQLLPFWFGCIAGIVPWIGLLFVYGIIISLFLLFNSFALVQERAYIVLSLVAKSSALA"
        print("\nUsing hardcoded sequences")
        print(f"Query (MO_024): {len(query_seq)} aa")
        print(f"Target (6SU3): {len(target_seq)} aa")
    
    # MMseqs2 alignment info from actual output
    query_start = 226
    query_end = 285
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
    print("\nDense Alignment Format:")
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
    
    # Show alignment positions
    print(f"\nMMseqs2 alignment coordinates:")
    print(f"Query: {query_start}-{query_end} (aligned region: {query_aln[:20]}...)")
    print(f"Target: {target_start}-{target_end} (aligned region: {target_aln[:20]}...)")


def process_mmseqs_results(alignment_file, query_sequences=None, ref_sequences=None):
    """
    Process MMseqs2 alignment results and create dense alignments.
    
    Args:
        alignment_file: Path to MMseqs2 alignment TSV file with qaln,taln columns
        query_sequences: Dict of query sequences (optional)
        ref_sequences: Dict of reference sequences (optional)
        
    Returns:
        List of alignment dictionaries with dense alignments
    """
    # Load alignments
    column_names = ['query_id', 'target_id', 'sequence_identity', 'alignment_length', 
                    'mismatches', 'gap_opens', 'query_start', 'query_end', 
                    'target_start', 'target_end', 'e_value', 'bit_score', 'qaln', 'taln']
    
    df = pd.read_csv(alignment_file, sep='\t', header=None, names=column_names)
    print(f"Loaded {len(df)} alignments")
    
    # Load sequences if not provided
    if query_sequences is None or ref_sequences is None:
        print("Loading sequences...")
        query_sequences, ref_sequences = load_sequences()
    
    # Process alignments
    dense_alignments = []
    
    for _, row in df.iterrows():
        query_id = row['query_id']
        target_id = row['target_id']
        
        # Get full sequences
        query_seq = query_sequences.get(query_id, '')
        target_seq = ref_sequences.get(target_id, '')
        
        if query_seq and target_seq:
            # Convert to dense alignment
            dense_aln = mmseqs_to_dense_alignment(
                query_seq, target_seq,
                row['query_start'], row['query_end'],
                row['target_start'], row['target_end'],
                row['qaln'], row['taln']
            )
            
            dense_aln['query_id'] = query_id
            dense_aln['target_id'] = target_id
            dense_aln['identity'] = row['sequence_identity']
            dense_aln['e_value'] = row['e_value']
            
            dense_alignments.append(dense_aln)
    
    return dense_alignments


if __name__ == "__main__":
    # Test with real sequences
    print("Testing MMseqs2 to dense alignment conversion")
    print("-" * 80)
    test_alignment_with_protos()
    
    # Process alignment file if it exists
    alignment_file = Path("temp/alignment_results.tsv")
    if alignment_file.exists():
        print(f"\n\nProcessing {alignment_file}...")
        alignments = process_mmseqs_results(alignment_file)
        print(f"Created {len(alignments)} dense alignments")
        
        # Show first alignment
        if alignments:
            aln = alignments[0]
            print(f"\nFirst alignment: {aln['query_id']} vs {aln['target_id']}")
            formatted = format_alignment_display(aln, aln['query_id'], aln['target_id'])
            print(formatted[:500] + "..." if len(formatted) > 500 else formatted)