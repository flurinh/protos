"""
CLI command for assigning GRNs to sequences using the entity system.

This command now integrates with the entity registry for tracking sequences
and their GRN annotations.
"""

import argparse
import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple, Optional

from protos.processing.grn import GRNProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.grn.grn_table_utils import (
    GRNConfigManager, init_grn_intervals, init_row_from_alignment,
    expand_annotation, is_sequential
)
from protos.processing.grn.grn_utils import sort_grns_str
from protos.processing.sequence.seq_alignment import (
    init_aligner, align_blosum62, format_alignment, mmseqs2_align2
)
from protos.io.data_access import GlobalRegistry


def get_pairwise_alignment(query_seq_dict, ref_seq_dict, best_matches, aligner=None):
    """
    Get pairwise alignment between query sequences and reference sequences.
    
    Args:
        query_seq_dict: Dictionary of query sequences
        ref_seq_dict: Dictionary of reference sequences
        best_matches: List of [query_id, ref_id] pairs to align
        aligner: Alignment object (created if None)
    
    Returns:
        Dictionary of alignments keyed by query_id
    """
    if aligner is None:
        aligner = init_aligner(open_gap_score=-22)
    
    alignment_dict = {}
    for query_id, ref_id in best_matches:
        print(f"Aligning {query_id} to {ref_id}")
        try:
            query_seq = query_seq_dict[query_id].replace('\n', '')
            ref_seq = ref_seq_dict[ref_id].replace('\n', '')
            alm = align_blosum62(query_seq, ref_seq, aligner)
            alignment_dict[query_id] = format_alignment(alm)
        except Exception as e:
            print(f"Error aligning {query_id}: {e}")
    
    return alignment_dict


def get_aligned_grns(grnp, alignments, best_matches, grns_str_strict):
    """
    Get aligned GRNs based on alignments.
    
    Args:
        grnp: GRN processor instance
        alignments: Dictionary of alignments
        best_matches: List of [query_id, ref_id] pairs
        grns_str_strict: List of strict GRNs to filter by
    
    Returns:
        Dictionary of aligned GRNs keyed by query_id
    """
    new_rows = {}
    for query_id, ref_id in best_matches:
        try:
            alignment = alignments[query_id]
            ref_row = grnp.data.loc[ref_id, :]
            ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
            seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
            new_row = init_row_from_alignment(alignment, seq_pos2grn)
            
            # Filter for single TM GRNs (both x and dot notation)
            new_row_grns = []
            for grn in grns_str_strict:
                if 'x' in grn:
                    tm_part = grn.split('x')[0]
                elif '.' in grn:
                    tm_part = grn.split('.')[0]
                else:
                    continue
                    
                if len(tm_part) < 2 and grn in new_row.index.tolist():
                    new_row_grns.append(grn)
            
            new_row_strict = new_row.loc[new_row_grns]
            if len(new_row_strict) == 0:
                print(f"{query_id}: No strict GRNs found in alignment.")
            new_rows[query_id] = new_row_strict
        except Exception as e:
            print(f"Error initializing row of sequence {query_id}: {e}")
    
    return new_rows


def annotate_sequence(query_id, query_seq, new_row, protein_family):
    """
    Annotate a sequence with GRNs.
    
    Args:
        query_id: Query sequence ID
        query_seq: Query sequence
        new_row: New row with GRN annotations
        protein_family: Protein family name
    
    Returns:
        Tuple of (query_id, GRN dictionary) or (None, None) if error
    """
    try:
        new_row_seq = ''.join([x[0] for x in new_row]).replace('-', '')
        
        alignment = align_blosum62(query_seq, new_row_seq, init_aligner(), verbose=0)
        alignment = format_alignment(alignment)
        grns, rns, missing = expand_annotation(
            new_row, query_seq, alignment,
            protein_family=protein_family,
            max_alignment_gap=1, verbose=0
        )
        
        if len(missing) == 0:
            grn_dict = dict(zip(rns, grns))
            if is_sequential(grn_dict):
                return query_id, dict(zip(grns, rns))
        else:
            print(f"Missing GRNs in {query_id}: {missing}")
            return query_id, dict(zip(grns, rns))
    except Exception as e:
        print(f"Error in processing {query_id}: {str(e)}")
        return None, None


def parallel_annotate_sequences(query_dict, new_rows, protein_family, num_cores=8):
    """
    Execute the annotation process in parallel.
    
    Args:
        query_dict: Dictionary of query sequences
        new_rows: Dictionary of new rows with aligned GRNs
        protein_family: Protein family name
        num_cores: Number of CPU cores to use
    
    Returns:
        Dictionary of GRN annotations keyed by query_id
    """
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [
            executor.submit(annotate_sequence, query_id, query_dict[query_id], 
                          new_rows[query_id], protein_family)
            for query_id in query_dict if query_id in new_rows
        ]
        
        results = {}
        for future in futures:
            try:
                query_id, grn_dict = future.result()
                if grn_dict is not None:
                    results[query_id] = grn_dict
            except Exception as e:
                print(f"Error in processing: {str(e)}")
        
        return results


def assign_grns_to_entities(
    protein_family: str,
    sequence_ids: Optional[List[str]] = None,
    reference_table: str = None,
    num_cores: int = 8,
    output_name: Optional[str] = None
) -> Optional[pd.DataFrame]:
    """
    Assign GRNs to sequences using the entity system.
    
    Args:
        protein_family: Protein family name ('gpcr_a', 'microbial_opsins', etc.)
        sequence_ids: Specific sequence IDs to process (None = all registered sequences)
        reference_table: Reference GRN table to use (None = use default for family)
        num_cores: Number of CPU cores for parallel processing
        output_name: Name for output GRN table (None = auto-generate)
    
    Returns:
        DataFrame with GRN annotations or None if failed
    """
    # Initialize processors
    seq_proc = SeqProcessor(name="grn_assignment")
    grn_proc = GRNProcessor(name="grn_assignment")
    
    # Get sequences to process
    if sequence_ids:
        # Specific sequences requested
        query_dict = {}
        for seq_id in sequence_ids:
            sequence = seq_proc.load_sequence_entity(seq_id)
            if sequence:
                query_dict[seq_id] = sequence
            else:
                print(f"Warning: Sequence {seq_id} not found")
    else:
        # Get all registered sequences
        registry = GlobalRegistry()
        query_dict = {}
        
        for entity_id, entity_data in registry.entity_registry.entities.items():
            if 'sequence' in entity_data.get('formats', {}):
                original_id = entity_data['original_id']
                sequence = seq_proc.load_sequence_entity(original_id)
                if sequence:
                    query_dict[original_id] = sequence
    
    if not query_dict:
        print("No sequences found to process")
        return None
    
    print(f"Processing {len(query_dict)} sequences")
    
    # Load reference GRN table
    if reference_table is None:
        # Use default reference for protein family
        if protein_family == 'gpcr_a':
            reference_table = 'ref'
        else:
            reference_table = 'mo_ref'
    
    try:
        grn_proc.load_grn_table(reference_table)
    except Exception as e:
        print(f"Error loading reference table '{reference_table}': {e}")
        return None
    
    # Get reference sequences
    ref_dict = {k: v.replace('-', '') for k, v in grn_proc.get_seq_dict().items()}
    
    # Clean query sequences
    query_dict = {k: v.replace('-', '') for k, v in query_dict.items()}
    
    # Get GRN configuration
    config = GRNConfigManager(protein_family=protein_family)
    grn_config_strict = config.get_config(strict=True)
    grns_str_strict = init_grn_intervals(grn_config_strict)
    
    # Find best matches using MMseqs2
    print("Finding best matches...")
    output = mmseqs2_align2(query_seqs=query_dict, ref_seqs=ref_dict)
    best_matches = output.loc[output.groupby('query_id')['e_value'].idxmin()][
        ['query_id', 'target_id']
    ].values.tolist()
    
    # Get pairwise alignments
    print("Performing pairwise alignments...")
    alignments = get_pairwise_alignment(query_dict, ref_dict, best_matches=best_matches)
    alignments = {k: v for k, v in alignments.items() if len(v) > 0}
    print(f"Number of sequences with alignments: {len(alignments)}")
    
    # Get aligned GRNs
    best_matches = [[x, y] for [x, y] in best_matches if x in alignments.keys()]
    new_rows = get_aligned_grns(grn_proc, alignments, best_matches, grns_str_strict)
    print(f"Number of sequences with aligned GRNs: {len(new_rows)}")
    
    # Filter query dict to only sequences with aligned GRNs
    query_dict = {k: v for k, v in query_dict.items() if k in new_rows.keys()}
    
    # Annotate sequences in parallel
    print(f"Annotating sequences using {num_cores} cores...")
    final_results = parallel_annotate_sequences(
        query_dict, new_rows, 
        protein_family=protein_family, 
        num_cores=num_cores
    )
    
    print(f"Number of sequences with GRNs: {len(final_results)}")
    
    if not final_results:
        print("No sequences were successfully processed.")
        return None
    
    # Create DataFrame
    df = pd.DataFrame.from_dict(final_results, orient='index')
    df = df.loc[:, sort_grns_str(df.columns.tolist())]
    
    # Filter based on Schiff base (crude check)
    if protein_family == 'microbial_opsins':
        # Check for lysine at position 7.50 (or 7x50 for old notation)
        if '7.50' in df.columns:
            df = df[df['7.50'].str.contains('K')]
        elif '7x50' in df.columns:
            df = df[df['7x50'].str.contains('K')]
    elif protein_family == 'gpcr_a':
        # Check for position 7.43 (or 7x43)
        if '7.43' in df.columns:
            df = df[df['7.43'].str.contains('K')]
        elif '7x43' in df.columns:
            df = df[df['7x43'].str.contains('K')]
    
    # Save the GRN table with entity IDs
    if output_name is None:
        output_name = f"grn_{protein_family}_{len(df)}_sequences"
    
    # Reset the processor data and save
    grn_proc.data = df
    grn_proc.ids = df.index.tolist()
    grn_proc.grns = df.columns.tolist()
    
    saved_path = grn_proc.save_grn_table(output_name, include_entity_ids=True)
    print(f"Saved GRN table to: {saved_path}")
    
    return df


def main():
    """Command-line entry point for GRN assignment."""
    parser = argparse.ArgumentParser(
        description='Assign GRN annotations to sequences using the entity system',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Assign GRNs to all registered microbial opsin sequences
  protos grn assign -p microbial_opsins
  
  # Assign GRNs to specific sequences
  protos grn assign -p gpcr_a -s P12345 P67890
  
  # Use custom reference table and output name
  protos grn assign -p microbial_opsins --reference my_ref_table -o my_results
  
  # Use more cores for parallel processing
  protos grn assign -p microbial_opsins -n 16
        """
    )
    
    parser.add_argument(
        '-p', '--protein-family',
        required=True,
        choices=['gpcr_a', 'microbial_opsins', 'gpcr_b', 'gpcr_c'],
        help='Protein family'
    )
    
    parser.add_argument(
        '-s', '--sequences',
        nargs='+',
        help='Specific sequence IDs to process (default: all registered sequences)'
    )
    
    parser.add_argument(
        '-r', '--reference',
        help='Reference GRN table name (default: family-specific default)'
    )
    
    parser.add_argument(
        '-n', '--num-cores',
        type=int,
        default=8,
        help='Number of cores for parallel processing (default: 8)'
    )
    
    parser.add_argument(
        '-o', '--output',
        help='Output table name (default: auto-generated)'
    )
    
    args = parser.parse_args()
    
    # Assign GRNs
    result = assign_grns_to_entities(
        protein_family=args.protein_family,
        sequence_ids=args.sequences,
        reference_table=args.reference,
        num_cores=args.num_cores,
        output_name=args.output
    )
    
    if result is not None:
        print(f"\nSuccess! Assigned GRNs to {len(result)} sequences")
        print(f"Columns: {', '.join(result.columns[:10])}{'...' if len(result.columns) > 10 else ''}")
    else:
        print("\nGRN assignment failed.")
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())