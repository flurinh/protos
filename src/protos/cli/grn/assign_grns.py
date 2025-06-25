import argparse
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from pathlib import Path
import os
from protos.processing.grn.grn_table_utils import *
from protos.io.fasta_utils import *
from protos.processing.sequence.seq_alignment import init_aligner, align_blosum62, format_alignment, mmseqs2_align2
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.io.paths.path_config import ProtosPaths


def get_pairwise_alignment(query_seq_dict, ref_seq_dict, best_matches, aligner=None):
    """
    Get pairwise alignment between query sequences and reference sequences.
    
    Args:
        query_seq_dict (dict): Dictionary of query sequences
        ref_seq_dict (dict): Dictionary of reference sequences
        best_matches (list): List of [query_id, ref_id] pairs to align
        aligner (object, optional): Alignment object. Defaults to None.
    
    Returns:
        dict: Dictionary of alignments keyed by query_id
    """
    if aligner is None:
        aligner = init_aligner(open_gap_score=-22)
    alignment_dict = {}
    for query_id, ref_id in best_matches:
        print(query_id, ref_id)
        try:
            query_seq = query_seq_dict[query_id].replace('\n', '')
            ref_seq = ref_seq_dict[ref_id].replace('\n', '')
            alm = align_blosum62(query_seq, ref_seq, aligner)
            alignment_dict[query_id] = format_alignment(alm)
        except Exception as e:
            print(f"Error aligning {query_id}: {query_seq_dict[query_id]}: {e}")
    return alignment_dict


def get_aligned_grns(grnp, alignments, best_matches, grns_str_strict):
    """
    Get aligned GRNs based on alignments.
    
    Args:
        grnp (GRNProcessor): GRN processor instance
        alignments (dict): Dictionary of alignments
        best_matches (list): List of [query_id, ref_id] pairs
        grns_str_strict (list): List of strict GRNs to filter by
    
    Returns:
        dict: Dictionary of aligned GRNs keyed by query_id
    """
    new_rows = {}
    for query_id, ref_id in best_matches:
        try:
            alignment = alignments[query_id]
            ref_row = grnp.data.loc[ref_id, :]
            ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
            seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
            new_row = init_row_from_alignment(alignment, seq_pos2grn)
            new_row_grns = [grn for grn in grns_str_strict if (len(grn.split('x')[0]) < 2) &
                               (grn in new_row.index.tolist())]
            new_row_strict = new_row.loc[new_row_grns]
            if len(new_row_strict) == 0:
                print(query_id, "No strict GRNs found in alignment.")
            new_rows[query_id] = new_row_strict
        except Exception as e:
            print(f"Error initializing row of sequence {query_id}: {e}")
    return new_rows


def annotate_sequence(query_id, query_seq, new_row, protein_family):
    """
    Annotate a sequence with GRNs.
    
    Args:
        query_id (str): Query sequence ID
        query_seq (str): Query sequence
        new_row (Series): New row with GRN annotations
        protein_family (str): Protein family name
    
    Returns:
        tuple: (query_id, GRN dictionary) or (None, None) if error
    """
    try:
        new_row_seq = ''.join([x[0] for x in new_row]).replace('-', '')

        alignment = align_blosum62(query_seq, new_row_seq, init_aligner(), verbose=0)
        alignment = format_alignment(alignment)
        grns, rns, missing = expand_annotation(new_row, query_seq, alignment, protein_family=protein_family,
                                           max_alignment_gap=1, verbose=0)
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


def main_parallel_execution(query_dict, new_rows, protein_family, num_cores=8):
    """
    Execute the annotation process in parallel.
    
    Args:
        query_dict (dict): Dictionary of query sequences
        new_rows (dict): Dictionary of new rows with aligned GRNs
        protein_family (str): Protein family name
        num_cores (int, optional): Number of CPU cores to use. Defaults to 8.
    
    Returns:
        dict: Dictionary of GRN annotations keyed by query_id
    """
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(annotate_sequence, query_id, query_dict[query_id], new_rows[query_id], protein_family)
                  for query_id in query_dict]
        results = {}
        for future in futures:
            try:
                query_id, grn_dict = future.result()
                if grn_dict is not None:
                    results[query_id] = grn_dict
                else:
                    print("why is the result in parallel execution None?")
            except Exception as e:
                print(f"Error in processing: {str(e)}")
        return results


def assign_grns(protein_family, dataset, num_cores=8, data_root=None):
    """
    Assign GRNs to sequences in a dataset.
    
    Args:
        protein_family (str): Protein family name ('gpcr_a', 'microbial_opsins', etc.)
        dataset (str): Dataset name
        num_cores (int, optional): Number of CPU cores to use. Defaults to 8.
        data_root (str, optional): Root data directory. If None, uses PROTOS_DATA_ROOT env var.
    
    Returns:
        DataFrame: DataFrame with GRN annotations
    """
    # Set up paths
    if data_root is None:
        data_root = os.environ.get('PROTOS_DATA_ROOT', 'data')
    
    data_root = Path(data_root).absolute()
    
    config = GRNConfigManager(protein_family=protein_family)
    grn_config_strict = config.get_config(strict=True)
    grns_str_strict = init_grn_intervals(grn_config_strict)

    # Initialize GRNBaseProcessor with proper paths
    if protein_family == 'gpcr_a':
        ref_dataset_id = 'ref/ref'
    else:
        ref_dataset_id = 'ref/mo_ref'
    
    grnp = GRNBaseProcessor(
        name=f'{protein_family}_processor',
        data_root=str(data_root),
        processor_data_dir='grn',
        dataset=ref_dataset_id,
        preload=True
    )

    ref_dict = {k: v.replace('-', '') for k, v in grnp.get_seq_dict().items()}
    
    # Read query sequences from fasta
    fasta_path = data_root / 'sequence' / 'fasta' / 'processed' / f'{dataset}.fasta'
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    
    query_dict = read_fasta(str(fasta_path))
    query_dict = {k: v.replace('-', '') for k, v in query_dict.items()}
    print(f"Loaded {len(query_dict)} sequences for processing")
    
    output = mmseqs2_align2(query_seqs=query_dict, ref_seqs=ref_dict)
    best_matches = output.loc[output.groupby('query_id')['e_value'].idxmin()][['query_id', 'target_id']].values.tolist()

    alignments = get_pairwise_alignment(query_dict, ref_dict, best_matches=best_matches)
    alignments = {k: v for k, v in alignments.items() if len(v) > 0}
    print(f"Number of sequences with alignments: {len(alignments)}")

    best_matches = [[x, y] for [x, y] in best_matches if x in alignments.keys()]
    print(f"Number of best matches: {len(best_matches)}")

    new_rows = get_aligned_grns(grnp, alignments, best_matches, grns_str_strict)
    print(f"Number of sequences with aligned GRNs: {len(new_rows)}")
    query_dict = {k: v for k, v in query_dict.items() if k in new_rows.keys()}

    final_results = main_parallel_execution(query_dict, new_rows, protein_family=protein_family, num_cores=num_cores)

    print(f"Number of sequences with GRNs: {len(final_results)}")
    
    if not final_results:
        print("No sequences were successfully processed.")
        return None

    df = pd.DataFrame.from_dict(final_results, orient='index')
    df = df.loc[:, sort_grns_str(df.columns.tolist())]

    # Crude check to see if the schiffbase is correct
    # Handle both x and dot notation
    if protein_family == 'microbial_opsins':
        if '7x50' in df.columns:
            df = df[df['7x50'].str.contains('K')]
        elif '7.50' in df.columns:
            df = df[df['7.50'].str.contains('K')]
    elif protein_family == 'gpcr_a':
        if '7x43' in df.columns:
            df = df[df['7x43'].str.contains('K')]
        elif '7.43' in df.columns:
            df = df[df['7.43'].str.contains('K')]

    # Save using path configuration
    output_path = data_root / 'grn' / 'datasets' / f'{dataset}.csv'
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=True)
    print(f"Saved GRN table to {output_path}")
    return df


def main():
    """Command-line entry point for GRN assignment."""
    parser = argparse.ArgumentParser(
        description='Process GRN annotations for a dataset',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Assign GRNs to microbial opsin sequences
  python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_dataset
  
  # Assign GRNs to GPCR sequences with custom data root
  python -m protos.cli.grn.assign_grns -p gpcr_a -d gpcr_dataset --data-root /path/to/data
  
  # Use more cores for parallel processing
  python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_dataset -n 16
        """
    )
    parser.add_argument('-p', '--protein_family', required=True, 
                      help='Protein family (e.g., gpcr_a, microbial_opsins)')
    parser.add_argument('-d', '--dataset', required=True, 
                      help='Dataset name (e.g., Bacteriorhodopsins)')
    parser.add_argument('-n', '--num_cores', type=int, default=8, 
                      help='Number of cores for parallel processing')
    parser.add_argument('--data-root', type=str, default=None,
                      help='Root data directory (default: uses PROTOS_DATA_ROOT env var or "data")')

    args = parser.parse_args()
    
    result = assign_grns(
        protein_family=args.protein_family,
        dataset=args.dataset,
        num_cores=args.num_cores,
        data_root=args.data_root
    )
    
    if result is not None:
        print("Done! GRN assignment completed successfully.")
        print(f"Assigned GRNs to {len(result)} sequences")
    else:
        print("GRN assignment failed.")


if __name__ == '__main__':
    main()