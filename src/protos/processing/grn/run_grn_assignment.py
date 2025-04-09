from concurrent.futures import ProcessPoolExecutor
from protos.processing.grn.grn_table_utils import *
from protos.processing.sequence.fasta_utils import *
from protos.processing.grn.seq_alignment import init_aligner, align_blosum62, format_alignment
import pandas as pd
import argparse


def get_pairwise_alignment(query_seq_dict, ref_seq_dict, best_matches, aligner=None):
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
                pass
    return results



if __name__ == '__main__':
    # Add argument parsing
    parser = argparse.ArgumentParser(description='Process GRN annotations for a dataset')
    parser.add_argument('-p', '--protein_family', required=True, help='Protein family (e.g., gpcr_a, microbial_opsins)')
    parser.add_argument('-d', '--dataset', required=True, help='Dataset name (e.g., Bacteriorhodopsins)')
    parser.add_argument('-n', '--num_cores', type=int, default=8, help='Number of cores for parallel processing')

    args = parser.parse_args()

    protein_family = args.protein_family
    dataset = args.dataset
    num_cores = args.num_cores

    config = GRNConfigManager(protein_family=protein_family)
    grn_config_strict = config.get_config(strict=True)
    grns_str_strict = init_grn_intervals(grn_config_strict)

    if protein_family == 'gpcr_a':
        grnp = GRNProcessor('ref', 'data/grn/ref/')
    else:
        grnp = GRNProcessor('mo_ref', 'data/grn/ref/')

    ref_dict = {k: v.replace('-', '') for k, v in grnp.get_seq_dict().items()}
    query_dict = read_fasta(f"data/fasta/processed/{dataset}.fasta")

    query_dict = {k: v.replace('-', '') for k, v in query_dict.items()}
    print(query_dict.keys())
    output = mmseqs2_align2(query_seqs=query_dict, ref_seqs=ref_dict)
    best_matches = output.loc[output.groupby('query_id')['e_value'].idxmin()][['query_id', 'target_id']].values.tolist()

    alignments = get_pairwise_alignment(query_dict, ref_dict, best_matches=best_matches)
    alignments = {k: v for k, v in alignments.items() if len(v) > 0}
    print("Number of sequences with alignments:", len(alignments))

    best_matches = [[x, y] for [x, y] in best_matches if x in alignments.keys()]
    print("Number of best matches:", len(best_matches))

    new_rows = get_aligned_grns(grnp, alignments, best_matches, grns_str_strict)
    print("Number of sequences with aligned GRNs:", len(new_rows))
    query_dict = {k: v for k, v in query_dict.items() if k in new_rows.keys()}

    final_results = main_parallel_execution(query_dict, new_rows, protein_family=protein_family, num_cores=num_cores)

    print("Number of sequences with GRNs:", len(final_results))
    print("done!")

    df = pd.DataFrame.from_dict(final_results, orient='index')
    df = df.loc[:, sort_grns_str(df.columns.tolist())]
    df.to_csv(f"data/grn/datasets/{dataset}.csv", index=True)