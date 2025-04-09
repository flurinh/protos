import time
import logging
import os
from protos.processing.grn.grn_table_utils import *
from protos.loaders.uniprot_utils import get_uniprot


def get_ref_table(path=r'data/gpcrdb/residue_table.csv', interval=None):
    table = pd.read_csv(path, index_col='GPCRdb(A)').T
    if len(interval) > 0:
        table = table.loc[:,interval]
    return table


def init_default_grn_table(ref_table, error_dict={}, max_seq_len=3000, min_score=.2, max_alignment_gap=1,
                           load=True, load_file='data/grn/all/grn_tableV1.csv',
                           out_file='data/grn/all/grn_tableV1.csv',
                           overwrite=False, limit=None, save=True,
                           organism='HUMAN',
                           sleep=2):
    # Initialize utilities
    logging.basicConfig(filename='error_log.log', level=logging.ERROR,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s - [%(filename)s:%(lineno)d]')

    aligner = init_aligner()

    results = {}
    if load & os.path.isfile(load_file):
        print('Loading file...', load_file)
        processed_table = pd.read_csv(load_file, index_col=[0])
        processed_genes = processed_table.index.tolist()
        print('Found file with {} processed genes!'.format(len(processed_genes)))
        all_cols = list(processed_table.columns)
    else:
        print("Initializing new output grn table!", out_file)
        processed_table = pd.DataFrame([], columns=STD_GRNS)
        processed_genes = []
        all_cols = STD_GRNS

    all_rows = []
    genes = ref_table.index.tolist()

    # processed_genes = [re.sub('[^A-Za-z0-9]+', ' ', g) for g in processed_genes]
    genes = [g for g in genes if g not in processed_genes]
    error_list = []

    if limit == None:
        limit = len(genes)

    done = 0
    for g, gene in enumerate(genes):
        print(g, gene)
        try:
            exists = False
            skip = False
            if (limit - done) > 0:
                if True:
                    valid_gene_name = re.sub('[^A-Za-z0-9]+', ' ', gene)

                    if gene in list(error_dict.keys()):
                        valid_gene_name = error_dict[gene]
                        print('Found query in error_dict -> new query:', valid_gene_name)
                    query_results = get_uniprot(query=valid_gene_name, batchsize=8)
                    time.sleep(sleep)
                    # MEMORY: problems with memory during alignment, remove sequences longer than max_seq_len (1000)
                    if len(query_results) > 0:
                        query_results = query_results[query_results['Sequence'].str.len(
                        ) <= max_seq_len].reset_index()
                        if len(query_results) > 0:
                            query_dict = dict(query_results[['Entry Name', 'Sequence']].set_index(
                                'Entry Name')['Sequence'].T)
                            if organism != None:
                                query_dict = {k: v for k, v in query_dict.items() if organism in k}

                            ref_dict = {}
                            ref_seq = get_seq(gene, ref_table)
                            ref_dict.update({gene: ref_seq})
                            # ALIGNMENTS
                            msa = msa_blosum62(seqs_query=query_dict, seqs_ref=ref_dict, aligner=aligner)  # aligns our reference to ALL queries

                            name, alignment = get_best_alignment(msa)

                            query_seq_dotted = alignment[0]
                            miss_match_nan = alignment[1]
                            ref_seq_dotted = alignment[2]

                            index_std_start = find_idx(miss_match_nan)
                            index_std_end = find_idx(miss_match_nan[::-1])

                            miss_match_nan = miss_match_nan[index_std_start:-index_std_end]

                            n_mm = miss_match_nan.count('.')
                            n_gap = miss_match_nan.count('-')

                            score = 1 - ((n_mm + n_gap) / (len(miss_match_nan)))
                            print('Best score:', round(score, 6))

                            if score < min_score:
                                print('Did not find good enough sequence alignment', gene)
                                done += 1
                                error_list.append(gene)
                                skip =True

                            # check if this already exists
                            elif not overwrite:
                                if name in processed_table.index.tolist():
                                    exists = True
                                    print("Exists in table and will not be overwritten.")

                            if ((overwrite) or (not exists)) and (not skip):
                                if exists:
                                    processed_table.drop(name, inplace=True)
                                print("Does not exists and / or will be overwritten in table.")
                                results.update(
                                    {gene: [round(score, 3), query_seq_dotted, miss_match_nan, ref_seq_dotted]})
                                print("Best reference gene:", name)

                                # GRN ANNOTATION
                                query_seq = query_dict[name]
                                row = ref_table.T[gene]
                                grn_list, rn_list, missing = expand_annotation(row, query_seq, alignment, max_alignment_gap)

                                all_cols_ = list(set(grn_list + all_cols))
                                if len(all_cols_) > len(all_cols):
                                    all_cols_float = sort_grns([parse_grn_str2float(x) for x in all_cols_])
                                    all_cols = [parse_grn_float2str(x) for x in all_cols_float]
                                dense_cells = []
                                for col in all_cols:
                                    if col not in grn_list:
                                        dense_cells.append('-')
                                    else:
                                        dense_cells.append(rn_list[grn_list.index(col)])
                                dense_row = pd.DataFrame([dense_cells], columns=all_cols)
                                dense_row.rename(index={0: name}, inplace=True)
                                all_rows.append(dense_row)

                                processed_table = pd.concat([processed_table, dense_row])
                                processed_table = processed_table[all_cols]
                                processed_table.fillna('-', inplace=True)
                                if len(missing) > 0:
                                    logging.error(f"Missing assignments for gene '{gene}': '{missing}")
                                    error_list.append(gene)
                                if save & overwrite:
                                    processed_table.to_csv(out_file)
                                done += 1
                                print('\n\n')
        except Exception as e:
            logging.error(f"Error processing gene '{gene}': {e}", exc_info=True)
            error_list.append(gene)
            done += 1
            print('\n\n')
            continue
    processed_table = processed_table.fillna('-')
    if save:
        processed_table.to_csv(out_file, index=True)
    return processed_table, error_list