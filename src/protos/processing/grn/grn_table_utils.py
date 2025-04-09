from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_assignment import *
from protos.processing.grn.grn_utils import *
from protos.processing.sequence.seq_alignment import *


def expand_annotation(new_row, query_seq, alignment, max_alignment_gap=1, protein_family: str = 'gpcr_a',
                      verbose=0):
    # Create reference_grn_dict from the input row
    aligned_grns = {v: k for (k, v) in new_row.to_dict().items() if (v != '-') & ('x' in k)}
    # Calculate the length of the query sequence
    query_gene_len = len(query_seq)

    # Create all_query_gene_numbers using the assign_gene_nr() function
    if isinstance(query_seq, str):
        all_query_gene_numbers = assign_gene_nr(query_seq)
    elif isinstance(query_seq, list):
        all_query_gene_numbers = query_seq
    else:
        print("Type of queryseq is not implemented, use ")

    # Initialize GRN intervals to determine missing grns
    config = GRNConfigManager(protein_family=protein_family)
    grn_config_std = config.get_config(strict=False)
    grns_str_std = init_grn_intervals(grn_config_std)
    grns_str_std = [grn for grn in grns_str_std if len(grn.split('x')[0]) < 2]
    grns_float_std = sort_grns([round(parse_grn_str2float(x), len(x) - 2) for x in grns_str_std])

    if verbose > 0:
        print("aligned_grns allowed in grn configuration", aligned_grns)

    # annotate missing standard grns
    aligned_grns = get_correctly_aligned_grns(all_query_gene_numbers, aligned_grns, alignment, max_alignment_gap)

    if verbose > 0:
        print("extended aligned grns", aligned_grns)

        # Calculate missing_query_gene_numbers
    missing_gene_numbers = calculate_missing_gene_numbers(all_query_gene_numbers, aligned_grns)

    # Annotate N-tail
    n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers,
                                                                      grns_float_std)
    if verbose > 0:
        print("after annotation: n_tail_list =", n_tail_list)

    # Annotate C-tail
    c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len,
                                                                     grns_float_std)

    if verbose > 0:
        print("after annotation: c_tail_list =", c_tail_list)
    # Combine the results and return the expanded GRN list, residue number list, and missing residue numbers
    expanded_grn_list = n_tail_list + list(aligned_grns.items()) + c_tail_list

    if verbose > 0:
        print("n-tail, tm, c-tail:", expanded_grn_list)

    missing_gene_numbers = calculate_missing_gene_numbers(all_query_gene_numbers, expanded_grn_list)

    if verbose > 0:
        print("missing_gene_numbers", missing_gene_numbers)

    # Annotate (gaps) missing standard GRNs
    missing_gene_numbers_int = [int(mgnr[1:]) for mgnr in missing_gene_numbers]
    present_seq_nr_grn_list = expanded_grn_list
    present_grns = [g[1] for g in present_seq_nr_grn_list]
    missing_std_grns = [grn for grn in grns_str_std if grn not in present_grns]

    # Missing SEQNRs
    missing = [x for x in missing_gene_numbers_int if (x > first_gene_number_int) & (x < last_gene_number_int)]

    grns, missing = assign_missing_std_grns(missing_std_grns, present_seq_nr_grn_list, query_seq, missing, grns_str_std)

    expanded_grn_list += grns

    # Annotate loops and gaps between transmembrane helices
    nloop, gaps, cloop = [], [], []
    if len(missing) > 0:
        nloop, gaps, cloop = annotate_gaps_and_loops(expanded_grn_list, missing, query_seq, grn_config_std, grns_str_std)

    if verbose > 0:
        print("nloop", nloop)
        print("gaps", gaps)
        print("cloop", cloop)

    expanded_grn_list = expanded_grn_list + nloop + gaps + cloop
    grn_list = [x[1] for x in expanded_grn_list]
    rn_list = [x[0] for x in expanded_grn_list]

    # Sort and complete the GRN/RN pairs
    grn_rn_pairs = list(zip(grn_list, rn_list))
    sorted_grn_f_list = sort_grns([parse_grn_str2float(x) for x in grn_list])
    sorted_grn_list = [parse_grn_float2str(x) for x in sorted_grn_f_list]

    grn_rn_pairs = sort_grn_rn_pairs(sorted_grn_list, grn_rn_pairs)
    grn_list, rn_list = zip(*grn_rn_pairs)
    grn_list = list(grn_list)
    rn_list = list(rn_list)

    seq_ids = [int(x[1:]) for x in rn_list]
    missing = [x + 1 for x in range(len(query_seq)) if (x + 1) not in seq_ids]

    return grn_list, rn_list, missing


def init_row_from_alignment(alignment, seq_pos2grn):
    qidx = 0
    ridx = 0
    new_row = []
    alignment[0] = alignment[0].replace('\x00', '')
    alignment[1] = alignment[1].replace('\x00', '')
    alignment[2] = alignment[2].replace('\x00', '')
    for idx in range(len(alignment[0])):
        qi = alignment[0][idx]
        if qi != '-':
            qidx += 1
        mm = alignment[1][idx]
        ri = alignment[2][idx]
        if ri != '-':
            ridx += 1
            grni = seq_pos2grn[ridx]
        if mm != '-':
            new_row.append((grni, qi + str(qidx)))
    new_row = pd.Series(dict(new_row))
    return new_row


def annotate_gpcr(grnp: GRNProcessor, query_name: str, query_seq: str,
                  add_to_GRNP: bool = False, verbose=0, protein_family='microbial_opsins',
                  reload=True):
    # Initialize aligner
    aligner = init_aligner()

    # Rset the grnp, we need to have loops in the alignment, otherwise we have 'gap-jumps', leading to errenous assignments
    if reload:
        grnp.reset_data()

    # Find the best match in the sequence database
    hits = mmseqs2_align(query_seq.replace('-', ''), grnp.get_seq_dict())
    best_match = hits['target_id'].iloc[0]

    if verbose > 0:
        print("best hit with mmseqs:", best_match)

    # Get reference sequence
    ref_seq = get_seq(best_match, grnp.data)

    # Perform MSA and remove gaps
    query_dict_no_gaps = remove_gaps_from_sequences({query_name: query_seq})
    ref_dict = {best_match: ref_seq}
    msa = msa_blosum62(query_dict_no_gaps, ref_dict, aligner)

    # Construct the initial annotation based on alignment
    entry = msa[query_name]
    name = entry[0]

    alignment = entry[2]
    if verbose > 0:
        print('alignment>>>')
        print(alignment[0])
        print(alignment[1])
        print(alignment[2])

    ref_row = grnp.data.loc[name, :]
    ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
    seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
    new_row = init_row_from_alignment(alignment, seq_pos2grn)

    # Filter by strict grns
    config = GRNConfigManager(protein_family=protein_family)
    grn_config_strict = config.get_config(strict=True)
    grns_str_strict = init_grn_intervals(grn_config_strict)
    grns_str_strict = [grn for grn in grns_str_strict if (len(grn.split('x')[0]) < 2) & (grn in new_row.index.tolist())]
    new_row = new_row[grns_str_strict]

    # Perform sequence expansion on our new row!
    new_row_seq = ''.join([x[0] for x in new_row.tolist()]).replace('-', '')
    alignment = msa_blosum62(query_dict_no_gaps, {'incomplete_new_row': new_row_seq}, aligner)

    if verbose > 0:
        print("new row:", new_row)
        print("alignment:", alignment)
    grn, rn, missing = expand_annotation(new_row, query_seq.replace('-', ''), alignment[query_name][2],
                                         max_alignment_gap=1, protein_family=protein_family,
                                         verbose=verbose)

    if len(missing) > 0:
        print("Missing residues:", missing)

    if verbose > 0:
        print("GRN:", grn)
        print("RN:", rn)

    # Create a row in the GRNP-format
    row = pd.Series(dict(zip(grn, rn)))

    if add_to_GRNP:
        merged_df = add_row_to_table(query_name, grn, rn, grnp.data, save=False)
        return merged_df
    else:
        return row

def add_row_to_table(name, grn_list, rn_list, table, save=False, outfile=''):
    all_cols = table.columns.tolist()
    all_cols_ = list(set(list(grn_list) + all_cols))
    if len(all_cols_) > len(all_cols):
        all_cols_float = sort_grns([parse_grn_str2float(x) for x in all_cols_])
        all_cols = [parse_grn_float2str(x) for x in all_cols_float]
    all_rows = []
    dense_cells = []
    for col in all_cols:
        if col not in grn_list:
            dense_cells.append('-')
        else:
            dense_cells.append(rn_list[grn_list.index(col)])
    dense_row = pd.DataFrame([dense_cells], columns=all_cols)
    dense_row.rename(index={0: name}, inplace=True)
    all_rows.append(dense_row)
    table = pd.concat([table, dense_row])
    table = table[all_cols]
    table.fillna('-', inplace=True)
    if save:
        table.to_csv(outfile)
    return table


def annotate_sequence(name: str, sequence: str, grnp: GRNProcessor, min_score: float = 0.25, add=True):
    aligner = init_aligner()
    query_dict = {name: sequence}
    best_msa = msa_blosum62(seqs_query=query_dict,
                            seqs_ref=grnp.get_seq_dict(),
                            aligner=aligner)  # aligns our reference to ALL queries

    ref_name = best_msa[name][0]
    alignment = best_msa[name][2]

    miss_match_nan = alignment[1]
    index_std_start = find_idx(miss_match_nan)
    index_std_end = find_idx(miss_match_nan[::-1])

    miss_match_nan = miss_match_nan[index_std_start:-index_std_end]
    n_mm = miss_match_nan.count('.')
    n_gap = miss_match_nan.count('-')

    score = 1 - ((n_mm + n_gap) / (len(miss_match_nan)))

    if score >= min_score:
        # GRN ANNOTATION => update existing residue annotations, no overwriting!
        row = grnp.data.T[ref_name]

        try:
            grn, rn, missing = expand_annotation(row, sequence, alignment)
        except:
            grn, rn, missing = None, None, None
    if add:
        merged_df = add_row_to_table(name, grn, rn, grnp.data, save=False, outfile=grnp.path + 'annotate_seq.csv')
        return merged_df
    elif grn != None:
        return dict(zip(grn, rn)), missing
    else:
        return None, None


def is_sequential(result):
    current = 1
    for value in result:
        if int(value[1:]) != current:
            return False
        current += 1
    return True