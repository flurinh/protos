from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_assignment import *
from protos.processing.grn.grn_utils import *
from protos.processing.sequence.seq_alignment import *


def expand_annotation(new_row, query_seq, alignment, max_alignment_gap=1, protein_family: str = 'gpcr_a',
                      verbose=0):

    # Create reference_grn_dict from the input row
    # Only check for dot notation (x notation is deprecated)
    aligned_grns = {v: k for (k, v) in new_row.to_dict().items() if (v != '-') & ('.' in k)}
    
    if verbose >= 1:
        print(f"\n=== Expand Annotation Started ===")
        print(f"Query length: {len(query_seq)}, Initial GRNs: {len(aligned_grns)}")
    
    if verbose >= 2:
        print(f"Input new_row has {len(new_row)} items")
        print(f"Alignment format: {len(alignment)} parts")
        if len(alignment) >= 3:
            print(f"  Query:  {alignment[0][:50]}...")
            print(f"  Match:  {alignment[1][:50]}...")
            print(f"  Ref:    {alignment[2][:50]}...")

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
    config = GRNConfigManager()
    grn_config_std = config.get_config(protein_family=protein_family, strict=False)
    
    # Generate all standard GRNs from the configuration
    # First, collect all GRNs for each TM region using get_grn_interval
    grns_str_std = []
    if grn_config_std:
        for region_name, (start_grn, end_grn) in grn_config_std.items():
            # Generate GRNs for this interval (auto-generate since we don't have a predefined list)
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_std.extend(region_grns)
    
    # Remove duplicates and sort
    grns_str_std = list(set(grns_str_std))
    grns_str_std = sort_grns_str(grns_str_std)
    
    if verbose >= 2:
        print(f"\nGRN Config for {protein_family}:")
        print(f"  Config keys: {list(grn_config_std.keys()) if grn_config_std else 'None'}")
        print(f"  Total standard GRNs: {len(grns_str_std)}")

    # Filter to single TM GRNs (exclude H8, etc.)
    grns_str_std_filtered = []
    for grn in grns_str_std:
        if '.' in grn:
            tm_part = grn.split('.')[0]
        else:
            continue
        # Only include single TM regions (not loops)
        if len(tm_part) == 1 and tm_part.isdigit():
            grns_str_std_filtered.append(grn)

    grns_str_std = grns_str_std_filtered
    grns_float_std = sort_grns([parse_grn_str2float(x) for x in grns_str_std])
    
    if verbose >= 2:
        print(f"  Filtered to single TM GRNs: {len(grns_str_std)}")
        print(f"  First few GRNs: {grns_str_std[:10] if grns_str_std else 'None'}")
        print("aligned_grns allowed in grn configuration", aligned_grns)

    # annotate missing standard grns
    # NOTE: get_correctly_aligned_grns is only needed when we have a new alignment
    # In expand_annotation, aligned_grns already contains query->grn mapping which is what we need
    # Skip get_correctly_aligned_grns since we already have the alignments from init_row_from_alignment

    if verbose >= 2:
        print("extended aligned grns", aligned_grns)

        # Calculate missing_query_gene_numbers
    missing_gene_numbers = calculate_missing_gene_numbers(all_query_gene_numbers, aligned_grns)

    # Annotate N-tail
    if verbose >= 1:
        print(f"\n1. Processing N-terminal region...")
        
    if verbose >= 2:
        print(f"  First aligned position: {list(aligned_grns.keys())[0] if aligned_grns else 'None'}")
        print(f"  First GRN: {list(aligned_grns.values())[0] if aligned_grns else 'None'}")
        print(f"  Missing positions before first: {[x for x in missing_gene_numbers if int(x[1:]) < (int(list(aligned_grns.keys())[0][1:]) if aligned_grns else 999)]}")
    
    n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers,
                                                                      grns_float_std)
    if verbose >= 1:
        print(f"   N-tail: {len(n_tail_list)} positions assigned")
        
    if verbose >= 2:
        print("   n_tail_list =", n_tail_list)
        print(f"   first_gene_number_int: {first_gene_number_int}")

    # Annotate C-tail
    if verbose >= 1:
        print(f"\n2. Processing C-terminal region...")
        
    if verbose >= 2:
        print(f"  Last aligned position: {list(aligned_grns.keys())[-1] if aligned_grns else 'None'}")
        print(f"  Last GRN: {list(aligned_grns.values())[-1] if aligned_grns else 'None'}")
        print(f"  Query length: {query_gene_len}")
        print(f"  Missing positions after last: {[x for x in missing_gene_numbers if int(x[1:]) > (int(list(aligned_grns.keys())[-1][1:]) if aligned_grns else 0)]}")
    
    c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len,
                                                                     grns_float_std)

    if verbose >= 1:
        print(f"   C-tail: {len(c_tail_list)} positions assigned")
        
    if verbose >= 2:
        print("   c_tail_list =", c_tail_list)
        print(f"   last_gene_number_int: {last_gene_number_int}")
    # Combine the results and return the expanded GRN list, residue number list, and missing residue numbers
    expanded_grn_list = n_tail_list + list(aligned_grns.items()) + c_tail_list

    if verbose >= 2:
        print(f"\n   Combined N-tail + aligned + C-tail: {len(expanded_grn_list)} total positions")

    missing_gene_numbers = calculate_missing_gene_numbers(all_query_gene_numbers, expanded_grn_list)

    if verbose >= 2:
        print(f"   Still missing: {len(missing_gene_numbers)} positions")

    # Annotate (gaps) missing standard GRNs
    if verbose >= 1:
        print(f"\n3. Processing missing standard GRNs...")
    
    missing_gene_numbers_int = [int(mgnr[1:]) for mgnr in missing_gene_numbers]
    present_seq_nr_grn_list = expanded_grn_list
    present_grns = [g[1] for g in present_seq_nr_grn_list]
    missing_std_grns = [grn for grn in grns_str_std if grn not in present_grns]
    
    if verbose >= 2:
        print(f"   Standard GRNs not yet assigned: {len(missing_std_grns)}")
        if missing_std_grns:
            print(f"   First few missing: {missing_std_grns[:5]}...")

    # Missing SEQNRs
    if first_gene_number_int is not None and last_gene_number_int is not None:
        missing = [x for x in missing_gene_numbers_int if (x > first_gene_number_int) & (x < last_gene_number_int)]
    else:
        missing = []
        
    if verbose >= 2:
        print(f"   Internal missing positions (between first and last): {len(missing)}")

    grns, missing = assign_missing_std_grns(missing_std_grns, present_seq_nr_grn_list, query_seq, missing, grns_str_std)
    
    if verbose >= 1 and grns:
        print(f"   Assigned {len(grns)} standard GRNs")

    expanded_grn_list += grns

    # Annotate loops and gaps between transmembrane helices
    if verbose >= 1:
        print(f"\n4. Processing loops and gaps...")
        
    nloop, gaps, cloop = [], [], []
    if len(missing) > 0:
        if verbose >= 2:
            print(f"   Processing {len(missing)} remaining missing positions")
            
        nloop, gaps, cloop = annotate_gaps_and_loops(expanded_grn_list, missing, query_seq, grn_config_std, grns_str_std)
        
        if verbose >= 1:
            total_annotations = len(nloop) + len(gaps) + len(cloop)
            if total_annotations > 0:
                print(f"   Annotated {total_annotations} positions:")
                if nloop:
                    print(f"     - N-side loops: {len(nloop)}")
                if gaps:
                    print(f"     - Gaps: {len(gaps)}")
                if cloop:
                    print(f"     - C-side loops: {len(cloop)}")

    if verbose >= 2:
        if nloop:
            print(f"   nloop details: {nloop[:3]}..." if len(nloop) > 3 else f"   nloop details: {nloop}")
        if gaps:
            print(f"   gaps details: {gaps[:3]}..." if len(gaps) > 3 else f"   gaps details: {gaps}")
        if cloop:
            print(f"   cloop details: {cloop[:3]}..." if len(cloop) > 3 else f"   cloop details: {cloop}")

    expanded_grn_list = expanded_grn_list + nloop + gaps + cloop
    
    if verbose >= 1:
        print(f"\n5. Finalizing annotation...")
        print(f"   Total positions annotated: {len(expanded_grn_list)}")
    
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
    
    if verbose >= 1:
        if missing:
            print(f"   WARNING: {len(missing)} positions could not be annotated")
        else:
            print(f"   SUCCESS: All {len(query_seq)} positions annotated")
            
    if verbose >= 1:
        print(f"\n=== Expand Annotation Completed ===\n")

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


def annotate_grnp(grnp: GRNProcessor, query_name: str, query_seq: str,
                  add_to_GRNP: bool = False, verbose=0, protein_family='microbial_opsins',
                  reload=True):
    """Annotate a sequence with GRN numbers using the GRN processor.
    
    Args:
        grnp: GRNProcessor instance with reference data
        query_name: Name for the query sequence
        query_seq: Query sequence to annotate
        add_to_GRNP: Whether to add the result to the GRN processor
        verbose: Verbosity level (0=silent, 1=major steps, 2=detailed)
        protein_family: Protein family for GRN configuration
        reload: Whether to reload the GRN processor data
    """
    if verbose >= 1:
        print(f"\n=== GRN Annotation Started ===")
        print(f"Query: {query_name}, Length: {len(query_seq)}")
        print(f"Protein family: {protein_family}")
        
    # Initialize aligner
    aligner = init_aligner()

    # Reset the grnp, we need to have loops in the alignment, otherwise we have 'gap-jumps', leading to erroneous assignments
    if reload:
        if verbose >= 2:
            print("\nReloading GRN processor data...")
        grnp.reset_data()

    # Find the best match in the sequence database
    if verbose >= 1:
        print(f"\n1. Finding best reference match...")
        
    hits = mmseqs2_align(query_seq.replace('-', ''), grnp.get_seq_dict())
    best_match = hits['target_id'].iloc[0]

    if verbose >= 1:
        print(f"   Best hit: {best_match}")
        if verbose >= 2:
            print(f"   Score: {hits.iloc[0]['score']:.3f}")
            print(f"   E-value: {hits.iloc[0]['evalue']:.2e}")

    # Get reference sequence
    ref_seq = get_seq(best_match, grnp.data)
    
    if verbose >= 2:
        print(f"\n   Reference length: {len(ref_seq)}")

    # Perform MSA and remove gaps
    if verbose >= 1:
        print(f"\n2. Performing sequence alignment...")
        
    query_dict_no_gaps = remove_gaps_from_sequences({query_name: query_seq})
    ref_dict = {best_match: ref_seq}
    msa = msa_blosum62(query_dict_no_gaps, ref_dict, aligner)

    # Construct the initial annotation based on alignment
    entry = msa[query_name]
    name = entry[0]

    alignment = entry[2]
    if verbose >= 2:
        print("\n   Alignment:")
        print(f"   Query:  {alignment[0][:50]}..." if len(alignment[0]) > 50 else f"   Query:  {alignment[0]}")
        print(f"   Match:  {alignment[1][:50]}..." if len(alignment[1]) > 50 else f"   Match:  {alignment[1]}")
        print(f"   Ref:    {alignment[2][:50]}..." if len(alignment[2]) > 50 else f"   Ref:    {alignment[2]}")

    ref_row = grnp.data.loc[name, :]
    ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
    seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])
    new_row = init_row_from_alignment(alignment, seq_pos2grn)
    
    if verbose >= 2:
        print(f"\n   Initial GRN mapping: {len(new_row)} positions")

    # Filter by strict grns
    if verbose >= 1:
        print(f"\n3. Filtering to strict GRN positions...")
        
    config = GRNConfigManager()
    grn_config_strict = config.get_config(protein_family=protein_family, strict=True)
    
    # Generate all strict GRNs from the configuration
    grns_str_strict = []
    if grn_config_strict:
        for region_name, (start_grn, end_grn) in grn_config_strict.items():
            # Generate GRNs for this interval
            region_grns = get_grn_interval(start_grn, end_grn, grns_str=None)
            grns_str_strict.extend(region_grns)
    
    # Remove duplicates and sort
    grns_str_strict = list(set(grns_str_strict))
    grns_str_strict = sort_grns_str(grns_str_strict)

    # Filter to single TM GRNs that are in new_row
    grns_str_strict_filtered = []
    for grn in grns_str_strict:
        if '.' in grn:
            tm_part = grn.split('.')[0]
        else:
            continue
        # Only include single TM regions that exist in new_row
        if len(tm_part) == 1 and tm_part.isdigit() and grn in new_row.index.tolist():
            grns_str_strict_filtered.append(grn)
    grns_str_strict = grns_str_strict_filtered
    new_row = new_row[grns_str_strict]
    
    if verbose >= 1:
        print(f"   Strict positions: {len(new_row)}")
        if verbose >= 2 and len(new_row) > 0:
            print(f"   Positions: {list(new_row.index)[:10]}..." if len(new_row) > 10 else f"   Positions: {list(new_row.index)}")

    # Perform sequence expansion on our new row!
    if verbose >= 1:
        print(f"\n4. Expanding annotation to full sequence...")
        
    new_row_seq = ''.join([x[0] for x in new_row.tolist()]).replace('-', '')
    alignment = msa_blosum62(query_dict_no_gaps, {'incomplete_new_row': new_row_seq}, aligner)

    if verbose >= 2:
        print(f"   Re-aligned with {len(new_row)} strict positions")
        
    grn, rn, missing = expand_annotation(new_row, query_seq.replace('-', ''), alignment[query_name][2],
                                         max_alignment_gap=1, protein_family=protein_family,
                                         verbose=verbose)

    if verbose >= 1:
        print(f"\n=== GRN Annotation Results ===")
        print(f"Total positions annotated: {len(grn)}")
        if missing:
            print(f"Missing positions: {len(missing)}")
            if verbose >= 2:
                print(f"  Missing: {missing[:10]}..." if len(missing) > 10 else f"  Missing: {missing}")
                
    if verbose >= 2:
        print(f"\nFirst 10 GRN assignments:")
        for i in range(min(10, len(grn))):
            print(f"  {rn[i]} -> {grn[i]}")
        if len(grn) > 10:
            print(f"  ... and {len(grn) - 10} more")

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


def is_sequential(rn_list):
    current = 1
    for value in rn_list:
        if int(value[1:]) != current:
            return False
        current += 1
    return True