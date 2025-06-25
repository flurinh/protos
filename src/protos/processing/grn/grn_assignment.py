"""
Functions for GRN assignment using sequence alignment.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Any, Optional, Union

def parse_grn_str2float(grn_str: str) -> float:
    """Convert GRN string to float."""
    if 'x' in grn_str:
        parts = grn_str.split('x')
        return float(parts[0]) + float(parts[1]) / 100
    elif '.' in grn_str:
        return float(grn_str)
    else:
        return float(grn_str)

def parse_grn_float2str(grn_float: float) -> str:
    """Convert GRN float to string using dot notation."""
    if grn_float < 0:
        return f"-{abs(int(grn_float)):03d}"
    elif grn_float >= 100:
        return f"{int(grn_float):03d}"
    else:
        # Always format with 2 decimals for consistency
        return f"{grn_float:.2f}"

def get_grn_interval(start_grn: str, end_grn: str, grns_str: List[str]) -> List[str]:
    """Get GRN interval between start and end."""
    start_float = parse_grn_str2float(start_grn)
    end_float = parse_grn_str2float(end_grn)
    return [g for g in grns_str if start_float <= parse_grn_str2float(g) <= end_float]

def assign_gene_nr(seq: str) -> List[str]:
    """Assign gene numbers to sequence residues."""
    return [seq[x]+str(x+1) for x in range(len(seq))]

def calculate_missing_gene_numbers(all_gene_numbers, aligned_or_expanded_grns):
    """Calculate missing gene numbers from alignment."""
    if isinstance(aligned_or_expanded_grns, dict):
        return [gnr for gnr in all_gene_numbers if gnr not in list(aligned_or_expanded_grns.keys())]
    else:
        residue_number_list = [gnr for gnr, _ in aligned_or_expanded_grns]
        return [gnr for gnr in all_gene_numbers if gnr not in residue_number_list]

def calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers, grns_float):
    """Calculate N-terminal tail GRN assignments."""
    beginning_tm1_float = grns_float[0]
    first_grn = list(aligned_grns.values())[0]
    first_grn_float = parse_grn_str2float(first_grn)
    first_gene_number = list(aligned_grns.keys())[0]
    
    first_gene_number_int = int(first_gene_number[1:]) - 1
    
    missing_grns_tm1_ = int(100 * (first_grn_float - beginning_tm1_float))
    missing_grns_tm1 = min(first_gene_number_int, missing_grns_tm1_)
    missing_tm1 = min(missing_grns_tm1, missing_grns_tm1_)
    missing_ntail = max(0, first_gene_number_int - missing_tm1)
    
    n_tail_float = grns_float[:grns_float.index(first_grn_float)]
    n_tail_float += [-(i + 1) for i in range(missing_ntail)]
    n_tail_float = sorted(n_tail_float)
    
    n_tail_str = [parse_grn_float2str(x) for x in n_tail_float]
    n_tail_list = list(zip(missing_gene_numbers[:first_gene_number_int][::-1], n_tail_str[::-1]))[::-1]
    return n_tail_list, first_gene_number_int

def calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len, grns_float):
    """Calculate C-terminal tail GRN assignments."""
    if not aligned_grns:
        return [], None
    
    ending_tmx_float = grns_float[-1]
    last_grn = list(aligned_grns.values())[-1]
    last_grn_float = parse_grn_str2float(last_grn)
    last_gene_number = list(aligned_grns.keys())[-1]
    last_gene_number_int = int(last_gene_number[1:])
    
    missing_grns_last_section_ = int(100 * (ending_tmx_float - last_grn_float))
    missing_grns_last_section = min(query_gene_len - last_gene_number_int, missing_grns_last_section_)
    missing_last_section = max(0, query_gene_len - last_gene_number_int - missing_grns_last_section)
    
    c_tail_float = grns_float[grns_float.index(last_grn_float) + 1:]
    c_tail_float += [100 + i + 1 for i in range(missing_last_section)]
    c_tail_float = sorted(c_tail_float)
    
    c_tail_str = [parse_grn_float2str(x) for x in c_tail_float]
    c_tail_list = list(zip(missing_gene_numbers[-(query_gene_len - last_gene_number_int):], c_tail_str))
    
    return c_tail_list, last_gene_number_int

def valid_jump(prev_ref_grn, curr_ref_grn, prev_query_key, curr_query_key, max_alignment_gap):
    """Check if GRN assignment jump is valid."""
    if prev_ref_grn is None:
        return True
    
    # Handle both 'x' and '.' notation
    if 'x' in prev_ref_grn:
        prev_grn_tm = int(prev_ref_grn.split('x')[0])
    else:
        prev_grn_tm = int(float(prev_ref_grn))
    
    if 'x' in curr_ref_grn:
        curr_grn_tm = int(curr_ref_grn.split('x')[0])
    else:
        curr_grn_tm = int(float(curr_ref_grn))
    
    if prev_grn_tm != curr_grn_tm:
        return True
    
    prev_grn_float = parse_grn_str2float(prev_ref_grn)
    curr_grn_float = parse_grn_str2float(curr_ref_grn)
    
    prev_query_num = int(prev_query_key[1:])
    curr_query_num = int(curr_query_key[1:])
    
    grn_diff = abs(curr_grn_float - prev_grn_float)
    rn_diff = abs(curr_query_num - prev_query_num)
    
    if (grn_diff <= .1) and (rn_diff == max_alignment_gap):
        return True
    
    return False

def get_correctly_aligned_grns(all_query_gene_numbers, reference_grn_dict, alignment, max_alignment_gap=1):
    """Extract correctly aligned GRNs from alignment."""
    query_seq, match_line, ref_seq = alignment
    pointer_B = 0
    pointer_A = 0
    ref_keys = list(reference_grn_dict.keys())
    result_dict = {}
    
    prev_pair = None
    prev_query_key = None
    
    for i in range(len(match_line)):
        if (pointer_B < len(all_query_gene_numbers)) and (pointer_A < len(ref_keys)):
            query_key = all_query_gene_numbers[pointer_B]
            ref_key = ref_keys[pointer_A]
            
            if match_line[i] == '|' or match_line[i] == '.':
                if ref_seq[i] != '-' and query_seq[i] != '-':
                    curr_pair = reference_grn_dict[ref_key]
                    
                    if valid_jump(prev_pair, curr_pair, prev_query_key, query_key, max_alignment_gap):
                        result_dict[query_key] = curr_pair
                        prev_pair = curr_pair
                        prev_query_key = query_key
            
            if query_seq[i] != '-':
                pointer_B += 1
            if ref_seq[i] != '-':
                pointer_A += 1
    
    return result_dict

def _remove_loop_grns(grn_str_list):
    """Remove loop GRNs from list."""
    grn_str_list = [x for x in grn_str_list if ((not 'c' in x) and (not 'n' in x))]
    grn_float_list = [parse_grn_str2float(x) for x in grn_str_list]
    grn_float_list = [x for x in grn_float_list if x < 10]
    return [parse_grn_float2str(x) for x in grn_float_list]

def _get_seq_nr_intervals(seq_nrs: list):
    """Group consecutive sequence numbers into intervals."""
    seq_nr_intervals = []
    seq_nr_stack = []
    for seq_nr in seq_nrs:
        if len(seq_nr_stack) == 0:
            seq_nr_stack.append(seq_nr)
        else:
            if (seq_nr_stack[-1] + 1) == seq_nr:
                seq_nr_stack.append(seq_nr)
            else:
                seq_nr_intervals.append(seq_nr_stack)
                seq_nr_stack = []
                seq_nr_stack.append(seq_nr)
    seq_nr_intervals.append(seq_nr_stack)
    return seq_nr_intervals

def _is_valid_gap(missing_grn: str, present_seq_nr_grn_list: list, grns_str):
    """Check if gap is valid for GRN assignment."""
    c_closest_seq_nr=None
    n_closest_seq_nr=None
    closest_type = None
    
    region = missing_grn[0]
    present_grns = [x[1] for x in present_seq_nr_grn_list]
    std_grns = [x for x in grns_str if x in present_grns]
    
    right_of_region = [x for x in std_grns if (region == x[0]) & 
                       (parse_grn_str2float(missing_grn) <= parse_grn_str2float(x))]
    left_of_region = [x for x in std_grns if (region == x[0]) & 
                      (parse_grn_str2float(missing_grn) >= parse_grn_str2float(x))]
    
    if len(right_of_region) == 0:
        return True, False, 0
    elif len(left_of_region) == 0:
        return False, True, 0
    else:
        present_seq_nr_grn_list_c = [g for g in present_seq_nr_grn_list if g[1] in right_of_region]
        present_seq_nr_grn_list_n = [g for g in present_seq_nr_grn_list if g[1] in left_of_region]
        
        c_grn_dists = sorted(
            [(abs(int(round(100 * (parse_grn_str2float(grn) - parse_grn_str2float(missing_grn)), 3))), grn)
             for grn in right_of_region])
        n_grn_dists = sorted(
            [(abs(int(round(100 * (parse_grn_str2float(grn) - parse_grn_str2float(missing_grn)), 3))), grn)
             for grn in left_of_region])
        
        all_grn_dists = sorted(c_grn_dists + n_grn_dists)
        
        if len(c_grn_dists) > 0:
            c_closest_grn = c_grn_dists[0][1]
            c_closest_seq_nr = [int(g[0][1:]) for g in present_seq_nr_grn_list
                                if g[1] == c_closest_grn][0]
            closest_type = 'c'
        
        if len(n_grn_dists) > 0:
            n_closest_grn = n_grn_dists[0][1]
            closest_type = 'n'
            n_closest_seq_nr = [int(g[0][1:]) for g in present_seq_nr_grn_list
                                if g[1] == n_closest_grn][0]
        
        if len(all_grn_dists) > 0:
            closest_grn = all_grn_dists[0][1]
            closest_type = 'c' if parse_grn_str2float(closest_grn) > parse_grn_str2float(missing_grn) else 'n'
        
        if abs(c_closest_seq_nr - n_closest_seq_nr) == 1:
            new_seqnr_grn = []
            return False, False, new_seqnr_grn
        else:
            if closest_type == 'c':
                new_seqnr_grn = (str(c_closest_seq_nr - 1), missing_grn)
            else:
                new_seqnr_grn = (str(n_closest_seq_nr + 1), missing_grn)
            return True, True, new_seqnr_grn

def _get_closest_present_seqnr(missing_seqnr: int, present_seq_nr_grn_list: list, loop_side='n'):
    """Find closest present sequence number."""
    min_dist = 1000
    closest = None
    for g in present_seq_nr_grn_list:
        seq_nr = int(g[0][1:])
        dist = seq_nr - missing_seqnr
        
        if loop_side == 'n' and dist >= 0:
            continue
        if loop_side == 'c' and dist <= 0:
            continue
        
        abs_dist = abs(dist)
        if abs_dist < abs(min_dist):
            min_dist = abs_dist
            closest = g
    return closest, min_dist

def _get_closest_present_grn(missing_grn: str, present_seq_nr_grn_list: list):
    """Find closest present GRN."""
    region = missing_grn[0]
    missing_grn_float = parse_grn_str2float(missing_grn)
    present_grns_float = [parse_grn_str2float(g[1]) for g in present_seq_nr_grn_list if g[1][0] == region]
    grn_dists = sorted([(round(abs(missing_grn_float - grn), 3), grn) for grn in present_grns_float])
    closest_float = grn_dists[0][1]
    return parse_grn_float2str(closest_float)

def _check_interval_is_gap(interval, present_seq_nr_grn_list, grn_config, grns_str):
    """Check if interval is a gap vs loop."""
    seqnr = interval[0]
    for name, tm in grn_config.items():
        grn_std_interval = get_grn_interval(tm[0], tm[1], grns_str=grns_str)
        closest_n, min_dist_n = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, 'n')
        closest_c, min_dist_c = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, 'c')
        if (closest_n[1] in grn_std_interval) and (closest_c[1] in grn_std_interval):
            if closest_n != closest_c:
                return True
    return False

def _check_for_duplicate_grns(seq_nr_grn_list):
    """Check for duplicate GRN assignments."""
    known_grns = [g[1] for g in seq_nr_grn_list]
    known_wo_dup = list(set(known_grns))
    if len(known_grns) != len(known_wo_dup):
        for kwd in known_wo_dup:
            known_grns.remove(kwd)
        return False, known_grns
    else:
        return True, []

def _annotate_missing_rns(interval, present_seq_nr_grn_list, query_seq, grn_config, grns_str):
    """Annotate missing residue numbers as gaps or loops."""
    known_grns = [g[1] for g in present_seq_nr_grn_list]
    
    gaps = []
    nloop = []
    cloop = []
    
    if _check_interval_is_gap(interval, present_seq_nr_grn_list, grn_config=grn_config, grns_str=grns_str):
        # TREATING A GAP!
        if len(interval) == 1:
            seqnr = interval[0]
            closest, min_dist = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, 'n')
            grn = parse_grn_float2str(parse_grn_str2float(closest[1]) + .001)
            seq_id = query_seq[seqnr - 1] + str(seqnr)
            gaps = [(seq_id, grn)]
        else:
            for seqnr in interval:
                closest, min_dist = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, 'n')
                grn = parse_grn_float2str(parse_grn_str2float(closest[1]) + .001)
                seq_id = query_seq[seqnr - 1] + str(seqnr)
                present_seq_nr_grn_list += [(seq_id, grn)]
                gaps.append((seq_id, grn))
    else:
        # TREATING A LOOP
        if len(interval) == 1:
            n_interval = interval
            c_interval = []
        else:
            c_len = int(len(interval) / 2)
            n_len = len(interval) - c_len
            assert c_len <= n_len, print("c length should be smaller or equal to n length, but...", c_len, n_len)
            if len(interval) % 2 == 0:
                c_interval = interval[-(n_len):]
            else:
                c_interval = interval[-(n_len - 1):]
            n_interval = interval[:n_len]
        
        for seqnr in n_interval:
            closest, min_dist = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, loop_side='n')
            region = closest[1][0]
            loop = region + str(int(region) + 1)
            min_dist = round(int(min_dist) * .01, 3)
            grn = loop + 'x' + str(min_dist).replace('.', '')
            if grn not in known_grns:
                seq_id = query_seq[seqnr - 1] + str(seqnr)
                nloop.append((seq_id, grn))
        
        not_found_error = True
        for seqnr in c_interval:
            closest, min_dist = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, loop_side='c')
            region = closest[1][0]
            loop = region + str(int(region) - 1)
            if loop != '10':
                grn = loop + 'x' + str(round(min_dist * .01, 3)).replace('.', '')
                if grn not in known_grns:
                    seq_id = query_seq[seqnr - 1] + str(seqnr)
                    cloop.append((seq_id, grn))
                if not_found_error:
                    not_found_error, b = _check_for_duplicate_grns(cloop)
                    if not not_found_error:
                        print("NOT FOUND ERROR")
                        print(seqnr, closest, min_dist, b)
                        print('cloop', cloop)
            else:
                print('10x000 ERROR INVALID LOOP DETECTED!')
                return gaps, nloop, cloop
    return gaps, nloop, cloop

def sort_grn_rn_pairs(sorted_grns, zippy):
    """Sort GRN-residue number pairs."""
    sorted_pairs = []
    zippy_dict = {p[0]: p for p in zippy}
    
    for value in sorted_grns:
        if value in zippy_dict:
            sorted_pairs.append(zippy_dict[value])
    
    return sorted_pairs

def assign_missing_std_grns(missing_std_grns, present_seq_nr_grn_list, query_seq, missing, grns_str):
    """Assign missing standard GRNs using pivot-based approach."""
    grns = []
    for missing_std_grn in missing_std_grns:
        if 'x' in missing_std_grn:
            c_pivot, n_pivot, new_seqnr_grn = _is_valid_gap(missing_std_grn, present_seq_nr_grn_list, grns_str)
            if c_pivot & n_pivot:
                if len(new_seqnr_grn) > 0:
                    seqnr = new_seqnr_grn[0][0]
                    seqnr = query_seq[int(seqnr) - 1] + str(seqnr)
                    grns += [(seqnr, new_seqnr_grn[1])]
            
            elif c_pivot:
                closest_grn = _get_closest_present_grn(missing_std_grn, present_seq_nr_grn_list)
                dist_seq_nr = int((float(closest_grn.split('x')[1]) - float(missing_std_grn.split('x')[1])))
                closest_seqnr = [int(x[0][1:]) for x in present_seq_nr_grn_list if x[1] == closest_grn][0]
                if (closest_seqnr - dist_seq_nr) in missing:
                    seq_nr = closest_seqnr - dist_seq_nr
                    grns += [(query_seq[seq_nr - 1] + str(seq_nr), missing_std_grn)]
                    missing.remove(closest_seqnr - dist_seq_nr)
            
            elif n_pivot:
                closest_grn = _get_closest_present_grn(missing_std_grn, present_seq_nr_grn_list)
                dist_seq_nr = int((float(missing_std_grn.split('x')[1]) - float(closest_grn.split('x')[1])))
                closest_seqnr = [int(x[0][1:]) for x in present_seq_nr_grn_list if x[1] == closest_grn][0]
                if (closest_seqnr + dist_seq_nr) in missing:
                    seq_nr = closest_seqnr + dist_seq_nr
                    grns += [(query_seq[seq_nr - 1] + str(seq_nr), missing_std_grn)]
                    missing.remove(closest_seqnr + dist_seq_nr)
    return grns, missing

def annotate_gaps_and_loops(present_seq_nr_grn_list, missing, query_seq, grn_config, grns_str):
    """Annotate gaps and loops in the sequence."""
    known_seqnrs = [int(g[0][1:]) for g in present_seq_nr_grn_list]
    missing_nrs = missing
    missing_ = [x for x in missing_nrs if x not in known_seqnrs]
    
    missing_seq_nr_intervals = _get_seq_nr_intervals(missing_)
    
    cloop = []
    nloop = []
    gaps = []
    
    for interval in missing_seq_nr_intervals:
        gaps_, nloop_, cloop_ = _annotate_missing_rns(interval, present_seq_nr_grn_list, query_seq, grn_config,
                                                      grns_str)
        cloop += cloop_
        nloop += nloop_
        gaps += gaps_
    return nloop, gaps, cloop