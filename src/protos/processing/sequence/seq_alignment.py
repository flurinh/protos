from Bio.Align import substitution_matrices
import numpy as np
from Bio import Align, AlignIO
from tqdm import tqdm
import pandas as pd
import subprocess
import os
from dotenv import load_dotenv
from typing import Dict

# Load environment variables from .env file
load_dotenv()


def format_alignment(alm):
    # If alm is a tuple of (seq1, seq2, score), return as a list
    if isinstance(alm, tuple) and len(alm) == 3:
        return list(alm)
    
    # Existing logic for formatting alignment objects
    text = format(alm)
    lines = [l for l in text.split('\n') if l != '']
    
    # Check if we have enough lines before accessing them
    if len(lines) < 3:
        # Return a default format if there are not enough lines
        return ['', '', 0]
        
    if ('-' in lines[-2]) or ('.' in lines[-2]) or ('|' in lines[-2]):
        l_3 = [l for l in lines[-3].split(' ') if l != ''][-2]
        l_2 = [l for l in lines[-2].split(' ') if l != ''][-2]
        l_1 = [l for l in lines[-1].split(' ') if l != ''][-2]
    else:
        l_1, l_2, l_3 = '', '', ''

    lines = lines[:-3]
    lines_query = ''.join([line.split(' ')[-1] for line in lines if 'query' in line]) + l_1
    lines_matching = ''.join([line.split(' ')[-1] for line in lines if (not ('query' in line) and
                                                                      not ('target' in line))]) + l_2
    lines_target = ''.join([line.split(' ')[-1] for line in lines if 'target' in line]) + l_3
    return [lines_target, lines_matching, lines_query]


def align_blosum62(a, b, aligner, verbose=0):
    alignments = list(aligner.align(a, b))
    best_score = 0
    best_alignment = None
    if verbose > 0:
        for alignment in sorted(alignments):
            print("Score = %.1f:" % alignment.score)
            print(alignment)
    for alignment in alignments:
        if alignment.score > best_score:
            best_score = alignment.score
            best_alignment = alignment
    return best_alignment


def contains_forbidden_residues(sequence, forbidden_residues):
    for residue in sequence:
        if residue in forbidden_residues:
            return True
    return False


def contains_only_valid_residues(sequence):
    valid_residues = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                      'Y', 'X', '*', 'J'}
    return all(residue in valid_residues for residue in sequence)


def msa_blosum62(seqs_query: Dict[str, str], seqs_ref: Dict[str, str], aligner):
    # assign each query sequence to a reference
    best_alignments = {}
    for query in tqdm(seqs_query):
        best_align_query = 0
        best_score_query = 0
        for ref in seqs_ref:
            # Skip sequences with forbidden residues and remove empty spaces (they can be present)
            if not contains_only_valid_residues(seqs_query[query].replace(' ', '')):
                print("invalid query sequence:", seqs_query[query])
                continue
            if not contains_only_valid_residues(seqs_ref[ref].replace(' ', '')):
                print("invalid reference sequence:", seqs_ref[ref])
                continue
            best = align_blosum62(seqs_query[query].replace(' ', ''), seqs_ref[ref].replace(' ', ''),
                                  aligner, verbose=False)
            if best_score_query < best.score:
                best_align_query = best
                best_ref = ref
                best_score_query = best.score
        alignment = format_alignment(best_align_query)
        best_alignments.update({query: [best_ref, round(best_score_query, 3), alignment]})
    return best_alignments


def calc_alignment_score_restricted_area(alignment):
    try:
        miss_match_nan = alignment[1]
        def find_idx(seq):
            index = 0
            prev = '-'
            for i, mm in enumerate(seq):
                if mm != prev:
                    index = i
                    break
                else:
                    prev = mm
            return index
        if miss_match_nan[0] == '|':
            index_std_start = 0
        else:
            index_std_start = find_idx(miss_match_nan)
        if miss_match_nan[-1] == '|':
            miss_match_nan = miss_match_nan[index_std_start:]
        else:
            index_std_end = find_idx(miss_match_nan[::-1])
            miss_match_nan = miss_match_nan[index_std_start:-index_std_end]
        n_mm = miss_match_nan.count('.')
        n_gap = miss_match_nan.count('-')
        return 1 - ((n_mm + n_gap) / (len(miss_match_nan)))
    except:
        return 0


def get_best_alignment(msa, score_type='restricted'):
    """
    :param msa: blosum62 MSA
    """
    best = 0
    alignment = None
    best_name = None
    for name in msa:
        align = msa[name]
        if score_type == 'default':
            score = align[1]
        else:
            score = calc_alignment_score_restricted_area(align[2])
        if score > best:
            best_name = name
            best = score
            alignment = align[2]
    return best_name, alignment


def init_aligner(open_gap_score=-10):
    # Load the original BLOSUM62 matrix
    blosum62 = substitution_matrices.load("BLOSUM62")

    # Convert BLOSUM62 to a regular dictionary for easier manipulation
    original_alphabet = blosum62.alphabet
    extended_alphabet = original_alphabet + 'J'  # Explicitly add 'J' to the alphabet (for Leucine/Isoleucine)

    # Initialize an extended matrix with default values (for example, mismatches)
    matrix_size = len(extended_alphabet)
    extended_matrix = np.full((matrix_size, matrix_size), -1)  # Default mismatch score

    # Copy existing scores from BLOSUM62
    for i, res1 in enumerate(original_alphabet):
        for j, res2 in enumerate(original_alphabet):
            extended_matrix[i, j] = blosum62[res1, res2]

    # Set scores for 'J' based on average scores of 'L' and 'I'
    j_index = extended_alphabet.index('J')
    l_index = original_alphabet.index('L')
    i_index = original_alphabet.index('I')
    for i in range(matrix_size - 1):  # Exclude 'J' itself
        l_score = extended_matrix[l_index, i]
        i_score = extended_matrix[i_index, i]
        j_score = (l_score + i_score) / 2
        extended_matrix[j_index, i] = j_score
        extended_matrix[i, j_index] = j_score
    extended_matrix[j_index, j_index] = (extended_matrix[l_index, l_index] + extended_matrix[i_index, i_index]) / 2  # 'J' with 'J'

    # Create the updated substitution matrix
    updated_matrix = substitution_matrices.Array(alphabet=extended_alphabet, data=extended_matrix)

    # Initialize the aligner with the updated matrix
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -0.5
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = -0.1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0
    aligner.substitution_matrix = updated_matrix
    return aligner


def find_idx(seq):
    index = 0
    prev = '-'
    for i, mm in enumerate(seq):
        if mm != prev:
            index = i
            break
        else:
            prev = mm
    return index


def get_score(alignment):
    miss_match_nan = alignment[1]
    index_std_start = find_idx(miss_match_nan)
    index_std_end = find_idx(miss_match_nan[::-1])
    # Check if there's at least one non-gap character in the string
    if index_std_start < len(miss_match_nan) - index_std_end:
        miss_match_nan = miss_match_nan[index_std_start:(len(miss_match_nan)-index_std_end)]
        n_mm = miss_match_nan.count('.')
        n_gap = miss_match_nan.count('-')
        score = 1 - ((n_mm + n_gap) / (len(miss_match_nan)))
        print('Best score:', score)
    else:
        print('No non-gap characters found in the alignment string')
    return score


def load_alignment_file(file_path):
    column_names = ['query_id', 'target_id', 'sequence_identity', 'alignment_length', 'mismatches', 'gap_openings',
                    'query_start', 'query_end', 'target_start', 'target_end', 'e_value', 'bit_score']

    df = pd.read_csv(file_path, sep='\t', header=None, names=column_names)
    return df


def mmseqs2_align(query_seq, seqs, temp_folder='temp'):
    """
    Calculates Alignment scores using MMseqs2 for a query sequence and a set of reference sequences
    :param query_seq:
    :param seqs:
    :param temp_folder:
    :return:
    """
    def write_fasta_file(seqs, filename):
        with open(filename, 'w') as fasta_file:
            for key, value in seqs.items():
                fasta_file.write(f'>{key}\n{value}\n')

    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    if not os.path.exists(os.path.join(temp_folder, "mmseqs_tmp")):
        os.makedirs(os.path.join(temp_folder, "mmseqs_tmp"))

    path_mmseqs = os.getenv("MMSEQS_PATH")
    write_fasta_file(seqs, os.path.join(temp_folder, 'sequences.fasta'))

    with open(os.path.join(temp_folder, 'query.fasta'), 'w') as query_file:
        query_file.write(f'>query\n{query_seq}\n')

    subprocess.run(['wsl', path_mmseqs, 'createdb', f"{temp_folder}/sequences.fasta",
                    f"{temp_folder}/mmseqs_tmp/sequences_db"])
    subprocess.run(['wsl', path_mmseqs, 'createdb', f"{temp_folder}/query.fasta",
                    f"{temp_folder}/mmseqs_tmp/query_db"])
    subprocess.run(['wsl', path_mmseqs, 'search', f"{temp_folder}/mmseqs_tmp/query_db",
                    f"{temp_folder}/mmseqs_tmp/sequences_db", f"{temp_folder}/mmseqs_tmp/results",
                    f"{temp_folder}/mmseqs_tmp"])
    subprocess.run(['wsl', path_mmseqs, 'convertalis', f"{temp_folder}/mmseqs_tmp/query_db",
                    f"{temp_folder}/mmseqs_tmp/sequences_db", f"{temp_folder}/mmseqs_tmp/results",
                    f"{temp_folder}/alignment_results.tsv"])
    subprocess.run(['wsl', 'rm', '-rf', f"{temp_folder}/mmseqs_tmp"])

    # Load the first round alignment results into a dataframe
    alignment_df = load_alignment_file(os.path.join(temp_folder, 'alignment_results.tsv'))
    os.remove(os.path.join(temp_folder, 'alignment_results.tsv'))

    return alignment_df


def mmseqs2_align2(query_seqs: Dict[str, str], ref_seqs: Dict[str, str], temp_folder: str = 'temp'):
    def write_fasta_file(seqs, filename):
        with open(filename, 'w') as fasta_file:
            for key, value in seqs.items():
                fasta_file.write(f'>{key}\n{value}\n')

    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    if not os.path.exists(os.path.join(temp_folder, "mmseqs_tmp")):
        os.makedirs(os.path.join(temp_folder, "mmseqs_tmp"))

    # path_mmseqs = "~/MMseqs2/build/bin/mmseqs"
    path_mmseqs = os.getenv("MMSEQS_PATH")
    write_fasta_file(ref_seqs, os.path.join(temp_folder, 'ref_seqs.fasta'))
    write_fasta_file(query_seqs, os.path.join(temp_folder, 'query_seqs.fasta'))

    subprocess.run(['wsl', path_mmseqs, 'createdb', f"{temp_folder}/ref_seqs.fasta",
                    f"{temp_folder}/mmseqs_tmp/sequences_db"])
    subprocess.run(['wsl', path_mmseqs, 'createdb', f"{temp_folder}/query_seqs.fasta",
                    f"{temp_folder}/mmseqs_tmp/query_db"])
    subprocess.run(['wsl', path_mmseqs, 'search', f"{temp_folder}/mmseqs_tmp/query_db",
                    f"{temp_folder}/mmseqs_tmp/sequences_db", f"{temp_folder}/mmseqs_tmp/results",
                    f"{temp_folder}/mmseqs_tmp"])
    subprocess.run(['wsl', path_mmseqs, 'convertalis', f"{temp_folder}/mmseqs_tmp/query_db",
                    f"{temp_folder}/mmseqs_tmp/sequences_db", f"{temp_folder}/mmseqs_tmp/results",
                    f"{temp_folder}/alignment_results.tsv"])
    subprocess.run(['wsl', 'rm', '-rf', f"{temp_folder}/mmseqs_tmp"])

    # Load the first round alignment results into a dataframe
    alignment_df = load_alignment_file(os.path.join(temp_folder, 'alignment_results.tsv'))
    os.remove(os.path.join(temp_folder, 'alignment_results.tsv'))

    return alignment_df


def clean_pdb_seq(sequence):
    """
    Removes empty spaces that are inherent to chain-sequences in structures (e.g. HOH, ligands etc)
    :param sequence:
    :return:
    """
    return sequence.replace(' ', '')


def check_chain_similarity(chain_dict, ref_dict, min_score=.3):
    output = mmseqs2_align2(chain_dict, ref_dict)
    idx = output.groupby('query_id')['sequence_identity'].idxmax()
    max_alm = output.loc[idx]
    # Chains that meet the minimum score threshold
    gpcr_chains = max_alm[max_alm['sequence_identity'] >= min_score]['query_id'].unique()
    # Dictionary with pdb_id_chain_id as keys and True/False as values
    gpcr_dict = {key: key in gpcr_chains for key in chain_dict.keys()}
    return gpcr_dict


# MSA using ClustalW2
def run_clustalw2_alignment(seqs, temp_folder='temp', clustalw_path='clustalw2'):
    def write_fasta_file(seqs, filename):
        with open(filename, 'w') as fasta_file:
            for key, value in seqs.items():
                fasta_file.write(f'>{key}\n{value}\n')

    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)

    fasta_file_path = os.path.join(temp_folder, 'sequences.fasta')
    alignment_file_path = os.path.join(temp_folder, 'aligned_sequences.aln')

    # Write sequences to a temporary FASTA file
    write_fasta_file(seqs, fasta_file_path)

    # Run ClustalW2 alignment
    subprocess.run([clustalw_path, '-INFILE=' + fasta_file_path, '-OUTFILE=' + alignment_file_path, '-ALIGN'])

    # Load the alignment result into a BioPython alignment object
    alignment = AlignIO.read(alignment_file_path, 'clustal')

    # Clean up temporary files
    os.remove(fasta_file_path)
    os.remove(alignment_file_path)

    return alignment