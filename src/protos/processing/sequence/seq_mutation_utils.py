import re
import warnings
import itertools


def parse_mutation_str(mutation, zero_index=True):
    """Parse a mutation string in the format '<original_aa>?<residue_number><mutant_residue>'. Returns the components."""
    pattern = re.compile(r"([A-Za-z])?(\d+)([A-Za-z])")
    match = pattern.match(mutation)
    if not match:
        raise ValueError(f"Invalid mutation format: {mutation}")
    original_aa, position, new_aa = match.groups()
    position = int(position) - (1 if zero_index else 0)
    return original_aa, position, new_aa


def parse_grn_mutation_str(mutation):
    pattern = re.compile(r"([A-Za-z])?(\d+x\d+|n\.\d+|c\.\d+)([A-Za-z])")
    match = pattern.match(mutation)
    if not match:
        raise ValueError(f"Invalid mutation format: {mutation}")
    original_aa, position, new_aa = match.groups()
    # If zero indexing is desired, you can adjust the position parsing logic here as needed
    return original_aa, position, new_aa


def generate_mutation_combinations(positions: list, possible_amino_acids: list[str],
                                   sequence: str = None):
    # List of all 20 standard amino acids
    all_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                       'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

    # Generate all possible combinations for each position
    all_combinations = itertools.product(*[
        [(pos, aa) for aa in (all_amino_acids if 'X' in aas or '*' in aas else aas)]
        for pos, aas in zip(positions, possible_amino_acids)
    ])

    if sequence != None:
        # Format combinations into a list of tuples
        formatted_combinations = []
        for combination in all_combinations:
            mutation_list = []
            for pos, aa in combination:
                # Fetch the original amino acid if the sequence is provided and flag is True
                original_aa = sequence[pos - 1]
                mutation_list.append(f"{original_aa}{pos}{aa}")
            formatted_combinations.append(mutation_list)
    else:
        formatted_combinations = [
            [str(pos) + aa for pos, aa in combination]
            for combination in all_combinations
        ]

    return formatted_combinations


def generate_rn_site_mutations(mutation_sites_list, sequence, name="", zero_index=True):
    """
    Generate all possible mutant sequences based on specified lists of mutation sites.

    :param mutation_sites_list: List of lists, each sublist specifying mutations in the format '<original_aa>?<residue_number><mutant_residue>'.
    :param sequence: Original sequence to mutate.
    :param name: Optional name to append to each mutation key.
    :param zero_index: Boolean indicating if the residue numbers are zero-indexed.
    :return: Dictionary of mutant sequences with mutation descriptions as keys.
    """
    results = {}

    for mutation_sites in mutation_sites_list:
        mutations = {}
        for site in mutation_sites:
            original_aa, position, new_aa = parse_mutation_str(site, zero_index=zero_index)
            if position in mutations:
                mutations[position].add(new_aa)
            else:
                mutations[position] = {new_aa}

        # Generate all combinations of mutations
        all_combinations = itertools.product(*[[(pos, aa) for aa in aas]
                                               for pos, aas in mutations.items()])

        # Process each combination
        for combination in all_combinations:
            mutant_sequence = list(sequence)
            key_parts = []
            for pos, aa in combination:
                if sequence[pos] != aa:  # Apply mutation only if it changes the sequence
                    print('applying mutation')
                    mutant_sequence[pos] = aa
                    key_parts.append(f"{sequence[pos]}{pos + 1 if zero_index else pos}{aa}")

            mutation_key = '_'.join(key_parts)
            if name:
                mutation_key = f"{name}_{mutation_key}"

            results[mutation_key] = ''.join(mutant_sequence)

    return results


def get_grn_for_rn(grn_dict, rns):
    grn_list = []
    for grn, res_name__rn in grn_dict.items():
        if int(res_name__rn[1:]) in rns:
            grn_list.append(grn)
    return grn_list


def check_rn_mutation(mutation, sequence):
    if isinstance(mutation, str):
        original_aa, pos, new_aa = parse_mutation_str(mutation)
    return sequence[pos] == new_aa


def apply_mutations_to_seq(sequence, mutations, zero_index=True):
    # Convert the sequence to a list to enable mutation
    sequence_list = list(sequence)

    # Define a pattern to extract the components of each mutation
    # Make the original amino acid optional
    for mutation in mutations:
        original_aa, position, new_aa = parse_mutation_str(mutation, zero_index=zero_index)

        # If original_aa is given, check if it matches the sequence at the position
        if original_aa:
            if sequence_list[position].upper() != original_aa.upper():
                # Issue a warning if the expected amino acid does not match
                warnings.warn(f"Expected {original_aa} at position {position + 1}, but found {sequence_list[position]}",
                              UserWarning)
                continue

        # Apply the mutation
        sequence_list[position] = new_aa

    # Join the list back into a string to form the new sequence
    return ''.join(sequence_list)


def apply_mutations_to_dict(data, mutation_list):
    """
    Apply a list of mutations to a dataframe row (dict format).
    """
    mutation_combinations = itertools.product(*[
        [(position, aa) for aa in (new_aa,)]
        for _, position, new_aa in [parse_grn_mutation_str(mutation) for mutation in mutation_list]
    ])

    for combination in mutation_combinations:
        mutated_row = data.copy()
        for pos, aa in combination:
            rn = int(mutated_row[pos][1:])
            mutated_row[pos] = aa + str(rn)
    return mutated_row


def apply_combinations_to_dict(data, combinations, name='wt'):
    seq_dict = {}
    seq_dict[name] = data
    for mutation_sites in combinations:
        seq_dict[name + '_' + '_'.join(mutation_sites)] = apply_mutations_to_dict(data, mutation_sites)
    return seq_dict
