AA_DICT = {
    'A': ('ala', 'alanine'),
    'R': ('arg', 'arginine'),
    'N': ('asn', 'asparagine'),
    'D': ('asp', 'aspartic acid'),
    'C': ('cys', 'cysteine'),
    'Q': ('gln', 'glutamine'),
    'E': ('glu', 'glutamic acid'),
    'G': ('gly', 'glycine'),
    'H': ('his', 'histidine'),
    'I': ('ile', 'isoleucine'),
    'L': ('leu', 'leucine'),
    'K': ('lys', 'lysine'),
    'M': ('met', 'methionine'),
    'F': ('phe', 'phenylalanine'),
    'P': ('pro', 'proline'),
    'S': ('ser', 'serine'),
    'T': ('thr', 'threonine'),
    'W': ('trp', 'tryptophan'),
    'Y': ('tyr', 'tyrosine'),
    'V': ('val', 'valine'),
}

AA_CHARACTERISTICS = {'Aromatic': ['F', 'W', 'Y'],
                      'Aliphatic': ['A', 'G', 'I', 'L', 'P', 'V'],
                      'Acidic': ['D', 'E'],
                      'Basic': ['R', 'H', 'K'],
                      'Hydroxylic': ['S', 'T'],
                      'Sulphur-containing': ['C', 'M'],
                      'Amidic': ['N', 'Q'],
                      'Hydrophobic': ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W'],
                      'Hydrophilic': ['N', 'D', 'Q', 'E', 'K', 'R']}
