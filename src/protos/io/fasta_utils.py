import hashlib
import pandas as pd


def read_fasta(file_path):
    sequences = {}
    current_sequence_id = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # Extract only the first word as the sequence ID
                current_sequence_id = line[1:].split()[0]
                sequences[current_sequence_id] = ''
            else:
                sequences[current_sequence_id] += line
    return sequences


def write_fasta(sequences, file_path):
    with open(file_path, 'w') as file:
        for sequence_id, sequence in sequences.items():
            file.write('>' + sequence_id + '\n')
            # Write sequence with line breaks every 80 characters
            for i in range(0, len(sequence), 80):
                file.write(sequence[i:i+80] + '\n')


def generate_protein_id(descriptor, protein_family='O', descriptor_length=10):
    # Create a hash of the descriptor and return the first 8 characters
    return protein_family + hashlib.md5(descriptor.encode()).hexdigest()[:descriptor_length-1]


def create_fasta_id_map(fasta_dict, descriptor_length=10):
    df = pd.DataFrame([(generate_protein_id(k, descriptor_length=descriptor_length), k) for k in fasta_dict.keys()])
    df.columns = ['protein_id', 'description']
    return df


def gen_input_files(fasta_path, filename, out_dir='input'):
    seq_dict = read_fasta(fasta_path + filename + '.fasta')
    df = create_fasta_id_map(seq_dict)
    seq_dict_mapped_protein_id = {generate_protein_id(k): v for k, v in seq_dict.items()}
    df.to_csv(out_dir + '/' + filename+'_'+'ids.csv', index=False)
    write_fasta(seq_dict_mapped_protein_id, out_dir+'/'+filename + '.fasta')

def initialize_properties(dataset=None):
    # go to fasta/processed, use the _ids files, and set up an empty dataset that we can use to store the predictions in
    pass

def validate_hashkey_consistency(dataset):
    # check which exist: grnp, embp, pp
    # check if we have consistent keys in all of them
    pass


def clean_sequence(seq):
    """Clean a sequence by removing whitespace and invalid characters and converting to uppercase."""
    # Remove whitespace
    seq = ''.join(seq.split())
    
    # Keep only valid amino acid characters
    valid_chars = 'ACDEFGHIKLMNPQRSTVWXY'
    seq = ''.join(c for c in seq.upper() if c in valid_chars)
    
    return seq


def validate_fasta_format(content):
    """Validate that a string is in proper FASTA format."""
    if not content.strip():
        return False
        
    lines = content.strip().split('\n')
    
    # First line must start with >
    if not lines[0].startswith('>'):
        return False
    
    current_header = None
    current_sequence = ""
    
    for line in lines:
        if line.startswith('>'):
            # If we already had a header, check that it had a sequence
            if current_header is not None and not current_sequence:
                return False
            
            current_header = line[1:]
            current_sequence = ""
        else:
            current_sequence += line.strip()
    
    # Check the last sequence
    if not current_sequence:
        return False
    
    return True