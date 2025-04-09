import argparse
import os
from protos.io.fasta_utils import read_fasta, write_fasta, validate_fasta_format, clean_sequence


def preprocess_fasta(input_path, output_path, filter_invalid=True, clean=True):
    """
    Preprocess a FASTA file by filtering invalid entries and cleaning sequences.
    
    Args:
        input_path (str): Path to the input FASTA file
        output_path (str): Path to save the processed FASTA file
        filter_invalid (bool, optional): Whether to filter invalid entries. Defaults to True.
        clean (bool, optional): Whether to clean sequences. Defaults to True.
    
    Returns:
        dict: Dictionary of processed sequences
    """
    # Read the input FASTA file
    sequences = read_fasta(input_path)
    
    # Process the sequences
    processed_sequences = {}
    for seq_id, sequence in sequences.items():
        # Skip invalid entries if requested
        if filter_invalid and not validate_fasta_format(f">{seq_id}\n{sequence}"):
            print(f"Skipping invalid entry: {seq_id}")
            continue
        
        # Clean the sequence if requested
        if clean:
            sequence = clean_sequence(sequence)
            
        # Skip sequences that became empty after cleaning
        if not sequence:
            print(f"Skipping empty sequence after cleaning: {seq_id}")
            continue
            
        processed_sequences[seq_id] = sequence
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    
    # Write the processed sequences to file
    write_fasta(processed_sequences, output_path)
    
    print(f"Processed {len(sequences)} sequences, kept {len(processed_sequences)}")
    print(f"Saved to {output_path}")
    
    return processed_sequences


def postprocess_fasta(input_path, ref_path, output_path):
    """
    Post-process a FASTA file using a reference FASTA file.
    
    Args:
        input_path (str): Path to the input FASTA file
        ref_path (str): Path to the reference FASTA file
        output_path (str): Path to save the post-processed FASTA file
    
    Returns:
        dict: Dictionary of post-processed sequences
    """
    # Read the input and reference FASTA files
    sequences = read_fasta(input_path)
    ref_sequences = read_fasta(ref_path)
    
    # Keep only sequences that exist in the reference
    postprocessed_sequences = {seq_id: sequence for seq_id, sequence in sequences.items() 
                              if seq_id in ref_sequences}
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    
    # Write the post-processed sequences to file
    write_fasta(postprocessed_sequences, output_path)
    
    print(f"Post-processed {len(sequences)} sequences, kept {len(postprocessed_sequences)}")
    print(f"Saved to {output_path}")
    
    return postprocessed_sequences


def main_preprocess():
    """Command-line entry point for FASTA preprocessing."""
    parser = argparse.ArgumentParser(description="Preprocess a FASTA file.")
    parser.add_argument('-i', '--input', type=str, required=True, 
                      help='Path to the input FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, 
                      help='Path to save the processed FASTA file')
    parser.add_argument('--no-filter', action='store_true', 
                      help='Disable filtering of invalid entries')
    parser.add_argument('--no-clean', action='store_true', 
                      help='Disable cleaning of sequences')
    
    args = parser.parse_args()
    
    preprocess_fasta(
        input_path=args.input,
        output_path=args.output,
        filter_invalid=not args.no_filter,
        clean=not args.no_clean
    )


def main_postprocess():
    """Command-line entry point for FASTA postprocessing."""
    parser = argparse.ArgumentParser(description="Post-process a FASTA file using a reference.")
    parser.add_argument('-i', '--input', type=str, required=True, 
                      help='Path to the input FASTA file')
    parser.add_argument('-r', '--reference', type=str, required=True, 
                      help='Path to the reference FASTA file')
    parser.add_argument('-o', '--output', type=str, required=True, 
                      help='Path to save the post-processed FASTA file')
    
    args = parser.parse_args()
    
    postprocess_fasta(
        input_path=args.input,
        ref_path=args.reference,
        output_path=args.output
    )


if __name__ == "__main__":
    # This script can be used for either preprocessing or postprocessing
    # Determine which operation to perform based on the script name
    if os.path.basename(__file__).startswith('preprocess'):
        main_preprocess()
    elif os.path.basename(__file__).startswith('postprocess'):
        main_postprocess()
    else:
        print("Script name must start with 'preprocess' or 'postprocess'")
        print("Defaulting to preprocessing")
        main_preprocess()