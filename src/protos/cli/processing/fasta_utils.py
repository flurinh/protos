import argparse
import os
from pathlib import Path
from typing import Dict, Optional, List
from protos.io.fasta_utils import read_fasta, write_fasta, validate_fasta_format, clean_sequence
from protos.processing.sequence import SequenceProcessor
from protos.io.data_access import GlobalRegistry
from protos.commands.register_data import register_sequence_file


def preprocess_fasta(input_path, output_path, filter_invalid=True, clean=True, 
                    register_output=False, processor_name="preprocess"):
    """
    Preprocess a FASTA file by filtering invalid entries and cleaning sequences.
    
    Now supports the entity system - can optionally register the output file.
    
    Args:
        input_path (str): Path to the input FASTA file (can be a registered entity name)
        output_path (str): Path to save the processed FASTA file
        filter_invalid (bool, optional): Whether to filter invalid entries. Defaults to True.
        clean (bool, optional): Whether to clean sequences. Defaults to True.
        register_output (bool, optional): Whether to register the output in entity system
        processor_name (str, optional): Name for the sequence processor
    
    Returns:
        dict: Dictionary of processed sequences
    """
    # Try to resolve input as an entity first
    seq_proc = SeqProcessor(name=processor_name)
    sequences = {}
    
    # Check if input_path is a registered entity name
    registry = GlobalRegistry()
    entity_found = False
    
    for entity_id, entity_data in registry.entity_registry.entities.items():
        if entity_data.get('original_id') == input_path and 'sequence' in entity_data.get('formats', {}):
            # Load sequences from this entity
            file_path = entity_data['formats']['sequence'].get('file_path')
            if file_path and Path(file_path).exists():
                sequences = read_fasta(file_path)
                entity_found = True
                print(f"Loading sequences from registered entity: {input_path}")
                break
    
    # If not found as entity, treat as file path
    if not entity_found:
        if not Path(input_path).exists():
            raise FileNotFoundError(f"Input not found as entity or file: {input_path}")
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
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write the processed sequences to file
    write_fasta(processed_sequences, str(output_path))
    
    print(f"Processed {len(sequences)} sequences, kept {len(processed_sequences)}")
    print(f"Saved to {output_path}")
    
    # Register the output file if requested
    if register_output:
        registered = register_sequence_file(output_path)
        print(f"Registered {len(registered)} sequences in entity system")
    
    return processed_sequences


def postprocess_fasta(input_path, ref_path, output_path, register_output=False, 
                     use_entities=True, processor_name="postprocess"):
    """
    Post-process a FASTA file using a reference FASTA file.
    
    Now supports the entity system - inputs can be entity names or file paths.
    
    Args:
        input_path (str): Path to the input FASTA file (can be entity name)
        ref_path (str): Path to the reference FASTA file (can be entity name)
        output_path (str): Path to save the post-processed FASTA file
        register_output (bool, optional): Whether to register output in entity system
        use_entities (bool, optional): Whether to try resolving inputs as entities
        processor_name (str, optional): Name for the sequence processor
    
    Returns:
        dict: Dictionary of post-processed sequences
    """
    seq_proc = SeqProcessor(name=processor_name)
    registry = GlobalRegistry()
    
    # Helper function to load sequences from entity or file
    def load_sequences(path_or_entity):
        if use_entities:
            # Try to find as registered entity
            for entity_id, entity_data in registry.entity_registry.entities.items():
                if entity_data.get('original_id') == path_or_entity and 'sequence' in entity_data.get('formats', {}):
                    file_path = entity_data['formats']['sequence'].get('file_path')
                    if file_path and Path(file_path).exists():
                        print(f"Loading from registered entity: {path_or_entity}")
                        return read_fasta(file_path)
        
        # Otherwise treat as file path
        if not Path(path_or_entity).exists():
            raise FileNotFoundError(f"Not found as entity or file: {path_or_entity}")
        return read_fasta(path_or_entity)
    
    # Load input and reference sequences
    sequences = load_sequences(input_path)
    ref_sequences = load_sequences(ref_path)
    
    # Keep only sequences that exist in the reference
    postprocessed_sequences = {seq_id: sequence for seq_id, sequence in sequences.items() 
                              if seq_id in ref_sequences}
    
    # Create the output directory if it doesn't exist
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write the post-processed sequences to file
    write_fasta(postprocessed_sequences, str(output_path))
    
    print(f"Post-processed {len(sequences)} sequences, kept {len(postprocessed_sequences)}")
    print(f"Saved to {output_path}")
    
    # Register the output file if requested
    if register_output:
        registered = register_sequence_file(output_path)
        print(f"Registered {len(registered)} sequences in entity system")
    
    return postprocessed_sequences


def main_preprocess():
    """Command-line entry point for FASTA preprocessing."""
    parser = argparse.ArgumentParser(
        description="Preprocess a FASTA file with entity system support.",
        epilog="""
Examples:
  # Process a regular file
  protos preprocess-fasta -i input.fasta -o cleaned.fasta
  
  # Process a registered entity
  protos preprocess-fasta -i P12345 -o cleaned.fasta
  
  # Process and register the output
  protos preprocess-fasta -i input.fasta -o cleaned.fasta --register
        """
    )
    parser.add_argument('-i', '--input', type=str, required=True, 
                      help='Input FASTA file path or registered entity name')
    parser.add_argument('-o', '--output', type=str, required=True, 
                      help='Path to save the processed FASTA file')
    parser.add_argument('--no-filter', action='store_true', 
                      help='Disable filtering of invalid entries')
    parser.add_argument('--no-clean', action='store_true', 
                      help='Disable cleaning of sequences')
    parser.add_argument('--register', '-r', action='store_true',
                      help='Register the output file in the entity system')
    parser.add_argument('--name', type=str, default='preprocess',
                      help='Name for the sequence processor (default: preprocess)')
    
    args = parser.parse_args()
    
    preprocess_fasta(
        input_path=args.input,
        output_path=args.output,
        filter_invalid=not args.no_filter,
        clean=not args.no_clean,
        register_output=args.register,
        processor_name=args.name
    )


def main_postprocess():
    """Command-line entry point for FASTA postprocessing."""
    parser = argparse.ArgumentParser(
        description="Post-process a FASTA file using a reference with entity system support.",
        epilog="""
Examples:
  # Process using file paths
  protos postprocess-fasta -i input.fasta -r reference.fasta -o filtered.fasta
  
  # Process using registered entities
  protos postprocess-fasta -i my_sequences -r ref_sequences -o filtered.fasta
  
  # Process and register the output
  protos postprocess-fasta -i input.fasta -r ref.fasta -o filtered.fasta --register
        """
    )
    parser.add_argument('-i', '--input', type=str, required=True, 
                      help='Input FASTA file path or registered entity name')
    parser.add_argument('-r', '--reference', type=str, required=True, 
                      help='Reference FASTA file path or registered entity name')
    parser.add_argument('-o', '--output', type=str, required=True, 
                      help='Path to save the post-processed FASTA file')
    parser.add_argument('--register', action='store_true',
                      help='Register the output file in the entity system')
    parser.add_argument('--no-entities', action='store_true',
                      help='Disable entity name resolution (treat all inputs as file paths)')
    parser.add_argument('--name', type=str, default='postprocess',
                      help='Name for the sequence processor (default: postprocess)')
    
    args = parser.parse_args()
    
    postprocess_fasta(
        input_path=args.input,
        ref_path=args.reference,
        output_path=args.output,
        register_output=args.register,
        use_entities=not args.no_entities,
        processor_name=args.name
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