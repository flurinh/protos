#!/usr/bin/env python3
"""
Extract sequences from a GRN table using GRNProcessor.

This script demonstrates how to:
1. Initialize GRNProcessor with proper paths
2. Load a GRN table
3. Extract sequences using get_seq_dict()
"""

import os
import sys
from pathlib import Path

# Add the protos source to the path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from protos.processing.grn.grn_processor import GRNProcessor
from protos.io.paths.path_config import ProtosPaths


def extract_sequences_from_grn_table(table_name):
    """
    Extract sequences from a GRN table.
    
    Args:
        table_name: Name of the GRN table to load (without .csv extension)
    
    Returns:
        Dictionary mapping protein IDs to sequences
    """
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/protos_data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize paths
    paths = ProtosPaths(
        data_root=str(data_dir.absolute())
    )
    
    # Create GRN processor
    print(f"Initializing GRN processor...")
    grn_processor = GRNProcessor(
        name="sequence_extractor",
        paths=paths
    )
    
    # Load the GRN table
    print(f"Loading GRN table: {table_name}")
    try:
        grn_processor.load_grn_table(table_name)
        print(f"Successfully loaded table with {len(grn_processor.ids)} sequences")
    except FileNotFoundError:
        print(f"Error: GRN table '{table_name}' not found")
        print(f"Available tables in {paths.get_subdir_path('grn', 'table_dir')}:")
        table_dir = Path(paths.get_subdir_path('grn', 'table_dir'))
        if table_dir.exists():
            for table_file in table_dir.glob("*.csv"):
                print(f"  - {table_file.stem}")
        return None
    
    # Extract sequences
    print("Extracting sequences...")
    seq_dict = grn_processor.get_seq_dict()
    
    # Display some example sequences
    print(f"\nExtracted {len(seq_dict)} sequences")
    print("\nFirst 5 sequences:")
    for i, (protein_id, sequence) in enumerate(seq_dict.items()):
        if i >= 5:
            break
        print(f"\n{protein_id}:")
        print(f"  Length: {len(sequence)}")
        print(f"  Sequence: {sequence[:60]}{'...' if len(sequence) > 60 else ''}")
    
    return seq_dict


def main():
    """Main function to demonstrate sequence extraction."""
    # Example usage with different tables
    tables_to_try = [
        "reference_grn_table",  # Common reference table
        "mo_reference_table",   # Microbial opsin reference
        "gpcr_reference_table", # GPCR reference
    ]
    
    print("GRN Sequence Extractor")
    print("=" * 50)
    
    # Try to find and load a table
    seq_dict = None
    for table_name in tables_to_try:
        print(f"\nTrying table: {table_name}")
        seq_dict = extract_sequences_from_grn_table(table_name)
        if seq_dict:
            break
    
    if seq_dict:
        print(f"\nSuccessfully extracted sequences!")
        print(f"You can now use the seq_dict to access any protein sequence by ID")
        
        # Example: Access a specific sequence
        if seq_dict:
            first_id = list(seq_dict.keys())[0]
            print(f"\nExample - accessing sequence for '{first_id}':")
            print(f"seq_dict['{first_id}'] = {seq_dict[first_id][:80]}...")
    else:
        print("\nNo tables found. Please ensure GRN tables are available in the data directory.")


if __name__ == "__main__":
    main()