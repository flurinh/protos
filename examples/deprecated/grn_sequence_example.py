#!/usr/bin/env python3
"""
Example of using GRNProcessor to load sequences for annotation.

This shows the pattern used in the annotation scripts:
1. Load GRN table using GRNProcessor
2. Extract sequences with get_seq_dict()
3. Use sequences for alignment/annotation
"""

import os
import sys
from pathlib import Path

# Add the protos source to the path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from protos.processing.grn.grn_processor import GRNProcessor
from protos.io.paths.path_config import ProtosPaths


def load_reference_sequences(grn_table_name="mo_ref"):
    """
    Load reference sequences from a GRN table.
    
    This is the pattern used in annotation scripts to load reference sequences
    for alignment and GRN transfer.
    
    Args:
        grn_table_name: Name of the GRN table to load
        
    Returns:
        Tuple of (grn_processor, seq_dict)
    """
    # Set up paths - CRITICAL for proper data loading
    data_dir = Path("/mnt/c/Users/hidbe/protos_data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize paths
    paths = ProtosPaths(
        data_root=str(data_dir.absolute())
    )
    
    # Create GRN processor
    grn_processor = GRNProcessor(
        name="annotation_processor",
        paths=paths
    )
    
    # Load the GRN table
    grn_processor.load_grn_table(grn_table_name)
    
    # Get sequences
    seq_dict = grn_processor.get_seq_dict()
    
    return grn_processor, seq_dict


def demonstrate_usage():
    """Demonstrate how to use GRN sequences in annotation context."""
    
    print("Loading reference sequences from GRN table...")
    grn_processor, seq_dict = load_reference_sequences()
    
    print(f"\nLoaded {len(seq_dict)} reference sequences")
    print(f"GRN columns available: {len(grn_processor.grns)} positions")
    
    # Example: Access a specific reference sequence for alignment
    print("\nExample usage for annotation:")
    
    # Get first reference sequence
    ref_id = list(seq_dict.keys())[0]
    ref_seq = seq_dict[ref_id]
    
    print(f"\nReference: {ref_id}")
    print(f"Sequence: {ref_seq[:60]}...")
    print(f"Length: {len(ref_seq)}")
    
    # Get the GRN row for this reference
    ref_grn_row = grn_processor.data.loc[ref_id]
    
    # Show some GRN positions
    print(f"\nGRN annotations for {ref_id}:")
    grn_positions = list(grn_processor.grns)[:10]  # First 10 positions
    for grn in grn_positions:
        residue = ref_grn_row[grn]
        if residue != '-':
            print(f"  {grn}: {residue}")
    
    # Example query sequence (would come from input file in real annotation)
    query_seq = "MVLDAEFRRSADQLGVSRPLIVIGGGESSPPGQVNRLKRWIHGLENVIFVAHMDHCKGMR"
    
    print(f"\nQuery sequence to annotate:")
    print(f"Sequence: {query_seq}")
    print(f"Length: {len(query_seq)}")
    
    # In the annotation scripts, this is where you would:
    # 1. Use MMseqs2 to find best matching reference
    # 2. Perform alignment (with Biopython or MMseqs2)
    # 3. Transfer GRN annotations based on alignment
    
    return grn_processor, seq_dict


if __name__ == "__main__":
    demonstrate_usage()