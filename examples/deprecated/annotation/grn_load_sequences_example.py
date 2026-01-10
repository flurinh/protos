#!/usr/bin/env python3
"""
Complete example showing how to load GRN tables and extract sequences.

This demonstrates two methods:
1. Loading from the tables/ directory using load_grn_table()
2. Loading directly from ref/ directory (as done in annotation scripts)
"""

import os
import sys
import pandas as pd
from pathlib import Path

# Add the protos source to the path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from protos.processing.grn.grn_processor import GRNProcessor
from protos.io.paths.path_config import ProtosPaths


def method1_load_from_tables_directory():
    """Method 1: Load GRN table from tables/ directory using load_grn_table()."""
    print("\n" + "="*60)
    print("METHOD 1: Loading from tables/ directory")
    print("="*60)
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/protos_data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize paths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Create GRN processor
    grn_processor = GRNProcessor(name="method1", paths=paths)
    
    # Check available tables
    table_dir = Path(paths.get_subdir_path('grn', 'table_dir'))
    print(f"\nAvailable tables in {table_dir}:")
    if table_dir.exists():
        tables = [f.stem for f in table_dir.glob("*.csv") if f.stem != "ref"]
        for table in tables[:5]:  # Show first 5
            print(f"  - {table}")
    
    # Try to load a table
    table_name = "test_grn"  # Choose an existing table
    try:
        grn_processor.load_grn_table(table_name)
        print(f"\nSuccessfully loaded '{table_name}'")
        print(f"Number of sequences: {len(grn_processor.ids)}")
        print(f"Number of GRN positions: {len(grn_processor.grns)}")
        
        # Get sequences
        seq_dict = grn_processor.get_seq_dict()
        print(f"\nExtracted {len(seq_dict)} sequences")
        
        # Show example
        if seq_dict:
            first_id = list(seq_dict.keys())[0]
            print(f"\nExample - {first_id}:")
            print(f"  Sequence: {seq_dict[first_id][:60]}...")
            
    except FileNotFoundError:
        print(f"Table '{table_name}' not found. Please check available tables above.")


def method2_load_from_ref_directory():
    """Method 2: Load reference table directly from ref/ directory (as in annotation scripts)."""
    print("\n" + "="*60)
    print("METHOD 2: Loading from ref/ directory (annotation script style)")
    print("="*60)
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/protos_data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize paths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Create GRN processor
    grn_processor = GRNProcessor(name="method2", paths=paths)
    
    # Load reference table directly
    ref_file = grn_processor.path_ref_dir / "mo_ref.csv"
    
    if ref_file.exists():
        print(f"\nLoading reference from: {ref_file}")
        
        # Load and clean the reference table (same as annotation scripts)
        grn_processor.data = pd.read_csv(ref_file, index_col=0)
        grn_processor.data = grn_processor.data.fillna('-')  # Replace NaN with '-'
        grn_processor.ids = grn_processor.data.index.tolist()
        grn_processor.grns = grn_processor.data.columns.tolist()
        
        print(f"Loaded {len(grn_processor.data)} reference sequences")
        print(f"Number of GRN positions: {len(grn_processor.data.columns)}")
        
        # Get reference sequences
        ref_sequences = grn_processor.get_seq_dict()
        print(f"\nBuilt {len(ref_sequences)} reference sequences")
        
        # Show first few sequences
        print("\nFirst 3 reference sequences:")
        for i, (ref_id, ref_seq) in enumerate(ref_sequences.items()):
            if i >= 3:
                break
            print(f"\n{ref_id}:")
            print(f"  Length: {len(ref_seq)}")
            print(f"  Sequence: {ref_seq[:60]}...")
            
            # Show some GRN annotations
            ref_grn_row = grn_processor.data.loc[ref_id]
            grn_positions = [grn for grn in grn_processor.grns[:10] 
                           if ref_grn_row[grn] != '-']
            if grn_positions:
                print("  GRN annotations (first few):")
                for grn in grn_positions[:5]:
                    print(f"    {grn}: {ref_grn_row[grn]}")
    else:
        print(f"Reference file not found: {ref_file}")
        print("\nAvailable reference files:")
        ref_dir = grn_processor.path_ref_dir
        if ref_dir.exists():
            for ref_file in ref_dir.glob("*.csv"):
                print(f"  - {ref_file.name}")


def demonstrate_grn_transfer_setup():
    """Show how sequences are used for GRN annotation transfer."""
    print("\n" + "="*60)
    print("GRN ANNOTATION TRANSFER SETUP")
    print("="*60)
    
    # Load reference sequences using method 2
    data_dir = Path("/mnt/c/Users/hidbe/protos_data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    grn_processor = GRNProcessor(name="annotation", paths=paths)
    
    # Load reference
    ref_file = grn_processor.path_ref_dir / "mo_ref.csv"
    if ref_file.exists():
        grn_processor.data = pd.read_csv(ref_file, index_col=0)
        grn_processor.data = grn_processor.data.fillna('-')
        grn_processor.ids = grn_processor.data.index.tolist()
        grn_processor.grns = grn_processor.data.columns.tolist()
        
        ref_sequences = grn_processor.get_seq_dict()
        
        # Example query sequence
        query_seq = "MVLDAEFRRSADQLGVSRPLIVIGGGESSPPGQVNRLKRWIHGLENVIFVAHMDHCKGMR"
        
        print(f"\nQuery sequence to annotate:")
        print(f"  {query_seq}")
        print(f"  Length: {len(query_seq)}")
        
        print(f"\nReference sequences available: {len(ref_sequences)}")
        print("\nIn the annotation pipeline:")
        print("1. This query would be aligned against all references using MMseqs2")
        print("2. The best matching reference would be selected")
        print("3. Full alignment would be performed (Biopython or MMseqs2)")
        print("4. GRN annotations would be transferred based on the alignment")
        
        return grn_processor, ref_sequences
    
    return None, None


def main():
    """Main function demonstrating all methods."""
    print("GRN Sequence Loading Examples")
    print("This demonstrates how to load sequences from GRN tables\n")
    
    # Method 1: Using load_grn_table()
    method1_load_from_tables_directory()
    
    # Method 2: Direct CSV loading (annotation script style)
    method2_load_from_ref_directory()
    
    # Show how it's used for annotation
    demonstrate_grn_transfer_setup()


if __name__ == "__main__":
    main()