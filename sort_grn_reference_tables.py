#!/usr/bin/env python3
"""
Script to load and sort GRN reference tables with proper column ordering.

This script:
1. Loads reference GRN tables (mo_ref.csv, gpcrdb_ref.csv)
2. Sorts columns according to GRN conventions (N-term, helices with loops, C-term)
3. Saves the sorted tables back to the reference data directory
"""

import pandas as pd
from pathlib import Path
import sys

# Add src to path for imports
sys.path.insert(0, 'src')

from protos.processing.grn.grn_utils import sort_grns_str
from protos.io.paths import ProtosPaths


def sort_grn_reference_table(file_path: Path, output_path: Path = None):
    """
    Load and sort a GRN reference table.
    
    Args:
        file_path: Path to the input CSV file
        output_path: Path to save sorted table (defaults to overwriting input)
    """
    if output_path is None:
        output_path = file_path
        
    print(f"\nProcessing: {file_path}")
    
    # Load the table without converting any values to NaN
    df = pd.read_csv(file_path, index_col=0, keep_default_na=False)
    print(f"  Loaded table with shape: {df.shape}")
    
    # Get column names (excluding index)
    columns = df.columns.tolist()
    print(f"  Original columns ({len(columns)}): {columns[:5]}...{columns[-5:]}")
    
    # Identify GRN columns vs metadata columns
    metadata_cols = []
    grn_cols = []
    
    for col in columns:
        # Check if it's a GRN position
        if any([
            col.startswith('n.'),      # N-terminal
            col.startswith('c.'),      # C-terminal  
            col[0].isdigit() and ('.' in col or 'x' in col)  # Standard GRN
        ]):
            grn_cols.append(col)
        else:
            metadata_cols.append(col)
    
    print(f"  Found {len(grn_cols)} GRN columns and {len(metadata_cols)} metadata columns")
    
    # Sort GRN columns using the proper algorithm
    sorted_grn_cols = sort_grns_str(grn_cols)
    print(f"  Sorted GRN columns: {sorted_grn_cols[:5]}...{sorted_grn_cols[-5:]}")
    
    # Combine metadata columns first, then sorted GRN columns
    sorted_columns = metadata_cols + sorted_grn_cols
    
    # Reorder the dataframe
    df_sorted = df[sorted_columns]
    
    # Save with '-' for missing values
    df_sorted.to_csv(output_path, na_rep='-')
    print(f"  Saved sorted table to: {output_path}")
    
    return df_sorted


def main():
    """Main function to sort reference GRN tables."""
    
    # Set up paths to reference data
    ref_data_root = Path("src/protos/reference_data")
    grn_ref_dir = ref_data_root / "grn" / "ref"
    
    if not grn_ref_dir.exists():
        print(f"Error: Reference directory not found: {grn_ref_dir}")
        return 1
        
    # Create output directory
    output_dir = Path("sorted_grn_tables")
    output_dir.mkdir(exist_ok=True)
    
    # Process mo_ref.csv
    mo_ref_path = grn_ref_dir / "mo_ref.csv"
    if mo_ref_path.exists():
        mo_output = output_dir / "mo_ref_sorted.csv"
        df_mo = sort_grn_reference_table(mo_ref_path, mo_output)
        
        # Show some statistics
        print("\nmo_ref.csv statistics:")
        print(f"  Total proteins: {len(df_mo)}")
        print(f"  Sample proteins: {list(df_mo.index[:5])}")
        
    else:
        print(f"Warning: mo_ref.csv not found at {mo_ref_path}")
    
    # Process gpcrdb_ref.csv
    gpcrdb_ref_path = grn_ref_dir / "gpcrdb_ref.csv"
    if gpcrdb_ref_path.exists():
        gpcrdb_output = output_dir / "gpcrdb_ref_sorted.csv"
        df_gpcrdb = sort_grn_reference_table(gpcrdb_ref_path, gpcrdb_output)
        
        print("\ngpcrdb_ref.csv statistics:")
        print(f"  Total proteins: {len(df_gpcrdb)}")
        print(f"  Sample proteins: {list(df_gpcrdb.index[:5])}")
    else:
        print(f"Warning: gpcrdb_ref.csv not found at {gpcrdb_ref_path}")
    
    print("\nDone! Reference tables have been sorted.")
    
    # Also update test data copies if they exist
    test_data_root = Path("tests/test-data")
    if test_data_root.exists():
        test_grn_ref_dir = test_data_root / "grn" / "ref"
        if test_grn_ref_dir.exists():
            print("\nUpdating test data copies...")
            
            # Copy sorted mo_ref.csv to test data
            test_mo_ref = test_grn_ref_dir / "mo_ref.csv"
            if mo_ref_path.exists() and test_mo_ref.parent.exists():
                import shutil
                shutil.copy2(mo_ref_path, test_mo_ref)
                print(f"  Updated: {test_mo_ref}")
                
            # Copy sorted gpcrdb_ref.csv to test data
            test_gpcrdb_ref = test_grn_ref_dir / "gpcrdb_ref.csv"
            if gpcrdb_ref_path.exists() and test_gpcrdb_ref.parent.exists():
                import shutil
                shutil.copy2(gpcrdb_ref_path, test_gpcrdb_ref)
                print(f"  Updated: {test_gpcrdb_ref}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())