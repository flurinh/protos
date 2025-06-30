import argparse
import pandas as pd
import re
import os
from pathlib import Path
from protos.processing.grn.grn_base_processor import GRNBaseProcessor


def validate_and_clean_row(row):
    """
    Validate and clean a row of GRN data.
    
    Args:
        row (list): A row of cell values
    
    Returns:
        tuple: (cleaned_row, is_erroneous)
    """
    clean_row = []
    previous_number = 0
    sequence_restart_found = False
    erroneous_sequence = False

    for cell in row:
        if cell == '-' or not isinstance(cell, str):
            clean_row.append(cell)
            continue

        # Extract the residue number
        match = re.match(r'([A-Za-z]+)(\d+)', cell)
        if not match:
            clean_row.append(cell)
            continue

        residue, number = match.groups()
        number = int(number)

        if sequence_restart_found:
            clean_row.append('-')
            continue

        if number == 1 and previous_number != 0:
            # Sequence restart detected, truncate the sequence here
            sequence_restart_found = True
            clean_row.append('-')
            continue

        if number != previous_number + 1 and previous_number != 0:
            # Erroneous sequence detected
            erroneous_sequence = True

        clean_row.append(cell)
        previous_number = number

    return clean_row, erroneous_sequence


def process_table(input_path, output_path, use_processor=True, include_entity_ids=True):
    """
    Process a GRN table by validating and cleaning rows.
    
    Now supports the entity system - GRN tables include entity_id column.
    
    Args:
        input_path (str): Path to the input CSV file (or dataset ID if use_processor=True)
        output_path (str): Path to save the cleaned CSV file (or dataset ID if use_processor=True)
        use_processor (bool): Whether to use GRNBaseProcessor for loading/saving
        include_entity_ids (bool): Whether to preserve entity_id column in output
    
    Returns:
        dict: Report of erroneous sequences
    """
    if use_processor:
        # Use GRNBaseProcessor for loading
        processor = GRNBaseProcessor(
            name='clean_grn',
            processor_data_dir='grn'
        )
        processor.load_grn_table(input_path)
        df = processor.data
    else:
        # Load directly from file path
        df = pd.read_csv(input_path, index_col=0)

    # Check if entity_id column exists
    has_entity_ids = 'entity_id' in df.columns
    entity_ids = None
    if has_entity_ids:
        # Preserve entity IDs
        entity_ids = df['entity_id']
        # Remove entity_id column temporarily for processing
        processing_df = df.drop(columns=['entity_id'])
    else:
        processing_df = df

    erroneous_sequences_report = {}

    # Process each row
    for index, row in processing_df.iterrows():
        clean_row, is_erroneous = validate_and_clean_row(row)

        # Update the row in the dataframe
        processing_df.loc[index] = clean_row

        # Record any erroneous sequences
        if is_erroneous:
            erroneous_sequences_report[index] = clean_row

    # Restore entity_id column if present and requested
    if has_entity_ids and include_entity_ids:
        processing_df['entity_id'] = entity_ids

    # Save the cleaned dataframe
    if use_processor:
        # Update processor data and save
        processor.data = processing_df
        processor.save_grn_table(
            dataset_id=output_path,
            include_entity_ids=include_entity_ids and has_entity_ids
        )
    else:
        processing_df.to_csv(output_path)

    # Report erroneous sequences
    if erroneous_sequences_report:
        print("Erroneous sequences found in the following rows:")
        for index, sequence in erroneous_sequences_report.items():
            print(f"Row {index}: {sequence}")
    else:
        print("No erroneous sequences found.")
        
    return erroneous_sequences_report


def clean_grn_table(input_path, output_path):
    """
    Clean a GRN table by removing erroneous sequences.
    
    Args:
        input_path (str): Path to the input CSV file
        output_path (str): Path to save the cleaned CSV file
    
    Returns:
        dict: Report of erroneous sequences
    """
    return process_table(input_path, output_path)


def main():
    """Command-line entry point for GRN table cleaning."""
    parser = argparse.ArgumentParser(
        description="Clean and validate residue sequences in a GRN table with entity system support.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Clean a GRN table using dataset IDs
  python -m protos.cli.grn.clean_grn_table -i ref/mo_ref -o ref/mo_ref_cleaned
  
  # Clean a GRN table using file paths
  python -m protos.cli.grn.clean_grn_table -i input.csv -o output.csv --use-files
  
  # Clean and preserve entity IDs
  python -m protos.cli.grn.clean_grn_table -i grn_assigned -o grn_cleaned --preserve-entities
  
  # Clean with custom data root
  python -m protos.cli.grn.clean_grn_table -i my_table -o my_table_clean --data-root /path/to/data
        """
    )
    parser.add_argument('-i', '--input_path', type=str, required=True, 
                      help='Input dataset ID (or file path if --use-files)')
    parser.add_argument('-o', '--output_path', type=str, required=True, 
                      help='Output dataset ID (or file path if --use-files)')
    parser.add_argument('--use-files', action='store_true',
                      help='Use file paths instead of dataset IDs')
    parser.add_argument('--preserve-entities', action='store_true',
                      help='Preserve entity_id column if present')
    parser.add_argument('--data-root', type=str, default=None,
                      help='Root data directory (default: uses PROTOS_DATA_ROOT env var)')

    args = parser.parse_args()
    
    # Set environment variable if data root provided
    if args.data_root:
        os.environ['PROTOS_DATA_ROOT'] = str(Path(args.data_root).absolute())

    process_table(
        args.input_path, 
        args.output_path, 
        use_processor=not args.use_files,
        include_entity_ids=args.preserve_entities
    )
    print(f"Cleaned table saved to {args.output_path}")


if __name__ == "__main__":
    main()