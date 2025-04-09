import argparse
import pandas as pd
import re


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


def process_table(input_path, output_path):
    """
    Process a GRN table by validating and cleaning rows.
    
    Args:
        input_path (str): Path to the input CSV file
        output_path (str): Path to save the cleaned CSV file
    
    Returns:
        dict: Report of erroneous sequences
    """
    # Load the table
    df = pd.read_csv(input_path, index_col=0)

    erroneous_sequences_report = {}

    # Process each row
    for index, row in df.iterrows():
        clean_row, is_erroneous = validate_and_clean_row(row)

        # Update the row in the dataframe
        df.loc[index] = clean_row

        # Record any erroneous sequences
        if is_erroneous:
            erroneous_sequences_report[index] = clean_row

    # Save the cleaned dataframe
    df.to_csv(output_path)

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
    parser = argparse.ArgumentParser(description="Clean and validate residue sequences in a GRN table.")
    parser.add_argument('-i', '--input_path', type=str, required=True, 
                      help='Path to the input CSV file')
    parser.add_argument('-o', '--output_path', type=str, required=True, 
                      help='Path to save the cleaned CSV file')

    args = parser.parse_args()

    process_table(args.input_path, args.output_path)
    print(f"Cleaned table saved to {args.output_path}")


if __name__ == "__main__":
    main()