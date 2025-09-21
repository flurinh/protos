import pandas as pd
import re
import argparse

def validate_and_clean_row(row):
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

def main():
    parser = argparse.ArgumentParser(description="Clean and validate residue sequences in a GRN table.")
    parser.add_argument('-i', '--input_path', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('-o', '--output_path', type=str, required=True, help='Path to save the cleaned CSV file')

    args = parser.parse_args()

    process_table(args.input_path, args.output_path)


if __name__ == "__main__":
    main()
    # python - m  scripts.clean_grn_table -i input/mosquito.csv -o output/mosquito2.csv

