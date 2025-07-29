#!/usr/bin/env python
"""
Fix NaN values in the microbial opsin reference table.
This script replaces NaN values with '-' to ensure compatibility with get_seq function.
"""

import pandas as pd
from pathlib import Path

# Path to the reference table
ref_file = Path(__file__).parent / "data" / "grn" / "ref" / "mo_ref.csv"

print(f"Loading reference table from: {ref_file}")

# Load the CSV
df = pd.read_csv(ref_file, index_col=0)

# Count NaN values before
nan_count_before = df.isna().sum().sum()
print(f"Found {nan_count_before} NaN values in the table")

# Replace NaN values with '-'
df = df.fillna('-')

# Count '-' values after
dash_count_after = (df == '-').sum().sum()
print(f"Total '-' values after replacement: {dash_count_after}")

# Save the cleaned table back
print(f"Saving cleaned table back to: {ref_file}")
df.to_csv(ref_file)

print("Done! The mo_ref.csv file has been updated with NaN values replaced by '-'")

# Show a sample of the data
print("\nSample of the cleaned data (first 5 rows, first 10 columns):")
print(df.iloc[:5, :10])