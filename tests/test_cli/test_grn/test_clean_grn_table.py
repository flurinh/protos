"""
Tests for the GRN table cleaning utility in protos.cli.grn.clean_grn_table
"""

import pytest
import pandas as pd
from unittest.mock import patch, mock_open, MagicMock
import os
from io import StringIO
from protos.cli.grn.clean_grn_table import validate_and_clean_row, process_table, clean_grn_table


def test_validate_and_clean_row_normal_sequence():
    """Test validation and cleaning of a normal sequence."""
    # A normal sequence with sequential residue numbers
    row = ['A1', 'R2', 'N3', 'D4', 'C5']
    clean_row, is_erroneous = validate_and_clean_row(row)
    
    assert clean_row == row  # Row should not be modified
    assert not is_erroneous  # No errors should be detected


def test_validate_and_clean_row_with_gaps():
    """Test validation and cleaning of a sequence with gaps."""
    # A sequence with some gaps (represented by '-')
    row = ['A1', '-', 'N3', '-', 'C5']
    clean_row, is_erroneous = validate_and_clean_row(row)
    
    assert clean_row == row  # Gaps should be preserved
    assert is_erroneous  # Errors should be detected due to non-sequential numbers


def test_validate_and_clean_row_with_restart():
    """Test validation and cleaning of a sequence with a restart."""
    # A sequence that restarts numbering
    row = ['A1', 'R2', 'N3', 'A1', 'C2']
    clean_row, is_erroneous = validate_and_clean_row(row)
    
    # After detecting restart, all subsequent positions should be gaps
    assert clean_row == ['A1', 'R2', 'N3', '-', '-']
    assert not is_erroneous  # Restart is handled separately from errors


def test_validate_and_clean_row_erroneous_sequence():
    """Test validation and cleaning of an erroneous sequence."""
    # A sequence with non-sequential residue numbers
    row = ['A1', 'R2', 'N5', 'D6', 'C8']
    clean_row, is_erroneous = validate_and_clean_row(row)
    
    assert clean_row == row  # Row is preserved but marked as erroneous
    assert is_erroneous  # Should be detected as erroneous


@patch('pandas.read_csv')
@patch('pandas.DataFrame.to_csv')
def test_process_table(mock_to_csv, mock_read_csv):
    """Test processing a GRN table."""
    # Create a mock DataFrame with some test data
    df = pd.DataFrame({
        'col1': ['A1', 'A1', 'A1'],
        'col2': ['R2', 'R2', 'R3'],  # Third row is erroneous
        'col3': ['N3', 'A1', 'N4']   # Second row has restart
    }, index=['seq1', 'seq2', 'seq3'])
    
    mock_read_csv.return_value = df
    
    # Process the table
    erroneous_report = process_table('input.csv', 'output.csv')
    
    # Verify the DataFrame was read and saved
    mock_read_csv.assert_called_once_with('input.csv', index_col=0)
    mock_to_csv.assert_called_once()
    
    # Check that erroneous sequences were correctly identified
    assert 'seq3' in erroneous_report


@patch('protos.cli.grn.clean_grn_table.process_table')
def test_clean_grn_table(mock_process_table):
    """Test the main clean_grn_table function."""
    # Set up the mock
    mock_process_table.return_value = {'seq3': ['A1', 'R3', 'N4']}
    
    # Call the function
    result = clean_grn_table('input.csv', 'output.csv')
    
    # Verify process_table was called with correct args
    mock_process_table.assert_called_once_with('input.csv', 'output.csv')
    
    # Check that the result matches what process_table returned
    assert result == {'seq3': ['A1', 'R3', 'N4']}