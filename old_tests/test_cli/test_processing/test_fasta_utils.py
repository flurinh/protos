"""
Tests for the FASTA processing utilities in protos.cli.processing.fasta_utils
"""

import pytest
from unittest.mock import patch, mock_open, MagicMock
import os
from protos.cli.processing.fasta_utils import preprocess_fasta, postprocess_fasta


@pytest.fixture
def mock_fasta_content():
    """Sample FASTA content for testing."""
    return {
        "protein1": "ACDEFGHIKLMNPQRSTVWY",
        "protein2": "ACDEFGHI-KLMNPQRSTVWY",  # Invalid character
        "protein3": "ACDEFGHIKLMNPQRSTVWY",
    }


@pytest.fixture
def mock_ref_fasta_content():
    """Sample reference FASTA content for testing."""
    return {
        "protein1": "ACDEFGHIKLMNPQRSTVWY",
        "protein3": "ACDEFGHIKLMNPQRSTVWY",
        # protein2 is not in reference
    }


@patch('protos.cli.processing.fasta_utils.read_fasta')
@patch('protos.cli.processing.fasta_utils.write_fasta')
@patch('protos.cli.processing.fasta_utils.validate_fasta_format')
@patch('protos.cli.processing.fasta_utils.clean_sequence')
@patch('os.makedirs')
def test_preprocess_fasta(mock_makedirs, mock_clean, mock_validate, mock_write, mock_read, mock_fasta_content):
    """Test preprocessing a FASTA file."""
    # Configure mocks
    mock_read.return_value = mock_fasta_content
    mock_validate.side_effect = lambda x: 'protein2' not in x  # protein2 is invalid
    mock_clean.side_effect = lambda x: x  # Identity function
    
    # Call the function
    result = preprocess_fasta('input.fasta', 'output.fasta')
    
    # Verify read_fasta was called
    mock_read.assert_called_once_with('input.fasta')
    
    # Verify write_fasta was called with correct arguments
    output_dict = mock_write.call_args[0][0]
    assert 'protein1' in output_dict
    assert 'protein2' not in output_dict  # Should be filtered out
    assert 'protein3' in output_dict
    
    # Verify output directory was created
    mock_makedirs.assert_called_once()
    
    # Verify the result
    assert len(result) == 2
    assert 'protein1' in result
    assert 'protein3' in result


@patch('protos.cli.processing.fasta_utils.read_fasta')
@patch('protos.cli.processing.fasta_utils.write_fasta')
@patch('os.makedirs')
def test_postprocess_fasta(mock_makedirs, mock_write, mock_read, mock_fasta_content, mock_ref_fasta_content):
    """Test postprocessing a FASTA file."""
    # Configure mocks
    mock_read.side_effect = [mock_fasta_content, mock_ref_fasta_content]
    
    # Call the function
    result = postprocess_fasta('input.fasta', 'ref.fasta', 'output.fasta')
    
    # Verify read_fasta was called twice
    assert mock_read.call_count == 2
    mock_read.assert_any_call('input.fasta')
    mock_read.assert_any_call('ref.fasta')
    
    # Verify write_fasta was called with correct arguments
    output_dict = mock_write.call_args[0][0]
    assert 'protein1' in output_dict
    assert 'protein2' not in output_dict  # Not in reference
    assert 'protein3' in output_dict
    
    # Verify output directory was created
    mock_makedirs.assert_called_once()
    
    # Verify the result
    assert len(result) == 2
    assert 'protein1' in result
    assert 'protein3' in result