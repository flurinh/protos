"""
Tests for the sequence alignment functions in protos.processing.sequence.seq_alignment
"""

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from protos.processing.sequence.seq_alignment import (
    format_alignment,
    align_blosum62,
    init_aligner,
    find_idx
)

def test_init_aligner():
    """Test initialization of BLOSUM62 aligner."""
    # Initialize the aligner with default settings
    aligner = init_aligner()
    
    # Verify basic properties of the aligner
    assert hasattr(aligner, 'mode')
    assert aligner.mode == 'global'
    assert hasattr(aligner, 'match_score')
    assert hasattr(aligner, 'mismatch_score')
    assert hasattr(aligner, 'open_gap_score')
    assert hasattr(aligner, 'extend_gap_score')
    
    # Test with custom open gap score
    custom_aligner = init_aligner(open_gap_score=-15)
    assert custom_aligner.open_gap_score == -15

@patch('protos.processing.sequence.seq_alignment.Align.PairwiseAligner.align')
def test_align_blosum62(mock_align):
    """Test the BLOSUM62 alignment function."""
    # Create a mock alignment object
    mock_alignment = MagicMock()
    mock_alignment.score = 42
    mock_align.return_value = [mock_alignment]
    
    # Create a mock aligner
    mock_aligner = MagicMock()
    mock_aligner.align = mock_align
    
    # Call the function
    result = align_blosum62('ACDEFG', 'ACEFG', mock_aligner, verbose=0)
    
    # Check the result
    assert result == mock_alignment
    
    # Verify the mock was called with correct parameters
    mock_align.assert_called_once()
    args, kwargs = mock_align.call_args
    assert args[0] == 'ACDEFG'
    assert args[1] == 'ACEFG'

def test_find_idx():
    """Test finding index of first non-gap character."""
    # Basic test
    assert find_idx('---ABC') == 3
    
    # All gaps
    assert find_idx('------') == 0
    
    # No gaps
    assert find_idx('ABCDEF') == 0
    
    # Mixed case
    assert find_idx('---abc---') == 3

def test_format_alignment():
    """Test the alignment formatting function."""
    # Test basic alignment formatting
    seq1 = 'ACDEFG'
    seq2 = 'AC-EFG'
    score = 5
    
    formatted = format_alignment((seq1, seq2, score))
    assert len(formatted) == 3
    assert formatted[0] == seq1
    assert formatted[1] == seq2
    assert formatted[2] == score
    
    # Test with longer sequences
    seq1 = 'ACDEFGHIJKLMN'
    seq2 = 'AC-EFGH-JKLMN'
    score = 10
    
    formatted = format_alignment((seq1, seq2, score))
    assert len(formatted) == 3
    assert formatted[0] == seq1
    assert formatted[1] == seq2
    assert formatted[2] == score