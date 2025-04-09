"""
Tests for the updated GRN utilities functions.

Tests the enhanced GRN parsing to handle both 'x' and '.' notation.
"""

import pytest
import pandas as pd
import numpy as np

# Import the original utilities for comparison
from protos.processing.grn.grn_utils import parse_grn_str2float as original_parse
from protos.processing.grn.grn_utils import parse_grn_float2str, sort_grns_str

# Import updated utilities for testing
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../src/protos/processing/grn')))
from grn_utils_updated import parse_grn_str2float as updated_parse
from grn_utils_updated import check_str_grn_valid, sort_grns_str as updated_sort


class TestGRNUtilsUpdated:
    """Test the updated GRN utilities."""
    
    def test_parse_grn_str2float_original(self):
        """Test the original parser with 'x' notation."""
        # Test TM positions
        assert original_parse("1x50") == 1.5
        assert original_parse("3x50") == 3.5
        assert original_parse("7x50") == 7.5
        
        # Test N-terminal
        assert original_parse("n.10") == -10
        
        # Test C-terminal
        assert original_parse("c.15") == 115
    
    def test_parse_grn_str2float_updated(self):
        """Test the updated parser with both 'x' and '.' notation."""
        # Test 'x' notation (same as original)
        assert updated_parse("1x50") == 1.5
        assert updated_parse("3x50") == 3.5
        assert updated_parse("7x50") == 7.5
        
        # Test '.' notation (new)
        assert updated_parse("1.50") == 1.5
        assert updated_parse("3.50") == 3.5
        assert updated_parse("7.50") == 7.5
        
        # Test N-terminal and C-terminal (should be unchanged)
        assert updated_parse("n.10") == -10
        assert updated_parse("c.15") == 115
        
        # Test error cases
        assert updated_parse("invalid") == 0.0
        assert updated_parse("1.5.0") == 0.0  # Invalid dot format
    
    def test_check_str_grn_valid(self):
        """Test validation of GRN strings."""
        # Valid GRNs with 'x' notation
        assert check_str_grn_valid("1x50")
        assert check_str_grn_valid("3x50")
        assert check_str_grn_valid("7x50")
        
        # Valid GRNs with '.' notation
        assert check_str_grn_valid("1.50")
        assert check_str_grn_valid("3.50") 
        assert check_str_grn_valid("7.50")
        
        # Invalid GRNs
        assert not check_str_grn_valid("0x50")  # TM0 doesn't exist
        assert not check_str_grn_valid("9x50")  # TM9 doesn't exist
        assert not check_str_grn_valid("0.50")  # TM0 doesn't exist
        assert not check_str_grn_valid("9.50")  # TM9 doesn't exist
    
    def test_sort_grns_mixed_notation(self):
        """Test sorting GRNs with mixed notation."""
        # Create a list with mixed notation
        mixed_grns = ["1x50", "3.50", "2x50", "7.49", "n.10", "c.5"]
        
        # Sort using updated function
        sorted_grns = updated_sort(mixed_grns)
        
        # Expected order: n-terminal, then TM regions in order, then c-terminal
        assert sorted_grns[0] == "n.10"
        assert "1x50" in sorted_grns[1:4]
        assert "2x50" in sorted_grns[1:4]
        assert "3x50" in sorted_grns[1:4]  # Note: 3.50 gets converted to 3x50
        assert "7x49" in sorted_grns[4:5]  # Note: 7.49 gets converted to 7x49
        assert sorted_grns[-1] == "c.5"
    
    def test_equivalence_of_notations(self):
        """Test that dot and x notations are treated equivalently."""
        # Test pairs of equivalent positions
        pairs = [
            ("1.50", "1x50"),
            ("3.50", "3x50"),
            ("7.49", "7x49")
        ]
        
        for dot_notation, x_notation in pairs:
            # Both should parse to the same float value
            assert updated_parse(dot_notation) == updated_parse(x_notation)
            
            # Both should be valid
            assert check_str_grn_valid(dot_notation)
            assert check_str_grn_valid(x_notation)
            
            # When sorted, they should be considered equivalent
            mixed = [dot_notation, x_notation, "2x50"]
            sorted_mixed = updated_sort(mixed)
            # Should only have 2 items after sorting (as duplicates get combined)
            assert len(sorted_mixed) == 2