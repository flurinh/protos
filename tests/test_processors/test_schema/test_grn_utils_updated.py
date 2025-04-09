"""
Tests for the updated GRN utility functions with corrected loop formatting.

This module tests the GRN utility functions defined in the 
protos.processing.schema.grn_utils_updated module, focusing on
correct handling of loop region formatting.
"""

import pytest
import re
from protos.processing.schema.grn_utils_updated import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string,
    sort_grns
)


class TestGRNUtilitiesUpdated:
    """Tests for updated GRN utility functions."""
    
    def test_parse_grn_str2float_standard(self):
        """Test parsing standard GRN format strings to floats."""
        assert parse_grn_str2float('1x50') == 1.5
        assert parse_grn_str2float('2x50') == 2.5
        assert parse_grn_str2float('7x01') == 7.01
        
    def test_parse_grn_str2float_n_term(self):
        """Test parsing N-terminal GRN format strings to floats."""
        assert parse_grn_str2float('n.10') == -0.1
        assert parse_grn_str2float('n.5') == -0.05
        assert parse_grn_str2float('n.100') == -1.0
        
    def test_parse_grn_str2float_c_term(self):
        """Test parsing C-terminal GRN format strings to floats."""
        assert parse_grn_str2float('c.5') == 100.05
        assert parse_grn_str2float('c.10') == 100.1
        assert parse_grn_str2float('c.100') == 101.0
        
    def test_parse_grn_str2float_loop(self):
        """Test parsing loop GRN format strings to floats."""
        # Loop between helix 1 and 2, closer to 1, distance 3
        assert parse_grn_str2float('12.003') == 12.003
        # Loop between helix 5 and 6, closer to 6, distance 11
        assert parse_grn_str2float('65.011') == 56.011
        # Loop between helix 3 and 7, closer to 3, distance 0
        assert parse_grn_str2float('37.000') == 37.0
        
    def test_parse_grn_str2float_invalid(self):
        """Test parsing invalid GRN format strings."""
        assert parse_grn_str2float('invalid') == 0.0
        assert parse_grn_str2float('') == 0.0
        assert parse_grn_str2float('1.2.3') == 0.0
        
    def test_parse_grn_float2str_standard(self):
        """Test converting standard GRN float values to strings."""
        assert parse_grn_float2str(1.5) == '1x50'
        assert parse_grn_float2str(2.5) == '2x50'
        assert parse_grn_float2str(7.01) == '7x01'
        
    def test_parse_grn_float2str_n_term(self):
        """Test converting N-terminal GRN float values to strings."""
        assert parse_grn_float2str(-0.1) == 'n.10'
        assert parse_grn_float2str(-0.05) == 'n.5'
        assert parse_grn_float2str(-1.0) == 'n.100'
        
    def test_parse_grn_float2str_c_term(self):
        """Test converting C-terminal GRN float values to strings."""
        assert parse_grn_float2str(100.05) == 'c.5'
        assert parse_grn_float2str(100.1) == 'c.10'
        assert parse_grn_float2str(101.0) == 'c.100'
        
    def test_parse_grn_float2str_loop(self):
        """Test converting loop GRN float values to strings."""
        # Loop between helix 1 and 2, closer to 1, distance 3
        assert parse_grn_float2str(12.003) == '12.003'
        # Loop between helix 5 and 6, closer to 6, distance 11
        assert parse_grn_float2str(56.011) == '56.011'
        # Loop between helix 3 and 7, closer to 3, distance 0
        assert parse_grn_float2str(37.0) == '37.000'
        
    def test_normalize_grn_format(self):
        """Test normalizing GRN strings to standardized format."""
        # Already in standard format
        assert normalize_grn_format('1x50') == '1x50'
        assert normalize_grn_format('n.10') == 'n.10'
        assert normalize_grn_format('c.5') == 'c.5'
        assert normalize_grn_format('12.003') == '12.003'
        
        # Legacy formats
        assert normalize_grn_format('12x05') == '12.005'
        assert normalize_grn_format('12.5') == '12.005'
        assert normalize_grn_format('1.50') == '1x50'
        
        # Invalid formats
        assert normalize_grn_format('invalid') == 'invalid'
        
    def test_validate_grn_string(self):
        """Test validating GRN strings."""
        # Valid GRNs
        assert validate_grn_string('1x50')[0] is True
        assert validate_grn_string('n.10')[0] is True
        assert validate_grn_string('c.5')[0] is True
        assert validate_grn_string('12.003')[0] is True
        
        # Invalid GRNs
        assert validate_grn_string('')[0] is False
        assert validate_grn_string('invalid')[0] is False
        assert validate_grn_string('9x50')[0] is False  # Helix out of range
        assert validate_grn_string('1x100')[0] is False  # Position out of range
        assert validate_grn_string('12.1000')[0] is False  # Distance too large
        
        # GRNs that can be normalized and then validated
        is_valid, message = validate_grn_string('12x05')
        assert is_valid is True
        assert "after normalization" in message
        
    def test_sort_grns(self):
        """Test sorting GRNs in standard order."""
        # Test with string GRNs
        grn_strings = ['3x50', '1x50', 'n.10', 'c.5', '12.003']
        sorted_strings = sort_grns(grn_strings)
        assert sorted_strings[0] == 'n.10'  # N-terminal first
        assert sorted_strings[-1] == 'c.5'  # C-terminal last
        assert sorted_strings[1] == '1x50'  # TM helices in order
        assert sorted_strings[2] == '3x50'
        assert sorted_strings[3] == '12.003'  # Loop between TM and C-terminal
        
        # Test with float GRNs
        grn_floats = [3.5, 1.5, -0.1, 100.05, 12.003]
        sorted_floats = sort_grns(grn_floats)
        assert sorted_floats[0] == -0.1  # N-terminal first
        assert sorted_floats[-1] == 100.05  # C-terminal last
        assert sorted_floats[1] == 1.5  # TM helices in order
        assert sorted_floats[2] == 3.5
        assert sorted_floats[3] == 12.003  # Loop between TM and C-terminal
        
    def test_roundtrip_conversion(self):
        """Test roundtrip conversion between string and float GRNs."""
        grn_strings = ['1x50', 'n.10', 'c.5', '12.003', '65.011']
        
        for grn in grn_strings:
            # Convert to float
            grn_float = parse_grn_str2float(grn)
            # Convert back to string
            grn_str = parse_grn_float2str(grn_float)
            # Normalize both for comparison
            assert normalize_grn_format(grn) == normalize_grn_format(grn_str)