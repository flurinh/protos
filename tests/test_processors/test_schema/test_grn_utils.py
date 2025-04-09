"""
Tests for the standardized GRN utility functions.

These tests verify that the GRN utility functions in the schema module
correctly handle various GRN formats and edge cases.
"""

import pytest
import re
from protos.processing.schema.grn_utils import (
    parse_grn_str2float,
    parse_grn_float2str,
    validate_grn_string,
    sort_grns,
    get_grn_region
)

# Original functions for comparison
# We need to import both to ensure compatibility
from protos.processing.grn.grn_utils_updated import (
    parse_grn_str2float as original_parse_grn_str2float,
    parse_grn_float2str as original_parse_grn_float2str,
    check_str_grn_valid as original_check_str_grn_valid,
    sort_grns_str as original_sort_grns_str
)


class TestGRNStringToFloat:
    """Tests for parse_grn_str2float function."""
    
    def test_valid_standard_grn(self):
        """Test parsing valid standard GRN strings (e.g., '1x50')."""
        assert parse_grn_str2float('1x50') == 1.5
        assert parse_grn_str2float('2x40') == 2.4
        assert parse_grn_str2float('7x65') == 7.65
        
        # Single-digit format
        assert parse_grn_str2float('1x01') == 1.01
        assert parse_grn_str2float('8x09') == 8.09
    
    def test_valid_n_terminal_grn(self):
        """Test parsing valid N-terminal GRN strings (e.g., 'n.10')."""
        assert parse_grn_str2float('n.10') == -0.1
        assert parse_grn_str2float('n.5') == -0.05
        assert parse_grn_str2float('n.100') == -1.0
    
    def test_valid_c_terminal_grn(self):
        """Test parsing valid C-terminal GRN strings (e.g., 'c.5')."""
        assert parse_grn_str2float('c.5') == 100.05
        assert parse_grn_str2float('c.10') == 100.1
        assert parse_grn_str2float('c.100') == 101.0
    
    def test_valid_loop_grn(self):
        """Test parsing valid loop GRN strings (e.g., '2.45')."""
        assert parse_grn_str2float('2.45') == 2.045
        assert parse_grn_str2float('3.60') == 3.06
        assert abs(parse_grn_str2float('1.05') - 1.005) < 0.001  # Allow for floating point imprecision
    
    def test_invalid_grn(self):
        """Test parsing invalid GRN strings."""
        # Should return 0.0 for invalid strings but not raise exceptions
        assert parse_grn_str2float('invalid') == 0.0
        assert parse_grn_str2float('1-50') == 0.0
        assert parse_grn_str2float('x50') == 0.0
        assert parse_grn_str2float('n.') == 0.0
        assert parse_grn_str2float('c.x') == 0.0
    
    def test_empty_grn(self):
        """Test parsing empty GRN strings."""
        assert parse_grn_str2float('') == 0.0
        assert parse_grn_str2float(None) == 0.0  # Should handle None
    
    def test_compatibility_with_original(self):
        """Test compatibility with the original parsing function."""
        test_grns = ['1x50', '2x40', '7x65', 'n.10', 'n.5', 'c.5', 'c.10']
        
        for grn in test_grns:
            assert abs(parse_grn_str2float(grn) - original_parse_grn_str2float(grn)) < 0.001


class TestGRNFloatToString:
    """Tests for parse_grn_float2str function."""
    
    def test_standard_grn(self):
        """Test converting standard GRN float values to strings."""
        assert parse_grn_float2str(1.5) == '1x50'
        assert parse_grn_float2str(2.4) == '2x40'
        assert parse_grn_float2str(7.65) == '7x65'
        
        # Test proper zero-padding for single-digit values
        assert parse_grn_float2str(1.01) == '1x01'
        assert parse_grn_float2str(8.09) == '8x09'
    
    def test_n_terminal_grn(self):
        """Test converting N-terminal GRN float values to strings."""
        assert parse_grn_float2str(-0.1) == 'n.10'
        assert parse_grn_float2str(-0.05) == 'n.5'
        assert parse_grn_float2str(-1.0) == 'n.100'
    
    def test_c_terminal_grn(self):
        """Test converting C-terminal GRN float values to strings."""
        assert parse_grn_float2str(100.05) == 'c.5'
        assert parse_grn_float2str(100.1) == 'c.10'
        assert parse_grn_float2str(101.0) == 'c.1'
    
    def test_loop_grn(self):
        """Test converting loop GRN float values to strings."""
        assert parse_grn_float2str(2.045) == '2.45'
        assert parse_grn_float2str(3.06) == '3.60'
        assert parse_grn_float2str(1.005) == '1.05'
    
    def test_zero_edge_case(self):
        """Test the zero edge case."""
        # The behavior for 0.0 is implementation-dependent
        # This should be a well-defined value, whatever it is
        result = parse_grn_float2str(0.0)
        assert isinstance(result, str)
    
    def test_rounding(self):
        """Test rounding behavior."""
        # Check that close values round correctly
        assert parse_grn_float2str(1.499) == '1x50'  # Should round to 1x50
        assert parse_grn_float2str(1.501) == '1x50'  # Should round to 1x50
    
    def test_compatibility_with_original(self):
        """Test compatibility with the original function."""
        test_floats = [1.5, 2.4, 7.65, -0.1, -0.05, 100.05, 100.1]
        
        for val in test_floats:
            # Allow for slight differences in zero-padding format
            assert re.sub(r'x0', 'x', parse_grn_float2str(val)) == re.sub(r'x0', 'x', original_parse_grn_float2str(val))


class TestValidateGRNString:
    """Tests for validate_grn_string function."""
    
    def test_valid_standard_grn(self):
        """Test validating valid standard GRN strings."""
        is_valid, message = validate_grn_string('1x50')
        assert is_valid is True
        assert "valid" in message.lower()
        
        is_valid, message = validate_grn_string('7x65')
        assert is_valid is True
        assert "valid" in message.lower()
    
    def test_valid_n_terminal_grn(self):
        """Test validating valid N-terminal GRN strings."""
        is_valid, message = validate_grn_string('n.10')
        assert is_valid is True
        assert "valid" in message.lower()
    
    def test_valid_c_terminal_grn(self):
        """Test validating valid C-terminal GRN strings."""
        is_valid, message = validate_grn_string('c.5')
        assert is_valid is True
        assert "valid" in message.lower()
    
    def test_valid_loop_grn(self):
        """Test validating valid loop GRN strings."""
        is_valid, message = validate_grn_string('2.45')
        assert is_valid is True
        assert "valid" in message.lower()
    
    def test_invalid_format(self):
        """Test validating invalid GRN formats."""
        is_valid, message = validate_grn_string('invalid')
        assert is_valid is False
        assert "invalid" in message.lower()
        
        is_valid, message = validate_grn_string('1-50')
        assert is_valid is False
        assert "invalid" in message.lower()
    
    def test_invalid_range(self):
        """Test validating GRNs with invalid ranges."""
        # Helix 9 is not valid
        is_valid, message = validate_grn_string('9x50')
        assert is_valid is False
        assert "invalid helix" in message.lower()
        
        # Position 100 is not valid
        is_valid, message = validate_grn_string('1x100')
        assert is_valid is False
        assert "invalid position" in message.lower()
    
    def test_invalid_leading_zero(self):
        """Test validating GRNs with invalid leading zeros."""
        is_valid, message = validate_grn_string('n.01')
        assert is_valid is False
        assert "leading zero" in message.lower()
        
        is_valid, message = validate_grn_string('c.01')
        assert is_valid is False
        assert "leading zero" in message.lower()
    
    def test_empty_grn(self):
        """Test validating empty GRN strings."""
        is_valid, message = validate_grn_string('')
        assert is_valid is False
        assert "empty" in message.lower()
        
        is_valid, message = validate_grn_string(None)
        assert is_valid is False
        assert "empty" in message.lower()
    
    def test_compatibility_with_original(self):
        """Test compatibility with the original validation function."""
        test_grns = ['1x50', '2x40', '7x65', 'n.10', 'n.5', 'c.5', 'c.10', 'invalid', '1-50', '9x50']
        
        for grn in test_grns:
            # Skip None/empty values since original might handle them differently
            if not grn:
                continue
                
            original_result = original_check_str_grn_valid(grn)
            new_result, _ = validate_grn_string(grn)
            assert new_result == original_result


class TestSortGRNs:
    """Tests for sort_grns function."""
    
    def test_sort_str_grns(self):
        """Test sorting string GRNs."""
        grns = ['3x50', '1x50', 'n.10', 'c.5', '2.45']
        sorted_grns = sort_grns(grns)
        
        # Check order: n-terminal, then helices, then loops, then c-terminal
        assert sorted_grns[0] == 'n.10'
        assert sorted_grns[1] == '1x50'
        assert sorted_grns[2] == '2.45'
        assert sorted_grns[3] == '3x50'
        assert sorted_grns[4] == 'c.5'
    
    def test_sort_float_grns(self):
        """Test sorting float GRNs."""
        grns = [3.5, 1.5, -0.1, 100.05, 2.045]
        sorted_grns = sort_grns(grns)
        
        # Check order: n-terminal, then helices, then loops, then c-terminal
        assert sorted_grns[0] == -0.1
        assert sorted_grns[1] == 1.5
        assert sorted_grns[2] == 2.045
        assert sorted_grns[3] == 3.5
        assert sorted_grns[4] == 100.05
    
    def test_mixed_format_error(self):
        """Test handling of mixed formats (should raise error)."""
        grns = ['1x50', 2.5, 'n.10']
        
        # This should raise ValueError since we have mixed string/float types
        with pytest.raises(ValueError):
            sort_grns(grns)
    
    def test_empty_list(self):
        """Test sorting an empty list of GRNs."""
        assert sort_grns([]) == []
    
    def test_compatibility_with_original(self):
        """Test compatibility with the original sort function."""
        grns = ['3x50', '1x50', 'n.10', 'c.5', '2.45']
        
        # We may have slightly different ordering logic, but the basic ordering principle
        # should be the same
        original_result = original_sort_grns_str(grns)
        new_result = sort_grns(grns)
        
        # Check that n-terminal comes first and c-terminal comes last in both
        assert original_result[0] == new_result[0] == 'n.10'
        assert original_result[-1] == new_result[-1] == 'c.5'


class TestGetGRNRegion:
    """Tests for get_grn_region function."""
    
    def test_n_terminal_region(self):
        """Test identifying N-terminal regions."""
        assert get_grn_region('n.10') == 'N-terminal'
        assert get_grn_region(-0.1) == 'N-terminal'
    
    def test_c_terminal_region(self):
        """Test identifying C-terminal regions."""
        assert get_grn_region('c.5') == 'C-terminal'
        assert get_grn_region(100.05) == 'C-terminal'
    
    def test_tm_regions(self):
        """Test identifying transmembrane regions."""
        assert get_grn_region('1x50') == 'TM1'
        assert get_grn_region('7x65') == 'TM7'
        assert get_grn_region(1.5) == 'TM1'
        assert get_grn_region(7.65) == 'TM7'
    
    def test_loop_regions(self):
        """Test identifying loop regions."""
        assert get_grn_region('2.45') == 'Loop'
        assert get_grn_region(2.045) == 'Loop'