"""
Tests for GRN CLI utilities.

These tests focus on pure functions that don't need file access.
"""

import pytest
import pandas as pd

from protos.cli.grn.clean_grn_table import validate_and_clean_row
from protos.processing.grn.grn_assignment import (
    assign_gene_nr, 
    calculate_missing_gene_numbers,
    calculate_missing_ntail_grns,
    calculate_missing_ctail_grns
)


class TestGRNUtils:
    """Test GRN utility functions with synthetic sequences."""
    
    def test_assign_gene_nr(self):
        """Test gene number assignment on a sequence fragment."""
        # Use a synthetic sequence fragment
        seq = "TAAVGADLLGDGRPETLWL"
        gene_numbers = assign_gene_nr(seq)
        
        # Check correct assignment
        assert len(gene_numbers) == len(seq)
        assert gene_numbers[0] == 'T1'
        assert gene_numbers[1] == 'A2'
        assert gene_numbers[18] == 'L19'
        
    def test_calculate_missing_gene_numbers(self):
        """Test calculation of missing gene numbers."""
        all_gene_numbers = ['M1', 'L2', 'E3', 'L4', 'L5', 'P6']
        aligned_grns = {'L2': '1.50', 'L4': '2.50'}
        
        missing = calculate_missing_gene_numbers(all_gene_numbers, aligned_grns)
        
        assert set(missing) == {'M1', 'E3', 'L5', 'P6'}
        
    def test_calculate_ntail_grns(self):
        """Test N-terminal tail GRN calculation."""
        aligned_grns = {'L45': '1.50', 'A100': '2.50', 'R150': '3.50'}
        missing_gene_numbers = ['M1', 'L2', 'E3'] + [f'X{i}' for i in range(4, 45)]
        grns_float = [1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50]
        
        n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(
            aligned_grns, missing_gene_numbers, grns_float
        )
        
        # Should have some N-terminal assignments
        assert len(n_tail_list) > 0
        assert first_gene_number_int == 44  # Position before L45
        
    def test_calculate_ctail_grns(self):
        """Test C-terminal tail GRN calculation."""
        aligned_grns = {'L45': '1.50', 'A100': '2.50', 'K230': '7.50'}
        missing_gene_numbers = [f'X{i}' for i in range(231, 263)]  # Typical opsin length ~262
        grns_float = [7.50, 7.51, 7.52, 7.53, 7.54, 7.55, 7.56]
        query_gene_len = 262
        
        c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(
            aligned_grns, missing_gene_numbers, query_gene_len, grns_float
        )
        
        # Should have some C-terminal assignments
        assert len(c_tail_list) > 0
        assert last_gene_number_int == 230  # Position of K230
        
    def test_validate_and_clean_row_patterns(self):
        """Test validation on typical GRN table patterns."""
        # Normal sequence
        row = ['P42', 'I43', 'Y44', 'E45', 'T46']
        clean_row, is_erroneous = validate_and_clean_row(row)
        
        assert clean_row == row
        assert not is_erroneous
        
        # Sequence with gap (non-sequential numbering)
        row = ['P42', 'I43', '-', '-', 'V85', 'I86']
        clean_row, is_erroneous = validate_and_clean_row(row)
        
        assert clean_row == row
        assert is_erroneous  # Gap causes error due to non-sequential
        
        # Sequence with restart
        row = ['A1', 'R2', 'N3', 'A1', 'C2']
        clean_row, is_erroneous = validate_and_clean_row(row)
        
        assert clean_row == ['A1', 'R2', 'N3', '-', '-']
        assert not is_erroneous  # Restart is handled separately