"""
Tests for the GRN assignment utility functions with improved loop handling.

This module tests the GRN assignment utility functions defined in the 
protos.processing.schema.grn_assignment_utils module, focusing on
correct handling of loop region assignment.
"""

import pytest
from typing import Dict, List, Tuple

from protos.processing.schema.grn_assignment_utils import (
    assign_gene_nr,
    get_closest_present_seqnr,
    annotate_loop_region,
    calculate_missing_ntail_grns,
    calculate_missing_ctail_grns,
    valid_jump,
    get_correctly_aligned_grns
)

from protos.processing.schema.grn_utils_updated import parse_grn_str2float


class TestGRNAssignmentUtils:
    """Tests for GRN assignment utility functions."""
    
    def test_assign_gene_nr(self):
        """Test assigning gene numbers to a sequence."""
        seq = "MEAKL"
        expected = ["M1", "E2", "A3", "K4", "L5"]
        assert assign_gene_nr(seq) == expected
        
    def test_get_closest_present_seqnr(self):
        """Test finding the closest present sequence number."""
        # Create a list of present sequence number and GRN pairs
        present_pairs = [
            ("A10", "1x50"),
            ("L20", "2x50"),
            ("K30", "3x50"),
            ("M40", "4x50")
        ]
        
        # Test finding the closest pair on the N-terminal side
        missing_seqnr = 25
        closest, min_dist = get_closest_present_seqnr(missing_seqnr, present_pairs, loop_side='n')
        assert closest == ("L20", "2x50")
        assert min_dist == -5
        
        # Test finding the closest pair on the C-terminal side
        closest, min_dist = get_closest_present_seqnr(missing_seqnr, present_pairs, loop_side='c')
        assert closest == ("K30", "3x50")
        assert min_dist == 5
        
    def test_annotate_loop_region(self):
        """Test annotating a loop region between two helices."""
        # Create a list of present sequence number and GRN pairs
        present_pairs = [
            ("A10", "1x50"),
            ("L15", "1x55"),
            ("K30", "2x45"),
            ("M35", "2x50")
        ]
        
        # Define a loop region between the two helices
        loop_interval = [20, 21, 22]
        
        # Test annotating the loop
        query_seq = "AAAAAAAAAALLLLLAAAAAKKKKKMMMMM"
        loop_pairs = annotate_loop_region(loop_interval, present_pairs, query_seq)
        
        # Check that we got the expected number of pairs
        assert len(loop_pairs) == 3
        
        # Check that the format is correct: <closer helix><further helix>.<distance>
        for pair in loop_pairs:
            assert len(pair[1]) == 6  # Format: 12.003
            assert pair[1][2] == "."  # Must have decimal point
            assert len(pair[1].split(".")[1]) == 3  # Must have 3 digits after decimal
            
        # Check that the first and second characters are the helix numbers
        for pair in loop_pairs:
            assert pair[1][0] in "12"
            assert pair[1][1] in "12"
            
    def test_valid_jump(self):
        """Test the valid_jump function for alignment continuity."""
        # Valid jumps within the same helix
        assert valid_jump("1x50", "1x51", "A10", "A11", max_alignment_gap=1) is True
        assert valid_jump("1x50", "1x52", "A10", "A12", max_alignment_gap=2) is True
        
        # Invalid jumps within the same helix (gap too large)
        assert valid_jump("1x50", "1x52", "A10", "A12", max_alignment_gap=1) is False
        
        # Valid jumps between different helices
        assert valid_jump("1x50", "2x50", "A10", "L20", max_alignment_gap=10) is True
        
        # First alignment position (no previous)
        assert valid_jump(None, "1x50", None, "A10", max_alignment_gap=1) is True
        
    def test_get_correctly_aligned_grns(self):
        """Test getting correctly aligned GRNs from an alignment."""
        # Mock alignment
        query_seq = "MEAKLF-GKIP"
        match_line = "|||||.|.|||"
        ref_seq = "MEAKLFRGK-P"
        alignment = (query_seq, match_line, ref_seq)
        
        # Mock reference GRNs
        reference_grn_dict = {
            "M1": "1x50",
            "E2": "1x51",
            "A3": "1x52",
            "K4": "1x53",
            "L5": "1x54",
            "F6": "1x55",
            "R7": "12.005",  # Loop region
            "G8": "2x50",
            "K9": "2x51",
            "P10": "2x52"
        }
        
        # Mock query gene numbers
        all_query_gene_numbers = ["M1", "E2", "A3", "K4", "L5", "F6", "G7", "K8", "I9", "P10"]
        
        # Get correctly aligned GRNs
        result_dict = get_correctly_aligned_grns(all_query_gene_numbers, reference_grn_dict, alignment)
        
        # Check the results
        assert "M1" in result_dict and result_dict["M1"] == "1x50"
        assert "K4" in result_dict and result_dict["K4"] == "1x53"
        assert "G7" in result_dict and result_dict["G7"] == "2x50"
        assert "P10" in result_dict and result_dict["P10"] == "2x52"
        
        # The 'R7' in reference maps to loop GRN, but it's not aligned in query
        assert "F6" in result_dict
        assert result_dict["F6"] == "1x55"
        
    def test_calculate_missing_ntail_grns(self):
        """Test calculating missing N-terminal GRNs."""
        # Mock alignment
        aligned_grns = {
            "K10": "1x50",
            "L11": "1x51"
        }
        
        # Mock missing gene numbers
        missing_gene_numbers = ["M1", "E2", "A3", "T4", "G5", "V6", "L7", "I8", "P9"]
        
        # Mock GRN float values
        grns_float = [1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51]
        
        # Calculate missing N-terminal GRNs
        n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers, grns_float)
        
        # Check the results
        assert len(n_tail_list) > 0
        assert first_gene_number_int == 9
        
        # Check that the format is correct
        for pair in n_tail_list:
            # Either standard format or N-terminal format
            is_standard = "x" in pair[1]
            is_n_term = pair[1].startswith("n.")
            assert is_standard or is_n_term
            
    def test_calculate_missing_ctail_grns(self):
        """Test calculating missing C-terminal GRNs."""
        # Mock alignment
        aligned_grns = {
            "K10": "1x50",
            "L11": "1x51"
        }
        
        # Mock missing gene numbers
        missing_gene_numbers = ["D12", "Y13", "S14", "T15", "G16", "W17", "Q18", "F19", "R20"]
        
        # Mock query gene length
        query_gene_len = 20
        
        # Mock GRN float values
        grns_float = [1.5, 1.51, 1.52, 1.53, 1.54, 1.55]
        
        # Calculate missing C-terminal GRNs
        c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len, grns_float)
        
        # Check the results
        assert len(c_tail_list) > 0
        assert last_gene_number_int == 11
        
        # Check that the format is correct
        for pair in c_tail_list:
            # Either standard format or C-terminal format
            is_standard = "x" in pair[1]
            is_c_term = pair[1].startswith("c.")
            assert is_standard or is_c_term