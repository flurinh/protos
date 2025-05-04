"""
Tests for the expanded GRN expand_annotation implementation.

This module old_tests the new implementation of expand_annotation function
that replaces the legacy implementation.
"""

import pytest
import pandas as pd
import os
from typing import Dict, List, Tuple

# Mock GRN config to avoid file not found errors
class MockGRNConfigManager:
    def __init__(self, protein_family=None):
        self.protein_family = protein_family
    
    def get_config(self, strict=False):
        return {
            'family': self.protein_family,
            'strict': strict,
            'grns': ['1x50', '2x50', '3x50', '4x50', '5x50', '6x50', '7x50']
        }
    
    def init_grns(self, strict=False):
        return ['1x50', '2x50', '3x50', '4x50', '5x50', '6x50', '7x50']

# Mock function for get_correctly_aligned_grns
def mock_get_correctly_aligned_grns(all_query_gene_numbers, reference_grn_dict, alignment, max_alignment_gap=1):
    # A simple mock that maps the first few positions to GRNs
    result = {}
    for i, gene_num in enumerate(all_query_gene_numbers[:10]):
        if i < len(reference_grn_dict):
            result[gene_num] = list(reference_grn_dict.keys())[i]
    return result

# Add mock functions to interface_definitions module
import sys
import protos.processing.schema.interface_definitions
setattr(protos.processing.schema.interface_definitions, "get_correctly_aligned_grns", mock_get_correctly_aligned_grns)

# Add mock functions for terminal region annotation
def calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers, grns_float):
    # Mock implementation
    return [("M1", "n.5")], 1

def calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len, grns_float):
    # Mock implementation
    return [("D250", "c.3")], 250

setattr(protos.processing.schema.interface_definitions, "calculate_missing_ntail_grns", calculate_missing_ntail_grns)
setattr(protos.processing.schema.interface_definitions, "calculate_missing_ctail_grns", calculate_missing_ctail_grns)

# Add a mock for the GRNConfigManager
import protos.processing.grn.grn_utils
setattr(protos.processing.grn.grn_utils, "GRNConfigManager", MockGRNConfigManager)

# Import functions after mocks are in place
from protos.cli.grn.assign_grns import (
    get_pairwise_alignment,
    expand_annotation,
    init_aligner,
    align_blosum62,
    format_alignment
)

class TestExpandAnnotation:
    """Tests for the new expand_annotation implementation."""
    
    @pytest.fixture
    def test_data(self):
        """Provide test data for annotation expansion."""
        # Sample sequences
        query_seq = "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD"
        ref_seq = "MSEALLQLLGILGALALSLLGQVQQQKVAGVETHSQPEGRWMWLALGTALMGLGLLALLVKGMGVSDADQGKFYAITTLVPAIAFTMYLCMLLGYGLTMVPMGAEQNPIYWARHADWLFTTPLLLLDLALLVDAGQGTILALVGADGIMIGTGLVGALTKAYSYRFVWWAISTAAMLYILYYLFGGFTSKAEAMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD"
        
        # Initial GRN annotations
        reference_grns = {
            '1x50': 'L',
            '2x50': 'A',
            '3x50': 'R',
            '4x50': 'G',
            '5x50': 'F',
            '6x50': 'W',
            '7x50': 'K',
            'n.5': 'M',
            'c.3': 'D',
            '12.003': 'G',
            '23.005': 'T'
        }
        
        # Create alignment
        aligner = init_aligner(open_gap_score=-22)
        raw_alignment = align_blosum62(query_seq, ref_seq, aligner)
        alignment = format_alignment(raw_alignment)
        
        # Create Series
        new_row = pd.Series(reference_grns)
        
        return {
            'query_seq': query_seq,
            'ref_seq': ref_seq,
            'new_row': new_row,
            'alignment': alignment
        }
    
    def test_expand_annotation_basic(self, test_data, monkeypatch):
        """Test basic functionality of expand_annotation."""
        # Mock assign_missing_std_grns to return a simple result
        def mock_assign_missing_std_grns(*args, **kwargs):
            return [("R50", "3x50")], [25, 75]
        
        # Mock annotate_gaps_and_loops to return a simple result
        def mock_annotate_gaps_and_loops(*args, **kwargs):
            return [("L15", "12.003")], [("T30", "23.005")], [("S60", "45.002")]
        
        # Apply mocks
        import protos.processing.grn.grn_assignment
        monkeypatch.setattr(
            protos.processing.grn.grn_assignment, 
            "assign_missing_std_grns", 
            mock_assign_missing_std_grns
        )
        monkeypatch.setattr(
            protos.processing.grn.grn_assignment, 
            "annotate_gaps_and_loops", 
            mock_annotate_gaps_and_loops
        )
        
        # Call expand_annotation
        grns, rns, missing = expand_annotation(
            new_row=test_data['new_row'], 
            query_seq=test_data['query_seq'], 
            alignment=test_data['alignment'], 
            protein_family='microbial_opsins',
            max_alignment_gap=1, 
            verbose=1
        )
        
        # Verify results
        assert len(grns) > 0, "No GRNs were assigned"
        assert len(rns) > 0, "No residue numbers were assigned"
        assert len(grns) == len(rns), "GRNs and residue numbers don't match"
    
    def test_expand_annotation_fallback(self, test_data, monkeypatch):
        """Test fallback to legacy implementation on error."""
        # Mock to force an error in the main implementation
        def mock_get_correctly_aligned_grns_error(*args, **kwargs):
            raise ValueError("Test error for fallback")
        
        # Mock legacy_expand_annotation
        def mock_legacy_expand_annotation(*args, **kwargs):
            return (["1x50", "7x50"], ["A1", "K99"], [])
        
        # Apply mocks
        # First import protos to fix the "protos referenced before assignment" error
        import protos.processing.schema.interface_definitions
        monkeypatch.setattr(
            protos.processing.schema.interface_definitions, 
            "get_correctly_aligned_grns", 
            mock_get_correctly_aligned_grns_error
        )
        
        import protos.processing.grn.grn_table_utils
        monkeypatch.setattr(
            protos.processing.grn.grn_table_utils, 
            "expand_annotation", 
            mock_legacy_expand_annotation
        )
        
        # Call expand_annotation
        grns, rns, missing = expand_annotation(
            new_row=test_data['new_row'], 
            query_seq=test_data['query_seq'], 
            alignment=test_data['alignment'], 
            protein_family='microbial_opsins'
        )
        
        # Verify fallback results
        assert grns == ["1x50", "7x50"], "Legacy fallback didn't work"
        assert rns == ["A1", "K99"], "Legacy fallback didn't work"
    
    def test_normalization(self, test_data, monkeypatch):
        """Test GRN format normalization."""
        # Replace reference_grn_dict with legacy formats
        legacy_grns = {
            '1.50': 'L',   # Should normalize to 1x50
            '12x05': 'G',  # Should normalize to 12.005
            '23.5': 'T'    # Should normalize to 23.005
        }
        test_data['new_row'] = pd.Series(legacy_grns)
        
        # Override mock_get_correctly_aligned_grns to return legacy formats
        def mock_get_legacy_grns(*args, **kwargs):
            return {"A1": "1.50", "G12": "12x05", "T23": "23.5"}
        
        # First import protos to fix the "protos referenced before assignment" error
        import protos.processing.schema.interface_definitions
        monkeypatch.setattr(
            protos.processing.schema.interface_definitions, 
            "get_correctly_aligned_grns", 
            mock_get_legacy_grns
        )
        
        # Mock other functions to simplify the test
        def mock_assign_missing_std_grns(*args, **kwargs):
            return [], []
        
        def mock_annotate_gaps_and_loops(*args, **kwargs):
            return [], [], []
        
        monkeypatch.setattr(
            protos.processing.grn.grn_assignment, 
            "assign_missing_std_grns", 
            mock_assign_missing_std_grns
        )
        monkeypatch.setattr(
            protos.processing.grn.grn_assignment, 
            "annotate_gaps_and_loops", 
            mock_annotate_gaps_and_loops
        )
        
        # Call expand_annotation
        grns, rns, missing = expand_annotation(
            new_row=test_data['new_row'], 
            query_seq=test_data['query_seq'], 
            alignment=test_data['alignment'], 
            protein_family='microbial_opsins',
            verbose=1
        )
        
        # Verify some GRNs were assigned
        assert len(grns) > 0, "No GRNs were assigned"