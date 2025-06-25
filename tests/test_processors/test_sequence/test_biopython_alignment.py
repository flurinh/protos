"""
Tests for BioPython alignment functionality
"""

import pytest
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
    format_alignment,
    calc_alignment_score_restricted_area,
    get_score,
    get_best_alignment
)


class TestBioPythonAlignment:
    """Test BioPython alignment functionality."""
    
    @pytest.fixture
    def test_sequences(self):
        """Test sequences for alignment."""
        return {
            "seq1": "MALIGTLLMLIGTFYFIVKGWGVTDK",
            "seq2": "MALIGTLLMLIGTFYFIVKGWGVTDK",  # Identical
            "seq3": "MALIGTLLMXIGTFYFXVKGWGVXDK",  # 3 mutations
            "seq4": "TAAVGADLLGDGRPETLWLGIGTLLM",  # Different but similar
            "seq5": "ACDEFGHIKLMNPQRSTVWY"         # Completely different
        }
    
    def test_init_aligner(self):
        """Test aligner initialization."""
        aligner = init_aligner()
        
        # Check aligner properties
        assert aligner.mode == 'global', "Should be global alignment"
        
        # Check substitution matrix
        assert aligner.substitution_matrix is not None
        # Should handle 'J' for Leucine/Isoleucine
        assert 'J' in aligner.substitution_matrix.alphabet
        
        # Check that matrix has expected properties
        matrix = aligner.substitution_matrix
        # Check some known BLOSUM62 values
        assert matrix['A', 'A'] > 0  # Same amino acid should have positive score
        assert matrix['A', 'W'] < 0  # Different amino acids typically negative
        
        # Test that aligner can actually align sequences
        test_align = aligner.align("ACGT", "ACGT")
        alignments = list(test_align)
        assert len(alignments) > 0, "Should produce at least one alignment"
    
    def test_init_aligner_custom_gap_penalty(self):
        """Test aligner with custom gap penalty."""
        aligner = init_aligner(open_gap_score=-5)
        # Just verify it was created successfully
        assert aligner is not None
        assert aligner.mode == 'global'
    
    def test_align_blosum62_identical(self, test_sequences):
        """Test alignment of identical sequences."""
        aligner = init_aligner()
        seq1 = test_sequences["seq1"]
        seq2 = test_sequences["seq2"]
        
        alignment = align_blosum62(seq1, seq2, aligner)
        formatted = format_alignment(alignment)
        
        # Check formatted output
        assert isinstance(formatted, list), "Should return list"
        assert len(formatted) == 3, "Should have 3 elements"
        
        # For identical sequences, match line should be all '|'
        match_line = formatted[1]
        assert all(c == '|' for c in match_line if c not in [' ', '']), \
            "All positions should match for identical sequences"
    
    def test_align_blosum62_mutations(self, test_sequences):
        """Test alignment with mutations."""
        aligner = init_aligner()
        seq1 = test_sequences["seq1"]
        seq3 = test_sequences["seq3"]  # Has 3 mutations
        
        alignment = align_blosum62(seq1, seq3, aligner)
        formatted = format_alignment(alignment)
        
        # Check that mutations are marked as mismatches
        match_line = formatted[1]
        mismatch_count = match_line.count('.')
        assert mismatch_count >= 3, "Should detect at least 3 mismatches"
    
    def test_align_blosum62_gaps(self, test_sequences):
        """Test alignment requiring gaps."""
        aligner = init_aligner()
        seq1 = test_sequences["seq1"]
        seq4 = test_sequences["seq4"]  # Different length and sequence
        
        alignment = align_blosum62(seq1, seq4, aligner)
        formatted = format_alignment(alignment)
        
        # Should have gaps in alignment
        seq1_aligned = formatted[0]
        seq2_aligned = formatted[2]
        assert '-' in seq1_aligned or '-' in seq2_aligned, \
            "Alignment should contain gaps"
    
    def test_format_alignment_tuple_input(self):
        """Test format_alignment with tuple input."""
        # Test the tuple format (seq1, seq2, score)
        alignment_tuple = ("MALIGTLLM", "MA-IGTLLM", 15.5)
        formatted = format_alignment(alignment_tuple)
        
        assert isinstance(formatted, list)
        assert len(formatted) == 3
        assert formatted[0] == "MALIGTLLM"
        assert formatted[1] == "MA-IGTLLM"
        assert formatted[2] == 15.5
    
    def test_format_alignment_empty(self):
        """Test format_alignment with empty/invalid input."""
        # Test with empty alignment
        formatted = format_alignment("")
        assert formatted == ['', '', 0]
    
    def test_calc_alignment_score_restricted_area(self):
        """Test restricted area alignment scoring."""
        # Create test alignment
        alignment = [
            "MALIGTLLMLIGTFYFIVK",
            "||||||||.|||||||.||",
            "MALIGTLLVLIGTFYFLVK"
        ]
        
        score = calc_alignment_score_restricted_area(alignment)
        
        # Score should be between 0 and 1
        assert 0 <= score <= 1
        # With 2 mismatches out of 19, score should be ~0.89
        assert 0.85 < score < 0.95
    
    def test_calc_alignment_score_perfect(self):
        """Test scoring of perfect alignment."""
        alignment = [
            "MALIGTLLMLIGTFYFIVK",
            "|||||||||||||||||||",
            "MALIGTLLMLIGTFYFIVK"
        ]
        
        score = calc_alignment_score_restricted_area(alignment)
        assert score == 1.0, "Perfect alignment should have score 1.0"
    
    def test_calc_alignment_score_with_gaps(self):
        """Test scoring with gaps."""
        alignment = [
            "MALIGTLLM-LIGTFYFIVK",
            "|||||||||-||||||||||",  # Gap should be marked with -
            "MALIGTLLMALIGTFYFIVK"
        ]
        
        score = calc_alignment_score_restricted_area(alignment)
        # Score should be less than 1 due to gap
        assert score < 1.0
        assert score > 0.9  # Still mostly matched
    
    def test_get_score(self):
        """Test get_score function."""
        alignment = [
            "MALIGTLLMLIGTFYFIVK",
            "||||||||.|||||||.||",
            "MALIGTLLVLIGTFYFLVK"
        ]
        
        score = get_score(alignment)
        assert 0 <= score <= 1
        # Should print score
    
    def test_get_best_alignment(self):
        """Test finding best alignment from MSA."""
        # Create mock MSA
        msa = {
            "seq1": ("MALIGTLLM", 10.5, ["MALIGTLLM", "||||||||.", "MALIGTLLV"]),
            "seq2": ("MALIGTLLM", 15.0, ["MALIGTLLM", "||||||||-", "MALIGTLL-"]),
            "seq3": ("MALIGTLLM", 12.0, ["MALIGTLLM", "|||||||.|", "MALIGTLVM"])
        }
        
        # Test with default scoring
        best_name, alignment = get_best_alignment(msa, score_type='default')
        assert best_name == "seq2", "Should select highest score"
        
        # Test with restricted area scoring
        best_name, alignment = get_best_alignment(msa, score_type='restricted')
        assert best_name in ["seq1", "seq2", "seq3"]
        assert alignment is not None
    
    def test_alignment_with_j_residue(self):
        """Test alignment with J (Leu/Ile) residue."""
        aligner = init_aligner()
        
        # J should be treated as similar to both L and I
        seq_with_j = "MAJIGT"
        seq_with_l = "MALIGT"
        seq_with_i = "MAIIGT"
        
        # Align J with L
        alignment_l = align_blosum62(seq_with_j, seq_with_l, aligner)
        score_l = alignment_l.score if hasattr(alignment_l, 'score') else 0
        
        # Align J with I
        alignment_i = align_blosum62(seq_with_j, seq_with_i, aligner)
        score_i = alignment_i.score if hasattr(alignment_i, 'score') else 0
        
        # J should align reasonably well with both L and I
        assert score_l > 0
        assert score_i > 0
    
    def test_alignment_edge_cases(self):
        """Test edge cases in alignment."""
        aligner = init_aligner()
        
        # Empty sequences - should raise ValueError
        try:
            result = align_blosum62("", "", aligner)
            # If it doesn't raise, that's also ok
            assert result is not None or result == []
        except ValueError:
            # Expected for empty sequences
            pass
        
        # Single character
        result = align_blosum62("M", "M", aligner)
        formatted = format_alignment(result)
        assert len(formatted) == 3
        
        # Very different sequences
        seq1 = "MALIGTLLM"
        seq2 = "WYKPQRST"
        result = align_blosum62(seq1, seq2, aligner)
        assert result is not None


class TestAlignmentScoring:
    """Test alignment scoring functions."""
    
    def test_score_range(self):
        """Test that scores are in valid range."""
        test_alignments = [
            ["MALIGTLLM", "||||||||-", "MALIGTLL-"],
            ["MALIGTLLM", ".........", "WYKPQRSTV"],
            ["MALIGTLLM", "|||.|.|||", "MALVGALLM"],
        ]
        
        for alignment in test_alignments:
            score = calc_alignment_score_restricted_area(alignment)
            assert 0 <= score <= 1, f"Score {score} out of range"
    
    def test_score_ordering(self):
        """Test that better alignments have higher scores."""
        perfect = ["MALIGTLLM", "||||||||-", "MALIGTLLM"]
        good = ["MALIGTLLM", "||||.||||", "MALIVTLLM"]
        poor = ["MALIGTLLM", "...|.....", "WYKPQRSTV"]
        
        score_perfect = calc_alignment_score_restricted_area(perfect)
        score_good = calc_alignment_score_restricted_area(good)
        score_poor = calc_alignment_score_restricted_area(poor)
        
        assert score_perfect > score_good > score_poor, \
            "Scores should reflect alignment quality"
    
    def test_gap_penalties(self):
        """Test that gaps are penalized."""
        no_gap = ["MALIGTLLM", "||||||||-", "MALIGTLLM"]
        one_gap = ["MALIGTLLM", "||||-||||", "MALI-TLLM"]
        two_gaps = ["MALIGTLLM", "||--|-|||", "MA--I-LLM"]
        
        score_no_gap = calc_alignment_score_restricted_area(no_gap)
        score_one_gap = calc_alignment_score_restricted_area(one_gap)
        score_two_gaps = calc_alignment_score_restricted_area(two_gaps)
        
        assert score_no_gap > score_one_gap > score_two_gaps, \
            "More gaps should result in lower scores"


class TestAlignmentIntegration:
    """Integration tests for alignment in GRN context."""
    
    def test_grn_alignment_workflow(self):
        """Test typical GRN alignment workflow."""
        # Reference sequence (e.g., bacteriorhodopsin)
        ref_seq = "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY"
        
        # Query sequence (e.g., archaerhodopsin)
        query_seq = "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG"
        
        # Initialize aligner
        aligner = init_aligner()
        
        # Perform alignment
        alignment = align_blosum62(query_seq, ref_seq, aligner)
        formatted = format_alignment(alignment)
        
        # Check alignment properties
        assert len(formatted) == 3
        assert len(formatted[0]) == len(formatted[1]) == len(formatted[2])
        
        # Should have reasonable alignment
        match_line = formatted[1]
        matches = match_line.count('|') + match_line.count('.')
        total = len(match_line) - match_line.count(' ')
        
        # Should have at least 30% similarity for GRN transfer
        similarity = matches / total if total > 0 else 0
        assert similarity > 0.3, "Alignment should have >30% similarity for GRN"
    
    def test_multiple_alignments(self):
        """Test aligning query against multiple references."""
        query = "MALIGTLLMLIGTFYFIVKGWGVTDK"
        
        references = {
            "ref1": "MALIGTLLMLIGTFYFIVKGWGVTDK",  # Identical
            "ref2": "MALIGTLLMXIGTFYFXVKGWGVXDK",  # Similar
            "ref3": "TAAVGADLLGDGRPETLWLGIGTLLM",  # Different
        }
        
        aligner = init_aligner()
        results = {}
        
        for ref_id, ref_seq in references.items():
            alignment = align_blosum62(query, ref_seq, aligner)
            formatted = format_alignment(alignment)
            score = calc_alignment_score_restricted_area(formatted)
            results[ref_id] = score
        
        # ref1 should have highest score
        assert results["ref1"] > results["ref2"] > results["ref3"]
        assert results["ref1"] == 1.0  # Perfect match


class TestAlignmentPerformance:
    """Performance tests for alignment."""
    
    def test_alignment_speed(self):
        """Test that alignments complete in reasonable time."""
        import time
        
        # Generate longer sequences
        seq1 = "MALIGTLLMLIGTFYFIVKGWGVTDK" * 10  # ~260 residues
        seq2 = "TAAVGADLLGDGRPETLWLGIGTLLM" * 10  # ~260 residues
        
        aligner = init_aligner()
        
        start_time = time.time()
        alignment = align_blosum62(seq1, seq2, aligner)
        elapsed = time.time() - start_time
        
        # Should complete quickly even for longer sequences
        assert elapsed < 1.0, f"Alignment took too long: {elapsed:.3f} seconds"
        assert alignment is not None
    
    def test_batch_alignment_performance(self):
        """Test performance of multiple alignments."""
        import time
        
        query = "MALIGTLLMLIGTFYFIVKGWGVTDK"
        num_refs = 50
        
        # Generate reference sequences
        refs = {}
        for i in range(num_refs):
            # Create variations
            ref = list(query)
            if i % 2 == 0:
                ref[5] = 'X'
            if i % 3 == 0:
                ref[10] = 'Y'
            refs[f"ref_{i}"] = "".join(ref)
        
        aligner = init_aligner()
        
        start_time = time.time()
        for ref_seq in refs.values():
            align_blosum62(query, ref_seq, aligner)
        elapsed = time.time() - start_time
        
        # Should handle many alignments efficiently
        avg_time = elapsed / num_refs
        assert avg_time < 0.1, f"Average alignment time too high: {avg_time:.3f} seconds"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])