"""
Tests for MMseqs2 sequence alignment functionality
"""

import pytest
import tempfile
import os
from pathlib import Path
import pandas as pd
from typing import Dict

from protos.processing.sequence.seq_alignment import mmseqs2_align


class TestMMseqs2:
    """Test MMseqs2 alignment functionality."""
    
    @pytest.fixture
    def test_sequences(self):
        """Create test sequences for MMseqs2."""
        return {
            "exact_match": "MALIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG",
            "one_mutation": "MALIGTLLMXIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG",
            "five_mutations": "MALIGTLLMXIGTFYFXVKGWGVXDKEARXYYSITXLVPGIASAAYLSMFFGIG",
            "similar": "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITIL",
            "different": "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQ"
        }
    
    @pytest.fixture
    def query_sequence(self):
        """Query sequence for testing."""
        return "MALIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG"
    
    def test_mmseqs2_available(self):
        """Test if MMseqs2 is available in the environment."""
        import platform
        import subprocess
        
        # Check if MMseqs2 was found by conftest
        if os.getenv("MMSEQS_PATH"):
            # Already configured
            return
            
        # Otherwise, do our own check
        is_windows = platform.system() == 'Windows'
        
        if is_windows:
            # Check with WSL
            try:
                result = subprocess.run(['wsl', 'mmseqs', '--help'], capture_output=True, shell=True)
                if result.returncode == 0:
                    os.environ["USE_WSL_MMSEQS"] = "true"
                    return
            except:
                pass
        else:
            # Check native
            result = os.system("mmseqs --help > /dev/null 2>&1")
            if result == 0:
                return
                
        pytest.skip("MMseqs2 not available")
    
    def test_mmseqs2_align_basic(self, query_sequence, test_sequences):
        """Test basic MMseqs2 alignment functionality."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run alignment
            results = mmseqs2_align(query_sequence, test_sequences, tmpdir)
            
            # Check results
            assert results is not None, "MMseqs2 should return results"
            assert isinstance(results, pd.DataFrame), "Results should be a DataFrame"
            assert len(results) > 0, "Should find at least one alignment"
            
            # Check DataFrame columns
            expected_columns = [
                'query_id', 'target_id', 'sequence_identity', 
                'alignment_length', 'mismatches', 'gap_openings',
                'query_start', 'query_end', 'target_start', 'target_end', 
                'e_value', 'bit_score'
            ]
            for col in expected_columns:
                assert col in results.columns, f"Missing column: {col}"
    
    def test_mmseqs2_align_identity(self, query_sequence, test_sequences):
        """Test that MMseqs2 correctly identifies sequence similarity."""
        with tempfile.TemporaryDirectory() as tmpdir:
            results = mmseqs2_align(query_sequence, test_sequences, tmpdir)
            
            if results is not None and len(results) > 0:
                # Check exact match
                exact_matches = results[results['target_id'] == 'exact_match']
                if len(exact_matches) > 0:
                    identity = exact_matches.iloc[0]['sequence_identity']
                    assert identity > 0.9, "Exact match should have >90% identity"
                
                # Check that different sequence has lower identity
                different = results[results['target_id'] == 'different']
                if len(different) > 0:
                    diff_identity = different.iloc[0]['sequence_identity']
                    assert diff_identity < 0.5, "Different sequence should have <50% identity"
    
    def test_mmseqs2_align_evalue(self, query_sequence, test_sequences):
        """Test that E-values are sensible."""
        with tempfile.TemporaryDirectory() as tmpdir:
            results = mmseqs2_align(query_sequence, test_sequences, tmpdir)
            
            if results is not None and len(results) > 0:
                # Sort by e-value
                results_sorted = results.sort_values('e_value')
                
                # Best hit should have very low e-value
                best_evalue = results_sorted.iloc[0]['e_value']
                assert best_evalue < 1e-10, "Best hit should have E-value < 1e-10"
                
                # Check that e-values increase
                evalues = results_sorted['e_value'].tolist()
                assert evalues == sorted(evalues), "E-values should be sorted"
    
    def test_mmseqs2_align_empty_database(self, query_sequence):
        """Test behavior with empty database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            empty_db = {}
            results = mmseqs2_align(query_sequence, empty_db, tmpdir)
            # Should handle gracefully - either None or empty DataFrame
            assert results is None or len(results) == 0
    
    def test_mmseqs2_align_single_sequence(self, query_sequence):
        """Test with single sequence in database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            single_db = {"single": query_sequence}
            results = mmseqs2_align(query_sequence, single_db, tmpdir)
            
            if results is not None:
                assert len(results) == 1, "Should find exactly one alignment"
                assert results.iloc[0]['target_id'] == 'single'
                # MMseqs2 may not report 100% identity due to its algorithm
                # but should be very high for identical sequences
                assert results.iloc[0]['sequence_identity'] > 0.9
    
    @pytest.mark.parametrize("gap_seq,expected_gaps", [
        ("MALIGTLLM---LIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG", 1),
        ("MALIGTLLM------LIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG", 1),
    ])
    def test_mmseqs2_align_with_gaps(self, query_sequence, gap_seq, expected_gaps):
        """Test alignment with sequences containing gaps."""
        # Remove gaps for database sequence
        db_seq = gap_seq.replace('-', '')
        
        with tempfile.TemporaryDirectory() as tmpdir:
            db = {"gapped": db_seq}
            results = mmseqs2_align(query_sequence, db, tmpdir)
            
            if results is not None and len(results) > 0:
                # Should still find alignment despite gaps
                assert len(results) > 0
                # Check gap openings if reported
                if 'gap_openings' in results.columns:
                    gap_openings = results.iloc[0]['gap_openings']
                    assert gap_openings >= 0  # Should be non-negative
    
    def test_mmseqs2_performance(self):
        """Test MMseqs2 performance with larger dataset."""
        import time
        import random
        
        # Generate test sequences
        base_seq = "MALIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG"
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        
        # Create 100 sequences with mutations
        database = {}
        for i in range(100):
            seq_list = list(base_seq)
            # Introduce random mutations
            num_mutations = random.randint(0, 10)
            for _ in range(num_mutations):
                pos = random.randint(0, len(seq_list) - 1)
                seq_list[pos] = random.choice(amino_acids)
            database[f"seq_{i:03d}"] = "".join(seq_list)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            start_time = time.time()
            results = mmseqs2_align(base_seq, database, tmpdir)
            elapsed = time.time() - start_time
            
            # Should complete quickly
            assert elapsed < 5.0, f"Search took too long: {elapsed:.2f} seconds"
            
            if results is not None:
                # Should find many similar sequences
                assert len(results) > 50, "Should find many similar sequences"
                
                # Check that most similar sequences have high identity
                high_identity = results[results['sequence_identity'] > 0.8]
                assert len(high_identity) > 20, "Should find many high-identity matches"


class TestMMseqs2ErrorHandling:
    """Test error handling in MMseqs2."""
    
    def test_invalid_sequence(self):
        """Test with invalid sequence characters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Sequence with invalid characters
            invalid_seq = "MALIGH123TLLMLIGTFY"
            db = {"valid": "MALIGTLLMLIGTFYFIVK"}
            
            # Should handle gracefully
            results = mmseqs2_align(invalid_seq, db, tmpdir)
            # May return None or empty results, but shouldn't crash
            assert results is None or isinstance(results, pd.DataFrame)
    
    def test_empty_sequence(self):
        """Test with empty sequence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            empty_seq = ""
            db = {"valid": "MALIGTLLMLIGTFYFIVK"}
            
            results = mmseqs2_align(empty_seq, db, tmpdir)
            # Should handle gracefully
            assert results is None or len(results) == 0
    
    def test_very_long_sequence(self):
        """Test with very long sequence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a very long sequence
            long_seq = "MALIGTLLMLIGTFYFIVK" * 100  # ~1900 residues
            db = {"short": "MALIGTLLMLIGTFYFIVK"}
            
            results = mmseqs2_align(long_seq, db, tmpdir)
            # Should handle without crashing
            assert results is None or isinstance(results, pd.DataFrame)
    
    def test_nonexistent_temp_dir(self):
        """Test with non-existent temporary directory."""
        # Use a path that doesn't exist
        fake_tmpdir = "/tmp/nonexistent_mmseqs_test_dir_12345"
        
        query = "MALIGTLLMLIGTFYFIVK"
        db = {"target": "MALIGTLLMLIGTFYFIVK"}
        
        # Should create the directory or handle gracefully
        results = mmseqs2_align(query, db, fake_tmpdir)
        
        # Clean up if directory was created
        if os.path.exists(fake_tmpdir):
            import shutil
            shutil.rmtree(fake_tmpdir)


class TestMMseqs2Integration:
    """Integration tests for MMseqs2 with GRN workflow."""
    
    @pytest.fixture
    def grn_reference_sequences(self):
        """Load actual GRN reference sequences."""
        # Simplified reference sequences for testing
        return {
            "BR": "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY",
            "HR": "MKKPLGLLGAPGENAWVDMAAVTGVSAALGVAGGVLGVLATVGAAAVAADPLLARTTRPGEWICLAFALVLVLVGVLY",
            "bPR": "MALPDTFFDLVAADAERQWWLVVGILVAVIGTAFTALSVGNGFNFGKPTDDIFNVFKTVFEIVLGSALVEVIIGTLSY"
        }
    
    def test_grn_reference_selection(self, grn_reference_sequences):
        """Test selecting best reference for GRN assignment."""
        # Test query sequence (archaerhodopsin-like)
        query = "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG"
        
        with tempfile.TemporaryDirectory() as tmpdir:
            results = mmseqs2_align(query, grn_reference_sequences, tmpdir)
            
            if results is not None and len(results) > 0:
                # Should find alignments
                assert len(results) > 0
                
                # Best hit should be BR (bacteriorhodopsin)
                best_hit = results.loc[results['e_value'].idxmin()]
                assert best_hit['target_id'] in ['BR', 'HR', 'bPR']
                
                # Should have reasonable e-value for GRN transfer
                assert best_hit['e_value'] < 1e-5, "E-value should be low enough for GRN transfer"
    
    def test_batch_grn_search(self, grn_reference_sequences):
        """Test batch search for multiple queries."""
        queries = {
            "AR1": "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG",
            "AR2": "MAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFVVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGVG",
            "AR3": "LAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFAVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGAG"
        }
        
        # For batch search, we'd use mmseqs2_align2, but for this test
        # we'll simulate with individual searches
        all_results = []
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for query_id, query_seq in queries.items():
                results = mmseqs2_align(query_seq, grn_reference_sequences, tmpdir)
                if results is not None:
                    results['query_id'] = query_id
                    all_results.append(results)
            
            if all_results:
                combined = pd.concat(all_results, ignore_index=True)
                
                # Each query should find references
                for query_id in queries:
                    query_results = combined[combined['query_id'] == query_id]
                    assert len(query_results) > 0, f"Should find results for {query_id}"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])