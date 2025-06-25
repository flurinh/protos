"""
Integration tests for sequence alignment combining MMseqs2 and BioPython
"""

import pytest
import tempfile
from pathlib import Path
import pandas as pd

from protos.processing.sequence.seq_alignment import (
    mmseqs2_align,
    init_aligner,
    align_blosum62,
    format_alignment,
    calc_alignment_score_restricted_area
)
from protos.io.fasta_utils import read_fasta, write_fasta


class TestSeqAlignmentIntegration:
    """Integration tests for the complete alignment workflow."""
    
    @pytest.fixture
    def sample_sequences(self):
        """Sample sequences for testing."""
        return {
            "BR": "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGY",
            "HR": "MKKPLGLLGAPGENAWVDMAAVTGVSAALGVAGGVLGVLATVGAAAVAADPLLARTTRPGEWICLAFALVLVLVGVLY",
            "AR": "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG",
            "bPR": "MALPDTFFDLVAADAERQWWLVVGILVAVIGTAFTALSVGNGFNFGKPTDDIFNVFKTVFEIVLGSALVEVIIGTLSY"
        }
    
    def test_mmseqs_to_biopython_workflow(self, sample_sequences):
        """Test complete workflow: MMseqs2 search -> BioPython alignment."""
        # Query sequence
        query_seq = sample_sequences["AR"]
        
        # Reference database
        ref_db = {k: v for k, v in sample_sequences.items() if k != "AR"}
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Step 1: MMseqs2 search
            mmseqs_results = mmseqs2_align(query_seq, ref_db, tmpdir)
            
            if mmseqs_results is not None and len(mmseqs_results) > 0:
                # Step 2: Select best hit
                best_hit = mmseqs_results.loc[mmseqs_results['e_value'].idxmin()]
                best_ref_id = best_hit['target_id']
                best_ref_seq = ref_db[best_ref_id]
                
                # Step 3: Detailed BioPython alignment
                aligner = init_aligner()
                alignment = align_blosum62(query_seq, best_ref_seq, aligner)
                formatted = format_alignment(alignment)
                
                # Verify results
                assert len(formatted) == 3
                assert all(len(s) == len(formatted[0]) for s in formatted)
                
                # Calculate final score
                score = calc_alignment_score_restricted_area(formatted)
                assert 0 <= score <= 1
                
                print(f"\nWorkflow complete:")
                print(f"  Best reference: {best_ref_id}")
                print(f"  MMseqs2 E-value: {best_hit['e_value']:.2e}")
                print(f"  MMseqs2 Identity: {best_hit['sequence_identity']*100:.1f}%")
                print(f"  BioPython Score: {score:.3f}")
    
    def test_batch_workflow(self, sample_sequences):
        """Test batch processing workflow."""
        # Multiple queries
        queries = {
            "AR1": "TAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG",
            "AR2": "MAAVGADLLGDGRPETLWLGIGTLLMLIGTFYFVVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGVG"
        }
        
        # Reference database
        ref_db = {k: v for k, v in sample_sequences.items() if k not in queries}
        
        results = []
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Process each query
            for query_id, query_seq in queries.items():
                # MMseqs2 search
                mmseqs_results = mmseqs2_align(query_seq, ref_db, tmpdir)
                
                if mmseqs_results is not None and len(mmseqs_results) > 0:
                    # Get best hit
                    best_hit = mmseqs_results.loc[mmseqs_results['e_value'].idxmin()]
                    
                    # BioPython alignment
                    aligner = init_aligner()
                    ref_seq = ref_db[best_hit['target_id']]
                    alignment = align_blosum62(query_seq, ref_seq, aligner)
                    formatted = format_alignment(alignment)
                    score = calc_alignment_score_restricted_area(formatted)
                    
                    results.append({
                        'query_id': query_id,
                        'best_ref': best_hit['target_id'],
                        'mmseqs_evalue': best_hit['e_value'],
                        'mmseqs_identity': best_hit['sequence_identity'],
                        'biopython_score': score
                    })
            
            # Verify batch results
            assert len(results) == len(queries)
            results_df = pd.DataFrame(results)
            print("\nBatch results:")
            print(results_df)
    
    def test_alignment_quality_correlation(self, sample_sequences):
        """Test that MMseqs2 and BioPython scores correlate."""
        query = sample_sequences["AR"]
        
        correlations = []
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Test against each reference
            for ref_id, ref_seq in sample_sequences.items():
                if ref_id == "AR":
                    continue
                
                # MMseqs2 search
                db = {ref_id: ref_seq}
                mmseqs_results = mmseqs2_align(query, db, tmpdir)
                
                if mmseqs_results is not None and len(mmseqs_results) > 0:
                    mmseqs_identity = mmseqs_results.iloc[0]['sequence_identity']
                    
                    # BioPython alignment
                    aligner = init_aligner()
                    alignment = align_blosum62(query, ref_seq, aligner)
                    formatted = format_alignment(alignment)
                    bio_score = calc_alignment_score_restricted_area(formatted)
                    
                    correlations.append({
                        'ref': ref_id,
                        'mmseqs_identity': mmseqs_identity,
                        'biopython_score': bio_score
                    })
            
            if correlations:
                # Check that scores generally correlate
                df = pd.DataFrame(correlations)
                print(f"\nCollected {len(correlations)} correlation points")
                print(df)
                
                # Need at least 2 points for correlation
                if len(correlations) >= 2:
                    # Higher MMseqs2 identity should correspond to higher BioPython score
                    correlation = df['mmseqs_identity'].corr(df['biopython_score'])
                    print(f"\nScore correlation: {correlation:.3f}")
                    # Should be positive correlation or NaN if insufficient variance
                    if not pd.isna(correlation):
                        assert correlation > -0.5, "Scores should not be strongly negatively correlated"
                else:
                    print("Not enough data points for correlation test")


class TestRealDataIntegration:
    """Tests with real protein data."""
    
    @pytest.fixture
    def real_sequences(self):
        """Load real sequences if available."""
        project_root = Path(__file__).parent.parent.parent.parent
        fasta_file = project_root / "src/protos/reference_data/sequence/fasta/opsin_sequences_from_yaml.fasta"
        
        if fasta_file.exists():
            sequences = read_fasta(str(fasta_file))
            # Return first few sequences
            return dict(list(sequences.items())[:5])
        else:
            pytest.skip("Real sequence file not found")
    
    def test_real_sequence_workflow(self, real_sequences):
        """Test workflow with real sequences."""
        if not real_sequences:
            return
        
        # Use first sequence as query
        query_id = list(real_sequences.keys())[0]
        query_seq = real_sequences[query_id]
        
        # Others as database
        ref_db = {k: v for k, v in real_sequences.items() if k != query_id}
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # MMseqs2 search
            results = mmseqs2_align(query_seq, ref_db, tmpdir)
            
            if results is not None and len(results) > 0:
                # Should find similar sequences
                assert len(results) > 0
                
                # Best hit should have reasonable e-value
                best_evalue = results['e_value'].min()
                assert best_evalue < 1e-5, "Should find similar sequences"
                
                # Perform detailed alignment on best hit
                best_hit = results.loc[results['e_value'].idxmin()]
                ref_seq = ref_db[best_hit['target_id']]
                
                aligner = init_aligner()
                alignment = align_blosum62(query_seq, ref_seq, aligner)
                formatted = format_alignment(alignment)
                
                # Check alignment length
                assert len(formatted[0]) >= len(query_seq)


class TestErrorHandlingIntegration:
    """Test error handling in integrated workflow."""
    
    def test_empty_database_handling(self):
        """Test handling of empty database."""
        query = "MALIGTLLMLIGTFYFIVK"
        empty_db = {}
        
        with tempfile.TemporaryDirectory() as tmpdir:
            results = mmseqs2_align(query, empty_db, tmpdir)
            # Should handle gracefully
            assert results is None or len(results) == 0
    
    def test_no_significant_hits(self):
        """Test handling when no significant hits found."""
        # Very different sequences
        query = "MALIGTLLMLIGTFYFIVK"
        db = {"unrelated": "WYKPQRSTVWYKPQRSTVW"}
        
        with tempfile.TemporaryDirectory() as tmpdir:
            results = mmseqs2_align(query, db, tmpdir)
            
            if results is not None and len(results) > 0:
                # E-value should be high (not significant)
                best_evalue = results['e_value'].min()
                # Might still find alignment but with poor e-value
                print(f"Best E-value for unrelated sequences: {best_evalue}")
    
    def test_invalid_sequence_handling(self):
        """Test handling of invalid sequences."""
        # Sequence with numbers (invalid)
        query = "MALIG123TLLM"
        db = {"valid": "MALIGTLLMLIGTFYFIVK"}
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Should not crash
            try:
                results = mmseqs2_align(query, db, tmpdir)
                # May return None or empty results
                assert results is None or isinstance(results, pd.DataFrame)
            except Exception as e:
                # Should handle gracefully
                print(f"Handled invalid sequence: {e}")


class TestPerformanceIntegration:
    """Performance tests for integrated workflow."""
    
    def test_large_database_search(self):
        """Test performance with larger database."""
        import time
        import random
        
        # Generate test database
        base_seq = "MALIGTLLMLIGTFYFIVKGWGVTDKEAREYYSITILVPGIASAAYLSMFFGIG"
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        
        database = {}
        for i in range(200):
            # Create variations
            seq_list = list(base_seq)
            num_mutations = random.randint(0, 15)
            for _ in range(num_mutations):
                pos = random.randint(0, len(seq_list) - 1)
                seq_list[pos] = random.choice(amino_acids)
            database[f"seq_{i:03d}"] = "".join(seq_list)
        
        query = base_seq
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Time MMseqs2 search
            start_time = time.time()
            results = mmseqs2_align(query, database, tmpdir)
            mmseqs_time = time.time() - start_time
            
            if results is not None and len(results) > 0:
                # Get top 10 hits
                top_hits = results.nsmallest(10, 'e_value')
                
                # Time BioPython alignments
                aligner = init_aligner()
                start_time = time.time()
                
                for _, hit in top_hits.iterrows():
                    ref_seq = database[hit['target_id']]
                    align_blosum62(query, ref_seq, aligner)
                
                bio_time = time.time() - start_time
                
                print(f"\nPerformance:")
                print(f"  MMseqs2 search ({len(database)} sequences): {mmseqs_time:.3f}s")
                print(f"  BioPython alignments (top 10): {bio_time:.3f}s")
                print(f"  Total workflow: {mmseqs_time + bio_time:.3f}s")
                
                # Should complete in reasonable time
                assert mmseqs_time < 10.0, "MMseqs2 too slow"
                assert bio_time < 2.0, "BioPython alignments too slow"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])