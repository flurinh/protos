"""
Test to verify the exact output format of embeddings, especially per-residue embeddings.

This is critical for mapping embeddings to residue positions, GRNs, etc.
"""

import pytest
import torch

# Handle optional dependencies
torch = pytest.importorskip("torch", reason="PyTorch required for embedding tests")
transformers = pytest.importorskip("transformers", reason="Transformers required for embedding tests")

from protos.processing.embedding.embedding_processor import EmbeddingProcessor


class TestEmbeddingOutputFormat:
    """Test the exact output format of embeddings for residue mapping."""
    
    def test_per_residue_output_length_esm2(self):
        """Test ESM-2 per-residue output format and length."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        # Test sequences of various lengths
        test_cases = [
            ("ACDEF", 5),      # 5 residues
            ("ACDEFGHIKL", 10), # 10 residues
            ("A", 1),           # Single residue
        ]
        
        for sequence, expected_residues in test_cases:
            result = processor.embed_sequences(sequence, embedding_type="per_residue")
            
            print(f"\nSequence: {sequence}")
            print(f"Expected residues: {expected_residues}")
            print(f"Embedding shape: {result.shape}")
            print(f"Embedding length: {result.shape[0]}")
            
            # ESM-2 uses: <cls> + sequence + <eos>
            # So total length should be len(sequence) + 2
            assert result.shape[0] == expected_residues + 2, \
                f"Expected {expected_residues + 2} tokens, got {result.shape[0]}"
            
            # Check embedding dimension
            assert result.shape[1] == 320  # ESM2 t6 8m dimension
        
        processor.clear_cache()
    
    def test_residue_mapping_correspondence(self):
        """Test that embeddings can be correctly mapped to residue positions."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        sequence = "ACDEFGHIKL"
        
        # Get per-residue embeddings
        embeddings = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # For ESM-2: 
        # Index 0: <cls> token
        # Index 1-10: Residues A, C, D, E, F, G, H, I, K, L
        # Index 11: <eos> token
        
        # Extract residue embeddings (excluding special tokens)
        residue_embeddings = embeddings[1:-1]  # Skip CLS and EOS
        
        print(f"\nSequence length: {len(sequence)}")
        print(f"Total embeddings: {embeddings.shape[0]}")
        print(f"Residue embeddings: {residue_embeddings.shape[0]}")
        
        # Verify we have exactly one embedding per residue
        assert residue_embeddings.shape[0] == len(sequence), \
            f"Expected {len(sequence)} residue embeddings, got {residue_embeddings.shape[0]}"
        
        # Create mapping example
        residue_mapping = {}
        for i, aa in enumerate(sequence):
            residue_mapping[i+1] = {  # 1-indexed position
                'amino_acid': aa,
                'embedding_index': i+1,  # +1 because of CLS token
                'embedding': residue_embeddings[i]
            }
        
        # Verify mapping
        assert len(residue_mapping) == len(sequence)
        assert residue_mapping[1]['amino_acid'] == 'A'
        assert residue_mapping[10]['amino_acid'] == 'L'
        
        processor.clear_cache()
    
    def test_batch_per_residue_consistency(self):
        """Test that batched sequences maintain correct per-residue mapping."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m", 
            batch_size=2,
            device="cpu"
        )
        
        sequences = {
            "seq1": "ACDEF",      # 5 residues
            "seq2": "GHIKLMNP",   # 8 residues
            "seq3": "QRS",        # 3 residues
        }
        
        embeddings = processor.embed_sequences(sequences, embedding_type="per_residue")
        
        # Verify each sequence has correct length
        assert embeddings["seq1"].shape[0] == 7   # 5 + 2 special tokens
        assert embeddings["seq2"].shape[0] == 10  # 8 + 2 special tokens
        assert embeddings["seq3"].shape[0] == 5   # 3 + 2 special tokens
        
        # Verify all have same embedding dimension
        for emb in embeddings.values():
            assert emb.shape[1] == 320
        
        processor.clear_cache()
    
    def test_special_tokens_positions(self):
        """Verify the position of special tokens in the output."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        sequence = "ACDEFG"
        
        # Get different embedding types
        per_res = processor.embed_sequences(sequence, embedding_type="per_residue")
        cls_emb = processor.embed_sequences(sequence, embedding_type="cls")
        
        # CLS embedding should match first position of per-residue
        assert torch.allclose(cls_emb, per_res[0], rtol=1e-5), \
            "CLS embedding doesn't match first token"
        
        print(f"\nPer-residue shape: {per_res.shape}")
        print(f"CLS embedding shape: {cls_emb.shape}")
        print(f"CLS matches first token: {torch.allclose(cls_emb, per_res[0])}")
        
        processor.clear_cache()
    
    def test_create_residue_embedding_mapping(self):
        """Test creating a practical mapping for downstream use."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        # Example with a protein sequence
        sequence = "MLELLPTAVE"  # First 10 residues of bacteriorhodopsin
        
        # Get embeddings
        embeddings = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # Create practical mapping for downstream use
        # This is how you would use it with GRN assignments
        residue_data = []
        
        for i, aa in enumerate(sequence):
            residue_info = {
                'position': i + 1,  # 1-indexed
                'amino_acid': aa,
                'embedding': embeddings[i + 1].numpy(),  # +1 to skip CLS
                # 'grn': grn_assignments.get(i + 1),  # Would come from GRN processor
                # 'secondary_structure': ss_assignments.get(i + 1),  # From structure
            }
            residue_data.append(residue_info)
        
        # Verify mapping
        assert len(residue_data) == len(sequence)
        assert residue_data[0]['amino_acid'] == 'M'
        assert residue_data[0]['position'] == 1
        assert residue_data[-1]['amino_acid'] == 'E'
        assert residue_data[-1]['position'] == 10
        
        # All embeddings should have correct dimension
        for res in residue_data:
            assert res['embedding'].shape == (320,)
        
        print(f"\nCreated mapping for {len(residue_data)} residues")
        print(f"First residue: {residue_data[0]['amino_acid']} at position {residue_data[0]['position']}")
        print(f"Last residue: {residue_data[-1]['amino_acid']} at position {residue_data[-1]['position']}")
        
        processor.clear_cache()
    
    def test_compare_model_families_tokenization(self):
        """Compare how different model families handle tokenization."""
        sequence = "ACDEFG"
        
        # Test ESM-2
        esm_processor = EmbeddingProcessor(
            name="esm_test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        esm_embeddings = esm_processor.embed_sequences(sequence, embedding_type="per_residue")
        print(f"\nESM-2 output shape: {esm_embeddings.shape}")
        print(f"ESM-2 tokens: {esm_embeddings.shape[0]} (sequence length: {len(sequence)})")
        esm_processor.clear_cache()
        
        # Note: We skip testing Ankh here because it's based on T5 architecture
        # which handles tokenization differently than BERT-based models like ESM-2.
        # T5 models don't use CLS/SEP tokens in the same way.
        
        # ESM-2 should add special tokens
        assert esm_embeddings.shape[0] == len(sequence) + 2  # CLS + sequence + EOS


class TestPracticalUsage:
    """Test practical usage patterns for residue-level analysis."""
    
    def test_get_residue_embeddings_method(self):
        """Test the built-in get_residue_embeddings helper method."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        sequence = "ACDEFGHIKL"
        full_embeddings = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # Test extracting residue embeddings
        residue_embeddings = processor.get_residue_embeddings(full_embeddings)
        
        # Should have exactly one embedding per residue
        assert residue_embeddings.shape[0] == len(sequence)
        assert residue_embeddings.shape[1] == 320
        
        # Test with special tokens included
        full_with_special = processor.get_residue_embeddings(full_embeddings, include_special_tokens=True)
        assert torch.equal(full_with_special, full_embeddings)
        
        # Verify CLS and EOS are excluded in residue embeddings
        assert not torch.equal(residue_embeddings[0], full_embeddings[0])  # First residue != CLS
        assert not torch.equal(residue_embeddings[-1], full_embeddings[-1])  # Last residue != EOS
        
        processor.clear_cache()
    
    def test_embedding_extraction_helper(self):
        """Test a helper function for extracting residue embeddings."""
        
        def extract_residue_embeddings(embeddings, skip_special_tokens=True):
            """
            Extract residue embeddings from per-residue output.
            
            Args:
                embeddings: Per-residue embeddings from EmbeddingProcessor
                skip_special_tokens: If True, remove CLS and EOS tokens
                
            Returns:
                Tensor of shape (num_residues, embedding_dim)
            """
            if skip_special_tokens:
                # Skip first (CLS) and last (EOS) tokens
                return embeddings[1:-1]
            return embeddings
        
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        sequence = "ACDEFGHIKL"
        full_embeddings = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # Extract just residue embeddings
        residue_embeddings = extract_residue_embeddings(full_embeddings)
        
        assert residue_embeddings.shape[0] == len(sequence)
        assert residue_embeddings.shape[1] == 320
        
        # Test with special tokens included
        full_with_special = extract_residue_embeddings(full_embeddings, skip_special_tokens=False)
        assert full_with_special.shape[0] == len(sequence) + 2
        
        processor.clear_cache()
    
    def test_embedding_grn_alignment(self):
        """Test aligning embeddings with GRN positions (simulation)."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        # Simulate a sequence with GRN assignments
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        
        # Simulate GRN assignments (position -> GRN)
        grn_assignments = {
            1: "1.50", 2: "1.51", 3: "1.52", 4: "1.53", 5: "1.54",
            6: "1.55", 7: "1.56", 8: "1.57", 9: "1.58", 10: "1.59",
            11: "2.50", 12: "2.51", 13: "2.52", 14: "2.53", 15: "2.54",
            16: "2.55", 17: "2.56", 18: "2.57", 19: "2.58", 20: "2.59"
        }
        
        # Get embeddings
        embeddings = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # Create GRN -> embedding mapping
        grn_embedding_map = {}
        for pos, grn in grn_assignments.items():
            # pos is 1-indexed, embedding index needs +1 for CLS token
            embedding_idx = pos  # Because we skip index 0 (CLS)
            if embedding_idx < embeddings.shape[0] - 1:  # Don't include EOS
                grn_embedding_map[grn] = embeddings[embedding_idx]
        
        # Verify mapping
        assert len(grn_embedding_map) == len(grn_assignments)
        assert "1.50" in grn_embedding_map
        assert "2.59" in grn_embedding_map
        
        # All embeddings should have correct shape
        for grn, emb in grn_embedding_map.items():
            assert emb.shape == (320,)
        
        print(f"\nCreated GRN embedding map with {len(grn_embedding_map)} positions")
        print(f"Sample GRNs: {list(grn_embedding_map.keys())[:5]}")
        
        processor.clear_cache()