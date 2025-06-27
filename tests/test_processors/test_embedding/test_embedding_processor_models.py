"""
Tests for EmbeddingProcessor models with different output formats and special sequences.

This module tests representative models from each family (ESM-2 and Ankh) using real models
with minimal memory usage, focusing on output format validation and special sequence handling.
"""

import pytest
import torch
import numpy as np
import os

# Handle optional dependencies
torch = pytest.importorskip("torch", reason="PyTorch required for embedding tests")
transformers = pytest.importorskip("transformers", reason="Transformers required for embedding tests")

from protos.processing.embedding.embedding_processor import EmbeddingProcessor


# Test only small models to save memory and time
# One model per family is sufficient to test output formats
TEST_MODELS = {
    "esm2_t6_8m": {"family": "esm2", "dim": 320},   # Smallest ESM-2 model
    "ankh_base": {"family": "ankh", "dim": 768}     # Smaller Ankh model
}


@pytest.fixture(scope="module")
def device():
    """Get the appropriate device for testing."""
    if torch.cuda.is_available() and torch.cuda.get_device_properties(0).total_memory > 2e9:  # 2GB
        return "cuda"
    return "cpu"


class TestModelOutputFormats:
    """Test output formats for representative models from each family."""
    
    @pytest.mark.parametrize("model_name,model_info", TEST_MODELS.items())
    def test_single_sequence_all_embedding_types(self, model_name, model_info, device):
        """Test single sequence embedding with all output types."""
        # Use GPU for ESM2 tiny model if available, CPU for others
        test_device = device if model_name == "esm2_t6_8m" else "cpu"
        
        processor = EmbeddingProcessor(
            name="test", 
            model_name=model_name,
            device=test_device
        )
        
        # Very short sequence to minimize memory
        sequence = "ACDEFG"
        
        # Test mean embedding
        mean_emb = processor.embed_sequences(sequence, embedding_type="mean")
        assert isinstance(mean_emb, torch.Tensor)
        assert mean_emb.shape == (model_info["dim"],)
        assert mean_emb.device.type == "cpu"  # Should be moved to CPU after generation
        
        # Test CLS embedding
        cls_emb = processor.embed_sequences(sequence, embedding_type="cls")
        assert isinstance(cls_emb, torch.Tensor)
        assert cls_emb.shape == (model_info["dim"],)
        
        # Test per-residue embedding
        per_res_emb = processor.embed_sequences(sequence, embedding_type="per_residue")
        assert isinstance(per_res_emb, torch.Tensor)
        assert per_res_emb.dim() == 2
        assert per_res_emb.shape[1] == model_info["dim"]
        # Should have embeddings for each token (including special tokens)
        assert per_res_emb.shape[0] >= len(sequence)
        
        # Clear cache to free memory
        processor.clear_cache()
    
    def test_dict_batching_behavior(self, device):
        """Test dictionary input with batching using smallest model."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name="esm2_t6_8m",
            batch_size=2,
            device=device
        )
        
        # Three very short sequences to test batching
        sequences = {
            "seq1": "ACD",
            "seq2": "FGH", 
            "seq3": "KLM"
        }
        
        embeddings = processor.embed_sequences(sequences, embedding_type="mean")
        
        # Check all embeddings are present and correct shape
        assert len(embeddings) == 3
        for seq_id, emb in embeddings.items():
            assert isinstance(emb, torch.Tensor)
            assert emb.shape == (320,)  # esm2_t6_8m dimension
            assert emb.device.type == "cpu"
        
        processor.clear_cache()
    
    @pytest.mark.parametrize("model_name", ["esm2_t6_8m"])  # Test only with smallest model
    def test_special_sequences(self, model_name):
        """Test handling of special sequences with X, U, and masked positions."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name=model_name,
            device="cpu"  # Use CPU to save memory
        )
        
        special_sequences = {
            "normal": "ACDEFG",
            "with_X": "ACDXFG",
            "with_U": "ACDUFG",
            "with_mask": "ACD<mask>FG",
            "mixed": "AXDUFG",
            "all_unknown": "XXX"
        }
        
        embeddings = processor.embed_sequences(special_sequences, embedding_type="mean")
        
        # All sequences should produce valid embeddings
        assert len(embeddings) == len(special_sequences)
        for seq_id, emb in embeddings.items():
            assert isinstance(emb, torch.Tensor)
            assert emb.shape == (320,)
            assert not torch.isnan(emb).any(), f"NaN in embedding for {seq_id}"
            assert not torch.isinf(emb).any(), f"Inf in embedding for {seq_id}"
        
        processor.clear_cache()
    
    def test_edge_cases(self):
        """Test edge cases like very short sequences using smallest model."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        edge_cases = {
            "single_aa": "A",
            "two_aa": "AC",
            "repeated": "AAA"
        }
        
        embeddings = processor.embed_sequences(edge_cases, embedding_type="mean")
        
        assert len(embeddings) == len(edge_cases)
        for seq_id, emb in embeddings.items():
            assert isinstance(emb, torch.Tensor)
            assert emb.shape == (320,)
        
        processor.clear_cache()


class TestBatchProcessing:
    """Test batch processing with real models."""
    
    def test_variable_length_sequences_in_batch(self):
        """Test batching with sequences of different lengths."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            batch_size=2,
            device="cpu"
        )
        
        sequences = {
            "short": "AC",
            "medium": "ACDEFG",
            "long": "ACDEFGHIKLM"
        }
        
        embeddings = processor.embed_sequences(sequences, embedding_type="mean")
        
        # All should produce same dimension embeddings despite different lengths
        assert all(emb.shape == (320,) for emb in embeddings.values())
        
        processor.clear_cache()
    
    def test_per_residue_length_handling(self):
        """Test that per-residue embeddings handle sequence length correctly."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            max_seq_length=20,  # Small for testing
            device="cpu"
        )
        
        sequences = {
            "short": "ACDE",           # 4 residues
            "medium": "ACDEFGHIKL",    # 10 residues
        }
        
        embeddings = processor.embed_sequences(sequences, embedding_type="per_residue")
        
        # Check that each sequence has appropriate length
        # Note: actual length includes special tokens (CLS, EOS)
        assert embeddings["short"].shape[0] >= 4  # At least sequence length
        assert embeddings["medium"].shape[0] >= 10
        
        # All should have correct embedding dimension
        assert all(emb.shape[1] == 320 for emb in embeddings.values())
        
        processor.clear_cache()


class TestMemoryEfficiency:
    """Test memory-efficient processing."""
    
    def test_clear_cache_functionality(self):
        """Test that cache clearing properly frees memory."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name="esm2_t6_8m",
            device="cpu"
        )
        
        # Generate some embeddings to load model
        processor.embed_sequences("ACDEFG")
        
        # Model should be loaded
        assert processor._model is not None
        assert processor._tokenizer is not None
        
        # Clear cache
        processor.clear_cache()
        
        # Model should be cleared
        assert processor._model is None
        assert processor._tokenizer is None
    
    def test_cpu_processing(self):
        """Test that embeddings work correctly on CPU."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cpu"  # Force CPU
        )
        
        embeddings = processor.embed_sequences(
            {"seq1": "ACDEFG", "seq2": "KLMNPQ"}, 
            embedding_type="mean"
        )
        
        # All embeddings should be on CPU
        for emb in embeddings.values():
            assert emb.device.type == 'cpu'
        
        processor.clear_cache()


@pytest.mark.skipif(not torch.cuda.is_available(), reason="CUDA not available")
class TestGPUProcessing:
    """Test GPU processing (only runs if CUDA available)."""
    
    def test_gpu_processing(self):
        """Test that model runs on GPU and embeddings are transferred to CPU."""
        if torch.cuda.get_device_properties(0).total_memory < 2e9:  # Skip if less than 2GB
            pytest.skip("Insufficient GPU memory")
        
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            device="cuda"
        )
        
        # Check model is on GPU
        assert processor.device == "cuda"
        
        embeddings = processor.embed_sequences("ACDEFG", embedding_type="mean")
        
        # Embeddings should be on CPU for storage
        assert embeddings.device.type == 'cpu'
        
        processor.clear_cache()


class TestModelComparison:
    """Test that different model families produce expected output formats."""
    
    def test_esm2_vs_ankh_output_format(self):
        """Compare output formats between ESM-2 and Ankh models."""
        sequences = {"test": "ACDEFG"}
        
        # Test ESM-2
        esm_processor = EmbeddingProcessor(
            name="esm_test",
            model_name="esm2_t6_8m",
            device="cpu"
        )
        esm_embeddings = esm_processor.embed_sequences(sequences, embedding_type="mean")
        assert esm_embeddings["test"].shape == (320,)
        esm_processor.clear_cache()
        
        # Test Ankh  
        ankh_processor = EmbeddingProcessor(
            name="ankh_test",
            model_name="ankh_base",
            device="cpu"
        )
        ankh_embeddings = ankh_processor.embed_sequences(sequences, embedding_type="mean")
        assert ankh_embeddings["test"].shape == (768,)
        ankh_processor.clear_cache()
        
        # Both should produce tensors but with different dimensions
        assert isinstance(esm_embeddings["test"], torch.Tensor)
        assert isinstance(ankh_embeddings["test"], torch.Tensor)
        assert esm_embeddings["test"].shape != ankh_embeddings["test"].shape