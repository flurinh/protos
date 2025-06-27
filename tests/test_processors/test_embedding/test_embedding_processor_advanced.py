"""
Advanced tests for EmbeddingProcessor functionality using real models.

This module contains comprehensive tests for the EmbeddingProcessor class,
covering edge cases, error handling, and integration scenarios with real models.
"""

import pytest
import json
import os
from pathlib import Path
import tempfile

# Handle optional dependencies
torch = pytest.importorskip("torch", reason="PyTorch required for embedding tests")
transformers = pytest.importorskip("transformers", reason="Transformers required for embedding tests")

from protos.processing.embedding.embedding_processor import EmbeddingProcessor


class TestEmbeddingProcessorDeviceHandling:
    """Test device selection and management."""
    
    def test_device_auto_selection(self):
        """Test automatic device selection."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device=None)
        
        if torch.cuda.is_available():
            assert processor.device == 'cuda'
        else:
            assert processor.device == 'cpu'
    
    def test_explicit_cpu_selection(self):
        """Test explicit CPU device selection."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device='cpu')
        assert processor.device == 'cpu'
    
    def test_cuda_memory_clearing(self):
        """Test GPU memory clearing in clear_cache."""
        if not torch.cuda.is_available():
            pytest.skip("CUDA not available")
            
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device='cuda')
        
        # Generate embedding to load model
        processor.embed_sequences("ACDEFG")
        
        # Clear cache should work without error
        processor.clear_cache()
        assert processor._model is None
        assert processor._tokenizer is None


class TestEmbeddingProcessorBatchProcessing:
    """Test batch processing functionality with real models."""
    
    def test_batch_processing_exact_multiple(self):
        """Test processing when sequence count is exact multiple of batch size."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name="esm2_t6_8m",
            batch_size=2,
            device="cpu"
        )
        
        sequences = {f"seq_{i}": "ACDEFG" for i in range(4)}
        
        embeddings = processor.embed_sequences(sequences, embedding_type="mean")
        
        assert len(embeddings) == 4
        assert all(isinstance(emb, torch.Tensor) for emb in embeddings.values())
        assert all(emb.shape == (320,) for emb in embeddings.values())
        
        processor.clear_cache()
    
    def test_batch_processing_with_remainder(self):
        """Test processing when sequence count is not multiple of batch size."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m", 
            batch_size=2,
            device="cpu"
        )
        
        sequences = {f"seq_{i}": "ACDEFG" for i in range(5)}
        
        embeddings = processor.embed_sequences(sequences, embedding_type="mean")
        
        assert len(embeddings) == 5
        assert all(emb.shape == (320,) for emb in embeddings.values())
        
        processor.clear_cache()
    
    def test_single_sequence_batching(self):
        """Test that single sequences work with batching."""
        processor = EmbeddingProcessor(
            name="test",
            model_name="esm2_t6_8m",
            batch_size=10,  # Larger than number of sequences
            device="cpu"
        )
        
        embedding = processor.embed_sequences("ACDEFG", embedding_type="mean")
        
        assert isinstance(embedding, torch.Tensor)
        assert embedding.shape == (320,)
        
        processor.clear_cache()


class TestEmbeddingProcessorFASTAProcessing:
    """Test FASTA file processing functionality."""
    
    def test_embed_from_fasta_success(self, tmp_path):
        """Test successful embedding generation from FASTA file."""
        # Create test FASTA file
        fasta_content = ">protein1\nACDEFGHIKL\n>protein2\nMNPQRSTVWY\n"
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)
        
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
        
        embeddings = processor.embed_from_fasta(str(fasta_file))
        
        assert len(embeddings) == 2
        assert 'protein1' in embeddings
        assert 'protein2' in embeddings
        assert all(emb.shape == (320,) for emb in embeddings.values())
        
        processor.clear_cache()
    
    def test_embed_from_fasta_with_save_dataset(self, tmp_path):
        """Test FASTA embedding with dataset saving."""
        # Create test FASTA file
        fasta_content = ">protein1\nACDEFGHIKL\n"
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(fasta_content)
        
        # Set up processor with temp data path
        with pytest.MonkeyPatch.context() as m:
            m.setenv("PROTOS_DATA_ROOT", str(tmp_path))
            processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
            processor.data_path = str(tmp_path)
            
            embeddings = processor.embed_from_fasta(
                str(fasta_file),
                save_dataset="test_dataset"
            )
            
            # Check dataset was saved
            dataset_path = tmp_path / "datasets" / "test_dataset"
            assert dataset_path.exists()
            assert (dataset_path / "embeddings.pt").exists()
            assert (dataset_path / "sequences.json").exists()
            assert (dataset_path / "metadata.json").exists()
            
            processor.clear_cache()


class TestEmbeddingProcessorErrorHandling:
    """Test error handling and edge cases."""
    
    def test_empty_sequence_handling(self):
        """Test handling of empty sequences."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
        
        # Empty string should still work (tokenizer will handle it)
        result = processor.embed_sequences("")
        assert isinstance(result, torch.Tensor)
        
        processor.clear_cache()
    
    def test_invalid_model_name(self):
        """Test initialization with invalid model name."""
        with pytest.raises(ValueError, match="Unknown model"):
            EmbeddingProcessor(name="test", model_name="invalid_model")
    
    def test_very_long_sequence(self):
        """Test handling of sequences longer than max_seq_length."""
        processor = EmbeddingProcessor(
            name="test", 
            model_name="esm2_t6_8m",
            max_seq_length=10,  # Very small for testing
            device="cpu"
        )
        
        long_sequence = "A" * 100  # Much longer than max_seq_length
        
        # Should truncate and still work
        result = processor.embed_sequences(long_sequence)
        assert isinstance(result, torch.Tensor)
        assert result.shape == (320,)
        
        processor.clear_cache()


class TestEmbeddingProcessorDatasetManagement:
    """Test dataset save/load functionality."""
    
    def test_save_and_load_embeddings(self, tmp_path):
        """Test saving and loading embeddings."""
        with pytest.MonkeyPatch.context() as m:
            m.setenv("PROTOS_DATA_ROOT", str(tmp_path))
            
            # Create and save embeddings
            processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
            processor.data_path = str(tmp_path)
            
            sequences = {'seq1': 'ACDEFG', 'seq2': 'KLMNPQ'}
            processor.embed_sequences(sequences, save_dataset="test_dataset")
            
            # Load embeddings
            loaded = processor.load_embeddings("test_dataset")
            
            assert len(loaded) == 2
            assert 'seq1' in loaded
            assert 'seq2' in loaded
            assert all(emb.shape == (320,) for emb in loaded.values())
            
            processor.clear_cache()
    
    def test_dataset_metadata(self, tmp_path):
        """Test that dataset metadata is properly saved."""
        with pytest.MonkeyPatch.context() as m:
            m.setenv("PROTOS_DATA_ROOT", str(tmp_path))
            
            processor = EmbeddingProcessor(
                name="test", 
                model_name="esm2_t6_8m",
                device="cpu"
            )
            processor.data_path = str(tmp_path)
            
            sequences = {'seq1': 'ACDEFG'}
            processor.embed_sequences(sequences, save_dataset="test_dataset")
            
            # Check metadata
            metadata_file = tmp_path / "datasets" / "test_dataset" / "metadata.json"
            with open(metadata_file, 'r') as f:
                metadata = json.load(f)
            
            assert metadata['model_name'] == "esm2_t6_8m"
            assert metadata['embedding_type'] == "mean"
            assert metadata['num_sequences'] == 1
            assert metadata['embedding_dim'] == 320
            
            processor.clear_cache()


class TestEmbeddingProcessorPerResidue:
    """Test per-residue embedding functionality."""
    
    def test_per_residue_embeddings(self):
        """Test generation of per-residue embeddings."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
        
        sequence = "ACDEFGHIKL"
        result = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # Should return 2D tensor
        assert isinstance(result, torch.Tensor)
        assert result.dim() == 2
        assert result.shape[1] == 320  # Embedding dimension
        
        # Length should include special tokens
        assert result.shape[0] > len(sequence)  # Due to CLS/EOS tokens
        
        processor.clear_cache()
    
    def test_per_residue_dict_input(self):
        """Test per-residue embeddings with dictionary input."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
        
        sequences = {
            "seq1": "ACDEF",
            "seq2": "GHIKL"
        }
        
        results = processor.embed_sequences(sequences, embedding_type="per_residue")
        
        assert len(results) == 2
        for seq_id, emb in results.items():
            assert emb.dim() == 2
            assert emb.shape[1] == 320
            assert emb.shape[0] > len(sequences[seq_id])
        
        processor.clear_cache()


class TestEmbeddingProcessorIntegration:
    """Integration tests with real models and longer sequences."""
    
    @pytest.mark.slow
    def test_real_protein_sequence(self):
        """Test with a real protein sequence."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
        
        # Bacteriorhodopsin sequence (truncated for testing)
        sequence = "MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITTLVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGTILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEVASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSRAIFGEAEAPEPSAGDGAAATSD"[:100]
        
        embedding = processor.embed_sequences(sequence)
        
        assert isinstance(embedding, torch.Tensor)
        assert embedding.shape == (320,)
        assert not torch.isnan(embedding).any()
        assert not torch.isinf(embedding).any()
        
        processor.clear_cache()
    
    @pytest.mark.slow  
    def test_multiple_embedding_types(self):
        """Test all embedding types with same sequence."""
        processor = EmbeddingProcessor(name="test", model_name="esm2_t6_8m", device="cpu")
        
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        
        mean_emb = processor.embed_sequences(sequence, embedding_type="mean")
        cls_emb = processor.embed_sequences(sequence, embedding_type="cls")
        per_res_emb = processor.embed_sequences(sequence, embedding_type="per_residue")
        
        # Check dimensions
        assert mean_emb.shape == (320,)
        assert cls_emb.shape == (320,)
        assert per_res_emb.shape[0] > len(sequence)
        assert per_res_emb.shape[1] == 320
        
        # Mean and CLS should be different
        assert not torch.allclose(mean_emb, cls_emb)
        
        processor.clear_cache()