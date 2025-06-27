"""
Tests for EmbeddingProcessor functionality.

These tests check the basic functionality of EmbeddingProcessor
without requiring torch/transformers to be installed.
"""

import pytest
import os
import json
from pathlib import Path

from protos.processing.embedding import EmbeddingProcessor


class TestEmbeddingProcessorBasic:
    """Test basic EmbeddingProcessor functionality without dependencies."""
    
    def test_initialization_without_dependencies(self):
        """Test that EmbeddingProcessor can be initialized without torch/transformers."""
        # Should work even without dependencies
        processor = EmbeddingProcessor(name="test_proc")
        
        assert processor.name == "test_proc"
        assert processor.model_name == "esm2_t12_35m"
        assert processor.batch_size == 8
        assert processor.max_seq_length == 1022
    
    def test_check_dependencies(self):
        """Test dependency checking."""
        processor = EmbeddingProcessor()
        deps = processor.check_dependencies()
        
        assert isinstance(deps, dict)
        assert 'torch' in deps
        assert 'transformers' in deps
        assert 'ready' in deps
        assert deps['ready'] == (deps['torch'] and deps['transformers'])
    
    def test_list_available_models(self):
        """Test listing available models."""
        processor = EmbeddingProcessor()
        models = processor.list_available_models()
        
        assert isinstance(models, dict)
        assert 'esm2_t6_8m' in models
        assert 'esm2_t33_650m' in models
        assert 'ankh_base' in models
        
        # Check model info structure
        for name, info in models.items():
            assert 'embedding_dim' in info
            assert 'description' in info
            assert isinstance(info['embedding_dim'], int)
    
    def test_get_embedding_dim(self):
        """Test getting embedding dimension."""
        processor = EmbeddingProcessor(model_name='esm2_t6_8m')
        assert processor.get_embedding_dim() == 320
        
        processor = EmbeddingProcessor(model_name='esm2_t33_650m')
        assert processor.get_embedding_dim() == 1280
    
    def test_invalid_model_name(self):
        """Test initialization with invalid model name."""
        with pytest.raises(ValueError, match="Unknown model"):
            EmbeddingProcessor(model_name="invalid_model")
    
    def test_metadata(self):
        """Test processor metadata."""
        processor = EmbeddingProcessor(
            name="test_meta",
            model_name="ankh_base",
            batch_size=16
        )
        
        assert processor.metadata['processor_type'] == 'EmbeddingProcessor'
        assert processor.metadata['name'] == 'test_meta'
        assert processor.metadata['model_name'] == 'ankh_base'
        assert processor.metadata['batch_size'] == 16
        assert processor.metadata['embedding_dim'] == 768
        assert 'dependencies_available' in processor.metadata
    
    def test_embed_sequences_without_dependencies(self):
        """Test that embed_sequences raises appropriate error without dependencies."""
        processor = EmbeddingProcessor()
        
        # Check if dependencies are available
        if not processor.dependencies_available:
            with pytest.raises(RuntimeError, match="missing dependencies"):
                processor.embed_sequences("ACDEFG")
    
    def test_model_property_without_dependencies(self):
        """Test that model property raises appropriate error without dependencies."""
        processor = EmbeddingProcessor()
        
        # Check if dependencies are available
        if not processor.dependencies_available:
            with pytest.raises(RuntimeError, match="missing dependencies"):
                _ = processor.model
    
    def test_clear_cache_without_dependencies(self):
        """Test that clear_cache works even without dependencies."""
        processor = EmbeddingProcessor()
        # Should not raise error even without dependencies
        processor.clear_cache()


class TestEmbeddingProcessorWithMocks:
    """Test EmbeddingProcessor with mocked dependencies."""
    
    @pytest.fixture
    def mock_torch(self, monkeypatch):
        """Mock torch module."""
        class MockTensor:
            def __init__(self, data):
                self.data = data
                self.shape = (len(data),) if isinstance(data, list) else data.shape
            
            def cpu(self):
                return self
            
            def sum(self, dim=None, keepdim=False):
                return MockTensor([1])
            
            def item(self):
                return 1
            
            def unsqueeze(self, dim):
                return self
            
            def __getitem__(self, idx):
                return self
            
            def __mul__(self, other):
                return self
        
        class MockTorch:
            @staticmethod
            def no_grad():
                class NoGradContext:
                    def __enter__(self):
                        return self
                    def __exit__(self, *args):
                        pass
                return NoGradContext()
            
            @staticmethod
            def sum(tensor, dim=None, keepdim=False):
                return MockTensor([1])
            
            class cuda:
                @staticmethod
                def is_available():
                    return False
                
                @staticmethod
                def empty_cache():
                    pass
            
            @staticmethod
            def save(obj, path):
                # Save as JSON instead
                import json
                with open(path.replace('.pt', '.json'), 'w') as f:
                    json.dump({"mock": "data"}, f)
            
            @staticmethod
            def load(path, map_location=None):
                # Load from JSON
                import json
                json_path = path.replace('.pt', '.json')
                if os.path.exists(json_path):
                    with open(json_path, 'r') as f:
                        return json.load(f)
                return {"mock": MockTensor([1, 2, 3])}
        
        monkeypatch.setattr('protos.processing.embedding.embedding_processor.torch', MockTorch())
        monkeypatch.setattr('protos.processing.embedding.embedding_processor._TORCH_AVAILABLE', True)
        return MockTensor
    
    @pytest.fixture
    def mock_transformers(self, monkeypatch):
        """Mock transformers module."""
        class MockTokenizer:
            def __call__(self, sequences, **kwargs):
                class TokenizerOutput:
                    def __init__(self):
                        self.data = {
                            'input_ids': [[1, 2, 3]],
                            'attention_mask': [[1, 1, 1]]
                        }
                    
                    def __getitem__(self, key):
                        return self.data[key]
                    
                    def to(self, device):
                        return self
                
                return TokenizerOutput()
        
        class MockModel:
            def __init__(self):
                pass
            
            def to(self, device):
                return self
            
            def eval(self):
                pass
            
            def __call__(self, **kwargs):
                class ModelOutput:
                    def __init__(self):
                        from protos.processing.embedding.embedding_processor import torch
                        self.last_hidden_state = torch.MockTensor([[1, 2, 3]])
                
                return ModelOutput()
        
        class MockAuto:
            @staticmethod
            def from_pretrained(name):
                return MockTokenizer() if 'tokenizer' in name else MockModel()
        
        monkeypatch.setattr('protos.processing.embedding.embedding_processor.AutoTokenizer', type('AutoTokenizer', (), {'from_pretrained': MockAuto.from_pretrained}))
        monkeypatch.setattr('protos.processing.embedding.embedding_processor.AutoModel', type('AutoModel', (), {'from_pretrained': MockAuto.from_pretrained}))
        monkeypatch.setattr('protos.processing.embedding.embedding_processor._TRANSFORMERS_AVAILABLE', True)
    
    def test_embed_sequences_with_mocks(self, mock_torch, mock_transformers):
        """Test embedding generation with mocked dependencies."""
        processor = EmbeddingProcessor(name="test_mock")
        
        # Test single sequence
        result = processor.embed_sequences("ACDEFG")
        assert result is not None
        
        # Test list of sequences
        results = processor.embed_sequences(["ACDEFG", "HIJKLM"])
        assert isinstance(results, dict)
        assert len(results) == 2
        
        # Test dict of sequences
        seq_dict = {"seq1": "ACDEFG", "seq2": "HIJKLM"}
        results = processor.embed_sequences(seq_dict)
        assert isinstance(results, dict)
        assert "seq1" in results
        assert "seq2" in results


# Skip these tests if dependencies are not available
try:
    import torch
    import transformers
    HAS_DEPENDENCIES = True
except ImportError:
    HAS_DEPENDENCIES = False


@pytest.mark.skipif(not HAS_DEPENDENCIES, reason="Requires torch and transformers")
class TestEmbeddingProcessorIntegration:
    """Integration tests that require actual dependencies."""
    
    def test_real_embedding_generation(self):
        """Test actual embedding generation with a small model."""
        processor = EmbeddingProcessor(
            name="test_real",
            model_name="esm2_t6_8m",  # Smallest model
            batch_size=2
        )
        
        # Test single sequence
        seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDL"
        embedding = processor.embed_sequences(seq, embedding_type="mean")
        
        assert isinstance(embedding, torch.Tensor)
        assert embedding.shape == (320,)  # esm2_t6_8m has 320 dim
        
        # Test multiple sequences
        seqs = {
            "seq1": "ACDEFGHIKLMNPQRSTVWY",
            "seq2": "WYTVRSQPNMLKIHGFEDCA"
        }
        embeddings = processor.embed_sequences(seqs, embedding_type="mean")
        
        assert len(embeddings) == 2
        assert all(isinstance(emb, torch.Tensor) for emb in embeddings.values())
        assert all(emb.shape == (320,) for emb in embeddings.values())
    
    def test_different_embedding_types(self):
        """Test different embedding extraction types."""
        processor = EmbeddingProcessor(model_name="esm2_t6_8m")
        seq = "ACDEFGHIKLMNPQRSTVWY"
        
        # Mean embedding
        mean_emb = processor.embed_sequences(seq, embedding_type="mean")
        assert mean_emb.shape == (320,)
        
        # CLS embedding
        cls_emb = processor.embed_sequences(seq, embedding_type="cls")
        assert cls_emb.shape == (320,)
        
        # Per-residue embedding
        residue_emb = processor.embed_sequences(seq, embedding_type="per_residue")
        assert residue_emb.shape[0] == len(seq)  # One per residue
        assert residue_emb.shape[1] == 320  # Embedding dim
    
    def test_save_and_load_embeddings(self, tmp_path):
        """Test saving and loading embeddings."""
        processor = EmbeddingProcessor(
            name="test_save",
            model_name="esm2_t6_8m"
        )
        
        # Set custom data path for testing
        processor.data_path = str(tmp_path)
        
        seqs = {
            "protein1": "ACDEFGHIKLMNPQRSTVWY",
            "protein2": "WYTVRSQPNMLKIHGFEDCA"
        }
        
        # Generate and save embeddings
        embeddings = processor.embed_sequences(
            seqs,
            embedding_type="mean",
            save_dataset="test_dataset"
        )
        
        # Check dataset was registered
        assert "test_dataset" in processor.list_datasets()
        
        # Load embeddings
        loaded = processor.load_embeddings("test_dataset")
        
        assert len(loaded) == 2
        assert "protein1" in loaded
        assert "protein2" in loaded
        
        # Check embeddings match
        for key in embeddings:
            assert torch.allclose(embeddings[key], loaded[key])