"""
Embedding Processor with BaseProcessor integration.

This module provides an EmbeddingProcessor class that extends BaseProcessor
to provide standardized data management capabilities for protein embeddings.

Note: This processor requires additional dependencies (torch, transformers) that
are not installed by default. Install them with:
    
For GPU support (recommended):
    pip install --no-build-isolation -e ".[gpu]"
    
For CPU-only support:
    pip install -e ".[embedding]"
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Union, Optional, Any, Literal
import warnings

from protos.io.core.base_processor import BaseProcessor
from protos.io.data_access import generate_entity_id

# Check for optional dependencies
_TORCH_AVAILABLE = False
_TRANSFORMERS_AVAILABLE = False
_DEPENDENCIES_ERROR = None

try:
    import torch
    _TORCH_AVAILABLE = True
except ImportError as e:
    _DEPENDENCIES_ERROR = "PyTorch not installed. Install with: pip install torch"

try:
    from transformers import AutoTokenizer, AutoModel
    _TRANSFORMERS_AVAILABLE = True
except ImportError as e:
    if _DEPENDENCIES_ERROR:
        _DEPENDENCIES_ERROR += "\nTransformers not installed. Install with: pip install transformers"
    else:
        _DEPENDENCIES_ERROR = "Transformers not installed. Install with: pip install transformers"

# Import sequence processor for loading sequences
try:
    from protos.processing.sequence import SequenceProcessor
except ImportError:
    SequenceProcessor = None

# Embedding type options
EmbeddingType = Literal["mean", "cls", "per_residue"]


class EmbeddingProcessor(BaseProcessor):
    """
    Processor for protein sequence embeddings.
    
    Handles generation, storage, and retrieval of protein embeddings
    using various transformer models (ESM-2, Ankh, etc.).
    
    Note: Requires torch and transformers libraries to be installed.
    """
    
    # Model registry with configurations
    MODEL_REGISTRY = {
        "esm2_t6_8m": {
            "hub_name": "facebook/esm2_t6_8M_UR50D",
            "embedding_dim": 320,
            "description": "ESM-2 tiny model (8M parameters)"
        },
        "esm2_t12_35m": {
            "hub_name": "facebook/esm2_t12_35M_UR50D", 
            "embedding_dim": 480,
            "description": "ESM-2 small model (35M parameters)"
        },
        "esm2_t30_150m": {
            "hub_name": "facebook/esm2_t30_150M_UR50D",
            "embedding_dim": 640,
            "description": "ESM-2 medium model (150M parameters)"
        },
        "esm2_t33_650m": {
            "hub_name": "facebook/esm2_t33_650M_UR50D",
            "embedding_dim": 1280,
            "description": "ESM-2 large model (650M parameters)"
        },
        "esm2_t36_3b": {
            "hub_name": "facebook/esm2_t36_3B_UR50D",
            "embedding_dim": 2560,
            "description": "ESM-2 extra large model (3B parameters)"
        },
        "esm2_t48_15b": {
            "hub_name": "facebook/esm2_t48_15B_UR50D",
            "embedding_dim": 5120,
            "description": "ESM-2 huge model (15B parameters)"
        },
        "ankh_base": {
            "hub_name": "ElnaggarLab/ankh-base",
            "embedding_dim": 768,
            "description": "Ankh base model"
        },
        "ankh_large": {
            "hub_name": "ElnaggarLab/ankh-large", 
            "embedding_dim": 1536,
            "description": "Ankh large model"
        }
    }
    
    def __init__(self,
                 name: str = "embedding_processor",
                 model_name: str = "esm2_t12_35m",
                 device: Optional[str] = None,
                 batch_size: int = 8,
                 max_seq_length: int = 1022):
        """
        Initialize the embedding processor.
        
        Args:
            name: Processor instance name
            model_name: Name of the model from MODEL_REGISTRY
            device: Device to use ('cuda', 'cpu', or None for auto)
            batch_size: Number of sequences to process at once
            max_seq_length: Maximum sequence length for tokenization
            processor_data_dir: Directory for embedding data
        """
        super().__init__(name=name)
        
        # Check dependencies
        self._check_dependencies()
        
        # Model configuration
        if model_name not in self.MODEL_REGISTRY:
            raise ValueError(f"Unknown model: {model_name}. Available: {list(self.MODEL_REGISTRY.keys())}")
        
        self.model_name = model_name
        self.model_config = self.MODEL_REGISTRY[model_name]
        self.batch_size = batch_size
        self.max_seq_length = max_seq_length
        
        # Set device only if torch is available
        if _TORCH_AVAILABLE:
            self.device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
        else:
            self.device = 'cpu'
        
        # Model and tokenizer (loaded on demand)
        self._model = None
        self._tokenizer = None
        
        # Update metadata
        self.metadata.update({
            "model_name": model_name,
            "device": self.device,
            "batch_size": batch_size,
            "max_seq_length": max_seq_length,
            "embedding_dim": self.model_config["embedding_dim"],
            "dependencies_available": _TORCH_AVAILABLE and _TRANSFORMERS_AVAILABLE
        })
        
        self.logger.info(f"Initialized EmbeddingProcessor with model {model_name}")
        if not (_TORCH_AVAILABLE and _TRANSFORMERS_AVAILABLE):
            self.logger.warning("Dependencies not available. Install with: pip install torch transformers")
    
    def _check_dependencies(self):
        """Check if required dependencies are available."""
        if not (_TORCH_AVAILABLE and _TRANSFORMERS_AVAILABLE):
            self.logger.warning(
                "EmbeddingProcessor requires additional dependencies.\n"
                "Install with one of:\n"
                "  pip install --no-build-isolation -e '.[gpu]'  # GPU support\n"
                "  pip install -e '.[embedding]'                  # CPU-only\n"
                "  pip install torch transformers                 # Manual install"
            )
    
    @property
    def dependencies_available(self) -> bool:
        """Check if all required dependencies are available."""
        return _TORCH_AVAILABLE and _TRANSFORMERS_AVAILABLE
    
    @property
    def model(self):
        """Lazy load the model."""
        if not self.dependencies_available:
            raise RuntimeError(
                "Cannot load model - missing dependencies.\n"
                f"{_DEPENDENCIES_ERROR}\n"
                "Install with: pip install -e '.[gpu]' (recommended) or pip install -e '.[embedding]'"
            )
        if self._model is None:
            self._load_model()
        return self._model
    
    @property 
    def tokenizer(self):
        """Lazy load the tokenizer."""
        if not self.dependencies_available:
            raise RuntimeError(
                "Cannot load tokenizer - missing dependencies.\n"
                f"{_DEPENDENCIES_ERROR}\n"
                "Install with: pip install -e '.[gpu]' (recommended) or pip install -e '.[embedding]'"
            )
        if self._tokenizer is None:
            self._load_model()
        return self._tokenizer
    
    def _load_model(self):
        """Load the model and tokenizer."""
        if not self.dependencies_available:
            raise RuntimeError("Cannot load model without torch and transformers installed.")
        
        self.logger.info(f"Loading model {self.model_name} from {self.model_config['hub_name']}")
        
        self._tokenizer = AutoTokenizer.from_pretrained(self.model_config['hub_name'])
        self._model = AutoModel.from_pretrained(self.model_config['hub_name'])
        self._model.to(self.device)
        self._model.eval()
        
        self.logger.info(f"Model loaded successfully on {self.device}")
    
    def embed_sequences(self,
                       sequences: Union[str, List[str], Dict[str, str]],
                       embedding_type: EmbeddingType = "mean",
                       save_dataset: Optional[str] = None,
                       register_entities: bool = True) -> Union["torch.Tensor", Dict[str, "torch.Tensor"]]:
        """
        Generate embeddings for protein sequences with entity support.
        
        Args:
            sequences: Single sequence, list of sequences, or dict mapping IDs to sequences
            embedding_type: Type of embedding ('mean', 'cls', 'per_residue')
            save_dataset: If provided, save embeddings as a dataset
            register_entities: Whether to register embeddings as entities
            
        Returns:
            Single tensor or dict of tensors depending on input
            
        Raises:
            RuntimeError: If dependencies are not installed
        """
        if not self.dependencies_available:
            raise RuntimeError(
                "Cannot generate embeddings - missing dependencies.\n"
                f"{_DEPENDENCIES_ERROR}\n"
                "Install with: pip install -e '.[gpu]' (recommended) or pip install -e '.[embedding]'"
            )
        
        # Standardize input
        is_single = isinstance(sequences, str)
        if is_single:
            seq_dict = {"seq_0": sequences}
        elif isinstance(sequences, list):
            seq_dict = {f"seq_{i}": seq for i, seq in enumerate(sequences)}
        else:
            seq_dict = sequences
        
        # Generate embeddings
        embeddings = self._generate_embeddings(seq_dict, embedding_type)
        
        # Save if requested
        if save_dataset:
            self._save_embeddings_dataset(embeddings, seq_dict, save_dataset, embedding_type, 
                                        register_entities=register_entities)
        
        # Register entities if requested (even without saving dataset)
        elif register_entities:
            self._register_embedding_entities(embeddings, seq_dict, embedding_type)
        
        # Return appropriate format
        return embeddings["seq_0"] if is_single else embeddings
    
    def _generate_embeddings(self,
                           seq_dict: Dict[str, str],
                           embedding_type: EmbeddingType) -> Dict[str, "torch.Tensor"]:
        """Generate embeddings for a dictionary of sequences."""
        results = {}
        seq_ids = list(seq_dict.keys())
        sequences = list(seq_dict.values())
        
        self.logger.info(f"Generating {embedding_type} embeddings for {len(sequences)} sequences")
        
        with torch.no_grad():
            for i in range(0, len(sequences), self.batch_size):
                batch_ids = seq_ids[i:i + self.batch_size]
                batch_seqs = sequences[i:i + self.batch_size]
                
                # Tokenize
                inputs = self.tokenizer(
                    batch_seqs,
                    padding="longest",
                    truncation=True,
                    max_length=self.max_seq_length,
                    return_tensors="pt"
                ).to(self.device)
                
                # Get model output
                outputs = self.model(**inputs)
                
                # Extract embeddings
                batch_embeddings = self._extract_embeddings(
                    outputs.last_hidden_state,
                    inputs["attention_mask"],
                    embedding_type
                )
                
                # Store results
                for j, seq_id in enumerate(batch_ids):
                    if embedding_type == "per_residue":
                        # Trim padding for per-residue embeddings
                        seq_len = inputs["attention_mask"][j].sum().item()
                        results[seq_id] = batch_embeddings[j, :seq_len].cpu()
                    else:
                        results[seq_id] = batch_embeddings[j].cpu()
        
        return results
    
    def _extract_embeddings(self,
                          hidden_states: "torch.Tensor",
                          attention_mask: "torch.Tensor",
                          embedding_type: EmbeddingType) -> "torch.Tensor":
        """Extract embeddings from model output."""
        if embedding_type == "mean":
            # Mean pooling with attention mask
            masked = hidden_states * attention_mask.unsqueeze(-1)
            summed = torch.sum(masked, dim=1)
            lengths = torch.sum(attention_mask, dim=1, keepdim=True)
            return summed / lengths
        elif embedding_type == "cls":
            # CLS token (first position)
            return hidden_states[:, 0, :]
        elif embedding_type == "per_residue":
            # Return all hidden states
            return hidden_states
        else:
            raise ValueError(f"Unknown embedding type: {embedding_type}")
    
    def _save_embeddings_dataset(self,
                                embeddings: Dict[str, "torch.Tensor"],
                                sequences: Dict[str, str],
                                dataset_name: str,
                                embedding_type: str,
                                register_entities: bool = True):
        """Save embeddings as a dataset with entity support."""
        # Create dataset directory
        dataset_path = os.path.join(self.data_path, "datasets", dataset_name)
        os.makedirs(dataset_path, exist_ok=True)
        
        # Save embeddings
        embeddings_file = os.path.join(dataset_path, "embeddings.pt")
        torch.save(embeddings, embeddings_file)
        
        # Save sequences
        sequences_file = os.path.join(dataset_path, "sequences.json")
        with open(sequences_file, 'w') as f:
            json.dump(sequences, f, indent=2)
        
        # Save metadata
        metadata = {
            "model_name": self.model_name,
            "embedding_type": embedding_type,
            "num_sequences": len(sequences),
            "embedding_dim": self.model_config["embedding_dim"],
            "created_at": pd.Timestamp.now().isoformat()
        }
        metadata_file = os.path.join(dataset_path, "metadata.json")
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        # Register dataset
        self._register_dataset(dataset_name, {
            "type": "embeddings",
            "model": self.model_name,
            "embedding_type": embedding_type,
            "num_sequences": len(sequences),
            "path": dataset_path
        })
        
        # Register entities if requested
        if register_entities:
            self._register_embedding_entities(embeddings, sequences, embedding_type, dataset_name)
        
        self.logger.info(f"Saved embeddings dataset: {dataset_name}")
    
    def load_embeddings(self, dataset_name: str) -> Dict[str, "torch.Tensor"]:
        """
        Load embeddings from a saved dataset.
        
        Note: Requires torch to be installed to load the tensors.
        """
        if not _TORCH_AVAILABLE:
            raise RuntimeError(
                "Cannot load embeddings - PyTorch not installed.\n"
                "Install with: pip install torch"
            )
        
        dataset_info = self.get_dataset_info(dataset_name)
        if not dataset_info:
            raise ValueError(f"Dataset not found: {dataset_name}")
        
        dataset_path = dataset_info["path"]
        embeddings_file = os.path.join(dataset_path, "embeddings.pt")
        
        if not os.path.exists(embeddings_file):
            raise FileNotFoundError(f"Embeddings file not found: {embeddings_file}")
        
        embeddings = torch.load(embeddings_file, map_location='cpu')
        self.logger.info(f"Loaded {len(embeddings)} embeddings from {dataset_name}")
        
        return embeddings
    
    def embed_from_fasta(self,
                        fasta_path: str,
                        embedding_type: EmbeddingType = "mean",
                        save_dataset: Optional[str] = None) -> Dict[str, "torch.Tensor"]:
        """
        Generate embeddings from a FASTA file.
        
        Args:
            fasta_path: Path to FASTA file
            embedding_type: Type of embedding to generate
            save_dataset: Optional dataset name to save results
            
        Returns:
            Dictionary mapping sequence IDs to embeddings
        """
        if SeqProcessor is None:
            raise ImportError("SeqProcessor not available. Cannot load FASTA files.")
        
        # Load sequences using SeqProcessor
        seq_proc = SeqProcessor(name="temp_seq_loader")
        sequences = seq_proc.load_sequences(fasta_path)
        
        # Generate embeddings
        return self.embed_sequences(sequences, embedding_type, save_dataset)
    
    def get_embedding_dim(self) -> int:
        """Get the embedding dimension for the current model."""
        return self.model_config["embedding_dim"]
    
    def list_available_models(self) -> Dict[str, Dict[str, Any]]:
        """List all available models with their configurations."""
        return {
            name: {
                "embedding_dim": config["embedding_dim"],
                "description": config.get("description", "")
            }
            for name, config in self.MODEL_REGISTRY.items()
        }
    
    def check_dependencies(self) -> Dict[str, bool]:
        """Check which dependencies are available."""
        return {
            "torch": _TORCH_AVAILABLE,
            "transformers": _TRANSFORMERS_AVAILABLE,
            "ready": self.dependencies_available
        }
    
    def get_residue_embeddings(self, embeddings: "torch.Tensor", include_special_tokens: bool = False) -> "torch.Tensor":
        """
        Extract residue embeddings from per-residue output.
        
        For ESM-2 and similar models, per-residue embeddings include special tokens:
        - Index 0: CLS token
        - Index 1 to N: Residue embeddings
        - Index N+1: EOS token
        
        Args:
            embeddings: Per-residue embeddings from embed_sequences with embedding_type="per_residue"
            include_special_tokens: If False (default), return only residue embeddings.
                                  If True, return all embeddings including CLS and EOS.
                                  
        Returns:
            Tensor of shape (num_residues, embedding_dim) if include_special_tokens=False
            Tensor of shape (num_residues+2, embedding_dim) if include_special_tokens=True
            
        Example:
            >>> processor = EmbeddingProcessor(name="test")
            >>> sequence = "ACDEFGHIKL"
            >>> full_embeddings = processor.embed_sequences(sequence, embedding_type="per_residue")
            >>> # full_embeddings has shape (12, 320) - includes CLS and EOS
            >>> 
            >>> residue_embeddings = processor.get_residue_embeddings(full_embeddings)
            >>> # residue_embeddings has shape (10, 320) - one per residue
            >>> 
            >>> # Map to positions
            >>> for i, aa in enumerate(sequence):
            >>>     position = i + 1  # 1-indexed
            >>>     embedding = residue_embeddings[i]
            >>>     print(f"Position {position}: {aa} -> embedding shape {embedding.shape}")
        """
        if include_special_tokens:
            return embeddings
        
        # Skip first (CLS) and last (EOS) tokens
        return embeddings[1:-1]
    
    def clear_cache(self):
        """Clear model from memory."""
        if self._model is not None:
            del self._model
            self._model = None
        if self._tokenizer is not None:
            del self._tokenizer
            self._tokenizer = None
        if _TORCH_AVAILABLE and torch.cuda.is_available():
            torch.cuda.empty_cache()
        self.logger.info("Cleared model cache")
    
    # Entity support methods
    def _register_embedding_entities(self, 
                                   embeddings: Dict[str, "torch.Tensor"],
                                   sequences: Dict[str, str],
                                   embedding_type: str,
                                   dataset_name: Optional[str] = None):
        """Register embeddings as entities."""
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            
            for seq_id in embeddings.keys():
                # Use same entity ID as the source sequence
                entity_id = generate_entity_id(seq_id)
                
                # Get embedding info
                embedding = embeddings[seq_id]
                if _TORCH_AVAILABLE:
                    embedding_shape = list(embedding.shape)
                else:
                    embedding_shape = []
                
                # Register entity
                global_registry.entity_registry.register_entity(
                    entity_id=entity_id,
                    entity_type="embedding",
                    original_id=seq_id,
                    file_path=None,  # Embeddings are in datasets/memory
                    metadata={
                        "model": self.model_name,
                        "embedding_type": embedding_type,
                        "shape": embedding_shape,
                        "dataset": dataset_name,
                        "sequence_length": len(sequences.get(seq_id, ""))
                    },
                    datasets=[dataset_name] if dataset_name else []
                )
            
            self.logger.info(f"Registered {len(embeddings)} embedding entities")
        except Exception as e:
            self.logger.warning(f"Could not register embedding entities: {e}")
    
    def load_embedding_entity(self, identifier: str) -> Optional["torch.Tensor"]:
        """
        Load a single embedding entity.
        
        Args:
            identifier: Sequence identifier (name or entity hash)
            
        Returns:
            Embedding tensor or None if not found
        """
        if not _TORCH_AVAILABLE:
            raise RuntimeError("Cannot load embeddings without torch installed")
        
        # Resolve identifier
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            entity_id = global_registry.entity_registry.resolve_identifier(identifier, format_type="embedding")
            
            # Get original ID
            original_id = global_registry.entity_registry.get_original_id(entity_id)
            if not original_id:
                original_id = identifier
        except:
            original_id = identifier
        
        # Check if we have embeddings loaded
        if hasattr(self, '_loaded_embeddings'):
            for dataset_embeddings in self._loaded_embeddings.values():
                if original_id in dataset_embeddings:
                    return dataset_embeddings[original_id]
        
        # Try to find which dataset contains this embedding
        try:
            global_registry = GlobalRegistry()
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info and 'embedding' in entity_info.get('formats', {}):
                dataset = entity_info['formats']['embedding']['metadata'].get('dataset')
                if dataset:
                    # Load the dataset
                    embeddings = self.load_embeddings(dataset)
                    if original_id in embeddings:
                        return embeddings[original_id]
        except:
            pass
        
        self.logger.warning(f"Embedding entity not found: {identifier}")
        return None
    
    def list_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all embedding entities.
        
        Args:
            dataset: Optional dataset to filter by
            
        Returns:
            List of sequence IDs (not hash IDs!)
        """
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            entity_ids = global_registry.entity_registry.list_entities(
                format_type="embedding",
                dataset=dataset
            )
            
            # Convert to original IDs
            original_ids = []
            for entity_id in entity_ids:
                original_id = global_registry.entity_registry.get_original_id(entity_id)
                if original_id:
                    original_ids.append(original_id)
            return original_ids
        except:
            return []
    
    def list_embedding_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all embedding entities (backward compatibility).
        
        Deprecated: Use list_entities() instead.
        
        Args:
            dataset: Optional dataset to filter by
            
        Returns:
            List of sequence IDs
        """
        return self.list_entities(dataset=dataset)
    
    def list_datasets(self) -> List[Dict[str, Any]]:
        """
        List available embedding datasets.
        
        Embedding datasets are collections of embeddings, typically organized by:
        - Model type (ESM-2, ProtTrans, etc.)
        - Embedding purpose (structural, functional, etc.)
        - Source dataset (sequences from specific experiments)
        
        Returns:
            List of dataset information dictionaries
        """
        datasets = []
        
        # Check dataset manager first
        if self.dataset_manager:
            return self.dataset_manager.list_datasets()
        
        # Otherwise, scan for embedding collections
        embeddings_dir = self.get_subdirectory_path('embeddings')
        
        # Look for organized embedding directories
        if embeddings_dir.exists():
            # Check for model-based organization
            for model_dir in embeddings_dir.iterdir():
                if model_dir.is_dir():
                    embedding_files = list(model_dir.glob("*.pkl")) + list(model_dir.glob("*.pt"))
                    if embedding_files:
                        datasets.append({
                            'id': model_dir.name,
                            'type': 'embedding_dataset',
                            'model': model_dir.name,
                            'embedding_count': len(embedding_files),
                            'path': str(model_dir)
                        })
        
        # Look for dataset JSON definitions
        dataset_file = embeddings_dir / 'datasets.json'
        if dataset_file.exists():
            try:
                with open(dataset_file, 'r') as f:
                    dataset_defs = json.load(f)
                    for dataset_id, dataset_info in dataset_defs.items():
                        datasets.append({
                            'id': dataset_id,
                            'type': 'embedding_dataset',
                            **dataset_info
                        })
            except:
                pass
        
        return datasets
    
    def load_entity(self, entity_id: str, format_type: str = "embedding") -> Optional[Any]:
        """
        Load an embedding entity by ID.
        
        Args:
            entity_id: Human-readable entity ID (e.g., 'P12345')
            format_type: Format type (default: 'embedding')
            
        Returns:
            Embedding tensor/array or None if not found
        """
        if format_type != "embedding":
            return None
            
        # Look for embedding file
        embeddings_dir = self.get_subdirectory_path('embeddings')
        
        # Try different file formats
        for ext in ['.npy', '.pkl', '.pt']:
            embedding_file = embeddings_dir / f"{entity_id}{ext}"
            if embedding_file.exists():
                if ext == '.npy':
                    return np.load(embedding_file)
                elif ext == '.pkl':
                    import pickle
                    with open(embedding_file, 'rb') as f:
                        return pickle.load(f)
                elif ext == '.pt' and _TORCH_AVAILABLE:
                    return torch.load(embedding_file)
        
        return None
    
    def save_entity(self, entity_id: str, data: Any, format_type: str = "embedding") -> bool:
        """
        Save an embedding entity.
        
        Args:
            entity_id: Human-readable entity ID
            data: Embedding data (numpy array or torch tensor)
            format_type: Format type (default: 'embedding')
            
        Returns:
            True if saved successfully
        """
        if format_type != "embedding":
            return False
            
        # Ensure embeddings directory exists
        embeddings_dir = self.get_subdirectory_path('embeddings')
        embeddings_dir.mkdir(parents=True, exist_ok=True)
        
        # Convert to numpy array if torch tensor
        if _TORCH_AVAILABLE and torch.is_tensor(data):
            data = data.cpu().numpy()
        
        # Save as numpy array
        embedding_file = embeddings_dir / f"{entity_id}.npy"
        np.save(embedding_file, data)
        
        # Register entity if registry available
        if self.entity_registry:
            self.entity_registry.register_entity(
                original_id=entity_id,
                format_type="embedding",
                file_path=embedding_file,
                metadata={
                    "shape": data.shape,
                    "dtype": str(data.dtype),
                    "model": self.model_name
                }
            )
        
        return True