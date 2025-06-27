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

from protos.core.base_processor import BaseProcessor

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
    from protos.processing.sequence.seq_processor import SeqProcessor
except ImportError:
    SeqProcessor = None

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
                 max_seq_length: int = 1022,
                 processor_data_dir: str = "embedding"):
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
        super().__init__(
            name=name,
            processor_data_dir=processor_data_dir
        )
        
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
                       save_dataset: Optional[str] = None) -> Union["torch.Tensor", Dict[str, "torch.Tensor"]]:
        """
        Generate embeddings for protein sequences.
        
        Args:
            sequences: Single sequence, list of sequences, or dict mapping IDs to sequences
            embedding_type: Type of embedding ('mean', 'cls', 'per_residue')
            save_dataset: If provided, save embeddings as a dataset
            
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
            self._save_embeddings_dataset(embeddings, seq_dict, save_dataset, embedding_type)
        
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
                                embedding_type: str):
        """Save embeddings as a dataset."""
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
        sequences = seq_proc.load_sequences_from_file(fasta_path)
        
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