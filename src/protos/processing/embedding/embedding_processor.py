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

import json
import os
import pickle
import warnings
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.io.formats.fasta_utils import read_fasta

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

# Embedding type options
EmbeddingType = Literal["mean", "sum", "cls", "per_residue"]


class EmbeddingProcessor(BaseProcessor):
    """
    Processor for protein sequence embeddings.
    
    Handles generation, storage, and retrieval of protein embeddings
    using various transformer models (ESM-2, Ankh, etc.).
    
    Note: Requires torch and transformers libraries to be installed.
    """

    processor_type = "embedding"
    
    # Model registry with configurations
    MODEL_REGISTRY = {
        "esm2_t6_8m": {
            "hub_name": "facebook/esm2_t6_8M_UR50D",
            "embedding_dim": 320,
            "description": "ESM-2 tiny models (8M parameters)"
        },
        "esm2_t12_35m": {
            "hub_name": "facebook/esm2_t12_35M_UR50D", 
            "embedding_dim": 480,
            "description": "ESM-2 small models (35M parameters)"
        },
        "esm2_t30_150m": {
            "hub_name": "facebook/esm2_t30_150M_UR50D",
            "embedding_dim": 640,
            "description": "ESM-2 medium models (150M parameters)"
        },
        "esm2_t33_650m": {
            "hub_name": "facebook/esm2_t33_650M_UR50D",
            "embedding_dim": 1280,
            "description": "ESM-2 large models (650M parameters)"
        },
        "esm2_t36_3b": {
            "hub_name": "facebook/esm2_t36_3B_UR50D",
            "embedding_dim": 2560,
            "description": "ESM-2 extra large models (3B parameters)"
        },
        "esm2_t48_15b": {
            "hub_name": "facebook/esm2_t48_15B_UR50D",
            "embedding_dim": 5120,
            "description": "ESM-2 huge models (15B parameters)"
        },
        "ankh_base": {
            "hub_name": "ElnaggarLab/ankh-base",
            "embedding_dim": 768,
            "description": "Ankh base models"
        },
        "ankh_large": {
            "hub_name": "ElnaggarLab/ankh-large", 
            "embedding_dim": 1536,
            "description": "Ankh large models"
        }
    }

    @classmethod
    def available_models(cls) -> Dict[str, Dict[str, Any]]:
        """Return the registry of available transformer models."""

        return {
            name: {
                "embedding_dim": config["embedding_dim"],
                "description": config.get("description", "")
            }
            for name, config in cls.MODEL_REGISTRY.items()
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
            model_name: Name of the models from MODEL_REGISTRY
            device: Device to use ('cuda', 'cpu', or None for auto)
            batch_size: Number of sequences to process at once
            max_seq_length: Maximum sequence length for tokenization
            processor_data_dir: Directory for embedding data
        """
        super().__init__(name=name)

        # Ensure processor directories exist
        self.embeddings_dir = self.get_subdirectory_path('embeddings_dir')
        self.datasets_dir = self.get_subdirectory_path('datasets_dir')
        self.embeddings_dir.mkdir(parents=True, exist_ok=True)
        self.datasets_dir.mkdir(parents=True, exist_ok=True)

        # Cache for loaded datasets
        self._loaded_embeddings: Dict[str, Dict[str, Any]] = {}
        self.metadata: Dict[str, Any] = {}
        
        # Check dependencies
        self._check_dependencies()
        
        # Model configuration
        if model_name not in self.MODEL_REGISTRY:
            raise ValueError(f"Unknown models: {model_name}. Available: {list(self.MODEL_REGISTRY.keys())}")
        
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
        self._is_encoder_decoder = False
        self._encoder = None
        
        # Update metadata
        self.metadata.update({
            "model_name": model_name,
            "device": self.device,
            "batch_size": batch_size,
            "max_seq_length": max_seq_length,
            "embedding_dim": self.model_config["embedding_dim"],
            "dependencies_available": _TORCH_AVAILABLE and _TRANSFORMERS_AVAILABLE,
            "is_encoder_decoder": self._is_encoder_decoder,
        })
        
        self.logger.info(f"Initialized EmbeddingProcessor with models {model_name}")
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
        """Lazy load the models."""
        if not self.dependencies_available:
            raise RuntimeError(
                "Cannot load models - missing dependencies.\n"
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
        """Load the models and tokenizer."""
        if not self.dependencies_available:
            raise RuntimeError("Cannot load models without torch and transformers installed.")
        
        self.logger.info(f"Loading models {self.model_name} from {self.model_config['hub_name']}")
        
        self._tokenizer = AutoTokenizer.from_pretrained(self.model_config['hub_name'])
        # Some tokenizers (e.g. T5) ship without an explicit pad token; default to EOS if needed
        if getattr(self._tokenizer, "pad_token", None) is None and getattr(self._tokenizer, "eos_token", None) is not None:
            self._tokenizer.pad_token = self._tokenizer.eos_token

        self._model = AutoModel.from_pretrained(self.model_config['hub_name'])
        self._model.to(self.device)
        self._model.eval()

        # Track whether we are working with an encoder-decoder architecture
        self._is_encoder_decoder = bool(getattr(self._model.config, "is_encoder_decoder", False))
        self._encoder = None
        if self._is_encoder_decoder:
            # Prefer the public getter when available so tied weights stay in sync
            if hasattr(self._model, "get_encoder"):
                self._encoder = self._model.get_encoder()
            elif hasattr(self._model, "encoder"):
                self._encoder = self._model.encoder
            else:
                self.logger.warning(
                    "Model %s reported encoder-decoder architecture but exposes no encoder; falling back to full forward pass",
                    self.model_name,
                )
                self._is_encoder_decoder = False
        self.metadata["is_encoder_decoder"] = self._is_encoder_decoder

        self.logger.info(f"Model loaded successfully on {self.device}")

    def embed_sequences(self,
                       sequences: Union[str, List[str], Dict[str, str]],
                       embedding_type: EmbeddingType = "mean",
                       save_dataset: Optional[str] = None,
                       register_entities: bool = True,
                       stream_to_disk: bool = False) -> Union["torch.Tensor", Dict[str, "torch.Tensor"]]:
        """
        Generate embeddings for protein sequences with entity support.
        
        Args:
            sequences: Single sequence, list of sequences, or dict mapping IDs to sequences
            embedding_type: Type of embedding ('mean', 'sum', 'cls', 'per_residue')
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
        
        # Persist as dataset if requested
        if save_dataset:
            self._save_embeddings_dataset(
                embeddings,
                seq_dict,
                save_dataset,
                embedding_type,
                register_entities=True,
            )
        elif register_entities:
            self._persist_embeddings(
                embeddings,
                seq_dict,
                embedding_type,
                dataset_name=None,
            )
        
        # Return appropriate format
        return embeddings["seq_0"] if is_single else embeddings
    
    def _generate_embeddings(
        self,
        seq_dict: Dict[str, str],
        embedding_type: EmbeddingType,
        stream_to_disk: bool = False,
        dataset_name: Optional[str] = None,
    ) -> Dict[str, "torch.Tensor"]:
        """Generate embeddings for a dictionary of sequences.

        Args:
            seq_dict: Mapping of sequence IDs to sequences.
            embedding_type: Type of embedding to generate.
            stream_to_disk: If True, write each batch immediately to disk to save memory.
            dataset_name: Dataset name for organizing streamed output files.

        Returns:
            Dictionary mapping sequence IDs to embedding tensors (empty if stream_to_disk=True).
        """
        from tqdm import tqdm

        results = {}
        seq_ids = list(seq_dict.keys())
        sequences = list(seq_dict.values())

        self.logger.info(f"Generating {embedding_type} embeddings for {len(sequences)} sequences")

        num_batches = (len(sequences) + self.batch_size - 1) // self.batch_size

        # Setup streaming directory if needed
        if stream_to_disk:
            target_dir = self.embeddings_dir / self.model_name / (dataset_name or "_stream")
            target_dir.mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Streaming embeddings to: {target_dir}")

        with torch.no_grad():
            for i in tqdm(
                range(0, len(sequences), self.batch_size),
                total=num_batches,
                desc=f"Embedding ({embedding_type})",
                unit="batch",
            ):
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

                # Forward pass through encoder (skip decoder for seq2seq models)
                hidden_states = self._forward_encoder(inputs)

                attention_mask = inputs.get("attention_mask")
                if attention_mask is None:
                    attention_mask = torch.ones(
                        hidden_states.shape[:2],
                        dtype=torch.long,
                        device=hidden_states.device,
                    )

                # Extract embeddings
                batch_embeddings = self._extract_embeddings(
                    hidden_states,
                    attention_mask,
                    embedding_type
                )

                # Process each sequence in the batch
                for j, seq_id in enumerate(batch_ids):
                    if embedding_type == "per_residue":
                        # Trim padding for per-residue embeddings
                        seq_len = attention_mask[j].sum().item()
                        emb_tensor = batch_embeddings[j, :seq_len].cpu()
                    else:
                        emb_tensor = batch_embeddings[j].cpu()

                    if stream_to_disk:
                        # Write immediately to disk
                        entity_name = self._build_entity_name(seq_id, embedding_type, dataset_name)
                        file_path = target_dir / f"{entity_name}.pkl"
                        array = self._to_numpy(emb_tensor)

                        payload = {
                            "embedding": array,
                            "models": self.model_name,
                            "embedding_type": embedding_type,
                            "source_sequence": seq_id,
                            "dataset": dataset_name,
                        }
                        with open(file_path, "wb") as handle:
                            pickle.dump(payload, handle)

                        # Register entity
                        relative_path = str(file_path.relative_to(self.paths.data_root))
                        metadata = {
                            "models": self.model_name,
                            "embedding_type": embedding_type,
                            "shape": list(array.shape),
                            "dtype": str(array.dtype),
                            "source_sequence": seq_id,
                        }
                        if dataset_name:
                            metadata["dataset"] = dataset_name

                        self.entity_registry.register_entity(
                            name=entity_name,
                            format_type=self.processor_type,
                            file_path=relative_path,
                            metadata=metadata,
                        )
                        # Don't keep in memory
                        del emb_tensor
                    else:
                        results[seq_id] = emb_tensor

                # Clear GPU memory after each batch
                if stream_to_disk:
                    del batch_embeddings, hidden_states, inputs
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()

        return results

    def _forward_encoder(self, inputs: "torch.Tensor") -> "torch.Tensor":
        """Run the transformer encoder and return hidden states."""

        model = self.model  # Ensure lazy loading takes place

        if self._is_encoder_decoder and self._encoder is not None:
            encoder_kwargs = self._prepare_encoder_kwargs(inputs)
            encoder_outputs = self._encoder(**encoder_kwargs)
            if isinstance(encoder_outputs, tuple):
                hidden_states = encoder_outputs[0]
            else:
                hidden_states = encoder_outputs.last_hidden_state
        else:
            outputs = model(**inputs)
            hidden_states = outputs.last_hidden_state if hasattr(outputs, "last_hidden_state") else outputs[0]

        return hidden_states

    def _prepare_encoder_kwargs(self, inputs: "torch.Tensor") -> Dict[str, Any]:
        """Filter tokenized inputs down to what encoders expect."""

        allowed_keys = {
            "input_ids",
            "attention_mask",
            "inputs_embeds",
            "position_ids",
            "head_mask",
        }
        encoder_kwargs: Dict[str, Any] = {}
        for key in allowed_keys:
            if key in inputs:
                encoder_kwargs[key] = inputs[key]

        if "return_dict" not in encoder_kwargs:
            encoder_kwargs["return_dict"] = getattr(self._model.config, "use_return_dict", True)

        return encoder_kwargs
    
    def _extract_embeddings(self,
                          hidden_states: "torch.Tensor",
                          attention_mask: "torch.Tensor",
                          embedding_type: EmbeddingType) -> "torch.Tensor":
        """Extract embeddings from models output."""
        if embedding_type == "mean":
            # Mean pooling with attention mask
            masked = hidden_states * attention_mask.unsqueeze(-1)
            summed = torch.sum(masked, dim=1)
            lengths = torch.sum(attention_mask, dim=1, keepdim=True)
            return summed / lengths
        if embedding_type == "sum":
            masked = hidden_states * attention_mask.unsqueeze(-1)
            return torch.sum(masked, dim=1)
        elif embedding_type == "cls":
            # CLS token (first position)
            return hidden_states[:, 0, :]
        elif embedding_type == "per_residue":
            # Return all hidden states
            return hidden_states
        else:
            raise ValueError(f"Unknown embedding type: {embedding_type}")
    
    def _to_numpy(self, tensor: Any) -> np.ndarray:
        """Convert an embedding tensor to a NumPy array."""

        if _TORCH_AVAILABLE and torch.is_tensor(tensor):
            return tensor.detach().cpu().numpy()
        if isinstance(tensor, np.ndarray):
            return tensor
        return np.asarray(tensor)

    def _build_entity_name(
        self,
        sequence_id: str,
        embedding_type: str,
        dataset_name: Optional[str] = None,
    ) -> str:
        """Construct a unique, human-readable embedding entity name."""

        parts = [sequence_id, self.model_name, embedding_type]
        if dataset_name:
            parts.append(dataset_name)
        raw_name = "__".join(parts)
        return self._sanitize_filename(raw_name.replace(os.sep, "_"))

    def _persist_embeddings(
        self,
        embeddings: Dict[str, "torch.Tensor"],
        sequences: Dict[str, str],
        embedding_type: str,
        dataset_name: Optional[str] = None,
        entity_name_map: Optional[Dict[str, str]] = None,
    ) -> Tuple[Path, List[str]]:
        """Persist embeddings to disk and register entities."""

        target_base = self.embeddings_dir / self.model_name
        target_dir = target_base / (dataset_name or "_entities")
        target_dir.mkdir(parents=True, exist_ok=True)

        entity_names: List[str] = []

        for seq_id, tensor in embeddings.items():
            if entity_name_map and seq_id in entity_name_map:
                entity_name = self._sanitize_filename(entity_name_map[seq_id])
            else:
                entity_name = self._build_entity_name(seq_id, embedding_type, dataset_name)
            file_path = target_dir / f"{entity_name}.pkl"
            array = self._to_numpy(tensor)

            payload = {
                "embedding": array,
                "models": self.model_name,
                "embedding_type": embedding_type,
                "source_sequence": seq_id,
                "dataset": dataset_name,
            }
            with open(file_path, "wb") as handle:
                pickle.dump(payload, handle)

            relative_path = str(file_path.relative_to(self.paths.data_root))
            metadata = {
                "models": self.model_name,
                "embedding_type": embedding_type,
                "shape": list(array.shape),
                "dtype": str(array.dtype),
                "source_sequence": seq_id,
            }
            if dataset_name:
                metadata["dataset"] = dataset_name

            self.entity_registry.register_entity(
                name=entity_name,
                format_type=self.processor_type,
                file_path=relative_path,
                metadata=metadata,
            )

            try:
                self.entity_registry.add_relationship(
                    source_name=entity_name,
                    target_name=seq_id,
                    rel_type="derived_from",
                    metadata={
                        "models": self.model_name,
                        "embedding_type": embedding_type,
                    },
                )
            except ValueError:
                # Source sequence not registered; skip relationship silently
                pass

            entity_names.append(entity_name)

        return target_dir, entity_names

    def _save_embeddings_dataset(
        self,
        embeddings: Dict[str, "torch.Tensor"],
        sequences: Dict[str, str],
        dataset_name: str,
        embedding_type: str,
        register_entities: bool = True,
    ) -> None:
        """Persist embeddings and register a dataset via the DatasetManager."""

        dataset_dir, entity_names = self._persist_embeddings(
            embeddings,
            sequences,
            embedding_type,
            dataset_name=dataset_name,
        )

        dataset_metadata = {
            "models": self.model_name,
            "embedding_type": embedding_type,
            "entity_count": len(entity_names),
            "embedding_dim": self.model_config["embedding_dim"],
            "artifact_path": str(dataset_dir.relative_to(self.paths.data_root)),
            "sequence_ids": list(sequences.keys()),
        }

        if self.dataset_manager.dataset_exists(dataset_name):
            self.dataset_manager.delete_dataset(dataset_name)

        self.create_dataset(dataset_name, entity_names, dataset_metadata)

        if register_entities:
            # Entities already registered during persistence. Nothing more to do.
            pass

        cached = {}
        for seq_id, tensor in embeddings.items():
            if _TORCH_AVAILABLE and torch.is_tensor(tensor):
                cached[seq_id] = tensor.detach().cpu()
            else:
                cached[seq_id] = tensor

        self._loaded_embeddings[dataset_name] = cached

        self.logger.info(
            "Saved embeddings dataset '%s' with %d entities", dataset_name, len(entity_names)
        )
    
    def load_embeddings(self, dataset_name: str) -> Dict[str, "torch.Tensor"]:
        """
        Load embeddings from a saved dataset.

        Note: Requires torch to be installed to load the tensors.
        """
        dataset_info = self.get_dataset_info(dataset_name)
        if not dataset_info:
            raise ValueError(f"Dataset not found: {dataset_name}")

        embeddings: Dict[str, "torch.Tensor"] = {}

        for entity in dataset_info.get("entities", []):
            entity_name = entity.get("name")
            if not entity_name:
                continue
            tensor = self.load_entity(entity_name)
            if tensor is None:
                continue

            source_sequence = entity.get("metadata", {}).get("source_sequence")
            # metadata isn't attached in dataset_info; fallback to registry lookup
            if not source_sequence:
                entity_info = self.entity_registry.find_entity(entity_name, self.processor_type)
                if entity_info:
                    source_sequence = entity_info.metadata.get("source_sequence")

            key = source_sequence or entity_name
            embeddings[key] = tensor

        self._loaded_embeddings[dataset_name] = embeddings
        self.logger.info("Loaded %d embeddings from dataset '%s'", len(embeddings), dataset_name)

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
        if not Path(fasta_path).exists():
            raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

        sequences = read_fasta(fasta_path)
        
        # Generate embeddings
        return self.embed_sequences(sequences, embedding_type, save_dataset)
    
    def get_embedding_dim(self) -> int:
        """Get the embedding dimension for the current models."""
        return self.model_config["embedding_dim"]
    
    def list_available_models(self) -> Dict[str, Dict[str, Any]]:
        """List all available models with their configurations."""
        return self.available_models()
    
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

    def collapse_per_residue(
        self,
        per_residue_embedding: Any,
        reduction: str = "mean",
    ) -> Any:
        """Aggregate per-residue embeddings into a single vector."""

        array = self._to_numpy(per_residue_embedding)

        if reduction == "mean":
            aggregated = array.mean(axis=0)
        elif reduction == "sum":
            aggregated = array.sum(axis=0)
        else:
            raise ValueError(f"Unsupported reduction: {reduction}")

        if _TORCH_AVAILABLE:
            return torch.from_numpy(aggregated)
        return aggregated

    def clear_cache(self):
        """Clear models from memory."""
        if self._model is not None:
            del self._model
            self._model = None
        if self._tokenizer is not None:
            del self._tokenizer
            self._tokenizer = None
        if _TORCH_AVAILABLE and torch.cuda.is_available():
            torch.cuda.empty_cache()
        self.logger.info("Cleared models cache")
    
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

        entity_info = self.entity_registry.find_entity(identifier, self.processor_type)

        if entity_info is None:
            # Try resolving by source sequence metadata
            for name in self.entity_registry.list_entities(self.processor_type):
                info = self.entity_registry.find_entity(name, self.processor_type)
                if info and info.metadata.get("source_sequence") == identifier:
                    entity_info = info
                    break

        if entity_info is None:
            self.logger.warning("Embedding entity not found: %s", identifier)
            return None

        file_path = Path(self.paths.data_root) / entity_info.file_path
        if not file_path.exists():
            self.logger.warning("Embedding file missing on disk: %s", file_path)
            return None

        with open(file_path, "rb") as handle:
            payload = pickle.load(handle)

        array = payload.get("embedding") if isinstance(payload, dict) else payload
        if _TORCH_AVAILABLE:
            return torch.from_numpy(np.asarray(array))
        return np.asarray(array)
    
    def list_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all embedding entities.

        Args:
            dataset: Optional dataset to filter by

        Returns:
            List of registered embedding entity names
        """
        if dataset:
            return self.dataset_manager.get_dataset_entities(dataset)
        return self.entity_registry.list_entities(self.processor_type)
    
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
    
    def load_entity(self, entity_id: str, format_type: str = "embedding") -> Optional[Any]:
        """
        Load an embedding entity by ID.

        Args:
            entity_id: Human-readable entity ID (e.g., 'P12345')
            format_type: Format type (default: 'embedding')
            
        Returns:
            Embedding tensor/array or None if not found
        """
        if format_type != self.processor_type:
            return None

        entity_info = self.entity_registry.find_entity(entity_id, self.processor_type)
        if entity_info is None:
            return None

        file_path = Path(self.paths.data_root) / entity_info.file_path
        if not file_path.exists():
            self.logger.warning("Embedding file missing on disk: %s", file_path)
            return None

        with open(file_path, "rb") as handle:
            payload = pickle.load(handle)

        array = payload.get("embedding") if isinstance(payload, dict) else payload
        array = np.asarray(array)

        if _TORCH_AVAILABLE:
            return torch.from_numpy(array)
        return array
    
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
        if format_type != self.processor_type:
            return False

        tensor_dict = {entity_id: data}
        sequence_stub = {entity_id: ""}
        self._persist_embeddings(
            tensor_dict,
            sequence_stub,
            embedding_type="custom",
            dataset_name=None,
            entity_name_map={entity_id: entity_id},
        )
        return True

    # --- Import compact embedding artifacts from external tools ---
    def ingest_from_artifact(
        self,
        dataset_name: str,
        artifact_path: Union[str, Path],
        *,
        model_name: Optional[str] = None,
        embedding_type: str = "per_residue",
        sequences: Optional[Dict[str, str]] = None,
        extra_metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        """Create an embedding dataset from a compact artifact (e.g., NPZ).

        Reads the artifact, converts each array into a standard per-entity PKL
        under embedding/<model>/<dataset>/, registers entities and finally the
        dataset JSON with metadata. Returns the dataset name.
        """
        path = Path(artifact_path)
        if not path.exists():
            raise FileNotFoundError(f"Embedding artifact not found: {path}")

        # Allow overriding model on the processor instance
        if model_name and model_name != self.model_name:
            # Create a new processor for the requested model to keep layout tidy
            other = EmbeddingProcessor(model_name=model_name)
            return other.ingest_from_artifact(
                dataset_name,
                path,
                model_name=model_name,
                embedding_type=embedding_type,
                sequences=sequences,
                extra_metadata=extra_metadata,
            )

        loaded = np.load(path, allow_pickle=False)
        # Build a mapping of id -> np.ndarray
        arrays: Dict[str, np.ndarray] = {k: np.asarray(loaded[k]) for k in loaded.files}

        # Persist using existing helpers to standardize entity layout
        persist_sequences = sequences or {k: "" for k in arrays.keys()}
        # Fake tensors: they can be numpy arrays; _persist_embeddings will handle conversion
        saved_dir, entity_names = self._persist_embeddings(
            arrays,  # type: ignore[arg-type]
            persist_sequences,
            embedding_type,
            dataset_name=dataset_name,
        )

        ds_meta = {
            "models": self.model_name,
            "embedding_type": embedding_type,
            "entity_count": len(entity_names),
            "embedding_dim": self.model_config["embedding_dim"],
            "artifact_path": str(saved_dir.relative_to(self.paths.data_root)),
            "source_artifact": str(path.relative_to(self.paths.data_root))
            if str(path).startswith(self.paths.data_root)
            else str(path),
        }
        if extra_metadata:
            ds_meta.update(extra_metadata)

        if self.dataset_manager.dataset_exists(dataset_name):
            self.dataset_manager.delete_dataset(dataset_name)
        self.create_dataset(dataset_name, entity_names, ds_meta)

        return dataset_name
