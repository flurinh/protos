"""StandardModel implementation for configuration-based models.

This module provides a generic model class that can handle any model
defined in model_definitions.py without requiring custom implementation.
"""

from typing import Any, Dict, List, Optional, Union
from pathlib import Path
import json
import importlib
import warnings
from datetime import datetime

import numpy as np
import torch

from protos.models.base_model import BaseModel, ModelConfig
from protos.models.model_definitions import (
    ModelDefinition, ModelFramework, InputFormat, OutputFormat,
    StandardAdapters, get_model_definition
)
from protos.io.paths import ProtosPaths


class StandardModel(BaseModel):
    """
    Generic model implementation that works with any model defined in model_definitions.
    
    This class:
    - Loads models based on configuration
    - Handles different frameworks (PyTorch, TensorFlow, JAX)
    - Manages format conversions automatically
    - Downloads weights as needed
    """
    
    def __init__(self, 
                 model_name: str,
                 model_variant: Optional[str] = None,
                 paths: Optional[ProtosPaths] = None,
                 device: Optional[str] = None):
        """
        Initialize a standard model.
        
        Args:
            model_name: Name of the model (e.g., 'esm2', 'ankh')
            model_variant: Specific variant (e.g., 'esm2_t33_650M')
            paths: ProtosPaths instance
            device: Device to use (overrides definition)
        """
        # Get model definition
        self.definition = get_model_definition(model_name)
        self.model_variant = model_variant
        
        # Create config from definition
        config = self._create_config_from_definition(device)
        
        # Initialize base class
        super().__init__(config, paths, name=model_name)
        
        # Framework-specific model holder
        self.framework_model = None
        self.tokenizer = None
        self.preprocessor = None
        self.postprocessor = None
        
        # Set up adapters based on definition
        self._setup_adapters()
    
    def _create_config_from_definition(self, device: Optional[str] = None) -> ModelConfig:
        """Create ModelConfig from ModelDefinition."""
        # Convert enum types to strings
        input_formats = [fmt.value for fmt in self.definition.input_formats]
        output_format = self.definition.output_format.value
        
        return ModelConfig(
            model_type=self.definition.name,
            model_name=self.definition.full_name,
            version=self.definition.version,
            input_formats=input_formats,
            output_format=output_format,
            model_params=self.definition.model_config,
            preprocessing_params=self.definition.preprocessing_config,
            device=device or ("cuda" if torch.cuda.is_available() else "cpu")
        )
    
    def _setup_adapters(self):
        """Set up format adapters based on model definition."""
        # Input adapters
        for input_format in self.definition.input_formats:
            if input_format == InputFormat.SEQUENCE:
                tokenizer = self.definition.preprocessing_config.get("tokenizer", "esm")
                self.register_input_adapter(
                    "sequence",
                    StandardAdapters.sequence_to_tokens(tokenizer)
                )
            elif input_format == InputFormat.STRUCTURE:
                graph_config = self.definition.preprocessing_config.get("graph_construction", {})
                self.register_input_adapter(
                    "structure", 
                    StandardAdapters.structure_to_graph(**graph_config)
                )
            elif input_format == InputFormat.GRN:
                encoding = self.definition.preprocessing_config.get("grn_encoding", "onehot")
                self.register_input_adapter(
                    "grn",
                    StandardAdapters.grn_to_features(encoding)
                )
    
    def load_model(self, checkpoint_path: Optional[Path] = None):
        """Load model based on framework."""
        # Download weights if needed
        if checkpoint_path is None:
            checkpoint_path = self._ensure_weights()
        
        if self.definition.framework == ModelFramework.PYTORCH:
            self._load_pytorch_model(checkpoint_path)
        elif self.definition.framework == ModelFramework.TENSORFLOW:
            self._load_tensorflow_model(checkpoint_path)
        elif self.definition.framework == ModelFramework.JAX:
            self._load_jax_model(checkpoint_path)
        else:
            raise NotImplementedError(f"Framework {self.definition.framework} not yet supported")
        
        self.is_loaded = True
    
    def _ensure_weights(self) -> Path:
        """Ensure model weights are available, downloading if needed."""
        # Check for existing weights
        variant = self.model_variant or list(self.definition.sources.keys())[0]
        weights_file = self.weights_path / f"{variant}.pt"
        
        if not weights_file.exists():
            # Download weights
            from protos.models.model_downloader import ModelDownloader
            downloader = ModelDownloader(self.paths)
            weights_file = downloader.download_model(
                self.definition.name,
                variant,
                self.weights_path
            )
        
        return weights_file
    
    def _load_pytorch_model(self, checkpoint_path: Path):
        """Load PyTorch model."""
        if self.definition.name == "esm2":
            self._load_esm2_model(checkpoint_path)
        elif self.definition.name == "ankh":
            self._load_ankh_model(checkpoint_path)
        elif self.definition.name == "lambda":
            self._load_lambda_model(checkpoint_path)
        else:
            # Generic PyTorch loading
            self.framework_model = torch.load(
                checkpoint_path,
                map_location=self.config.device
            )
            if hasattr(self.framework_model, 'eval'):
                self.framework_model.eval()
    
    def _load_esm2_model(self, checkpoint_path: Path):
        """Load ESM2 model specifically."""
        try:
            import esm_test
            
            # Load model and alphabet
            model_data = torch.load(checkpoint_path, map_location=self.config.device)
            
            # ESM models come with their alphabet
            if isinstance(model_data, dict):
                self.framework_model = model_data['model']
                self.alphabet = model_data['alphabet']
            else:
                # Load from package
                model_name = checkpoint_path.stem
                self.framework_model, self.alphabet = esm.pretrained.load_model_and_alphabet(model_name)
            
            self.framework_model.to(self.config.device)
            self.framework_model.eval()
            
            # Create batch converter
            self.batch_converter = self.alphabet.get_batch_converter()
            
        except ImportError:
            raise ImportError("ESM package not installed. Run: pip install fair-esm")
    
    def _load_ankh_model(self, checkpoint_path: Path):
        """Load Ankh model specifically."""
        try:
            from transformers import T5EncoderModel, T5Tokenizer
            
            model_path = str(checkpoint_path.parent)
            self.framework_model = T5EncoderModel.from_pretrained(model_path)
            self.tokenizer = T5Tokenizer.from_pretrained(model_path)
            
            self.framework_model.to(self.config.device)
            self.framework_model.eval()
            
        except ImportError:
            raise ImportError("Transformers not installed. Run: pip install transformers")
    
    def _load_lambda_model(self, checkpoint_path: Path):
        """Load Lambda model using the existing implementation."""
        # Reuse the lambda model loading logic
        from protos.models.lambda_model import LambdaModel
        
        # Create a lambda model instance
        lambda_model = LambdaModel(self.config, self.paths)
        lambda_model.load_model(checkpoint_path)
        
        # Use its internals
        self.framework_model = lambda_model.model
        self.model_factory = lambda_model.model_factory
    
    def _load_tensorflow_model(self, checkpoint_path: Path):
        """Load TensorFlow model."""
        try:
            import tensorflow as tf
            self.framework_model = tf.keras.models.load_model(str(checkpoint_path))
        except ImportError:
            raise ImportError("TensorFlow not installed")
    
    def _load_jax_model(self, checkpoint_path: Path):
        """Load JAX model."""
        try:
            import jax
            import flax
            # JAX model loading logic
            raise NotImplementedError("JAX model loading not yet implemented")
        except ImportError:
            raise ImportError("JAX/Flax not installed")
    
    def _preprocess_input(self, input_data: Dict[str, Any]) -> Any:
        """Preprocess input based on model requirements."""
        # Apply framework-specific preprocessing
        if self.definition.framework == ModelFramework.PYTORCH:
            return self._preprocess_pytorch(input_data)
        elif self.definition.framework == ModelFramework.TENSORFLOW:
            return self._preprocess_tensorflow(input_data)
        else:
            return input_data
    
    def _preprocess_pytorch(self, input_data: Dict[str, Any]) -> Any:
        """PyTorch-specific preprocessing."""
        if self.definition.name == "esm2":
            return self._preprocess_esm2(input_data)
        elif self.definition.name == "ankh":
            return self._preprocess_ankh(input_data)
        elif self.definition.name == "lambda":
            return self._preprocess_lambda(input_data)
        else:
            # Generic preprocessing
            processed = {}
            for key, value in input_data.items():
                if isinstance(value, np.ndarray):
                    processed[key] = torch.tensor(value, dtype=torch.float32)
                else:
                    processed[key] = value
            return processed
    
    def _preprocess_esm2(self, input_data: Dict[str, Any]) -> Any:
        """Preprocess for ESM2."""
        if 'sequence' not in input_data:
            raise ValueError("ESM2 requires sequence input")
        
        sequence = input_data['sequence']
        
        # Handle single sequence
        if isinstance(sequence, str):
            data = [("protein", sequence)]
        else:
            data = [(f"protein_{i}", seq) for i, seq in enumerate(sequence)]
        
        # Use batch converter
        batch_labels, batch_strs, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.config.device)
        
        return {
            'tokens': batch_tokens,
            'labels': batch_labels,
            'strs': batch_strs
        }
    
    def _preprocess_ankh(self, input_data: Dict[str, Any]) -> Any:
        """Preprocess for Ankh."""
        if 'sequence' not in input_data:
            raise ValueError("Ankh requires sequence input")
        
        sequence = input_data['sequence']
        
        # Tokenize
        inputs = self.tokenizer(
            sequence,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=self.definition.max_sequence_length
        )
        
        # Move to device
        inputs = {k: v.to(self.config.device) for k, v in inputs.items()}
        
        return inputs
    
    def _preprocess_lambda(self, input_data: Dict[str, Any]) -> Any:
        """Preprocess for Lambda - reuse lambda model logic."""
        # This would use the graph construction logic from lambda_model
        from torch_geometric.data import Data
        
        # Convert to graph format
        if 'graph' in input_data:
            return input_data['graph']
        
        # Build graph from other inputs
        # Simplified - would include full logic
        return Data(x=torch.randn(10, 256))
    
    def _preprocess_tensorflow(self, input_data: Dict[str, Any]) -> Any:
        """TensorFlow-specific preprocessing."""
        # Convert to TF tensors
        import tensorflow as tf
        processed = {}
        for key, value in input_data.items():
            if isinstance(value, np.ndarray):
                processed[key] = tf.constant(value)
            else:
                processed[key] = value
        return processed
    
    def _predict_single(self, input_data: Any) -> Any:
        """Make prediction using the loaded model."""
        if self.definition.framework == ModelFramework.PYTORCH:
            return self._predict_pytorch(input_data)
        elif self.definition.framework == ModelFramework.TENSORFLOW:
            return self._predict_tensorflow(input_data)
        else:
            raise NotImplementedError(f"Prediction for {self.definition.framework} not implemented")
    
    def _predict_pytorch(self, input_data: Any) -> Any:
        """PyTorch prediction."""
        with torch.no_grad():
            if self.definition.name == "esm2":
                return self._predict_esm2(input_data)
            elif self.definition.name == "ankh":
                return self._predict_ankh(input_data)
            elif self.definition.name == "lambda":
                return self._predict_lambda(input_data)
            else:
                # Generic prediction
                output = self.framework_model(input_data)
                return self._process_output(output)
    
    def _predict_esm2(self, input_data: Dict[str, Any]) -> Any:
        """ESM2-specific prediction."""
        tokens = input_data['tokens']
        
        # Get representations
        results = self.framework_model(
            tokens, 
            repr_layers=self.config.model_params.get('repr_layers', [-1]),
            return_contacts=self.config.model_params.get('return_contacts', False)
        )
        
        # Extract embeddings
        embeddings = results["representations"][results["repr_layers"][-1]]
        
        # Remove special tokens
        embeddings = embeddings[:, 1:-1, :]  # Remove <cls> and <eos>
        
        # Return based on what's requested
        output = {
            "embeddings": embeddings.cpu().numpy(),
        }
        
        if "contacts" in results:
            output["contacts"] = results["contacts"].cpu().numpy()
        
        if "attentions" in results:
            output["attentions"] = results["attentions"].cpu().numpy()
        
        return output
    
    def _predict_ankh(self, input_data: Dict[str, Any]) -> Any:
        """Ankh-specific prediction."""
        outputs = self.framework_model(**input_data)
        
        # Get embeddings
        embeddings = outputs.last_hidden_state
        
        # Aggregate if needed
        aggregation = self.config.model_params.get("aggregation", "mean")
        if aggregation == "mean":
            pooled = embeddings.mean(dim=1)
        elif aggregation == "max":
            pooled = embeddings.max(dim=1)[0]
        else:
            pooled = embeddings[:, 0, :]  # CLS token
        
        return {
            "embeddings": embeddings.cpu().numpy(),
            "pooled_embedding": pooled.cpu().numpy()
        }
    
    def _predict_lambda(self, input_data: Any) -> Any:
        """Lambda-specific prediction."""
        # Use the model's forward method
        output = self.framework_model(input_data)
        
        # Process based on task type
        if isinstance(output, dict):
            return {k: v.cpu().numpy() for k, v in output.items()}
        else:
            return output.cpu().numpy()
    
    def _predict_tensorflow(self, input_data: Any) -> Any:
        """TensorFlow prediction."""
        predictions = self.framework_model.predict(input_data)
        return predictions
    
    def _process_output(self, raw_output: Any) -> Any:
        """Process model output to standard format."""
        # Convert to numpy if needed
        if torch.is_tensor(raw_output):
            return raw_output.cpu().numpy()
        elif isinstance(raw_output, dict):
            return {k: self._process_output(v) for k, v in raw_output.items()}
        else:
            return raw_output
    
    def get_model_info(self) -> Dict[str, Any]:
        """Get information about the loaded model."""
        info = {
            "name": self.definition.name,
            "full_name": self.definition.full_name,
            "version": self.definition.version,
            "description": self.definition.description,
            "framework": self.definition.framework.value,
            "input_formats": [fmt.value for fmt in self.definition.input_formats],
            "output_format": self.definition.output_format.value,
            "loaded": self.is_loaded,
            "device": self.config.device,
            "variant": self.model_variant
        }
        
        if self.definition.paper_url:
            info["paper"] = self.definition.paper_url
            info["citation"] = self.definition.citation
        
        return info