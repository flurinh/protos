"""Base Model class for integrating AI models into Protos.

This module provides the foundation for integrating various AI models
(protein language models, structure predictors, property predictors, etc.)
into the Protos framework while maintaining consistency with the existing
processor architecture.
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, List, Optional, Union, Type
from pathlib import Path
import json
from datetime import datetime
from dataclasses import dataclass, field

from protos.io.core.base_processor import BaseProcessor
from protos.io.paths import ProtosPaths
from protos.io.core.entity_registry import EntityRegistry
from protos.io.core.dataset_manager import DatasetManager


@dataclass
class ModelConfig:
    """Configuration for a models."""
    model_type: str
    model_name: str
    version: str
    input_formats: List[str]  # e.g., ['sequence', 'structure', 'graph']
    output_format: str  # e.g., 'embedding', 'property', 'structure'
    model_params: Dict[str, Any] = field(default_factory=dict)
    preprocessing_params: Dict[str, Any] = field(default_factory=dict)
    device: str = "cpu"
    
    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'model_type': self.model_type,
            'model_name': self.model_name,
            'version': self.version,
            'input_formats': self.input_formats,
            'output_format': self.output_format,
            'model_params': self.model_params,
            'preprocessing_params': self.preprocessing_params,
            'device': self.device
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'ModelConfig':
        """Create from dictionary."""
        return cls(**data)


@dataclass
class ModelPrediction:
    """Container for models predictions."""
    entity_name: str
    model_name: str
    prediction: Any
    confidence: Optional[float] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class BaseModel(ABC):
    """
    Base class for all AI models in Protos.
    
    This class provides:
    - Consistent interface for models loading and prediction
    - Integration with Protos path system
    - Automatic entity tracking
    - Dataset support
    - Model registry integration
    
    Following Protos patterns:
    - Zero configuration by default
    - Human-readable names everywhere
    - Automatic path management through ProtosPaths
    """
    
    def __init__(self, 
                 config: ModelConfig,
                 paths: Optional[ProtosPaths] = None,
                 name: Optional[str] = None):
        """
        Initialize base models.
        
        Args:
            config: Model configuration
            paths: ProtosPaths instance (created if not provided)
            name: Human-readable name for this models instance
        """
        # Follow Protos pattern: automatic path management
        self.paths = paths or ProtosPaths()
        self.config = config
        self.name = name or f"{config.model_type}_{config.model_name}"
        
        # Model-specific paths
        self.models_path = Path(self.paths.get_processor_path("models"))
        self.model_path = self.models_path / config.model_type / config.model_name
        self.weights_path = self.model_path / "weights"
        self.cache_path = self.model_path / "cache"
        self.predictions_path = self.model_path / "predictions"
        
        # Ensure directories exist
        for path in [self.weights_path, self.cache_path, self.predictions_path]:
            path.mkdir(parents=True, exist_ok=True)
        
        # Entity tracking (follows Protos pattern)
        self.entity_registry = EntityRegistry(self.paths)
        self.dataset_manager = DatasetManager(self.paths)
        
        # Model registry for tracking available models
        self.registry_path = self.models_path / "model_registry.json"
        self._register_model()
        
        # Initialize the actual models
        self.model = None
        self.is_loaded = False
        
        # Adapters for input/output format conversion
        self.input_adapters = {}
        self.output_adapter = None
        
    def _register_model(self):
        """Register this models in the models registry."""
        registry = {}
        if self.registry_path.exists():
            with open(self.registry_path, 'r') as f:
                registry = json.load(f)
        
        model_key = f"{self.config.model_type}/{self.config.model_name}"
        registry[model_key] = {
            'config': self.config.to_dict(),
            'path': str(self.model_path),
            'registered': datetime.now().isoformat()
        }
        
        with open(self.registry_path, 'w') as f:
            json.dump(registry, f, indent=2)
    
    @abstractmethod
    def load_model(self, checkpoint_path: Optional[Path] = None):
        """
        Load the models weights.
        
        Args:
            checkpoint_path: Path to models checkpoint (uses default if not provided)
        """
        pass
    
    @abstractmethod
    def _predict_single(self, input_data: Any) -> Any:
        """
        Make prediction on single input.
        
        This is the core prediction method that subclasses must implement.
        
        Args:
            input_data: Preprocessed input data
            
        Returns:
            Model output
        """
        pass
    
    def predict(self, 
                entity_names: Union[str, List[str]], 
                save: bool = True,
                batch_size: Optional[int] = None) -> Union[ModelPrediction, List[ModelPrediction]]:
        """
        Make predictions on entities.
        
        Args:
            entity_names: Single entity name or list of names
            save: Whether to save predictions
            batch_size: Batch size for prediction (if supported)
            
        Returns:
            Single prediction or list of predictions
        """
        if not self.is_loaded:
            raise RuntimeError("Model not loaded. Call load_model() first.")
        
        # Handle single vs batch
        single_input = isinstance(entity_names, str)
        if single_input:
            entity_names = [entity_names]
        
        predictions = []
        
        for entity_name in entity_names:
            # Get input data from various sources
            input_data = self._prepare_input(entity_name)
            
            # Make prediction
            output = self._predict_single(input_data)
            
            # Create prediction object
            prediction = ModelPrediction(
                entity_name=entity_name,
                model_name=self.name,
                prediction=output,
                metadata={
                    'model_version': self.config.version,
                    'input_formats': self.config.input_formats
                }
            )
            
            predictions.append(prediction)
            
            # Save if requested
            if save:
                self._save_prediction(prediction)
        
        return predictions[0] if single_input else predictions
    
    def _prepare_input(self, entity_name: str) -> Any:
        """
        Prepare input data for the models.
        
        This method:
        1. Finds the entity in the registry
        2. Loads data from required formats
        3. Applies format adapters
        4. Preprocesses for models input
        """
        # Get entity info
        entity_info = self.entity_registry.find_entity(entity_name)
        if not entity_info:
            raise ValueError(f"Entity '{entity_name}' not found in registry")
        
        # Collect input data from different formats
        input_data = {}
        for format_type in self.config.input_formats:
            if format_type in entity_info.formats:
                # Load data using appropriate processor
                processor = self._get_processor_for_format(format_type)
                data = processor.load_entity(entity_name)
                
                # Apply format adapter if available
                if format_type in self.input_adapters:
                    data = self.input_adapters[format_type](data)
                
                input_data[format_type] = data
            else:
                raise ValueError(f"Entity '{entity_name}' missing required format: {format_type}")
        
        # Model-specific preprocessing
        return self._preprocess_input(input_data)
    
    @abstractmethod
    def _preprocess_input(self, input_data: Dict[str, Any]) -> Any:
        """
        Preprocess input data for models.
        
        Subclasses implement this to convert from Protos formats
        to models-specific input format.
        """
        pass
    
    def _save_prediction(self, prediction: ModelPrediction):
        """Save prediction to disk."""
        # Save to predictions directory with human-readable name
        pred_file = self.predictions_path / f"{prediction.entity_name}_{prediction.timestamp}.json"
        
        with open(pred_file, 'w') as f:
            json.dump({
                'entity_name': prediction.entity_name,
                'model_name': prediction.model_name,
                'prediction': self._serialize_prediction(prediction.prediction),
                'confidence': prediction.confidence,
                'metadata': prediction.metadata,
                'timestamp': prediction.timestamp
            }, f, indent=2)
        
        # Register prediction in entity registry
        self.entity_registry.register_entity(
            name=prediction.entity_name,
            format_type=f"model_{self.config.output_format}",
            file_path=str(pred_file),
            metadata={
                'models': self.name,
                'timestamp': prediction.timestamp
            }
        )
    
    def _serialize_prediction(self, prediction: Any) -> Any:
        """
        Serialize prediction for storage.
        
        Override this for complex prediction types.
        """
        return prediction
    
    def _get_processor_for_format(self, format_type: str) -> BaseProcessor:
        """Get the appropriate processor for a data format."""
        # Import here to avoid circular imports
        if format_type == "structure":
            from protos.processing.structure import CifBaseProcessor
            return CifBaseProcessor(paths=self.paths)
        elif format_type == "sequence":
            from protos.processing.sequence import SeqProcessor
            return SeqProcessor(paths=self.paths)
        elif format_type == "grn":
            from protos.processing.grn import GRNBaseProcessor
            return GRNBaseProcessor(paths=self.paths)
        elif format_type == "property":
            from protos.processing.property import PropertyProcessor
            return PropertyProcessor(paths=self.paths)
        elif format_type == "embedding":
            from protos.processing.embedding import EmbeddingProcessor
            return EmbeddingProcessor(paths=self.paths)
        else:
            raise ValueError(f"Unknown format type: {format_type}")
    
    def register_input_adapter(self, format_type: str, adapter_fn):
        """
        Register an adapter function for input format conversion.
        
        Args:
            format_type: The input format type (e.g., 'sequence', 'structure')
            adapter_fn: Function that converts from Protos format to models input
        """
        self.input_adapters[format_type] = adapter_fn
    
    def register_output_adapter(self, adapter_fn):
        """
        Register an adapter function for output format conversion.
        
        Args:
            adapter_fn: Function that converts models output to Protos format
        """
        self.output_adapter = adapter_fn
    
    def create_dataset(self, dataset_name: str, entity_names: List[str], 
                      metadata: Optional[dict] = None):
        """
        Create a dataset for this models.
        
        Args:
            dataset_name: Name for the dataset
            entity_names: List of entity names to include
            metadata: Optional metadata for the dataset
        """
        self.dataset_manager.create_dataset(
            name=dataset_name,
            entity_names=entity_names,
            processor_type=f"model_{self.config.model_type}",
            metadata=metadata or {}
        )
    
    def predict_dataset(self, dataset_name: str, 
                       save: bool = True,
                       batch_size: Optional[int] = None) -> List[ModelPrediction]:
        """
        Make predictions on an entire dataset.
        
        Args:
            dataset_name: Name of the dataset
            save: Whether to save predictions
            batch_size: Batch size for prediction
            
        Returns:
            List of predictions
        """
        # Get dataset info
        dataset = self.dataset_manager.get_dataset(
            dataset_name, 
            f"model_{self.config.model_type}"
        )
        
        if not dataset:
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        # Get entity names
        entity_names = [
            self.entity_registry.get_human_name(eid) 
            for eid in dataset.entity_ids
        ]
        
        # Predict on all entities
        return self.predict(entity_names, save=save, batch_size=batch_size)
    
    def list_predictions(self, entity_name: Optional[str] = None) -> List[str]:
        """
        List available predictions.
        
        Args:
            entity_name: Filter by entity name (returns all if None)
            
        Returns:
            List of prediction file paths
        """
        pattern = f"{entity_name}_*.json" if entity_name else "*.json"
        return [str(p) for p in self.predictions_path.glob(pattern)]
    
    def load_prediction(self, entity_name: str, 
                       timestamp: Optional[str] = None) -> ModelPrediction:
        """
        Load a saved prediction.
        
        Args:
            entity_name: Entity name
            timestamp: Specific timestamp (uses latest if None)
            
        Returns:
            Loaded prediction
        """
        # Find prediction files
        predictions = self.list_predictions(entity_name)
        if not predictions:
            raise ValueError(f"No predictions found for '{entity_name}'")
        
        # Use latest or specific timestamp
        if timestamp:
            pred_file = self.predictions_path / f"{entity_name}_{timestamp}.json"
        else:
            pred_file = sorted(Path(p) for p in predictions)[-1]  # Latest
        
        # Load prediction
        with open(pred_file, 'r') as f:
            data = json.load(f)
        
        return ModelPrediction(
            entity_name=data['entity_name'],
            model_name=data['model_name'],
            prediction=data['prediction'],
            confidence=data.get('confidence'),
            metadata=data.get('metadata', {}),
            timestamp=data['timestamp']
        )