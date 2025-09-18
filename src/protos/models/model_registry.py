"""Model Registry for tracking and managing AI models in Protos.

This module provides a centralized registry for discovering, loading,
and managing different AI models integrated into Protos.
"""

from typing import Dict, List, Optional, Type, Any
from pathlib import Path
import json
import importlib
from datetime import datetime

from protos.io.paths import ProtosPaths
from protos.models.base_model import BaseModel, ModelConfig


class ModelRegistry:
    """
    Central registry for AI models in Protos.
    
    The registry:
    - Tracks available models and their configurations
    - Provides model discovery and loading
    - Manages model versioning
    - Follows Protos patterns (zero config, human names)
    """
    
    # Built-in model types and their implementations
    BUILTIN_MODELS = {
        'lambda': 'protos.models.lambda_model.LambdaModel',
        'esm': 'protos.models.esm_model.ESMModel',
        'alphafold': 'protos.models.alphafold_model.AlphaFoldModel',
        'property_predictor': 'protos.models.property_model.PropertyModel',
        'graph_model': 'protos.models.graph_model.GraphModel'
    }
    
    def __init__(self, paths: Optional[ProtosPaths] = None):
        """
        Initialize model registry.
        
        Args:
            paths: ProtosPaths instance (created if not provided)
        """
        self.paths = paths or ProtosPaths()
        self.models_path = self.paths.get_processor_path("models")
        self.registry_file = self.models_path / "model_registry.json"
        
        # Ensure directory exists
        self.models_path.mkdir(parents=True, exist_ok=True)
        
        # Load existing registry
        self._registry = self._load_registry()
        
        # Model class cache
        self._model_classes: Dict[str, Type[BaseModel]] = {}
    
    def _load_registry(self) -> dict:
        """Load registry from disk."""
        if self.registry_file.exists():
            with open(self.registry_file, 'r') as f:
                return json.load(f)
        return {}
    
    def _save_registry(self):
        """Save registry to disk."""
        with open(self.registry_file, 'w') as f:
            json.dump(self._registry, f, indent=2)
    
    def register_model_class(self, model_type: str, model_class: Type[BaseModel]):
        """
        Register a model class for dynamic loading.
        
        Args:
            model_type: Type identifier for the model
            model_class: The model class (must inherit from BaseModel)
        """
        if not issubclass(model_class, BaseModel):
            raise TypeError(f"{model_class} must inherit from BaseModel")
        
        self._model_classes[model_type] = model_class
    
    def discover_models(self) -> Dict[str, Dict[str, Any]]:
        """
        Discover all available models.
        
        Returns:
            Dictionary of model info keyed by model path
        """
        models = {}
        
        # Check for model directories
        for model_type_dir in self.models_path.iterdir():
            if model_type_dir.is_dir() and not model_type_dir.name.startswith('_'):
                for model_name_dir in model_type_dir.iterdir():
                    if model_name_dir.is_dir():
                        model_key = f"{model_type_dir.name}/{model_name_dir.name}"
                        
                        # Check for config file
                        config_file = model_name_dir / "config.json"
                        if config_file.exists():
                            with open(config_file, 'r') as f:
                                config_data = json.load(f)
                                models[model_key] = {
                                    'path': str(model_name_dir),
                                    'config': config_data,
                                    'discovered': datetime.now().isoformat()
                                }
        
        # Merge with registry
        self._registry.update(models)
        self._save_registry()
        
        return models
    
    def list_models(self, model_type: Optional[str] = None) -> List[str]:
        """
        List available models.
        
        Args:
            model_type: Filter by model type (e.g., 'lambda', 'esm')
            
        Returns:
            List of model identifiers
        """
        models = list(self._registry.keys())
        
        if model_type:
            models = [m for m in models if m.startswith(f"{model_type}/")]
        
        return sorted(models)
    
    def get_model_info(self, model_identifier: str) -> Optional[dict]:
        """
        Get information about a specific model.
        
        Args:
            model_identifier: Model identifier (e.g., 'lambda/opsin_v1')
            
        Returns:
            Model information or None if not found
        """
        return self._registry.get(model_identifier)
    
    def load_model(self, model_identifier: str, 
                   device: Optional[str] = None) -> BaseModel:
        """
        Load a model by identifier.
        
        Args:
            model_identifier: Model identifier (e.g., 'lambda/opsin_v1')
            device: Device to load model on (overrides config)
            
        Returns:
            Loaded model instance
        """
        # Get model info
        model_info = self.get_model_info(model_identifier)
        if not model_info:
            raise ValueError(f"Model '{model_identifier}' not found")
        
        # Parse identifier
        model_type, model_name = model_identifier.split('/', 1)
        
        # Get model class
        model_class = self._get_model_class(model_type)
        
        # Create config
        config_data = model_info['config'].copy()
        if device:
            config_data['device'] = device
        
        config = ModelConfig(**config_data)
        
        # Create model instance
        model = model_class(config=config, paths=self.paths, name=model_name)
        
        # Load weights if available
        weights_dir = Path(model_info['path']) / 'weights'
        if weights_dir.exists():
            latest_checkpoint = self._find_latest_checkpoint(weights_dir)
            if latest_checkpoint:
                model.load_model(latest_checkpoint)
        
        return model
    
    def _get_model_class(self, model_type: str) -> Type[BaseModel]:
        """Get model class for a type."""
        # Check cache
        if model_type in self._model_classes:
            return self._model_classes[model_type]
        
        # Check built-in models
        if model_type in self.BUILTIN_MODELS:
            module_path = self.BUILTIN_MODELS[model_type]
            module_name, class_name = module_path.rsplit('.', 1)
            
            try:
                module = importlib.import_module(module_name)
                model_class = getattr(module, class_name)
                self._model_classes[model_type] = model_class
                return model_class
            except (ImportError, AttributeError) as e:
                raise ImportError(f"Could not load model class for '{model_type}': {e}")
        
        raise ValueError(f"Unknown model type: '{model_type}'")
    
    def _find_latest_checkpoint(self, weights_dir: Path) -> Optional[Path]:
        """Find the latest checkpoint in a weights directory."""
        checkpoints = list(weights_dir.glob("*.pth")) + list(weights_dir.glob("*.pt"))
        if checkpoints:
            # Sort by modification time
            return max(checkpoints, key=lambda p: p.stat().st_mtime)
        return None
    
    def create_model(self, model_type: str, model_name: str, 
                    config: ModelConfig) -> BaseModel:
        """
        Create and register a new model.
        
        Args:
            model_type: Type of model (e.g., 'lambda')
            model_name: Name for this model instance
            config: Model configuration
            
        Returns:
            Created model instance
        """
        # Get model class
        model_class = self._get_model_class(model_type)
        
        # Ensure config has correct model type
        config.model_type = model_type
        config.model_name = model_name
        
        # Create model
        model = model_class(config=config, paths=self.paths, name=model_name)
        
        # Register in registry
        model_key = f"{model_type}/{model_name}"
        self._registry[model_key] = {
            'config': config.to_dict(),
            'path': str(model.model_path),
            'created': datetime.now().isoformat()
        }
        self._save_registry()
        
        return model
    
    def delete_model(self, model_identifier: str, remove_files: bool = False):
        """
        Delete a model from the registry.
        
        Args:
            model_identifier: Model to delete
            remove_files: Whether to remove model files
        """
        if model_identifier not in self._registry:
            raise ValueError(f"Model '{model_identifier}' not found")
        
        model_info = self._registry[model_identifier]
        
        # Remove from registry
        del self._registry[model_identifier]
        self._save_registry()
        
        # Optionally remove files
        if remove_files:
            model_path = Path(model_info['path'])
            if model_path.exists():
                import shutil
                shutil.rmtree(model_path)
    
    def get_models_for_entity(self, entity_name: str, 
                             output_format: Optional[str] = None) -> List[str]:
        """
        Find models that can process a given entity.
        
        Args:
            entity_name: Entity to process
            output_format: Filter by output format
            
        Returns:
            List of compatible model identifiers
        """
        from protos.io.entity_registry import EntityRegistry
        
        # Get entity info
        registry = EntityRegistry(self.paths)
        entity_info = registry.find_entity(entity_name)
        
        if not entity_info:
            return []
        
        # Find compatible models
        compatible_models = []
        
        for model_id, model_info in self._registry.items():
            config = model_info['config']
            
            # Check if entity has required input formats
            required_formats = set(config.get('input_formats', []))
            entity_formats = set(entity_info.formats.keys())
            
            if required_formats.issubset(entity_formats):
                # Check output format if specified
                if output_format is None or config.get('output_format') == output_format:
                    compatible_models.append(model_id)
        
        return compatible_models