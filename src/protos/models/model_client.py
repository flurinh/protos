"""Unified models client that supports both local and containerized models.

This module provides a single interface for working with models regardless
of whether they run locally or in Docker containers.
"""

from typing import Any, Dict, List, Optional, Union
from pathlib import Path
from enum import Enum
import warnings

import numpy as np
import pandas as pd

from protos.models.base_model import BaseModel, ModelConfig
from protos.models.model_definitions import get_model_definition, InputFormat
from protos.models.model_service import (
    RemoteModelService, ModelServiceManager, MODEL_SERVICES
)
from protos.models.format_schemas import (
    SequenceFormat, StructureFormat, EmbeddingFormat,
    GRNFormat, PropertyFormat
)
from protos.models.format_validators import FormatValidator, FormatAdapter
from protos.io.paths import ProtosPaths
from protos.io.core.entity_registry import EntityRegistry


class ModelBackend(Enum):
    """Model execution backend."""
    LOCAL = "local"      # Run locally (requires dependencies)
    DOCKER = "docker"    # Run in Docker container
    AUTO = "auto"        # Auto-detect best option


class UnifiedModelClient:
    """
    Unified client for all Protos models.
    
    This client:
    - Automatically selects between local and Docker execution
    - Handles format validation and conversion
    - Manages entity integration
    - Provides consistent interface regardless of backend
    """
    
    def __init__(self, 
                 paths: Optional[ProtosPaths] = None,
                 backend: ModelBackend = ModelBackend.AUTO,
                 docker_host: Optional[str] = None):
        """
        Initialize the unified models client.
        
        Args:
            paths: ProtosPaths instance
            backend: Which backend to use for models
            docker_host: Docker host URL (if not local)
        """
        self.paths = paths or ProtosPaths()
        self.backend = backend
        self.docker_host = docker_host
        
        # Components
        self.entity_registry = EntityRegistry(self.paths)
        self.format_validator = FormatValidator()
        self.format_adapter = FormatAdapter()
        
        # Model instances
        self.local_models: Dict[str, BaseModel] = {}
        self.service_manager: Optional[ModelServiceManager] = None
        self.remote_services: Dict[str, RemoteModelService] = {}
        
        # Initialize Docker if needed
        if self.backend in [ModelBackend.DOCKER, ModelBackend.AUTO]:
            try:
                self.service_manager = ModelServiceManager(docker_host)
            except Exception as e:
                if self.backend == ModelBackend.DOCKER:
                    raise RuntimeError(f"Docker backend requested but unavailable: {e}")
                else:
                    warnings.warn(f"Docker not available: {e}. Using local backend.")
                    self.backend = ModelBackend.LOCAL
    
    def list_models(self) -> List[str]:
        """List all available models."""
        from protos.models.model_definitions import STANDARD_MODELS
        return list(STANDARD_MODELS.keys())
    
    def get_model_info(self, model_name: str) -> Dict[str, Any]:
        """Get information about a models."""
        definition = get_model_definition(model_name)
        
        info = {
            "name": definition.name,
            "full_name": definition.full_name,
            "description": definition.description,
            "framework": definition.framework.value,
            "input_formats": [fmt.value for fmt in definition.input_formats],
            "output_format": definition.output_format.value,
            "variants": list(definition.sources.keys()),
            "paper": definition.paper_url,
            "backend_available": []
        }
        
        # Check local availability
        try:
            self._check_local_requirements(model_name)
            info["backend_available"].append("local")
        except:
            pass
        
        # Check Docker availability
        if self.service_manager and model_name in MODEL_SERVICES:
            info["backend_available"].append("docker")
        
        return info
    
    def predict(self,
                model_name: str,
                entity_name: Optional[str] = None,
                inputs: Optional[Dict[str, Any]] = None,
                model_variant: Optional[str] = None,
                backend: Optional[ModelBackend] = None,
                **kwargs) -> Dict[str, Any]:
        """
        Make a prediction using a models.
        
        Args:
            model_name: Name of the models to use
            entity_name: Entity to load data for (optional)
            inputs: Direct input data (if not using entity)
            model_variant: Specific models variant
            backend: Override default backend
            **kwargs: Additional models-specific parameters
            
        Returns:
            Prediction results
        """
        # Determine backend
        backend = backend or self.backend
        if backend == ModelBackend.AUTO:
            backend = self._select_backend(model_name)
        
        # Get input data
        if entity_name:
            inputs = self._load_entity_data(entity_name, model_name)
        elif inputs is None:
            raise ValueError("Either entity_name or inputs must be provided")
        
        # Validate inputs
        definition = get_model_definition(model_name)
        self._validate_inputs(inputs, definition.input_formats)
        
        # Make prediction based on backend
        if backend == ModelBackend.DOCKER:
            return self._predict_docker(model_name, inputs, model_variant, **kwargs)
        else:
            return self._predict_local(model_name, inputs, model_variant, **kwargs)
    
    def _select_backend(self, model_name: str) -> ModelBackend:
        """Auto-select the best backend for a models."""
        # Check if Docker service is available
        if self.service_manager and model_name in MODEL_SERVICES:
            # Check if service is running
            services = self.service_manager.list_services()
            if services.get(model_name) and services[model_name].value == "ready":
                return ModelBackend.DOCKER
        
        # Check local requirements
        try:
            self._check_local_requirements(model_name)
            return ModelBackend.LOCAL
        except:
            pass
        
        # Try to start Docker service
        if self.service_manager and model_name in MODEL_SERVICES:
            return ModelBackend.DOCKER
        
        # Fallback to local and let it fail with clear error
        return ModelBackend.LOCAL
    
    def _check_local_requirements(self, model_name: str):
        """Check if local execution is possible."""
        definition = get_model_definition(model_name)
        
        # Check framework
        if definition.framework.value == "pytorch":
            import torch
        elif definition.framework.value == "tensorflow":
            import tensorflow
        elif definition.framework.value == "jax":
            import jax
        
        # Model-specific checks
        if model_name == "esm2":
            import esm_test
        elif model_name == "ankh":
            import transformers
        # Add other models checks as needed
    
    def _predict_local(self, 
                       model_name: str,
                       inputs: Dict[str, Any],
                       model_variant: Optional[str] = None,
                       **kwargs) -> Dict[str, Any]:
        """Make prediction using local models."""
        # Get or create models instance
        model_key = f"{model_name}_{model_variant or 'default'}"
        
        if model_key not in self.local_models:
            # Import dynamically to avoid dependency issues
            try:
                from protos.models.standard_model import StandardModel
                model = StandardModel(
                    model_name=model_name,
                    model_variant=model_variant,
                    paths=self.paths,
                    device=kwargs.get("device", "cpu")
                )
                model.load_model()
                self.local_models[model_key] = model
            except ImportError as e:
                raise RuntimeError(
                    f"Cannot run {model_name} locally. Missing dependencies: {e}\n"
                    f"Consider using Docker backend: backend=ModelBackend.DOCKER"
                )
        
        model = self.local_models[model_key]
        
        # Adapt inputs for models
        adapted_inputs = self._adapt_inputs(inputs, model_name, kwargs)
        
        # Make prediction
        result = model.predict(inputs=adapted_inputs)
        
        return result
    
    def _predict_docker(self,
                        model_name: str,
                        inputs: Dict[str, Any],
                        model_variant: Optional[str] = None,
                        **kwargs) -> Dict[str, Any]:
        """Make prediction using Docker service."""
        if not self.service_manager:
            raise RuntimeError("Docker backend not available")
        
        # Start service if needed
        if model_name not in self.remote_services:
            print(f"Starting {model_name} service...")
            service = self.service_manager.start_service(model_name)
            self.remote_services[model_name] = service
        
        service = self.remote_services[model_name]
        
        # Adapt inputs for models
        adapted_inputs = self._adapt_inputs(inputs, model_name, kwargs)
        
        # Add models-specific parameters
        request_data = {
            **adapted_inputs,
            "model_variant": model_variant
        }
        
        # Make prediction
        result = service.predict(request_data)
        
        return result
    
    def _load_entity_data(self, entity_name: str, model_name: str) -> Dict[str, Any]:
        """Load entity data in formats required by models."""
        definition = get_model_definition(model_name)
        required_formats = [fmt.value for fmt in definition.input_formats]
        
        inputs = {}
        
        # Load each required format
        for fmt in required_formats:
            if fmt == "sequence":
                # Try sequence processor
                try:
                    from protos.processing.sequence.sequence_processor import SequenceProcessor
                    proc = SequenceProcessor(paths=self.paths)
                    if entity_name in proc.list_entities():
                        sequence = proc.load_entity(entity_name)
                        inputs["sequence"] = SequenceFormat(
                            sequence=sequence,
                            sequence_id=entity_name
                        )
                except:
                    pass
            
            elif fmt == "structure":
                # Try structure processor
                try:
                    from protos.processing.structure.structure_processor import StructureProcessor
                    proc = StructureProcessor(paths=self.paths)
                    if entity_name in proc.list_entities():
                        structure_df = proc.load_entity(entity_name)
                        if structure_df is not None:
                            chains = structure_df['auth_chain_id'].unique()
                            inputs["structure"] = StructureFormat(
                                coordinates=structure_df,
                                pdb_id=entity_name,
                                chains=list(chains)
                            )
                except:
                    pass
            
            elif fmt == "grn":
                # Try GRN processor
                try:
                    from protos.processing.grn.grn_processor import GRNProcessor
                    proc = GRNProcessor(paths=self.paths)
                    if entity_name in proc.list_entities():
                        grn_data = proc.load_entity(entity_name)
                        if grn_data is not None:
                            inputs["grn"] = GRNFormat(
                                grn_positions=grn_data,
                                sequence=inputs.get("sequence", {}).get("sequence", ""),
                                grn_system="ballesteros_weinstein"
                            )
                except:
                    pass
        
        if not inputs:
            raise ValueError(f"No data found for entity {entity_name}")
        
        # Check we have all required formats
        missing = set(required_formats) - set(inputs.keys())
        if missing:
            raise ValueError(f"Entity {entity_name} missing required formats: {missing}")
        
        return inputs
    
    def _validate_inputs(self, inputs: Dict[str, Any], required_formats: List[InputFormat]):
        """Validate that inputs match required formats."""
        required = [fmt.value for fmt in required_formats]
        
        # Check all required formats present
        missing = set(required) - set(inputs.keys())
        if missing:
            raise ValueError(f"Missing required input formats: {missing}")
        
        # Validate each input
        for fmt_name, data in inputs.items():
            is_valid, error = self.format_validator.validate_input(data, fmt_name)
            if not is_valid:
                raise ValueError(f"Invalid {fmt_name} format: {error}")
    
    def _adapt_inputs(self, 
                      inputs: Dict[str, Any], 
                      model_name: str,
                      config: Dict[str, Any]) -> Dict[str, Any]:
        """Adapt inputs for specific models requirements."""
        return self.format_adapter.adapt_for_model(inputs, model_name, config)
    
    def shutdown(self):
        """Shutdown all services and cleanup."""
        # Stop Docker services
        if self.service_manager:
            for service_name in self.remote_services:
                self.service_manager.stop_service(service_name)
        
        # Clear local models
        self.local_models.clear()
        self.remote_services.clear()


# Convenience functions
def predict(model_name: str, 
            entity_name: Optional[str] = None,
            inputs: Optional[Dict[str, Any]] = None,
            backend: ModelBackend = ModelBackend.AUTO,
            **kwargs) -> Dict[str, Any]:
    """
    Quick prediction function.
    
    Examples:
        # Using entity
        result = predict("esm2", entity_name="BACR_HALSA")
        
        # Using direct input
        result = predict("esm2", inputs={"sequence": "MKTAYIAK"})
        
        # Force Docker backend
        result = predict("lambda", entity_name="protein1", backend=ModelBackend.DOCKER)
    """
    client = UnifiedModelClient(backend=backend)
    try:
        return client.predict(model_name, entity_name, inputs, **kwargs)
    finally:
        client.shutdown()