"""Protos Models - AI model integration for structural biology.

This package provides a unified interface for integrating various AI models
into the Protos framework, including:
- Protein language models (ESM, ProtTrans, etc.)
- Structure prediction models (AlphaFold, RoseTTAFold)
- Property prediction models
- Graph neural networks
- Custom models

The model system follows Protos design principles:
- Zero configuration by default
- Human-readable names everywhere
- Automatic path management
- Entity tracking and dataset support

Key Features:
- Docker-based model isolation (each model runs in its own container)
- Automatic dependency management (no version conflicts)
- Unified API regardless of backend (local or Docker)
- Format validation and conversion
"""

# Import lightweight components first (no heavy ML dependencies)
from protos.models.base_model import BaseModel, ModelConfig, ModelPrediction
from protos.models.model_definitions import (
    ModelDefinition, ModelFramework, InputFormat, OutputFormat,
    get_model_definition, list_available_models, 
    get_models_by_input, get_models_by_output
)
from protos.models.format_schemas import (
    SequenceFormat, StructureFormat, EmbeddingFormat,
    MSAFormat, GraphFormat, GRNFormat, PropertyFormat,
    StructureOutput, AttentionOutput, LogitsOutput,
    FormatConverter
)
from protos.models.format_validators import (
    FormatValidator, FormatAdapter,
    validate_model_compatibility, suggest_format_conversions
)

# Service components (Docker-based execution)
from protos.models.model_service import (
    ServiceStatus, ServiceConfig, RemoteModelService,
    ModelServiceManager, MODEL_SERVICES
)
from protos.models.model_client import (
    ModelBackend, UnifiedModelClient, predict
)

# Lazy imports for components with heavy dependencies
_LAZY_IMPORTS = {
    "StandardModel": "protos.models.standard_model",
    "ModelRegistry": "protos.models.model_registry",
    "ModelDownloader": "protos.models.model_downloader",
    "ModelInstaller": "protos.models.model_installer",
    "ModelManager": "protos.models.model_manager",
}

def __getattr__(name):
    """Lazy import for heavy components to avoid dependency issues."""
    if name in _LAZY_IMPORTS:
        import importlib
        module = importlib.import_module(_LAZY_IMPORTS[name])
        return getattr(module, name)
    raise AttributeError(f"module 'protos.models' has no attribute '{name}'")

__all__ = [
    # Base classes
    'BaseModel',
    'ModelConfig', 
    'ModelPrediction',
    'StandardModel',
    
    # Management
    'ModelRegistry',
    'ModelDownloader',
    'ModelInstaller',
    
    # Definitions
    'ModelDefinition',
    'ModelFramework',
    'InputFormat',
    'OutputFormat',
    'get_model_definition',
    'list_available_models',
    'get_models_by_input',
    'get_models_by_output',
    
    # Format schemas
    'SequenceFormat',
    'StructureFormat',
    'EmbeddingFormat',
    'MSAFormat',
    'GraphFormat',
    'GRNFormat',
    'PropertyFormat',
    'StructureOutput',
    'AttentionOutput',
    'LogitsOutput',
    'FormatConverter',
    
    # Validation
    'FormatValidator',
    'FormatAdapter',
    'validate_model_compatibility',
    'suggest_format_conversions',
    
    # Service components
    "ServiceStatus", "ServiceConfig", "RemoteModelService",
    "ModelServiceManager", "MODEL_SERVICES",
    "ModelBackend", "UnifiedModelClient", "predict",
    
    # Lazy-loaded components
    "StandardModel", "ModelRegistry", "ModelDownloader", "ModelInstaller",
    "ModelManager"
]

# Model type registry for dynamic loading - populate lazily
MODEL_TYPES = {}

def register_model_type(model_type: str, model_class):
    """Register a new model type."""
    MODEL_TYPES[model_type] = model_class

# Register standard model type lazily
def _ensure_standard_model():
    """Ensure StandardModel is registered."""
    if 'standard' not in MODEL_TYPES:
        try:
            from protos.models.standard_model import StandardModel
            MODEL_TYPES['standard'] = StandardModel
        except ImportError:
            pass