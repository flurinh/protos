"""Reusable model definitions and utilities for structural biology.

This package provides lightweight building blocks for integrating AI models
with ProtOS, including:
- Protein language models (ESM, ProtTrans, etc.)
- Structure prediction models (AlphaFold, RoseTTAFold)
- Property prediction models
- Graph neural networks
- Custom models

The models system follows Protos design principles:
- Zero configuration by default
- Human-readable names everywhere
- Automatic path management
- Entity tracking and dataset support

The package intentionally does not provide model execution orchestration.
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

# Lazy imports for components with heavy dependencies
_LAZY_IMPORTS = {
    "StandardModel": "protos.models.standard_model",
    "ModelRegistry": "protos.models.model_registry",
    "ModelDownloader": "protos.models.model_downloader",
    "ModelInstaller": "protos.models.model_installer",
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
    
    # Lazy-loaded components
    "StandardModel", "ModelRegistry", "ModelDownloader", "ModelInstaller"
]

# Model type registry for dynamic loading - populate lazily
MODEL_TYPES = {}

def register_model_type(model_type: str, model_class):
    """Register a new models type."""
    MODEL_TYPES[model_type] = model_class

# Register standard models type lazily
def _ensure_standard_model():
    """Ensure StandardModel is registered."""
    if 'standard' not in MODEL_TYPES:
        try:
            from protos.models.standard_model import StandardModel
            MODEL_TYPES['standard'] = StandardModel
        except ImportError:
            pass
