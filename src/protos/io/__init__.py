"""
Protos I/O module - handles all input/output operations.

This module is organized into:
- core/: Core data management infrastructure (EntityRegistry, DatasetManager, etc.)
- formats/: Format-specific utilities (CIF, FASTA, SDF handlers)
- paths/: Path management (ProtosPaths)
- schemas/: Data schemas and interface definitions
- utils/: General utilities (file operations, identifiers, data access)
"""

# Maintain backward compatibility by re-exporting from submodules
from .core import (
    BaseProcessor,
    EntityRegistry,
    DatasetManager,
    InputManager,
    RegistryHealthChecker,
    ConflictResolver
)

from .formats import (
    FormatRegistry,
    FileFormat,
    ProcessorType,
    FormatCategory,
    FormatHandler,
    CifHandler,
    FORMATS
)

# Re-export utilities for backward compatibility
from .formats import *  # All format utilities
from .schemas import *  # All schemas and interfaces
from .utils import *    # All general utilities

__all__ = [
    # Core infrastructure
    'BaseProcessor',
    'EntityRegistry',
    'DatasetManager',
    'InputManager',
    'RegistryHealthChecker',
    'ConflictResolver',
    # Formats
    'FormatRegistry',
    'FileFormat',
    'ProcessorType',
    'FormatCategory',
    'FormatHandler',
    'CifHandler',
    'FORMATS',
]