"""Format-specific utilities and handlers."""

from .format_registry import (
    FormatRegistry,
    FileFormat,
    ProcessorType,
    FormatCategory,
    FORMATS
)
from .formats import FormatHandler
# Import CifHandler later to avoid circular import
from .cif_utils import *  # Re-export all CIF utilities
from .fasta_utils import *  # Re-export all FASTA utilities
from .sdf_utils import *  # Re-export all SDF utilities

# Now import CifHandler after FormatHandler is available
from .cif_handler import CifHandler

__all__ = [
    'FormatRegistry',
    'FileFormat',
    'ProcessorType',
    'FormatCategory',
    'FormatHandler',
    'CifHandler',
    'FORMATS'
]