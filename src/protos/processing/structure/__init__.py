"""
Structure analysis modules for protein structure analysis.

This package provides tools for working with protein structures,
including loading, filtering, and analyzing structural data.
"""

# Version information
__version__ = "0.1.0"

# Import the main structure processor
from protos.processing.structure.structure_processor import StructureProcessor

# Export the main processor
__all__ = ['StructureProcessor']
