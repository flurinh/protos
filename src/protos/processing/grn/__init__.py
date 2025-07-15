"""
GRN processing modules for protein structure analysis.

This package provides tools for working with Generic Residue Numbers (GRNs),
including assignment, validation, and conversion between different formats.
"""

# Version information
__version__ = "0.1.0"

# Import the main GRN processor
from protos.processing.grn.grn_processor import GRNProcessor

# Export the main processor
__all__ = ['GRNProcessor']