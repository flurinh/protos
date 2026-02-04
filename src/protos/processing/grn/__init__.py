"""
GRN analysis modules for protein structure analysis.

This package provides tools for working with Generic Residue Numbers (GRNs),
including assignment, validation, and conversion between different formats.
"""

# Version information
__version__ = "0.1.0"

# Import the main GRN processor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.binding_domain import BindingDomainConfig
from protos.processing.grn.grn_utils import normalize_grn_format, sort_grns_str

# Export the main processor and utilities
__all__ = [
    'GRNProcessor',
    'BindingDomainConfig',
    'normalize_grn_format',
    'sort_grns_str',
]