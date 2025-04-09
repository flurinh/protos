"""
Structure utilities with standardized schemas.

This module provides refactored structure processing functions that
use the standardized schemas and interfaces defined in the schema module.
"""

from .structure_utils import (
    load_structure,
    normalize_structures,
    get_ca_ret_coords
)

__all__ = [
    'load_structure',
    'normalize_structures',
    'get_ca_ret_coords'
]