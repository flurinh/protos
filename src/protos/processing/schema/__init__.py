"""
Schema definitions for standardized data models and interfaces.

This package provides:
1. Standard data models and formats
2. Interface definitions
3. Validation utilities
4. Schema conversion tools
"""

from .schema_definitions import (
    # Constants
    GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL,
    
    # Patterns
    GRN_PATTERNS, GRN_FORMAT_DOCS,
    
    # Validation functions
    validate_structure_df, validate_grn_table
)

from .interface_definitions import (
    # Utility functions
    get_correctly_aligned_grns,
    calculate_missing_gene_numbers,
    calculate_missing_ntail_grns,
    calculate_missing_ctail_grns,
    
    # Interfaces
    StructureInterface,
    GRNInterface,
    SequenceInterface,
    
    # Decorators
    validate_structure_operation,
    validate_sequence_operation,
    validate_grn_operation
)

from .grn_utils_updated import (
    # GRN utilities
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    sort_grns,
    validate_grn_string
)

# Version information
__version__ = "0.1.0"