"""
Protos workflows module for cross-format operations.

This module provides high-level workflow functions that coordinate
operations across multiple processors and data formats.
"""

from .cross_format import (
    sequence_to_structure_workflow,
    structure_to_sequence_workflow,
    sequence_to_grn_workflow,
    any_format_to_embeddings_workflow,
    track_conversion_lineage
)

__all__ = [
    'sequence_to_structure_workflow',
    'structure_to_sequence_workflow', 
    'sequence_to_grn_workflow',
    'any_format_to_embeddings_workflow',
    'track_conversion_lineage'
]