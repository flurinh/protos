"""
Embedding analysis module for protein sequences.

This module provides functionality for generating embeddings from protein sequences
using various transformer models (ESM-2, Ankh, etc.).
"""

from .embedding_processor import EmbeddingProcessor
from .esmc_adapter import (
    ESMCArtifactError,
    ESMCBackend,
    ESMCBatchOutput,
    ESMCEmbeddingResult,
    ESMCLoadPolicy,
    ESMCModelProvenance,
    ESMCOutputProvenance,
    ESMCSequenceConflictError,
    ESMCShapeError,
    ESMCSourceLineage,
    ESMCTokenMapping,
    ESMCTokenPolicy,
    ESMCTokenPolicyError,
    ESMCTransformersProvenance,
    ESMCTruncationError,
    ESMCValidationError,
    HuggingFaceESMCBackend,
    normalize_sequence,
    sequence_sha256,
    validate_transformers_runtime,
)

__all__ = [
    "EmbeddingProcessor",
    "ESMCArtifactError",
    "ESMCBackend",
    "ESMCBatchOutput",
    "ESMCEmbeddingResult",
    "ESMCLoadPolicy",
    "ESMCModelProvenance",
    "ESMCOutputProvenance",
    "ESMCSequenceConflictError",
    "ESMCShapeError",
    "ESMCSourceLineage",
    "ESMCTokenMapping",
    "ESMCTokenPolicy",
    "ESMCTokenPolicyError",
    "ESMCTransformersProvenance",
    "ESMCTruncationError",
    "ESMCValidationError",
    "HuggingFaceESMCBackend",
    "normalize_sequence",
    "sequence_sha256",
    "validate_transformers_runtime",
]
