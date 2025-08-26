"""
Embedding processing module for protein sequences.

This module provides functionality for generating embeddings from protein sequences
using various transformer models (ESM-2, Ankh, etc.).
"""

from .embedding_processor import EmbeddingProcessor

__all__ = ['EmbeddingProcessor']