"""Sequence analysis module."""

# Import the main sequence processor
from .sequence_processor import SequenceProcessor

# Re-export alignment helpers from analysis package
from protos.analysis.sequence import (
    init_biopython_aligner,
    perform_pairwise_alignment,
    format_alignment_human,
    run_mmseqs_alignment,
    MMseqsUnavailableError,
)

__all__ = [
    'SequenceProcessor',
    'init_biopython_aligner',
    'perform_pairwise_alignment',
    'format_alignment_human',
    'run_mmseqs_alignment',
    'MMseqsUnavailableError',
]
