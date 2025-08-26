"""Sequence processing module."""

# Import the main sequence processor
from .sequence_processor import SequenceProcessor

from .seq_alignment import (
    init_aligner,
    align_blosum62,
    format_alignment,
    calc_alignment_score_restricted_area,
    get_best_alignment,
    get_score,
    mmseqs2_align,
    mmseqs2_align2,
    msa_blosum62,
    check_chain_similarity,
)

from .mmseqs_utils import (
    detect_mmseqs2,
    get_mmseqs_command,
    ensure_mmseqs2,
    windows_to_wsl_path,
)

__all__ = [
    # Main processor
    'SequenceProcessor',
    # Alignment functions
    'init_aligner',
    'align_blosum62', 
    'format_alignment',
    'calc_alignment_score_restricted_area',
    'get_best_alignment',
    'get_score',
    'mmseqs2_align',
    'mmseqs2_align2',
    'msa_blosum62',
    'check_chain_similarity',
    # MMseqs2 utilities
    'detect_mmseqs2',
    'get_mmseqs_command',
    'ensure_mmseqs2',
    'windows_to_wsl_path',
]