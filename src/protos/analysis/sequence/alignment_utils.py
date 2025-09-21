"""Biopython-based sequence alignment utilities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from protos.processing.sequence.seq_alignment import (
    init_aligner,
    align_blosum62,
)


@dataclass
class AlignmentResult:
    score: float
    alignment_lines: List[str]
    format: str = "human"


def init_biopython_aligner(open_gap_score: int = -10):
    return init_aligner(open_gap_score)


def perform_pairwise_alignment(
    seq1: str,
    seq2: str,
    aligner,
    seq1_id: Optional[str] = None,
    seq2_id: Optional[str] = None,
) -> AlignmentResult:
    alignment = align_blosum62(seq1, seq2, aligner)
    from protos.processing.sequence.seq_alignment import format_alignment

    formatted = format_alignment(alignment)
    return AlignmentResult(score=alignment.score, alignment_lines=formatted)


def format_alignment_human(result: AlignmentResult) -> str:
    return "\n".join(result.alignment_lines)


def format_alignment_table(result: AlignmentResult) -> str:
    rows = ["seq_id\talignment"]
    for idx, line in enumerate(result.alignment_lines, start=1):
        rows.append(f"line_{idx}\t{line}")
    return "\n".join(rows)
