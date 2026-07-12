"""High-level sequence alignment engine wrapping Biopython and MMseqs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .alignment_utils import (
    init_biopython_aligner,
    perform_pairwise_alignment,
    format_alignment_human,
)
from .mmseqs_interface import run_mmseqs_alignment, MMseqsUnavailableError


@dataclass
class SequenceAlignmentResult:
    seq1_id: str
    seq2_id: str
    score: float
    alignment: List[str]
    method: str
    format: str = "human"


class SequenceAlignmentEngine:
    def __init__(
        self,
        *,
        open_gap_score: float = -10.0,
        extend_gap_score: float = -0.05,
        end_gap_score: float = 0.0,
    ) -> None:
        self.gap_parameters = {
            "open_gap_score": float(open_gap_score),
            "extend_gap_score": float(extend_gap_score),
            "end_gap_score": float(end_gap_score),
        }
        self._aligner = init_biopython_aligner(**self.gap_parameters)

    def align_pairwise(
        self,
        seq1_id: str,
        seq1: str,
        seq2_id: str,
        seq2: str,
    ) -> SequenceAlignmentResult:
        result = perform_pairwise_alignment(seq1, seq2, self._aligner, seq1_id, seq2_id)
        return SequenceAlignmentResult(
            seq1_id=seq1_id,
            seq2_id=seq2_id,
            score=result.score,
            alignment=result.alignment_lines,
            method="biopython",
            format="human",
        )

    def align_pairwise_mmseqs(
        self,
        sequences: Dict[str, str],
    ) -> List[str]:
        try:
            return run_mmseqs_alignment(sequences)
        except MMseqsUnavailableError:
            raise

    @staticmethod
    def to_human(result: SequenceAlignmentResult) -> str:
        return format_alignment_human(result)
