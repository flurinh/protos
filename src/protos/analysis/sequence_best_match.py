"""Sequence best-match utilities (Biopython + optional MMseqs)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Optional, Tuple

from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine


@dataclass
class BestMatchResult:
    query_id: str
    best_id: Optional[str]
    best_score: float
    best_alignment: Iterable[str]
    method: str


def find_best_match(
    query_id: str,
    query_sequence: str,
    references: Dict[str, str],
    *,
    use_mmseqs: bool = True,
    fallback_to_biopython: bool = True,
) -> BestMatchResult:
    """Return the best-scoring reference for the given query sequence."""

    engine = SequenceAlignmentEngine()
    best_id: Optional[str] = None
    best_score = float("-inf")
    best_alignment = []
    method_used = ""

    if use_mmseqs:
        try:
            mmseqs_output = engine.align_pairwise_mmseqs(
                {
                    query_id: query_sequence,
                    **references,
                }
            )
            for line in mmseqs_output:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                ref_id = parts[1]
                if ref_id == query_id or ref_id not in references:
                    continue
                raw_score = float(parts[2])
                normalized = raw_score / max(len(query_sequence), 1)
                if normalized > best_score:
                    best_id = ref_id
                    best_score = normalized
                    best_alignment = []
                    method_used = "mmseqs"
        except Exception:
            if not fallback_to_biopython:
                raise

    if best_id is None and fallback_to_biopython:
        for ref_id, ref_seq in references.items():
            result = engine.align_pairwise(
                query_id,
                query_sequence,
                ref_id,
                ref_seq,
            )
            normalized = result.score / max(len(query_sequence), 1)
            if normalized > best_score:
                best_id = ref_id
                best_score = normalized
                best_alignment = result.alignment
                method_used = "biopython"

    return BestMatchResult(
        query_id=query_id,
        best_id=best_id,
        best_score=best_score,
        best_alignment=best_alignment,
        method=method_used or ("mmseqs" if use_mmseqs else "biopython"),
    )


def find_best_matches(
    sequences: Dict[str, str],
    references: Dict[str, str],
    *,
    use_mmseqs: bool = True,
    fallback_to_biopython: bool = True,
) -> Dict[str, BestMatchResult]:
    """Run best-match search for multiple sequences."""

    results: Dict[str, BestMatchResult] = {}
    for seq_id, sequence in sequences.items():
        results[seq_id] = find_best_match(
            seq_id,
            sequence,
            references,
            use_mmseqs=use_mmseqs,
            fallback_to_biopython=fallback_to_biopython,
        )
    return results
