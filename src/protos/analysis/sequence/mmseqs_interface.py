"""Wrapper helpers for MMseqs alignments."""

from __future__ import annotations

from typing import Dict, List

from protos.processing.sequence.seq_alignment import mmseqs2_align2


class MMseqsUnavailableError(RuntimeError):
    """Raised when MMseqs binaries are unavailable or execution fails."""


def run_mmseqs_alignment(sequences: Dict[str, str]) -> List[str]:
    try:
        return mmseqs2_align2(sequences, sequences)
    except FileNotFoundError as exc:
        raise MMseqsUnavailableError("MMseqs binary not found") from exc
    except Exception as exc:
        raise MMseqsUnavailableError(str(exc)) from exc
