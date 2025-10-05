#!/usr/bin/env python3
"""Sequence and GRN bootstrap workflow using real Protos datasets."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean
from typing import Dict

import protos
from protos.io.formats.fasta_utils import read_fasta
from protos.processing.sequence import SequenceProcessor

DATA_RELATIVE_ROOT = Path(__file__).resolve().parents[2] / "data"
SEQUENCE_DATASET = "gpcr_agonist_antagonist_sequences"
SEQUENCE_FASTA = DATA_RELATIVE_ROOT / "sequence" / "fasta" / f"{SEQUENCE_DATASET}.fasta"
REFERENCE_TABLE = "gpcrdb_ref"
PROTEIN_FAMILY = "gpcr_a"
OUTPUT_DIR = DATA_RELATIVE_ROOT / "reports"
SUMMARY_PATH = OUTPUT_DIR / "sequence_grn_bootstrap.json"


def ensure_data_root() -> Path:
    DATA_RELATIVE_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_RELATIVE_ROOT))
    return DATA_RELATIVE_ROOT


@dataclass
class SequenceStats:
    dataset: str
    sequence_count: int
    min_length: int
    max_length: int
    average_length: float


@dataclass
class AnnotationStats:
    reference_table: str
    protein_family: str
    annotated_sequences: int
    total_sequences: int
    mean_coverage: float
    min_coverage: float
    max_coverage: float


def summarize_sequences(sequences: Dict[str, str]) -> SequenceStats:
    lengths = [len(seq) for seq in sequences.values()]
    return SequenceStats(
        dataset=SEQUENCE_DATASET,
        sequence_count=len(sequences),
        min_length=min(lengths) if lengths else 0,
        max_length=max(lengths) if lengths else 0,
        average_length=mean(lengths) if lengths else 0.0,
    )


def annotate_with_grn(sequence_processor: SequenceProcessor, sequences: Dict[str, str]) -> AnnotationStats:
    dataframe, summary = sequence_processor.annotate_with_grn(
        sequences=sequences,
        reference_table=REFERENCE_TABLE,
        protein_family=PROTEIN_FAMILY,
        output_table="gpcr_sequence_grn_bootstrap",
        return_summary=True,
        allow_create=True,
        metadata={"source_fasta": str(SEQUENCE_FASTA.relative_to(DATA_RELATIVE_ROOT))},
    )

    coverages = [
        info.get("coverage", 0.0)
        for info in summary.get("per_sequence", {}).values()
    ]

    return AnnotationStats(
        reference_table=REFERENCE_TABLE,
        protein_family=PROTEIN_FAMILY,
        annotated_sequences=sum(1 for value in coverages if value > 0.0),
        total_sequences=len(coverages),
        mean_coverage=mean(coverages) if coverages else 0.0,
        min_coverage=min(coverages) if coverages else 0.0,
        max_coverage=max(coverages) if coverages else 0.0,
    )


def main() -> None:
    ensure_data_root()

    if not SEQUENCE_FASTA.exists():
        raise FileNotFoundError(f"Expected FASTA file missing: {SEQUENCE_FASTA}")

    sequences = read_fasta(str(SEQUENCE_FASTA))
    sequence_stats = summarize_sequences(sequences)

    sequence_processor = SequenceProcessor()
    annotation_stats = annotate_with_grn(sequence_processor, sequences)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    result = {
        "sequence": asdict(sequence_stats),
        "annotation": asdict(annotation_stats),
    }
    SUMMARY_PATH.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
