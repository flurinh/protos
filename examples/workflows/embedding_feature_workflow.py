#!/usr/bin/env python3
"""Embedding feature extraction workflow with dependency awareness."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import protos
from protos.io.formats.fasta_utils import read_fasta
from protos.processing.embedding import EmbeddingProcessor

DATA_RELATIVE_ROOT = Path(__file__).resolve().parents[2] / "data"
SEQUENCE_DATASET = "gpcr_agonist_antagonist_sequences"
SEQUENCE_FASTA = DATA_RELATIVE_ROOT / "sequence" / "fasta" / f"{SEQUENCE_DATASET}.fasta"
OUTPUT_DIR = DATA_RELATIVE_ROOT / "reports"
SUMMARY_PATH = OUTPUT_DIR / "embedding_feature_summary.json"
MODEL_NAME = "esm2_t6_8m"
MAX_SEQUENCES = 3


@dataclass
class EmbeddingSummary:
    model: str
    dependencies_available: bool
    sequence_count: int
    embedding_type: str
    vector_length: int


def ensure_data_root() -> Path:
    DATA_RELATIVE_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_RELATIVE_ROOT))
    return DATA_RELATIVE_ROOT


def load_sequences() -> Dict[str, str]:
    if not SEQUENCE_FASTA.exists():
        raise FileNotFoundError(f"Expected FASTA file missing: {SEQUENCE_FASTA}")
    sequences = read_fasta(str(SEQUENCE_FASTA))
    return dict(list(sequences.items())[:MAX_SEQUENCES])


def main() -> None:
    ensure_data_root()
    subset = load_sequences()

    embedding_processor = EmbeddingProcessor(model_name=MODEL_NAME, batch_size=MAX_SEQUENCES)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not embedding_processor.dependencies_available:
        payload: Dict[str, object] = {
            "model": MODEL_NAME,
            "dependencies_available": False,
            "message": "Install torch and transformers to enable embeddings.",
            "sequence_subset": list(subset.keys()),
        }
        SUMMARY_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(json.dumps(payload, indent=2))
        return

    embeddings = embedding_processor.embed_sequences(
        subset,
        embedding_type="mean",
        register_entities=False,
    )

    first_vector = next(iter(embeddings.values()))
    vector_length = int(first_vector.shape[-1]) if hasattr(first_vector, "shape") else len(first_vector)
    summary = EmbeddingSummary(
        model=MODEL_NAME,
        dependencies_available=True,
        sequence_count=len(embeddings),
        embedding_type="mean",
        vector_length=vector_length,
    )

    payload = {
        "embedding": summary.__dict__,
        "sequences": list(embeddings.keys()),
    }
    SUMMARY_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
