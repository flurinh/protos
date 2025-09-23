#!/usr/bin/env python3
"""Demonstrate embeddings on GPCR structure-derived sequence datasets."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.processing.sequence import SequenceProcessor


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def load_gpcr_structure_sequences(
    preferred_dataset: str = "gpcr_chains_real",
) -> Tuple[str, Dict[str, str]]:
    seq_proc = SequenceProcessor()

    try:
        sequences = seq_proc.load_dataset(preferred_dataset)
        return preferred_dataset, sequences
    except (FileNotFoundError, KeyError):
        pass

    available = seq_proc.list_datasets()
    chain_datasets = [name for name in available if name.startswith("gpcr_chain_dataset_")]

    sequences: Dict[str, str] = {}
    for dataset_name in chain_datasets:
        try:
            dataset_sequences = seq_proc.load_dataset(dataset_name)
            if dataset_sequences:
                sequences.update(dataset_sequences)
        except Exception as exc:  # noqa: BLE001
            print(f"Warning: failed to load dataset '{dataset_name}': {exc}")

    if sequences:
        dataset_label = preferred_dataset if preferred_dataset in available else "gpcr_chain_dataset"
        return dataset_label, sequences

    raise RuntimeError(
        "No GPCR chain datasets available. Run test_cross_processor_annotation.py first "
        "to register structure-derived sequences."
    )


MODEL_ALIASES = {
    "esm2": "esm2_t12_35m",
    "esm2_small": "esm2_t12_35m",
    "esm2_medium": "esm2_t30_150m",
    "esm2_large": "esm2_t36_3b",
}


def resolve_model_selection(
    requested: Iterable[str] | None,
    *,
    all_models: Dict[str, Dict[str, str]],
    run_all: bool,
) -> List[str]:
    if run_all:
        return list(all_models.keys())

    targets = list(requested or ("ankh_large", "esm2"))
    resolved: List[str] = []

    for name in targets:
        canonical = MODEL_ALIASES.get(name, name)
        if canonical not in all_models:
            print(f"Skipping unknown model '{name}' (resolved '{canonical}')")
            continue
        resolved.append(canonical)

    if not resolved:
        raise ValueError(
            "No valid models selected. Available models: " + ", ".join(sorted(all_models.keys()))
        )

    seen = set()
    unique_models: List[str] = []
    for model in resolved:
        if model not in seen:
            unique_models.append(model)
            seen.add(model)
    return unique_models


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run embedding demos on GPCR chain sequences")
    parser.add_argument(
        "--models",
        nargs="+",
        help="Embedding models to run (aliases: esm2 -> esm2_t12_35m)",
    )
    parser.add_argument(
        "--all-models",
        action="store_true",
        help="Run all available embedding models",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_data_root()

    from protos.processing.embedding import EmbeddingProcessor

    dataset_name, sequences = load_gpcr_structure_sequences()
    print(f"Loaded dataset '{dataset_name}' with {len(sequences)} sequences")

    models = EmbeddingProcessor.available_models()

    try:
        target_models = resolve_model_selection(args.models, all_models=models, run_all=args.all_models)
    except ValueError as exc:
        print(exc)
        return

    for model_name in target_models:
        try:
            emb_proc = EmbeddingProcessor(model_name=model_name)
        except Exception as exc:  # noqa: BLE001
            print(f"Skipping {model_name}: {exc}")
            continue

        if not emb_proc.dependencies_available:
            print(f"Skipping {model_name}: install torch/transformers for embeddings")
            emb_proc.clear_cache()
            continue

        print(f"\n=== Embedding with {model_name} ===")
        for embedding_type in ("per_residue", "mean", "sum"):
            dataset_tag = f"{dataset_name}__{model_name}__{embedding_type}"
            try:
                emb_proc.embed_sequences(
                    sequences,
                    embedding_type=embedding_type,
                    save_dataset=dataset_tag,
                    register_entities=True,
                )
                print(f"  • Stored {embedding_type} embeddings -> {dataset_tag}")
            except Exception as exc:  # noqa: BLE001
                print(f"  ! Failed to embed with {embedding_type}: {exc}")

        emb_proc.clear_cache()


if __name__ == "__main__":
    main()
