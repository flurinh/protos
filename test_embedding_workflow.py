#!/usr/bin/env python3
"""Demonstrate embedding generation across all registered models."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


REFERENCE_SEQUENCES = {
    "3sn6_chain_A": "MKTIIALSYIFCLVFADYKDDDDAAAFVVVLG",
    "5d5a_chain_A": "MNTSVYIFCLVFADVTDKDNRTLLGFFVASLL",
    "6b73_chain_A": "MKSVLIFCLVFADYKDDDAAGGMVLLVFVVIL",
}


def ensure_sequence_dataset(dataset_name: str = "gpcr_sequences") -> None:
    from protos.processing.sequence import SequenceProcessor
    from protos.io.paths import get_protos_paths
    from protos.io.ingest.sequence_loader import SequenceLoader

    seq_proc = SequenceProcessor()

    try:
        seq_proc.load_dataset(dataset_name)
        return
    except Exception:
        pass

    paths = get_protos_paths()
    input_dir = Path(paths.get_processor_path('input'))
    input_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = input_dir / f"{dataset_name}.fasta"
    with open(fasta_path, "w") as handle:
        for seq_id, seq in REFERENCE_SEQUENCES.items():
            handle.write(f">{seq_id}\n{seq}\n")

    loader = SequenceLoader(processor=seq_proc)
    loader.download_and_register(
        str(fasta_path),
        name=dataset_name,
        materialize_entities=True,
    )


def main() -> None:
    ensure_data_root()
    ensure_sequence_dataset()

    from protos.processing.sequence import SequenceProcessor
    from protos.processing.embedding import EmbeddingProcessor

    seq_proc = SequenceProcessor()
    dataset_name = "gpcr_sequences"
    sequences = seq_proc.load_dataset(dataset_name)

    models = EmbeddingProcessor.available_models()

    for model_name in models.keys():
        try:
            emb_proc = EmbeddingProcessor(model_name=model_name)
        except Exception as exc:
            print(f"Skipping {model_name}: {exc}")
            continue

        if not emb_proc.dependencies_available:
            print(f"Skipping {model_name}: install torch/transformers for embeddings")
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
            except Exception as exc:
                print(f"  ! Failed to embed with {embedding_type}: {exc}")

        emb_proc.clear_cache()


if __name__ == "__main__":
    main()

