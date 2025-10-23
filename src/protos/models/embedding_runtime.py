"""Runtime embedding execution for ModelManager embedding cards.

This module provides a simple callable that the ConfigurableRuntimeAdapter can
invoke to compute embeddings for a sequence dataset using EmbeddingProcessor.

The callable writes a compact artifact (NPZ + JSON metadata) into the provided
outputs directory and returns a structure compatible with RuntimeResult.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import json
import numpy as np

from protos.processing.embedding import EmbeddingProcessor
from protos.models.model_specs import ArtifactBundle, ArtifactSpec


def _infer_model_from_card(card_name: str) -> str:
    """Map a card name like 'embedding_esm2_t12_35m' -> 'esm2_t12_35m'."""
    if card_name.startswith("embedding_"):
        return card_name[len("embedding_") :]
    return card_name


def run_embedding(
    *,
    card,  # ModelCard
    request,  # ModelRequest
    inputs,  # List[ArtifactBundle]
    outputs_dir: str,  # Provided by ConfigurableRuntimeAdapter
    paths=None,
) -> Dict[str, Any]:
    debug = bool(getattr(request, "config", {}).get("debug", False))
    # Resolve sequences from the provided dataset artifact
    seq_bundle = next((b for b in inputs if b.spec.name == "sequence_dataset"), None)
    if seq_bundle is None:
        raise ValueError("Embedding runtime requires a 'sequence_dataset' input")

    sequences: Dict[str, str] = seq_bundle.metadata.get("sequences", {})
    dataset_name: str = seq_bundle.metadata.get("dataset") or "sequences"
    if not sequences:
        raise ValueError("No sequences resolved for embedding runtime")

    # Determine model and embedding type from card/config
    model_name = _infer_model_from_card(card.name)
    embedding_type = str(request.config.get("embedding_type", "per_residue"))

    if debug:
        print(f"[embed-runtime] dataset={dataset_name} sequences={len(sequences)}")

    # Try to initialize the embedding processor lazily
    try:
        emb = EmbeddingProcessor(model_name=model_name)
    except Exception as exc:  # noqa: BLE001
        # Surface a structured failure so tests can skip gracefully
        return {
            "outputs": {"status": "error", "error": str(exc)},
            "metadata": {
                "model_name": model_name,
                "embedding_type": embedding_type,
                "dataset": dataset_name,
            },
            "artifacts": [],
        }

    if not emb.dependencies_available:
        return {
            "outputs": {"status": "skipped", "reason": "missing_dependencies"},
            "metadata": {
                "model_name": model_name,
                "embedding_type": embedding_type,
                "dataset": dataset_name,
            },
            "artifacts": [],
        }

    # Compute embeddings in-memory without registering a dataset
    try:
        tensors = emb.embed_sequences(
            sequences,
            embedding_type=embedding_type,
            save_dataset=None,
            register_entities=False,
        )
    except Exception as exc:  # noqa: BLE001
        return {
            "outputs": {"status": "error", "error": str(exc)},
            "metadata": {
                "model_name": model_name,
                "embedding_type": embedding_type,
                "dataset": dataset_name,
            },
            "artifacts": [],
        }

    if debug:
        print(
            f"[embed-runtime] computed embeddings model={model_name} type={embedding_type} n={len(tensors)}"
        )

    # Convert to numpy and persist a compact artifact
    npz_path = Path(outputs_dir) / f"{dataset_name}__{model_name}__{embedding_type}.npz"
    arrays: Dict[str, np.ndarray] = {}
    shapes: Dict[str, List[int]] = {}
    for key, tensor in tensors.items():
        if hasattr(tensor, "detach"):
            arr = tensor.detach().cpu().numpy()
        else:
            arr = np.asarray(tensor)
        arrays[key] = arr
        shapes[key] = list(arr.shape)
    # Save NPZ
    np.savez_compressed(npz_path, **arrays)

    # Write sidecar metadata JSON
    meta = {
        "model_name": model_name,
        "embedding_type": embedding_type,
        "dataset": dataset_name,
        "entity_count": len(arrays),
        "sequence_ids": list(arrays.keys()),
        "shapes_preview": shapes,
        "artifact_path": str(npz_path.relative_to(paths.data_root)) if paths else str(npz_path),
    }
    with open(Path(outputs_dir) / "embedding_meta.json", "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=2)

    artifact = ArtifactBundle(
        spec=ArtifactSpec(
            name="embedding_artifact",
            kind="embedding",
            provider="embedding_runtime",
            format="npz",
        ),
        path=npz_path,
        metadata=meta,
    )

    return {
        "outputs": {
            "status": "ok",
            "entities": len(arrays),
        },
        "artifacts": [artifact],
        "metadata": meta,
    }
