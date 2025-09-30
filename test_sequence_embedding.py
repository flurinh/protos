"""Tests for EmbeddingProcessor covering encoder-only and encoder-decoder models."""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import pytest

import protos
from protos.processing.embedding import EmbeddingProcessor


@pytest.fixture(scope="session")
def data_root(tmp_path_factory: pytest.TempPathFactory) -> Path:
    path = tmp_path_factory.mktemp("protos_embedding")
    protos.set_data_path(str(path))
    return path


@pytest.fixture(scope="module")
def sample_sequences() -> Dict[str, str]:
    return {
        "opsin1": "MNKLLVIATAVCLLFSYYQYQGVG",
        "opsin2": "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEP"[:25],
    }


@pytest.mark.parametrize("model_name", ["esm2_t6_8m", "ankh_base"])
@pytest.mark.parametrize("embedding_type", ["mean", "sum", "cls", "per_residue"])
@pytest.mark.parametrize("register", [False, True])
def test_embedding_processor_outputs_expected_shapes(
    data_root: Path,
    sample_sequences: Dict[str, str],
    model_name: str,
    embedding_type: str,
    register: bool,
) -> None:
    processor = EmbeddingProcessor(model_name=model_name, batch_size=2)
    if not processor.dependencies_available:
        pytest.skip("Optional embedding dependencies (torch, transformers) missing.")

    try:
        outputs = processor.embed_sequences(
            sample_sequences,
            embedding_type=embedding_type,
            save_dataset=None,
            register_entities=register,
        )
    except OSError as error:
        pytest.skip(f"Model {model_name} unavailable locally: {error}")

    assert set(outputs) == set(sample_sequences)
    for seq_id, tensor in outputs.items():
        assert tensor is not None, f"No tensor returned for {seq_id}"
        if embedding_type == "per_residue":
            expected_len = len(sample_sequences[seq_id])
            assert tensor.shape[0] >= expected_len
        else:
            assert tensor.ndim == 1, f"Expected 1D tensor for {embedding_type}, got {tensor.ndim}D"

    processor.clear_cache()
