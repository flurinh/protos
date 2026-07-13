"""Offline contract tests for the pinned Biohub ESMC-6B adapter."""

from __future__ import annotations

import json
import os
import shutil
import urllib.request
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import numpy as np
import pytest

from protos.processing.embedding import (
    ESMCArtifactError,
    ESMCBatchOutput,
    ESMCLoadPolicy,
    ESMCModelProvenance,
    ESMCOutputProvenance,
    ESMCSequenceConflictError,
    ESMCShapeError,
    ESMCSourceLineage,
    ESMCTransformersProvenance,
    ESMCTokenMapping,
    ESMCTokenPolicy,
    ESMCTokenPolicyError,
    ESMCTruncationError,
    EmbeddingProcessor,
    sequence_sha256,
)
from protos.processing.embedding.esmc_adapter import ESMCEmbeddingAdapter


class FakeESMCBackend:
    """Return deterministic padded token embeddings without model downloads."""

    def __init__(self, mode: str = "valid", dimension: int = 2560) -> None:
        self.mode = mode
        self.dimension = dimension
        self.calls: list[tuple[str, ...]] = []
        self.load_policies: list[ESMCLoadPolicy] = []

    def encode(
        self,
        sequences: Sequence[str],
        *,
        model: ESMCModelProvenance,
        transformers: ESMCTransformersProvenance,
        load_policy: ESMCLoadPolicy,
        token_policy: ESMCTokenPolicy,
        device: str,
    ) -> ESMCBatchOutput:
        assert transformers == ESMCTransformersProvenance()
        del model, token_policy, device
        self.calls.append(tuple(sequences))
        self.load_policies.append(load_policy)
        width = max(len(sequence) + 2 for sequence in sequences)
        dimension = self.dimension - 1 if self.mode == "wrong_dimension" else self.dimension
        embeddings = np.zeros((len(sequences), width, dimension), dtype=np.float16)
        input_ids = np.ones((len(sequences), width), dtype=np.int64)
        attention = np.zeros((len(sequences), width), dtype=np.int64)
        special = np.ones((len(sequences), width), dtype=np.int64)

        for batch_index, sequence in enumerate(sequences):
            active = len(sequence) + 2
            input_ids[batch_index, 0] = 0
            input_ids[batch_index, 1 : active - 1] = np.arange(4, 4 + len(sequence))
            input_ids[batch_index, active - 1] = 2
            attention[batch_index, :active] = 1
            special[batch_index, 1 : active - 1] = 0
            base = sum(sequence.encode("ascii")) % 101
            for token_index in range(active):
                embeddings[batch_index, token_index, :] = base + token_index

        mapping = ESMCTokenMapping(0, 2, 1, 3)
        if self.mode == "missing_mapping":
            mapping = replace(mapping, pad_token_id=None)
        elif self.mode == "ambiguous_mapping":
            mapping = replace(mapping, eos_token_id=0)
        elif self.mode == "wrong_bos":
            input_ids[:, 0] = 7
        elif self.mode == "ambiguous_special_mask":
            special[:, 1] = 1
        elif self.mode == "wrong_mask_shape":
            special = special[:, :-1]
        elif self.mode == "non_float":
            embeddings = embeddings.astype(np.int32)

        return ESMCBatchOutput(
            embeddings=embeddings,
            input_ids=input_ids,
            attention_mask=attention,
            special_tokens_mask=special,
            token_mapping=mapping,
            truncated=tuple(self.mode == "truncated" for _ in sequences),
            output_field="embeddings" if self.mode == "wrong_output" else "last_hidden_state",
            final_layer=self.mode != "not_final",
        )


@pytest.fixture
def isolated_processor(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    import protos.io.core as core
    import protos.io.paths.path_config as path_config
    from protos.io.core.entity_registry import EntityRegistry

    monkeypatch.setenv("PROTOS_DATA_ROOT", str(tmp_path / "protos-data"))
    path_config._paths_instance = None
    path_config.ProtosPaths._instance = None
    path_config.ProtosPaths._initialized = False
    core._registry_instance = None
    EntityRegistry._instance = None

    def factory(*, batch_size: int = 2, device: str = "cpu") -> EmbeddingProcessor:
        return EmbeddingProcessor(
            name="esmc-test",
            model_name="esmc_6b",
            device=device,
            batch_size=batch_size,
        )

    yield factory

    path_config._paths_instance = None
    path_config.ProtosPaths._instance = None
    path_config.ProtosPaths._initialized = False
    core._registry_instance = None
    EntityRegistry._instance = None


def test_registry_and_two_sequence_padded_batch_contract(isolated_processor) -> None:
    processor = isolated_processor(batch_size=2)
    backend = FakeESMCBackend()
    lineage = {
        "long": ESMCSourceLineage(
            "fasta", "lambda-v2.fasta", "27ec586", "a" * 64
        ),
        "short": ESMCSourceLineage("fasta", "lambda-v2.fasta", "27ec586", "a" * 64),
    }

    results = processor.embed_esmc_sequences(
        [("long", "ACDE"), ("short", "m n")],
        backend=backend,
        source_lineage=lineage,
    )

    assert list(results) == ["long", "short"]
    assert backend.calls == [("ACDE", "MN")]
    assert backend.load_policies == [ESMCLoadPolicy()]
    assert results["long"].embedding.shape == (4, 2560)
    assert results["short"].embedding.shape == (2, 2560)
    assert results["long"].embedding.dtype == np.float32
    long_base = sum(b"ACDE") % 101
    np.testing.assert_array_equal(
        results["long"].embedding[:, 0],
        np.asarray([long_base + 1, long_base + 2, long_base + 3, long_base + 4]),
    )

    metadata = results["long"].metadata
    assert metadata["model_id"] == "biohub/ESMC-6B"
    assert metadata["model_revision"] == "45b0fa5d7fb06faefbd5e3b89bdcef35d564e79a"
    assert metadata["code_revision"] == "ba4d7124864eed323a93bf3cfefcd958f573b75a"
    assert metadata["model"]["layer_count"] == 80
    assert metadata["embedding_dimension"] == 2560
    assert metadata["dtype"] == "float32"
    assert metadata["inference_dtype"] == "bfloat16"
    assert metadata["load_policy"]["cuda_max_memory"] == "28GiB"
    assert metadata["token_policy"]["truncation"] is False
    assert metadata["source_lineage"] == {
        "source_kind": "fasta",
        "source_id": "lambda-v2.fasta",
        "source_revision": "27ec586",
        "source_sha256": "a" * 64,
    }
    assert metadata["sequence_sha256"] == sequence_sha256("ACDE")
    assert metadata["shape"] == [4, 2560]
    assert "esmc_6b" in results["long"].artifact_path.parts

    model = EmbeddingProcessor.available_models()["esmc_6b"]
    assert model["embedding_dim"] == 2560
    assert model["layers"] == 80
    assert model["inference_dtype"] == "bfloat16"
    assert model["cuda_device_map"] == "auto"
    assert model["cuda_max_memory"] == "28GiB"
    assert EmbeddingProcessor.available_models()["esm2_t6_8m"]["embedding_dim"] == 320
    assert EmbeddingProcessor.available_models()["ankh_large"]["embedding_dim"] == 1536


def test_resume_recompute_ordering_and_determinism(isolated_processor) -> None:
    processor = isolated_processor(batch_size=2)
    sequences = [("b", "MNP"), ("a", "AC")]
    first_backend = FakeESMCBackend()
    first = processor.embed_esmc_sequences(sequences, backend=first_backend)

    resume_backend = FakeESMCBackend(mode="wrong_dimension")
    resumed = processor.embed_esmc_sequences(sequences, backend=resume_backend)
    assert resume_backend.calls == []
    assert list(resumed) == ["b", "a"]
    assert all(result.reused for result in resumed.values())
    np.testing.assert_array_equal(first["a"].embedding, resumed["a"].embedding)

    processor.batch_size = 1
    recompute_backend = FakeESMCBackend()
    recomputed = processor.embed_esmc_sequences(
        sequences,
        backend=recompute_backend,
        force_recompute=True,
    )
    assert recompute_backend.calls == [("MNP",), ("AC",)]
    assert not any(result.reused for result in recomputed.values())
    assert [result.artifact_path for result in recomputed.values()] == [
        result.artifact_path for result in first.values()
    ]
    for entity_id in ("b", "a"):
        np.testing.assert_array_equal(
            first[entity_id].embedding, recomputed[entity_id].embedding
        )


def test_rejects_input_and_stored_sequence_conflicts(isolated_processor) -> None:
    processor = isolated_processor()
    backend = FakeESMCBackend()
    with pytest.raises(ESMCSequenceConflictError, match="conflicting sequences"):
        processor.embed_esmc_sequences(
            [("same", "ACD"), ("same", "ACE")], backend=backend
        )
    assert backend.calls == []

    processor.embed_esmc_sequences({"same": "ACD"}, backend=backend)
    with pytest.raises(ESMCSequenceConflictError, match="conflicting sequence"):
        processor.embed_esmc_sequences({"same": "ACE"}, backend=backend)


def test_rejects_policy_and_backend_truncation(isolated_processor) -> None:
    processor = isolated_processor()
    with pytest.raises(ESMCTruncationError, match="policy maximum"):
        processor.embed_esmc_sequences(
            {"too-long": "ACDE"},
            backend=FakeESMCBackend(),
            token_policy=ESMCTokenPolicy(max_residues=3),
        )
    with pytest.raises(ESMCTruncationError, match="reported truncated"):
        processor.embed_esmc_sequences(
            {"truncated": "ACD"}, backend=FakeESMCBackend("truncated")
        )


@pytest.mark.parametrize(
    "mode,exception",
    [
        ("missing_mapping", ESMCTokenPolicyError),
        ("ambiguous_mapping", ESMCTokenPolicyError),
        ("wrong_bos", ESMCTokenPolicyError),
        ("ambiguous_special_mask", ESMCTokenPolicyError),
        ("wrong_mask_shape", ESMCShapeError),
        ("wrong_dimension", ESMCShapeError),
        ("non_float", ESMCShapeError),
        ("wrong_output", ESMCShapeError),
        ("not_final", ESMCShapeError),
    ],
)
def test_rejects_token_and_output_contract_violations(
    isolated_processor, mode: str, exception: type[Exception]
) -> None:
    processor = isolated_processor()
    with pytest.raises(exception):
        processor.embed_esmc_sequences(
            {f"bad-{mode}": "ACD"}, backend=FakeESMCBackend(mode)
        )


def test_exact_reuse_rejects_tampered_array(isolated_processor) -> None:
    processor = isolated_processor()
    result = processor.embed_esmc_sequences(
        {"tamper": "ACD"}, backend=FakeESMCBackend()
    )["tamper"]
    with np.load(result.artifact_path, allow_pickle=False) as archive:
        array = np.asarray(archive["embedding"])
        metadata = str(archive["metadata"].item())
    array[0, 0] += 1
    np.savez_compressed(result.artifact_path, embedding=array, metadata=np.asarray(metadata))

    with pytest.raises(ESMCArtifactError, match="content hash mismatch"):
        processor.embed_esmc_sequences(
            {"tamper": "ACD"}, backend=FakeESMCBackend()
        )


def _direct_adapter(
    processor: EmbeddingProcessor,
    backend: FakeESMCBackend,
    *,
    model: ESMCModelProvenance = ESMCModelProvenance(),
    token_policy: ESMCTokenPolicy = ESMCTokenPolicy(),
    output: ESMCOutputProvenance = ESMCOutputProvenance(),
) -> ESMCEmbeddingAdapter:
    return ESMCEmbeddingAdapter(
        embeddings_dir=processor.embeddings_dir,
        data_root=Path(processor.paths.data_root),
        entity_registry=processor.entity_registry,
        batch_size=processor.batch_size,
        device=processor.device,
        backend=backend,
        model_provenance=model,
        token_policy=token_policy,
        output_provenance=output,
    )


def test_no_reuse_across_ankh_600m_revision_adapter_or_policy(
    isolated_processor,
) -> None:
    processor = isolated_processor()
    ankh_dir = processor.embeddings_dir / "ankh_large" / "_entities"
    ankh_dir.mkdir(parents=True)
    (ankh_dir / "protein.pkl").write_bytes(b"legacy-ankh-artifact")

    old_600m_model = replace(
        ESMCModelProvenance(),
        model_id="biohub/ESMC-600M",
        model_revision="a7e82012c83126b9eedb055fea9fa84b6c02f094",
        layer_count=36,
        embedding_dimension=1152,
    )
    with pytest.raises(ESMCShapeError):
        _direct_adapter(processor, FakeESMCBackend(dimension=1152), model=old_600m_model,
                        output=replace(ESMCOutputProvenance(), dimension=1152))

    production_backend = FakeESMCBackend()
    production = processor.embed_esmc_sequences(
        {"protein": "ACD"}, backend=production_backend
    )["protein"]
    assert production_backend.calls == [("ACD",)]
    assert "esmc_6b" in production.artifact_path.parts
    assert "ankh_large" not in production.artifact_path.parts

    variants = [
        (
            replace(
                ESMCModelProvenance(),
                model_revision="0000000000000000000000000000000000000000",
            ),
            ESMCTokenPolicy(),
        ),
        (
            replace(ESMCModelProvenance(), adapter_id="foreign.esmc.adapter.v1"),
            ESMCTokenPolicy(),
        ),
        (ESMCModelProvenance(), ESMCTokenPolicy(max_residues=2045)),
    ]
    for index, (model, policy) in enumerate(variants):
        backend = FakeESMCBackend()
        result = _direct_adapter(
            processor, backend, model=model, token_policy=policy
        ).embed({"protein": "ACD"})["protein"]
        assert backend.calls == [("ACD",)], index
        assert result.artifact_path != production.artifact_path

    cuda_processor = isolated_processor(device="cuda")
    cuda_result = cuda_processor.embed_esmc_sequences(
        {"protein": "ACD"}, backend=FakeESMCBackend()
    )["protein"]
    assert cuda_result.artifact_path != production.artifact_path
    assert cuda_result.metadata["execution_device"] == "cuda"
    assert production.metadata["execution_device"] == "cpu"


def test_generic_embedding_path_is_not_used_for_esmc(isolated_processor) -> None:
    processor = isolated_processor()
    with pytest.raises(ValueError, match="embed_esmc_sequences"):
        processor.embed_sequences("ACD", register_entities=False)


def test_backend_typeerror_propagates_without_provenance_retry(isolated_processor) -> None:
    processor = isolated_processor()

    class BrokenBackend(FakeESMCBackend):
        def __init__(self):
            super().__init__()
            self.invocations = 0

        def encode(self, *args, **kwargs):
            self.invocations += 1
            assert kwargs["transformers"] == ESMCTransformersProvenance()
            raise TypeError("backend bug")

    backend = BrokenBackend()
    with pytest.raises(TypeError, match="backend bug"):
        processor.embed_esmc_sequences({"protein": "ACD"}, backend=backend)
    assert backend.invocations == 1


def _available_memory_bytes() -> int:
    meminfo = Path("/proc/meminfo")
    if not meminfo.exists():
        return 0
    for line in meminfo.read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024
    return 0


@pytest.mark.integration
@pytest.mark.network
@pytest.mark.slow
@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_live_esmc_6b_two_sequence_smoke(isolated_processor, device: str) -> None:
    """Opt-in 25.4-GB download; never runs as an ordinary test."""

    if os.environ.get("PROTOS_RUN_ESMC_6B_LIVE") != "1":
        pytest.skip("set PROTOS_RUN_ESMC_6B_LIVE=1 to enable the 25.4-GB live smoke")
    try:
        import torch
        import transformers  # noqa: F401
        import accelerate  # noqa: F401
    except ImportError as exc:
        pytest.skip(f"live ESMC-6B dependency unavailable: {exc}")
    if device == "cuda" and not torch.cuda.is_available():
        pytest.skip("CUDA live smoke requested but torch.cuda.is_available() is false")
    if device == "cuda":
        free_bytes, _total_bytes = torch.cuda.mem_get_info()
        if free_bytes < 28 * 1024**3:
            pytest.skip(f"CUDA live smoke needs 28 GiB free; found {free_bytes / 1024**3:.1f} GiB")
    elif _available_memory_bytes() < 32 * 1024**3:
        pytest.skip("CPU live smoke needs at least 32 GiB available RAM")
    free_disk = shutil.disk_usage(Path.home()).free
    if free_disk < 35 * 1024**3:
        pytest.skip(f"live ESMC-6B cache needs 35 GiB free; found {free_disk / 1024**3:.1f} GiB")
    try:
        with urllib.request.urlopen(
            "https://huggingface.co/biohub/ESMC-6B/resolve/"
            "45b0fa5d7fb06faefbd5e3b89bdcef35d564e79a/config.json",
            timeout=10,
        ) as response:
            config = json.load(response)
    except Exception as exc:
        pytest.skip(f"pinned Hugging Face config is unreachable: {exc}")
    if config.get("d_model") != 2560 or config.get("n_layers") != 80:
        pytest.fail("pinned live config is not the approved 80-layer d_model=2560 model")

    processor = isolated_processor(batch_size=2, device=device)
    results = processor.embed_esmc_sequences(
        {"smoke-a": "ACDEFGHIK", "smoke-b": "MNPQRSTVWY"},
        source_lineage={
            key: ESMCSourceLineage("live_smoke", "RR-20260713-1704")
            for key in ("smoke-a", "smoke-b")
        },
        force_recompute=True,
    )
    assert results["smoke-a"].embedding.shape == (9, 2560)
    assert results["smoke-b"].embedding.shape == (10, 2560)
    assert all(result.embedding.dtype == np.float32 for result in results.values())
