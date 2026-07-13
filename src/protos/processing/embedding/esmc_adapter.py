"""Pinned Biohub ESMC residue embedding adapter.

The adapter is deliberately independent of model/job orchestration.  It owns
token validation, output validation, and provenance-complete artifact reuse for
the one ESMC model approved by :class:`EmbeddingProcessor`.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Protocol, Sequence, Tuple, Union

import numpy as np


ESMC_MODEL_ID = "biohub/ESMC-6B"
ESMC_MODEL_REVISION = "45b0fa5d7fb06faefbd5e3b89bdcef35d564e79a"
ESMC_CODE_REVISION = "ba4d7124864eed323a93bf3cfefcd958f573b75a"
ESMC_EMBEDDING_DIMENSION = 2560
ESMC_LAYER_COUNT = 80
ESMC_TRANSFORMERS_REPOSITORY = "https://github.com/Biohub/transformers.git"
ESMC_TRANSFORMERS_REVISION = "ef32577f55da19a4989cd7b22e004dc43a4998cb"


class ESMCValidationError(ValueError):
    """Base class for ESMC input, token, output, and artifact failures."""


class ESMCTruncationError(ESMCValidationError):
    """A sequence or token batch was truncated."""


class ESMCTokenPolicyError(ESMCValidationError):
    """Special-token identity or placement is missing or ambiguous."""


class ESMCShapeError(ESMCValidationError):
    """The backend output does not have the required residue shape."""


class ESMCSequenceConflictError(ESMCValidationError):
    """One entity ID refers to different normalized sequences."""


class ESMCArtifactError(ESMCValidationError):
    """A stored artifact fails exact provenance or content validation."""


@dataclass(frozen=True)
class ESMCModelProvenance:
    """Immutable model and implementation lineage."""

    adapter_id: str = "protos.embedding.esmc.v1"
    model_id: str = ESMC_MODEL_ID
    model_revision: str = ESMC_MODEL_REVISION
    code_revision: str = ESMC_CODE_REVISION
    layer_count: int = ESMC_LAYER_COUNT
    embedding_dimension: int = ESMC_EMBEDDING_DIMENSION

    def __post_init__(self) -> None:
        expected = {
            "adapter_id": "protos.embedding.esmc.v1",
            "model_id": ESMC_MODEL_ID,
            "model_revision": ESMC_MODEL_REVISION,
            "code_revision": ESMC_CODE_REVISION,
            "layer_count": ESMC_LAYER_COUNT,
            "embedding_dimension": ESMC_EMBEDDING_DIMENSION,
        }
        for field, value in expected.items():
            if getattr(self, field) != value:
                raise ESMCValidationError(
                    f"ESMC canonical {field} must remain {value!r}"
                )


@dataclass(frozen=True)
class ESMCTransformersProvenance:
    """Exact Biohub Transformers runtime required by the ESMC loader."""

    distribution: str = "transformers"
    repository: str = ESMC_TRANSFORMERS_REPOSITORY
    revision: str = ESMC_TRANSFORMERS_REVISION
    install_contract: str = "pep508-vcs-direct-url-v1"

    def __post_init__(self) -> None:
        if self.distribution != "transformers":
            raise ESMCValidationError("ESMC runtime distribution must be transformers")
        if self.repository != ESMC_TRANSFORMERS_REPOSITORY:
            raise ESMCValidationError("ESMC requires the pinned Biohub Transformers fork")
        if self.revision != ESMC_TRANSFORMERS_REVISION:
            raise ESMCValidationError("ESMC requires the exact pinned Transformers commit")
        if self.install_contract != "pep508-vcs-direct-url-v1":
            raise ESMCValidationError("Unsupported Transformers installation contract")


@dataclass(frozen=True)
class ESMCLoadPolicy:
    """Pinned inference precision and conservative device-placement policy."""

    policy_id: str = "esmc-6b-bfloat16-device-map-v1"
    inference_dtype: str = "bfloat16"
    cuda_device_map: str = "auto"
    cuda_max_memory: str = "28GiB"
    cpu_max_memory: str = "64GiB"
    cpu_offload: bool = True
    low_cpu_mem_usage: bool = True

    def __post_init__(self) -> None:
        if self.inference_dtype != "bfloat16":
            raise ESMCValidationError("Pinned ESMC-6B inference dtype is bfloat16")
        if self.cuda_device_map != "auto":
            raise ESMCValidationError("Pinned ESMC-6B CUDA device_map policy is 'auto'")


@dataclass(frozen=True)
class ESMCTokenPolicy:
    """Tokenization policy included in every artifact identity."""

    policy_id: str = "esmc-residue-bos-eos-padding-v1"
    add_special_tokens: bool = True
    padding: str = "longest"
    truncation: bool = False
    max_residues: int = 2046
    residue_selection: str = "attention_mask=1-and-special_tokens_mask=0"
    require_bos: bool = True
    require_eos: bool = True
    require_pad_mapping: bool = True
    reject_unknown_tokens: bool = True

    def __post_init__(self) -> None:
        if self.truncation:
            raise ESMCTokenPolicyError("ESMC truncation must remain disabled")
        if not self.add_special_tokens:
            raise ESMCTokenPolicyError("ESMC BOS/EOS tokens are required")
        if self.padding != "longest":
            raise ESMCTokenPolicyError("ESMC batches must use longest-sequence padding")
        if not (1 <= self.max_residues <= 2046):
            raise ESMCTokenPolicyError("max_residues must be between 1 and 2046")
        if not (self.require_bos and self.require_eos and self.require_pad_mapping):
            raise ESMCTokenPolicyError("BOS, EOS, and padding validation is mandatory")
        if not self.reject_unknown_tokens:
            raise ESMCTokenPolicyError("Unknown residue-token rejection is mandatory")
        if self.residue_selection != "attention_mask=1-and-special_tokens_mask=0":
            raise ESMCTokenPolicyError("The pinned residue-selection policy cannot change")


@dataclass(frozen=True)
class ESMCOutputProvenance:
    """Required final-layer output contract."""

    output_field: str = "last_hidden_state"
    layer: str = "final"
    embedding_type: str = "per_residue"
    dimension: int = ESMC_EMBEDDING_DIMENSION
    storage_dtype: str = "float32"

    def __post_init__(self) -> None:
        if self.output_field != "last_hidden_state" or self.layer != "final":
            raise ESMCShapeError("ESMC output must be final last_hidden_state")
        if self.embedding_type != "per_residue" or self.dimension != ESMC_EMBEDDING_DIMENSION:
            raise ESMCShapeError("ESMC output must be per-residue dimension 2560")
        if self.storage_dtype != "float32":
            raise ESMCShapeError("ESMC artifacts must use float32 storage")


@dataclass(frozen=True)
class ESMCSourceLineage:
    """Caller-supplied source identity retained in artifact provenance."""

    source_kind: str
    source_id: str
    source_revision: Optional[str] = None
    source_sha256: Optional[str] = None

    def __post_init__(self) -> None:
        if not self.source_kind or not self.source_id:
            raise ESMCValidationError("Source lineage kind and ID must be non-empty")
        if self.source_sha256 is not None and not re.fullmatch(
            r"[0-9a-f]{64}", self.source_sha256
        ):
            raise ESMCValidationError("Source lineage SHA-256 must be 64 lowercase hex digits")


@dataclass(frozen=True)
class ESMCTokenMapping:
    """Concrete tokenizer IDs observed for one backend invocation."""

    bos_token_id: Optional[int]
    eos_token_id: Optional[int]
    pad_token_id: Optional[int]
    unk_token_id: Optional[int]


@dataclass(frozen=True)
class ESMCBatchOutput:
    """Narrow injectable backend result used by production and fake backends."""

    embeddings: np.ndarray
    input_ids: np.ndarray
    attention_mask: np.ndarray
    special_tokens_mask: np.ndarray
    token_mapping: ESMCTokenMapping
    truncated: Tuple[bool, ...]
    output_field: str = "last_hidden_state"
    final_layer: bool = True


class ESMCBackend(Protocol):
    """Backend boundary that permits deterministic tests without model weights."""

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
        """Return padded final-layer token embeddings and token masks."""


@dataclass(frozen=True)
class ESMCEmbeddingResult:
    """One validated residue embedding and its complete stored provenance."""

    entity_id: str
    embedding: np.ndarray
    metadata: Mapping[str, Any]
    artifact_path: Path
    reused: bool


SequenceInput = Union[Mapping[str, str], Sequence[Tuple[str, str]]]


def normalize_sequence(sequence: str) -> str:
    """Return the stable sequence representation used for hashing and inference."""

    if not isinstance(sequence, str):
        raise TypeError("Protein sequences must be strings")
    normalized = "".join(sequence.split()).upper()
    if not normalized:
        raise ESMCValidationError("Protein sequences must not be empty")
    if not re.fullmatch(r"[A-Z]+", normalized):
        raise ESMCValidationError(
            "Protein sequences may contain ASCII letters and whitespace only"
        )
    return normalized


def sequence_sha256(sequence: str) -> str:
    """Hash an already normalized sequence with stable UTF-8 encoding."""

    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _canonical_json(value: Mapping[str, Any]) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _normalize_inputs(sequences: SequenceInput) -> List[Tuple[str, str]]:
    items = list(sequences.items()) if isinstance(sequences, Mapping) else list(sequences)
    if not items:
        raise ESMCValidationError("At least one sequence is required")

    normalized_by_id: Dict[str, str] = {}
    ordered: List[Tuple[str, str]] = []
    for entity_id, sequence in items:
        if not isinstance(entity_id, str) or not entity_id:
            raise ESMCValidationError("Sequence entity IDs must be non-empty strings")
        normalized = normalize_sequence(sequence)
        previous = normalized_by_id.get(entity_id)
        if previous is not None:
            if previous != normalized:
                raise ESMCSequenceConflictError(
                    f"Entity '{entity_id}' was supplied with conflicting sequences"
                )
            continue
        normalized_by_id[entity_id] = normalized
        ordered.append((entity_id, normalized))
    return ordered


def _validate_token_mapping(mapping: ESMCTokenMapping) -> Tuple[int, int, int, int]:
    values = (
        mapping.bos_token_id,
        mapping.eos_token_id,
        mapping.pad_token_id,
        mapping.unk_token_id,
    )
    if any(value is None or isinstance(value, bool) for value in values):
        raise ESMCTokenPolicyError(
            "Tokenizer must define unambiguous BOS, EOS, padding, and unknown token IDs"
        )
    integer_values = tuple(int(value) for value in values)
    if len(set(integer_values)) != len(integer_values):
        raise ESMCTokenPolicyError(
            "BOS, EOS, padding, and unknown token IDs must be distinct"
        )
    return integer_values  # type: ignore[return-value]


def validate_transformers_runtime(
    provenance: ESMCTransformersProvenance = ESMCTransformersProvenance(),
) -> str:
    """Validate the installed fork from its PEP 610 ``direct_url.json``.

    Returns the installed distribution version after proving that the code came
    from the exact pinned Biohub repository and resolved Git commit. A stock
    wheel, an editable/local checkout, or another fork revision is rejected.
    """

    try:
        distribution = importlib.metadata.distribution(provenance.distribution)
    except importlib.metadata.PackageNotFoundError as exc:
        raise RuntimeError("The pinned Biohub Transformers fork is not installed") from exc
    direct_url_text = distribution.read_text("direct_url.json")
    if not direct_url_text:
        raise RuntimeError(
            "ESMC requires a PEP 508 VCS installation of the pinned Biohub "
            "Transformers fork; stock Transformers is incompatible"
        )
    try:
        direct_url = json.loads(direct_url_text)
        vcs_info = direct_url["vcs_info"]
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise RuntimeError("Transformers direct_url.json has no usable VCS provenance") from exc
    repository = str(direct_url.get("url", "")).rstrip("/")
    expected_repository = provenance.repository.rstrip("/")
    if repository != expected_repository:
        raise RuntimeError(
            f"ESMC requires Transformers from {expected_repository}, found {repository or 'unknown'}"
        )
    if vcs_info.get("vcs") != "git" or vcs_info.get("commit_id") != provenance.revision:
        raise RuntimeError(
            "Installed Biohub Transformers commit does not match the pinned ESMC runtime"
        )
    return distribution.version


class HuggingFaceESMCBackend:
    """Direct pinned Hugging Face loader for Biohub ESMC-6B."""

    def __init__(self) -> None:
        self._model: Any = None
        self._tokenizer: Any = None

    @staticmethod
    def tokenizer_load_kwargs(model: ESMCModelProvenance) -> Dict[str, Any]:
        """Return the pinned, remote-code-disabled tokenizer load arguments."""

        return {
            "revision": model.model_revision,
            "trust_remote_code": False,
        }

    @staticmethod
    def model_load_kwargs(
        model: ESMCModelProvenance,
        load_policy: ESMCLoadPolicy,
        device: str,
        *,
        torch_dtype: Any,
    ) -> Dict[str, Any]:
        """Construct the exact production model-loader arguments without loading."""

        is_cuda = device.startswith("cuda")
        device_map: Any = load_policy.cuda_device_map if is_cuda else {"": "cpu"}
        max_memory: Dict[Union[int, str], str]
        if is_cuda:
            max_memory = {0: load_policy.cuda_max_memory}
            if load_policy.cpu_offload:
                max_memory["cpu"] = load_policy.cpu_max_memory
        else:
            max_memory = {"cpu": load_policy.cpu_max_memory}
        return {
            "revision": model.model_revision,
            "trust_remote_code": False,
            "torch_dtype": torch_dtype,
            "device_map": device_map,
            "max_memory": max_memory,
            "low_cpu_mem_usage": load_policy.low_cpu_mem_usage,
        }

    @staticmethod
    def validate_config(config: Any, model: ESMCModelProvenance) -> None:
        """Validate the pinned model architecture without allocating weights."""

        if getattr(config, "d_model", None) != model.embedding_dimension:
            raise ESMCShapeError("Pinned ESMC config does not declare d_model=2560")
        if getattr(config, "n_layers", None) != model.layer_count:
            raise ESMCShapeError("Pinned ESMC config does not declare 80 layers")
        if getattr(config, "model_type", None) != "esmc":
            raise ESMCValidationError("Pinned model config is not model_type='esmc'")

    def _load(
        self,
        model_provenance: ESMCModelProvenance,
        transformers_provenance: ESMCTransformersProvenance,
        load_policy: ESMCLoadPolicy,
        device: str,
    ) -> None:
        try:
            import torch
            from transformers import AutoModelForMaskedLM, AutoTokenizer
        except ImportError as exc:  # pragma: no cover - depends on optional extras
            raise RuntimeError(
                "ESMC inference requires PyTorch, accelerate, and the exact pinned "
                "Biohub Transformers fork"
            ) from exc

        validate_transformers_runtime(transformers_provenance)

        if self._model is not None:
            return
        self._tokenizer = AutoTokenizer.from_pretrained(
            model_provenance.model_id,
            **self.tokenizer_load_kwargs(model_provenance),
        )
        self._model = AutoModelForMaskedLM.from_pretrained(
            model_provenance.model_id,
            **self.model_load_kwargs(
                model_provenance,
                load_policy,
                device,
                torch_dtype=torch.bfloat16,
            ),
        )
        config = self._model.config
        self.validate_config(config, model_provenance)
        self._model.eval()

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
        self._load(model, transformers, load_policy, device)
        import torch

        encoded = self._tokenizer(
            list(sequences),
            add_special_tokens=token_policy.add_special_tokens,
            padding=token_policy.padding,
            truncation=False,
            return_special_tokens_mask=True,
            return_tensors="pt",
        )
        special_tokens_mask = encoded.pop("special_tokens_mask")
        input_device = next(self._model.parameters()).device
        model_inputs = {key: value.to(input_device) for key, value in encoded.items()}
        with torch.inference_mode():
            output = self._model(**model_inputs, return_dict=True)
        embeddings = getattr(output, "last_hidden_state", None)
        if embeddings is None:
            raise ESMCShapeError("ESMC output is missing the required 'last_hidden_state' field")

        bos_token_id = getattr(self._tokenizer, "bos_token_id", None)
        if bos_token_id is None:
            bos_token_id = getattr(self._tokenizer, "cls_token_id", None)
        mapping = ESMCTokenMapping(
            bos_token_id=bos_token_id,
            eos_token_id=getattr(self._tokenizer, "eos_token_id", None),
            pad_token_id=getattr(self._tokenizer, "pad_token_id", None),
            unk_token_id=getattr(self._tokenizer, "unk_token_id", None),
        )
        return ESMCBatchOutput(
            embeddings=embeddings.detach().float().cpu().numpy(),
            input_ids=encoded["input_ids"].cpu().numpy(),
            attention_mask=encoded["attention_mask"].cpu().numpy(),
            special_tokens_mask=special_tokens_mask.cpu().numpy(),
            token_mapping=mapping,
            truncated=tuple(False for _ in sequences),
            output_field="last_hidden_state",
            final_layer=True,
        )


class ESMCEmbeddingAdapter:
    """Validate, persist, and exactly resume pinned ESMC embeddings."""

    def __init__(
        self,
        *,
        embeddings_dir: Path,
        data_root: Path,
        entity_registry: Any,
        batch_size: int,
        device: str,
        backend: Optional[ESMCBackend] = None,
        model_provenance: ESMCModelProvenance = ESMCModelProvenance(),
        transformers_provenance: ESMCTransformersProvenance = ESMCTransformersProvenance(),
        load_policy: ESMCLoadPolicy = ESMCLoadPolicy(),
        token_policy: ESMCTokenPolicy = ESMCTokenPolicy(),
        output_provenance: ESMCOutputProvenance = ESMCOutputProvenance(),
    ) -> None:
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        self.embeddings_dir = Path(embeddings_dir)
        self.data_root = Path(data_root)
        self.entity_registry = entity_registry
        self.batch_size = batch_size
        self.device = device
        self.backend = backend or HuggingFaceESMCBackend()
        self.model_provenance = model_provenance
        self.transformers_provenance = transformers_provenance
        self.load_policy = load_policy
        self.token_policy = token_policy
        self.output_provenance = output_provenance
        self.execution_device = self._normalize_device(device)
        if (
            self.model_provenance.embedding_dimension
            != self.output_provenance.dimension
        ):
            raise ESMCShapeError(
                "Model and output provenance declare different embedding dimensions"
            )
        if (
            self.output_provenance.output_field != "last_hidden_state"
            or self.output_provenance.layer != "final"
            or self.output_provenance.embedding_type != "per_residue"
            or self.output_provenance.storage_dtype != "float32"
        ):
            raise ESMCShapeError("ESMC output provenance must remain final per-residue float32")

    @property
    def namespace_dir(self) -> Path:
        model_namespace = re.sub(
            r"[^a-z0-9]+", "_", self.model_provenance.model_id.rsplit("/", 1)[-1].lower()
        ).strip("_")
        namespace_hash = _sha256_bytes(
            _canonical_json(
                {
                    "adapter_id": self.model_provenance.adapter_id,
                    "model_id": self.model_provenance.model_id,
                }
            ).encode("utf-8")
        )[:12]
        return self.embeddings_dir / model_namespace / f"artifacts-v1-{namespace_hash}"

    @staticmethod
    def _normalize_device(device: str) -> str:
        normalized = str(device).strip().lower()
        if not normalized:
            raise ESMCValidationError("ESMC execution device must be non-empty")
        return normalized

    def embed(
        self,
        sequences: SequenceInput,
        *,
        source_lineage: Optional[Mapping[str, ESMCSourceLineage]] = None,
        resume: bool = True,
        force_recompute: bool = False,
    ) -> Dict[str, ESMCEmbeddingResult]:
        normalized_items = _normalize_inputs(sequences)
        for entity_id, sequence in normalized_items:
            if len(sequence) > self.token_policy.max_residues:
                raise ESMCTruncationError(
                    f"Entity '{entity_id}' has {len(sequence)} residues; policy maximum is "
                    f"{self.token_policy.max_residues}, and truncation is forbidden"
                )

        prepared: List[Tuple[str, str, ESMCSourceLineage, Dict[str, Any], Path]] = []
        results: Dict[str, ESMCEmbeddingResult] = {}
        for entity_id, sequence in normalized_items:
            lineage = (
                source_lineage.get(entity_id)
                if source_lineage and entity_id in source_lineage
                else ESMCSourceLineage("direct_sequence", entity_id)
            )
            if not isinstance(lineage, ESMCSourceLineage):
                raise TypeError("source_lineage values must be ESMCSourceLineage instances")
            identity = self._identity(entity_id, sequence, lineage)
            artifact_path = self._artifact_path(entity_id, identity["identity_sha256"])
            self._validate_entity_history(entity_id, sequence, artifact_path.parent)
            if artifact_path.exists() and resume and not force_recompute:
                results[entity_id] = self._load_artifact(artifact_path, identity, reused=True)
            else:
                prepared.append((entity_id, sequence, lineage, identity, artifact_path))

        for start in range(0, len(prepared), self.batch_size):
            batch = prepared[start : start + self.batch_size]
            backend_output = self.backend.encode(
                [item[1] for item in batch],
                model=self.model_provenance,
                transformers=self.transformers_provenance,
                load_policy=self.load_policy,
                token_policy=self.token_policy,
                device=self.device,
            )
            arrays, token_mapping = self._extract_residue_arrays(
                [item[1] for item in batch], backend_output
            )
            for item, array in zip(batch, arrays):
                entity_id, _sequence, _lineage, identity, artifact_path = item
                metadata = self._metadata(identity, array, token_mapping)
                self._save_artifact(artifact_path, array, metadata)
                self._register(entity_id, artifact_path, metadata)
                results[entity_id] = ESMCEmbeddingResult(
                    entity_id=entity_id,
                    embedding=array,
                    metadata=metadata,
                    artifact_path=artifact_path,
                    reused=False,
                )

        return {entity_id: results[entity_id] for entity_id, _ in normalized_items}

    def _identity(
        self, entity_id: str, sequence: str, lineage: ESMCSourceLineage
    ) -> Dict[str, Any]:
        provenance = {
            "schema_version": "protos-esmc-artifact-v1",
            "source_entity_id": entity_id,
            "normalized_sequence": sequence,
            "sequence_sha256": sequence_sha256(sequence),
            "sequence_length": len(sequence),
            "model": asdict(self.model_provenance),
            "transformers_runtime": asdict(self.transformers_provenance),
            "load_policy": asdict(self.load_policy),
            "token_policy": asdict(self.token_policy),
            "output": asdict(self.output_provenance),
            "source_lineage": asdict(lineage),
            "execution_device": self.execution_device,
        }
        return {
            **provenance,
            "identity_sha256": _sha256_bytes(_canonical_json(provenance).encode("utf-8")),
        }

    def _entity_dir(self, entity_id: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", entity_id).strip("._") or "entity"
        id_hash = _sha256_bytes(entity_id.encode("utf-8"))[:16]
        return self.namespace_dir / f"{safe[:80]}--{id_hash}"

    def _artifact_path(self, entity_id: str, identity_sha256: str) -> Path:
        return self._entity_dir(entity_id) / f"{identity_sha256}.npz"

    def _extract_residue_arrays(
        self, sequences: Sequence[str], output: ESMCBatchOutput
    ) -> Tuple[List[np.ndarray], ESMCTokenMapping]:
        bos_id, eos_id, pad_id, unk_id = _validate_token_mapping(output.token_mapping)
        if output.output_field != self.output_provenance.output_field or not output.final_layer:
            raise ESMCShapeError("Backend did not return final-layer output.last_hidden_state")
        if len(output.truncated) != len(sequences) or any(output.truncated):
            raise ESMCTruncationError("Backend reported truncated ESMC tokenization")

        embeddings = np.asarray(output.embeddings)
        input_ids = np.asarray(output.input_ids)
        attention = np.asarray(output.attention_mask)
        special = np.asarray(output.special_tokens_mask)
        if embeddings.ndim != 3:
            raise ESMCShapeError("ESMC backend embeddings must have shape [batch,tokens,dim]")
        if not np.issubdtype(embeddings.dtype, np.floating):
            raise ESMCShapeError("ESMC backend embeddings must have a floating dtype")
        if embeddings.shape[0] != len(sequences):
            raise ESMCShapeError("ESMC backend batch dimension does not match inputs")
        if embeddings.shape[2] != self.output_provenance.dimension:
            raise ESMCShapeError(
                f"ESMC embedding dimension must be {self.output_provenance.dimension}, "
                f"got {embeddings.shape[2]}"
            )
        expected_mask_shape = embeddings.shape[:2]
        for name, value in (
            ("input_ids", input_ids),
            ("attention_mask", attention),
            ("special_tokens_mask", special),
        ):
            if value.shape != expected_mask_shape:
                raise ESMCShapeError(
                    f"{name} shape {value.shape} does not match token shape {expected_mask_shape}"
                )

        arrays: List[np.ndarray] = []
        for index, sequence in enumerate(sequences):
            active = attention[index].astype(bool)
            active_count = int(active.sum())
            if active_count < len(sequence) + 2:
                raise ESMCTruncationError(
                    f"Sequence {index} retained only {active_count - 2} residue tokens"
                )
            if active_count != len(sequence) + 2:
                raise ESMCTokenPolicyError(
                    f"Sequence {index} token count is not exactly L+2"
                )
            if not np.array_equal(
                active, np.arange(expected_mask_shape[1]) < active_count
            ):
                raise ESMCTokenPolicyError("Attention mask must be contiguous then padded")
            if input_ids[index, 0] != bos_id or input_ids[index, active_count - 1] != eos_id:
                raise ESMCTokenPolicyError("BOS/EOS token placement is missing or ambiguous")
            if np.any(input_ids[index, active_count:] != pad_id):
                raise ESMCTokenPolicyError("Padding positions do not use the declared pad token")

            expected_special = np.ones(expected_mask_shape[1], dtype=bool)
            expected_special[1 : active_count - 1] = False
            if not np.array_equal(special[index].astype(bool), expected_special):
                raise ESMCTokenPolicyError(
                    "special_tokens_mask must identify exactly BOS, EOS, and padding"
                )
            residue_selector = active & ~special[index].astype(bool)
            residue_ids = input_ids[index, residue_selector]
            if self.token_policy.reject_unknown_tokens and np.any(residue_ids == unk_id):
                raise ESMCTokenPolicyError("Sequence produced an unknown residue token")
            if int(residue_selector.sum()) != len(sequence):
                raise ESMCTruncationError(
                    "Verified special-token removal did not retain exactly L residues"
                )
            array = np.ascontiguousarray(
                embeddings[index, residue_selector, :], dtype=np.float32
            )
            expected_shape = (len(sequence), self.output_provenance.dimension)
            if array.shape != expected_shape:
                raise ESMCShapeError(
                    f"Residue output shape must be {expected_shape}, got {array.shape}"
                )
            if not np.isfinite(array).all():
                raise ESMCShapeError("Residue output contains non-finite values")
            arrays.append(array)
        return arrays, output.token_mapping

    def _metadata(
        self,
        identity: Mapping[str, Any],
        array: np.ndarray,
        token_mapping: ESMCTokenMapping,
    ) -> Dict[str, Any]:
        return {
            **dict(identity),
            "adapter_id": self.model_provenance.adapter_id,
            "model_id": self.model_provenance.model_id,
            "model_revision": self.model_provenance.model_revision,
            "code_revision": self.model_provenance.code_revision,
            "transformers_repository": self.transformers_provenance.repository,
            "transformers_revision": self.transformers_provenance.revision,
            "execution_device": self.execution_device,
            "embedding_dimension": self.output_provenance.dimension,
            "inference_dtype": self.load_policy.inference_dtype,
            "dtype": self.output_provenance.storage_dtype,
            "shape": list(array.shape),
            "token_mapping": asdict(token_mapping),
            "embedding_sha256": _sha256_bytes(array.tobytes(order="C")),
        }

    def _save_artifact(
        self, artifact_path: Path, array: np.ndarray, metadata: Mapping[str, Any]
    ) -> None:
        artifact_path.parent.mkdir(parents=True, exist_ok=True)
        fd, temp_name = tempfile.mkstemp(
            dir=artifact_path.parent, prefix=".tmp-esmc-", suffix=".npz"
        )
        os.close(fd)
        try:
            np.savez_compressed(
                temp_name,
                embedding=array,
                metadata=np.asarray(_canonical_json(metadata)),
            )
            os.replace(temp_name, artifact_path)
        finally:
            if os.path.exists(temp_name):
                os.unlink(temp_name)

    def _load_artifact(
        self,
        artifact_path: Path,
        expected_identity: Mapping[str, Any],
        *,
        reused: bool,
    ) -> ESMCEmbeddingResult:
        try:
            with np.load(artifact_path, allow_pickle=False) as archive:
                if set(archive.files) != {"embedding", "metadata"}:
                    raise ESMCArtifactError("ESMC artifact has unexpected members")
                array = np.asarray(archive["embedding"])
                metadata = json.loads(str(archive["metadata"].item()))
        except ESMCArtifactError:
            raise
        except Exception as exc:
            raise ESMCArtifactError(f"Cannot read ESMC artifact {artifact_path}") from exc

        for key, expected in expected_identity.items():
            if metadata.get(key) != expected:
                raise ESMCArtifactError(
                    f"Stored ESMC provenance mismatch for '{key}' at {artifact_path}"
                )
        exact_fields = {
            "adapter_id": self.model_provenance.adapter_id,
            "model_id": self.model_provenance.model_id,
            "model_revision": self.model_provenance.model_revision,
            "code_revision": self.model_provenance.code_revision,
            "transformers_repository": self.transformers_provenance.repository,
            "transformers_revision": self.transformers_provenance.revision,
            "execution_device": self.execution_device,
            "embedding_dimension": self.output_provenance.dimension,
            "inference_dtype": self.load_policy.inference_dtype,
            "dtype": self.output_provenance.storage_dtype,
        }
        for key, expected in exact_fields.items():
            if metadata.get(key) != expected:
                raise ESMCArtifactError(
                    f"Stored ESMC metadata mismatch for '{key}' at {artifact_path}"
                )
        expected_shape = (
            expected_identity["sequence_length"],
            self.output_provenance.dimension,
        )
        if array.shape != expected_shape or array.dtype != np.dtype(np.float32):
            raise ESMCArtifactError(
                f"Stored ESMC array must be float32 {expected_shape}, "
                f"got {array.dtype} {array.shape}"
            )
        if metadata.get("shape") != list(array.shape) or metadata.get("dtype") != "float32":
            raise ESMCArtifactError("Stored ESMC shape/dtype metadata does not match its array")
        if metadata.get("embedding_sha256") != _sha256_bytes(array.tobytes(order="C")):
            raise ESMCArtifactError("Stored ESMC embedding content hash mismatch")
        mapping = ESMCTokenMapping(**metadata.get("token_mapping", {}))
        _validate_token_mapping(mapping)
        return ESMCEmbeddingResult(
            entity_id=str(expected_identity["source_entity_id"]),
            embedding=np.ascontiguousarray(array),
            metadata=metadata,
            artifact_path=artifact_path,
            reused=reused,
        )

    def _validate_entity_history(
        self, entity_id: str, sequence: str, entity_dir: Path
    ) -> None:
        if not entity_dir.exists():
            return
        expected_hash = sequence_sha256(sequence)
        for artifact_path in sorted(entity_dir.glob("*.npz")):
            try:
                with np.load(artifact_path, allow_pickle=False) as archive:
                    metadata = json.loads(str(archive["metadata"].item()))
            except Exception as exc:
                raise ESMCArtifactError(
                    f"Cannot validate stored history for entity '{entity_id}'"
                ) from exc
            if metadata.get("source_entity_id") != entity_id:
                raise ESMCArtifactError("ESMC entity namespace collision detected")
            if metadata.get("sequence_sha256") != expected_hash:
                raise ESMCSequenceConflictError(
                    f"Stored ESMC entity '{entity_id}' has a conflicting sequence"
                )

    def _register(
        self, entity_id: str, artifact_path: Path, metadata: Mapping[str, Any]
    ) -> None:
        model_namespace = re.sub(
            r"[^a-z0-9]+", "_", self.model_provenance.model_id.rsplit("/", 1)[-1].lower()
        ).strip("_")
        entity_name = (
            f"{entity_id}__{model_namespace}__{str(metadata['identity_sha256'])[:16]}"
        )
        try:
            relative_path = str(artifact_path.relative_to(self.data_root))
        except ValueError as exc:
            raise ESMCArtifactError("ESMC artifact is outside the Protos data root") from exc
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="embedding",
            file_path=relative_path,
            metadata=dict(metadata),
        )


__all__ = [
    "ESMCArtifactError",
    "ESMCBackend",
    "ESMCBatchOutput",
    "ESMCEmbeddingAdapter",
    "ESMCEmbeddingResult",
    "ESMCLoadPolicy",
    "ESMCModelProvenance",
    "ESMCOutputProvenance",
    "ESMCSequenceConflictError",
    "ESMCShapeError",
    "ESMCSourceLineage",
    "ESMCTokenMapping",
    "ESMCTokenPolicy",
    "ESMCTokenPolicyError",
    "ESMCTransformersProvenance",
    "ESMCTruncationError",
    "ESMCValidationError",
    "HuggingFaceESMCBackend",
    "normalize_sequence",
    "sequence_sha256",
    "validate_transformers_runtime",
]
