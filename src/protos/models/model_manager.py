"""Model orchestration utilities for Protos models.

The ModelManager coordinates artifact assembly through processors and delegates
execution to model-specific adapters. It supports both external configuration
jobs (e.g., Boltz) and in-process models (e.g., Lambda).
"""

from __future__ import annotations

import os
import json
import shutil
from abc import ABC, abstractmethod
from datetime import datetime
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, MutableMapping, Optional

import importlib

import numpy as np
import pandas as pd
import yaml
import shlex
import inspect

from protos.io.paths import ProtosPaths
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.grn import GRNProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor

from .model_specs import (
    ArtifactBundle,
    ArtifactSpec,
    ExecutionSpec,
    ModelCard,
    ModelInvocation,
    PreparedJob,
    RuntimeResult,
)

_lambda_runtime = importlib.import_module("protos.models.lambda.runtime_utils")
InMemoryEmbeddingAdapter = _lambda_runtime.InMemoryEmbeddingAdapter
InMemoryGRNAdapter = _lambda_runtime.InMemoryGRNAdapter
align_embeddings_to_grn = _lambda_runtime.align_embeddings_to_grn
build_grn_assignments = _lambda_runtime.build_grn_assignments
build_property_rows = _lambda_runtime.build_property_rows
copy_if_missing = _lambda_runtime.copy_if_missing
ensure_positional_map = _lambda_runtime.ensure_positional_map
prepare_predictor = _lambda_runtime.prepare_predictor
build_positional_map_from_binding = _lambda_runtime.build_positional_map_from_binding
normalize_binding_config = _lambda_runtime.normalize_binding_config


class ModelAdapterBase(ABC):
    """Base class for model adapters."""

    def __init__(self, manager: "ModelManager") -> None:
        self.manager = manager
        self.paths = manager.paths

    @staticmethod
    def _require_bundle(bundles: List[ArtifactBundle], name: str) -> ArtifactBundle:
        for bundle in bundles:
            if bundle.spec.name == name:
                return bundle
        raise KeyError(f"Required artifact '{name}' not available")

    def _ensure_pdb_file(self, source: Path, ctx: "ModelRunContext", *, purpose: str) -> Path:
        """Ensure a PDB file exists for a given source structure path.

        - If source is CIF/MMCIF, convert to PDB using gemmi and write into ctx.inputs_dir
        - If source is PDB, copy into ctx.inputs_dir for provenance
        Returns the path to the ensured PDB file.
        """
        pdb_src = Path(source)
        if pdb_src.suffix.lower() in {".cif", ".mmcif"}:
            try:
                import gemmi  # type: ignore
                st = gemmi.read_structure(str(pdb_src))
                pdb_dst = ctx.inputs_dir / (pdb_src.stem + ".pdb")
                pdb_dst.parent.mkdir(parents=True, exist_ok=True)
                st.write_pdb(str(pdb_dst))
                return pdb_dst
            except Exception as exc:
                raise RuntimeError(f"[{purpose}] CIF→PDB conversion failed: {exc}")
        else:
            pdb_dst = ctx.inputs_dir / pdb_src.name
            try:
                shutil.copy2(pdb_src, pdb_dst)
            except Exception:
                pass
            return pdb_dst

    def prepare(
        self,
        card: ModelCard,
        request: "ModelRequest",
    ) -> ModelInvocation:
        inputs = self.manager.assemble_inputs(card, request)
        invocation = self._prepare(card, request, inputs)
        invocation.model = card.name
        invocation.card = card
        invocation.inputs = inputs
        return invocation

    @abstractmethod
    def _prepare(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> ModelInvocation:
        """Create a model invocation using resolved input artifacts."""


class ExternalJobAdapter(ModelAdapterBase):
    """Adapter for models executed outside of the current process."""

    def _prepare(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> ModelInvocation:
        job = self.build_job(card, request, inputs)
        return ModelInvocation(model=card.name, card=card, job=job, inputs=inputs)

    @abstractmethod
    def build_job(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        """Create the executable payload for this model."""


class PlaceholderAdapter(ModelAdapterBase):
    """Non-executing adapter used while wiring up ModelCards.

    This adapter allows model cards to be listed and described before
    their full adapters and artifact providers are implemented. Any
    attempt to prepare a job will raise a clear error.
    """

    def _prepare(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> ModelInvocation:  # pragma: no cover - placeholder behavior
        raise NotImplementedError(
            f"Model adapter for '{card.name}' is not implemented yet."
        )


class RuntimeAdapter(ModelAdapterBase):
    """Adapter for models executed within the current process."""

    def _prepare(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> ModelInvocation:
        runtime = self.run_runtime(card, request, inputs)
        invocation = ModelInvocation(
            model=card.name,
            card=card,
            runtime=runtime,
            inputs=inputs,
            outputs=runtime.artifacts,
            metadata=dict(runtime.metadata),
        )
        return invocation

    @abstractmethod
    def run_runtime(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> RuntimeResult:
        """Execute the model locally and return results."""


class ModelRunContext:
    """Filesystem context for a single model run under data/models/<model>/.

    Creates a stable directory layout expected by our workflow:
    - work_dir:    data/models/<model>/<run_id>/
    - inputs_dir:  .../inputs/
    - outputs_dir: .../outputs/
    - config_path: optional, when a config file is created
    """

    def __init__(self, paths: ProtosPaths, card: ModelCard, run_prefix: str = "job") -> None:
        self.paths = paths
        self.card = card
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.work_dir = (
            Path(self.paths.data_root) / "models" / card.name / f"{run_prefix}_{self.run_id}"
        )
        self.inputs_dir = self.work_dir / "inputs"
        self.outputs_dir = self.work_dir / "outputs"
        self.config_path: Optional[Path] = None

    def create(self) -> None:
        self.inputs_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_dir.mkdir(parents=True, exist_ok=True)

    def package_inputs(self, bundles: List[ArtifactBundle]) -> Dict[str, str]:
        """Copy input artifacts into inputs/ and return mapping of names to paths."""
        packaged: Dict[str, str] = {}
        for bundle in bundles:
            src = Path(bundle.path)
            if not src.exists():
                continue
            dest_name = f"{bundle.spec.name}__{src.name}"
            dest = self.inputs_dir / dest_name
            try:
                if src.is_dir():
                    # If dest exists from a previous run, remove it first
                    if dest.exists():
                        shutil.rmtree(dest)
                    shutil.copytree(src, dest)
                else:
                    shutil.copy2(src, dest)
            except Exception:
                continue
            packaged[bundle.spec.name] = str(dest)
        return packaged

    def as_metadata(self) -> Dict[str, Optional[str]]:
        return {
            "work_dir": str(self.work_dir),
            "inputs_dir": str(self.inputs_dir),
            "outputs_dir": str(self.outputs_dir),
            "config_path": str(self.config_path) if self.config_path else None,
        }

class ConfigurableRuntimeAdapter(RuntimeAdapter):
    """Runtime adapter driven by ModelCard.execution.entrypoint.

    The entrypoint should be a dotted path "module.sub:callable". The callable
    will be invoked with a filtered subset of the following keyword args,
    depending on its signature: card, request, inputs, paths.

    The callable may return one of:
    - RuntimeResult
    - dict with optional keys 'outputs', 'artifacts', 'metadata'
    - any other object, which will be wrapped under outputs={'result': obj}
    """

    def run_runtime(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> RuntimeResult:
        entry = card.execution.entrypoint
        if not entry:
            raise ValueError(
                f"Model '{card.name}' missing execution.entrypoint for runtime"
            )

        module_name, _, attr = entry.partition(":")
        if not module_name or not attr:
            raise ValueError(
                "Runtime entrypoint must be in 'module.path:callable' form"
            )

        module = importlib.import_module(module_name)
        func = getattr(module, attr)

        kwargs = {"card": card, "request": request, "inputs": inputs, "paths": self.paths}
        sig = inspect.signature(func)
        filtered = {k: v for k, v in kwargs.items() if k in sig.parameters}

        # Prepare a run context so runtime callables can write outputs
        ctx = ModelRunContext(self.paths, card, run_prefix="run")
        ctx.create()
        # Pass context paths via metadata-style kwargs if accepted
        if "work_dir" in sig.parameters:
            filtered["work_dir"] = str(ctx.work_dir)
        if "outputs_dir" in sig.parameters:
            filtered["outputs_dir"] = str(ctx.outputs_dir)
        if "inputs_dir" in sig.parameters:
            filtered["inputs_dir"] = str(ctx.inputs_dir)

        result = func(**filtered)

        # Normalize to RuntimeResult
        if isinstance(result, RuntimeResult):  # type: ignore[unreachable]
            return result

        outputs: Dict[str, Any] = {}
        artifacts: List[ArtifactBundle] = []
        metadata: Dict[str, Any] = {}

        if isinstance(result, dict):
            if any(k in result for k in ("outputs", "artifacts", "metadata")):
                outputs = result.get("outputs", {})  # type: ignore[assignment]
                artifacts = result.get("artifacts", [])  # type: ignore[assignment]
                metadata = result.get("metadata", {})  # type: ignore[assignment]
            else:
                outputs = result
        else:
            outputs = {"result": result}

        # Always include run context paths in metadata
        meta = dict(metadata)
        meta.update({
            "work_dir": str(ctx.work_dir),
            "outputs_dir": str(ctx.outputs_dir),
        })
        return RuntimeResult(outputs=outputs, artifacts=artifacts, metadata=meta)


class ConfigurableExternalAdapter(ExternalJobAdapter):
    """External adapter that builds a PreparedJob from ExecutionSpec.

    - Uses 'entrypoint' as the base command (split with shlex)
    - If 'expected_config' is provided, writes assembled config (yaml/json)
      into a per-job working directory under data_root/models/<model>/<ts>
      and appends the config path to the command.
    """

    def build_job(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        entry = card.execution.entrypoint
        if not entry:
            raise ValueError(
                f"Model '{card.name}' missing execution.entrypoint for external job"
            )

        cmd = shlex.split(entry)
        ctx = ModelRunContext(self.paths, card)
        ctx.create()
        working_dir = ctx.work_dir

        # Package inputs alongside the job
        packaged_inputs = ctx.package_inputs(inputs)

        config_path = None
        if card.execution.expected_config:
            ext_hint = str(card.execution.expected_config).lower()
            if ext_hint.endswith(".json") or "json" in ext_hint:
                config_path = working_dir / "config.json"
                payload = self._build_config_payload(card, request, inputs, packaged_inputs)
                with open(config_path, "w", encoding="utf-8") as fh:
                    json.dump(payload, fh, indent=2)
            else:
                config_path = working_dir / "config.yaml"
                payload = self._build_config_payload(card, request, inputs, packaged_inputs)
                with open(config_path, "w", encoding="utf-8") as fh:
                    yaml.safe_dump(payload, fh, sort_keys=False)
            cmd.append(str(config_path))
        ctx.config_path = config_path

        job = PreparedJob(
            command=cmd,
            working_dir=working_dir,
            artifacts=list(inputs),
            metadata={
                "entrypoint": entry,
                "config_path": str(config_path) if config_path else None,
                "context": ctx.as_metadata(),
                "packaged_inputs": packaged_inputs,
            },
        )
        return job

    def _working_dir_for_job(self, card: ModelCard) -> Path:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        return Path(self.paths.data_root) / "models" / card.name / f"job_{ts}"

    @staticmethod
    def _build_config_payload(
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
        packaged: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        artifacts = {
            b.spec.name: {"path": str(b.path), "kind": b.spec.kind, "format": b.spec.format}
            for b in inputs
        }
        return {
            "model": {"name": card.name, "version": card.version},
            "execution": card.execution.mode,
            "inputs": artifacts,
            "packaged_inputs": packaged or {},
            "config": dict(request.config),
            "metadata": dict(request.metadata),
        }

class ModelRequest:
    """Simple container for invocation parameters."""

    def __init__(
        self,
        inputs: Optional[MutableMapping[str, Any]] = None,
        config: Optional[MutableMapping[str, Any]] = None,
        metadata: Optional[MutableMapping[str, Any]] = None,
    ) -> None:
        self.inputs: MutableMapping[str, Any] = inputs or {}
        self.config: MutableMapping[str, Any] = config or {}
        self.metadata: MutableMapping[str, Any] = metadata or {}

    def get_input(self, key: str, default: Any = None) -> Any:
        return self.inputs.get(key, default)


class BoltzYamlDumper(yaml.SafeDumper):
    """Custom YAML dumper to match Boltz formatting expectations."""


def _boltz_represent_list(dumper: BoltzYamlDumper, data: List[Any]):
    flow_style = all(not isinstance(item, dict) for item in data)
    return yaml.SafeDumper.represent_sequence(
        dumper,
        "tag:yaml.org,2002:seq",
        data,
        flow_style=flow_style,
    )


BoltzYamlDumper.add_representer(list, _boltz_represent_list)
BoltzYamlDumper.add_representer(tuple, _boltz_represent_list)
BoltzYamlDumper.add_representer(
    OrderedDict,
    lambda dumper, data: yaml.SafeDumper.represent_dict(dumper, data.items()),
)


class BoltzAdapter(ExternalJobAdapter):
    """Adapter for preparing Boltz configuration jobs."""

    MODEL_DIR = "boltz2"

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        sequence_bundle = self._require_bundle(inputs, "sequence_dataset")
        sequences = sequence_bundle.metadata.get("sequences", {})

        entity_name = request.get_input("entity")
        if not entity_name:
            raise ValueError("Boltz adapter requires 'entity' in inputs")
        if entity_name not in sequences:
            raise ValueError(f"Entity '{entity_name}' not found in dataset")

        mutations = request.config.get("mutations", [])
        entity_sequences = {entity_name: sequences[entity_name]}
        if mutations:
            entity_sequences = self._apply_mutations(entity_sequences, mutations)

        config_id = self._config_identifier(mutations) if mutations else "wild_type"

        yaml_data = self._generate_yaml(entity_sequences, request.config)
        input_dir = self._get_input_dir(entity_name, config_id)
        input_dir.mkdir(parents=True, exist_ok=True)

        yaml_path = input_dir / "config.yaml"
        fasta_path = input_dir / "sequences.fasta"

        with open(yaml_path, "w", encoding="utf-8") as fh:
            yaml.dump(
                yaml_data,
                fh,
                Dumper=BoltzYamlDumper,
                default_flow_style=False,
                sort_keys=False,
                indent=2,
            )

        with open(fasta_path, "w", encoding="utf-8") as fh:
            for seq_name, seq in entity_sequences.items():
                fh.write(f">{seq_name}\n{seq}\n")

        metadata_path = input_dir / "metadata.json"
        metadata_content = {
            "entity": entity_name,
            "config_id": config_id,
            "mutations": mutations,
            "dataset": sequence_bundle.metadata.get("dataset"),
        }
        with open(metadata_path, "w", encoding="utf-8") as fh:
            json.dump(metadata_content, fh, indent=2)

        job = PreparedJob(
            command=["boltz", "predict", str(yaml_path)],
            working_dir=input_dir,
            artifacts=[
                sequence_bundle,
                ArtifactBundle(
                    spec=ArtifactSpec(
                        name="boltz_config",
                        kind="config",
                        provider="boltz_adapter",
                        format="yaml",
                    ),
                    path=yaml_path,
                    metadata={"config_id": config_id},
                ),
                ArtifactBundle(
                    spec=ArtifactSpec(
                        name="boltz_sequences",
                        kind="sequence",
                        provider="boltz_adapter",
                        format="fasta",
                    ),
                    path=fasta_path,
                    metadata={"config_id": config_id},
                ),
            ],
            metadata={
                "entity": entity_name,
                "config_id": config_id,
                "mutations": mutations,
            },
        )
        return job

    def _config_identifier(self, mutations: Iterable[Dict[str, Any]]) -> str:
        labels = [mut.get("name") or f"{mut['position']}{mut['mutant']}" for mut in mutations]
        return "_".join(labels)

    def _get_input_dir(self, entity: str, config_id: str) -> Path:
        return (
            Path(self.paths.data_root)
            / "models"
            / self.MODEL_DIR
            / f"{entity}_{config_id}"
        )

    def _apply_mutations(
        self,
        sequences: Dict[str, str],
        mutations: Iterable[Dict[str, Any]],
    ) -> Dict[str, str]:
        mutated: Dict[str, str] = {}
        for name, seq in sequences.items():
            base = list(seq)
            labels: List[str] = []
            for mut in mutations:
                pos = mut.get("position")
                mutant = mut.get("mutant")
                if pos is None or mutant is None:
                    continue
                idx = pos - 1
                if idx < 0 or idx >= len(base):
                    continue
                original = mut.get("original")
                if original and base[idx] != original:
                    continue
                base[idx] = mutant
                labels.append(mut.get("name") or f"{original or seq[idx]}{pos}{mutant}")
            new_name = f"{name}_{'_'.join(labels)}" if labels else name
            mutated[new_name] = "".join(base)
        return mutated

    def _generate_yaml(
        self,
        sequences: Dict[str, str],
        config: MutableMapping[str, Any],
    ) -> OrderedDict:
        yaml_sections: OrderedDict[str, Any] = OrderedDict()
        yaml_sections["version"] = 1
        seq_entries = self._prepare_sequences_section(sequences, config)

        # Optional ligand handling
        ligand_cfg = config.get("ligand")
        ligand_id = None
        if isinstance(ligand_cfg, dict):
            ligand_id = ligand_cfg.get("id") or ligand_cfg.get("name") or "LIG"
            ligand_smiles = ligand_cfg.get("smiles")
        else:
            ligand_id = config.get("ligand_id")
            ligand_smiles = config.get("ligand_smiles")

        if ligand_smiles:
            seq_entries.append(
                {
                    "ligand": {
                        "id": ligand_id or "LIG",
                        "smiles": ligand_smiles,
                    }
                }
            )

        yaml_sections["sequences"] = seq_entries

        for section in ("constraints", "properties"):
            if section in config:
                yaml_sections[section] = deepcopy(config[section])

        additional = config.get("additional_sections", {})
        for key, value in additional.items():
            yaml_sections[key] = deepcopy(value)

        # Default properties: affinity when ligand present and no explicit properties provided
        if ligand_smiles and "properties" not in yaml_sections:
            yaml_sections["properties"] = [
                {"affinity": {"binder": (ligand_id or "LIG")}}
            ]

        return yaml_sections

    def _prepare_sequences_section(
        self,
        sequences: Dict[str, str],
        config: MutableMapping[str, Any],
    ) -> List[Dict[str, Any]]:
        if "sequences" in config:
            return deepcopy(config["sequences"])

        entries: List[Dict[str, Any]] = []
        overrides: Dict[str, Dict[str, Any]] = config.get("sequence_overrides", {})
        default_type: str = config.get("default_sequence_type", "protein")

        for index, (name, sequence_value) in enumerate(sequences.items(), start=1):
            override = overrides.get(name, {})
            seq_type: str = override.get("type", default_type)
            entry_body: Dict[str, Any] = OrderedDict()

            identifier = override.get("id") or self._default_identifier(index)
            entry_body["id"] = identifier if isinstance(identifier, list) else [identifier]

            override_sequence = override.get("sequence")
            if override_sequence is not None:
                entry_body["sequence"] = override_sequence
            elif seq_type == "protein":
                entry_body["sequence"] = sequence_value

            for key, value in override.get("fields", {}).items():
                entry_body[key] = value

            entries.append({seq_type: entry_body})

        for extra in config.get("extra_sequences", []):
            entries.append(deepcopy(extra))

        return entries

    @staticmethod
    def _default_identifier(index: int) -> str:
        base = (index - 1) % 26
        letter = chr(ord("A") + base)
        suffix = (index - 1) // 26 + 1
        return f"{letter}{suffix}"


class LambdaAdapter(RuntimeAdapter):
    """Adapter executing Lambda predictions within the current process."""

    DEFAULT_RUN_ID = "007061"
    DEFAULT_BINDING_CONFIG = Path("grn/configs/binding_domain2.json")
    DEFAULT_POSITIONAL_MAP = Path("grn/configs/final_mapping7.csv")

    def __init__(self, manager: "ModelManager") -> None:
        super().__init__(manager)
        self.lambda_root = Path(__file__).resolve().parent / "lambda"
        self._protos_data_root = Path(__file__).resolve().parents[3] / "data"

    def run_runtime(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> RuntimeResult:
        import logging
        logger = logging.getLogger("LambdaAdapter")

        # Minimal run context for optional debug artifacts
        ctx = ModelRunContext(self.manager.paths, card, run_prefix="run")
        ctx.create()

        logger.info("[lambda] Resolving inputs")
        sequence_bundle = self._require_bundle(inputs, "sequence_dataset")
        grn_bundle = self._require_bundle(inputs, "grn_table")
        embedding_bundle = next((b for b in inputs if b.spec.name == "embedding_dataset"), None)

        sequences: Dict[str, str] = sequence_bundle.metadata.get("sequences", {})
        grn_table: pd.DataFrame = grn_bundle.metadata.get("table")
        dataset_name: str = sequence_bundle.metadata.get("dataset") or "lambda_dataset"

        logger.info("[lambda] Sequences: %d | GRN rows: %d", len(sequences), int(len(grn_table) if grn_table is not None else 0))

        embeddings_map: Dict[str, np.ndarray] = {}
        embedding_dataset_name: Optional[str] = None
        embedding_model = str(request.config.get("embedding_model", "ankh_large"))
        embedding_type = str(request.config.get("embedding_type", "per_residue"))
        ingest_embeddings = bool(request.config.get("ingest_embeddings", False))

        if not sequences:
            raise ValueError("Lambda adapter requires sequence data")
        if grn_table is None or grn_table.empty:
            raise ValueError("Lambda adapter requires a populated GRN table")

        protein_family = (
            request.get_input("protein_family")
            or request.config.get("protein_family")
        )
        if not protein_family:
            raise ValueError("Lambda adapter requires 'protein_family' input")

        # Resolve embeddings: either from provided dataset or compute via embedding card
        if embedding_bundle is not None:
            logger.info("[lambda] Using provided embedding dataset '%s'", embedding_bundle.metadata.get("dataset"))
            embeddings_map = embedding_bundle.metadata.get("embeddings", {})
            embedding_dataset_name = embedding_bundle.metadata.get("dataset")
        else:
            logger.info(
                "[lambda] No embedding dataset provided; invoking embedding card model=%s type=%s",
                embedding_model,
                embedding_type,
            )
            try:
                emb_inv = self.manager.prepare(
                    f"embedding_{embedding_model}",
                    inputs={"sequence_dataset": dataset_name},
                    config={"embedding_type": embedding_type},
                )
            except Exception as exc:
                raise RuntimeError(f"Failed to prepare embedding runtime: {exc}")

            if not emb_inv.runtime:
                raise RuntimeError("Embedding runtime did not execute")

            status = emb_inv.runtime.outputs.get("status") if emb_inv.runtime.outputs else None
            if status == "skipped":
                reason = emb_inv.runtime.outputs.get("reason", "unknown")
                raise RuntimeError(f"Embedding dependencies missing ({reason})")
            if status == "error":
                err = emb_inv.runtime.outputs.get("error", "unknown error")
                raise RuntimeError(f"Embedding runtime error: {err}")

            # Load NPZ artifact
            art = None
            for b in emb_inv.runtime.artifacts:
                if b.spec.kind == "embedding":
                    art = b
                    break
            if art is None:
                raise RuntimeError("Embedding artifact not found in runtime outputs")

            try:
                npz = np.load(art.path, allow_pickle=False)
                for key in npz.files:
                    arr = np.asarray(npz[key])
                    embeddings_map[key] = arr
                logger.info("[lambda] Loaded embeddings: %d entities from NPZ", len(embeddings_map))
            except Exception as exc:
                raise RuntimeError(f"Failed to load embedding artifact: {exc}")

            # Optional ingestion as a formal dataset for traceability
            if ingest_embeddings:
                try:
                    from protos.processing.embedding import EmbeddingProcessor

                    ds_name = request.config.get(
                        "embedding_dataset_name",
                        f"{dataset_name}__{embedding_model}__{embedding_type}",
                    )
                    ep = EmbeddingProcessor(model_name=embedding_model)
                    ep.ingest_from_invocation(emb_inv, dataset_name=ds_name)
                    embedding_dataset_name = ds_name
                    logger.info("[lambda] Ingested embeddings as dataset '%s'", ds_name)
                except Exception as exc:
                    logger.warning("[lambda] Failed to ingest embeddings dataset: %s", exc)

        if not embeddings_map:
            raise ValueError("Lambda adapter requires embedding tensors")

        assignments = build_grn_assignments(grn_table)
        grn_dict, aligned_embeddings = align_embeddings_to_grn(assignments, embeddings_map)
        if not grn_dict:
            raise RuntimeError(
                "No embeddings aligned to GRN positions; aborting Lambda prediction"
            )

        logger.info(
            "[lambda] Aligned embeddings for %d proteins; first id: %s",
            len(grn_dict),
            next(iter(grn_dict.keys())) if grn_dict else "-",
        )

        available_ids = sorted(grn_dict.keys())
        filtered_table = grn_table.loc[available_ids]
        grn_adapter = InMemoryGRNAdapter(filtered_table, grn_dict)
        embedding_adapter = InMemoryEmbeddingAdapter(aligned_embeddings)

        run_id = str(request.config.get("run_id", self.DEFAULT_RUN_ID)).zfill(6)
        # Always prefer packaged Lambda configs (no global data install required)
        lambda_cfg_dir = self.lambda_root / "lmda" / "configs"
        binding_in = request.config.get("binding_config")
        if binding_in:
            cand = Path(binding_in)
            if not cand.is_absolute():
                # Resolve relative to packaged configs first
                cand2 = lambda_cfg_dir / cand.name
                binding_src = cand2 if cand2.exists() else cand
            else:
                binding_src = cand
        else:
            binding_src = lambda_cfg_dir / "binding_domain2.json"

        cfg_out_dir = ctx.work_dir / "configs"
        cfg_out_dir.mkdir(parents=True, exist_ok=True)
        binding_config_path = cfg_out_dir / "binding_domain2.json"
        try:
            if binding_src.exists():
                shutil.copy2(binding_src, binding_config_path)
            binding_config_path = normalize_binding_config(binding_config_path)
        except Exception:
            # Fall back to empty normalized file
            with open(binding_config_path, "w", encoding="utf-8") as fh:
                fh.write("{}")

        positional_in = request.config.get("positional_map")
        if positional_in:
            pcand = Path(positional_in)
            if not pcand.is_absolute():
                pcand2 = lambda_cfg_dir / pcand.name
                positional_src = pcand2 if pcand2.exists() else pcand
            else:
                positional_src = pcand
        else:
            positional_src = lambda_cfg_dir / "final_mapping7.csv"

        if positional_src.exists():
            positional_map_path = positional_src
        else:
            # Generate a default table scoped to the run context
            positional_map_path = cfg_out_dir / "final_mapping7.csv"
            ensure_positional_map(positional_map_path, protein_family=protein_family)

        if request.config.get("debug"):
            print(f"[lambda-debug] binding_config: {binding_config_path}")
            print(f"[lambda-debug] positional_map: {positional_map_path}")
        logger.info("[lambda] Prepared predictor | run_id=%s", run_id)

        predictor = prepare_predictor(
            run_id,
            self.lambda_root,
            positional_map=positional_map_path,
            binding_config=binding_config_path,
            verbose=bool(request.config.get("verbose", False)),
            prefer_checkpoint=bool(request.config.get("prefer_checkpoint", False)),
        )

        batch_size = int(request.config.get("batch_size", 8))
        collect_attention = bool(request.config.get("collect_attention", False))

        predictions_df, attention_data, _ = predictor.predict(
            grn_adapter,
            embedding_adapter,
            protein_family,
            batch_size=batch_size,
            collect_attention=collect_attention,
        )

        if predictions_df.empty:
            raise RuntimeError("Lambda predictor returned no predictions")

        property_rows = build_property_rows(
            predictions_df,
            protein_family=protein_family,
            run_id=predictor.run_id,
        )

        property_table = request.config.get("property_table") or (
            f"{dataset_name}_lambda_{predictor.run_id}"
        )

        prop_proc = PropertyProcessor()
        prop_proc.record_properties(
            property_table,
            property_rows,
            metadata={
                "model": "lambda",
                "run_id": predictor.run_id,
                "protein_family": protein_family,
                "source_sequence_dataset": dataset_name,
                "grn_table": grn_bundle.metadata.get("table_name"),
                "embedding_dataset": embedding_dataset_name,
                "embedding_model": embedding_model,
                "embedding_type": embedding_type,
            },
            allow_create=True,
            materialize_entries=False,
        )

        property_path = prop_proc.tables_dir / f"{property_table}.csv"
        property_spec = ArtifactSpec(
            name="property_table",
            kind="property",
            provider="property_processor",
            format="csv",
        )
        property_bundle = ArtifactBundle(
            spec=property_spec,
            path=property_path,
            metadata={
                "table": property_table,
                "row_count": len(property_rows),
            },
        )

        # Persist a lightweight run-config snapshot for debugging
        try:
            snapshot = {
                "run_id": predictor.run_id,
                "protein_family": protein_family,
                "sequence_dataset": dataset_name,
                "embedding_dataset": embedding_dataset_name,
                "embedding_model": embedding_model,
                "embedding_type": embedding_type,
                "batch_size": batch_size,
                "collect_attention": collect_attention,
            }
            snap_path = ctx.outputs_dir / "lambda_run_config.json"
            with open(snap_path, "w", encoding="utf-8") as fh:
                json.dump(snapshot, fh, indent=2)
        except Exception:
            pass

        runtime_metadata = {
            "property_table": property_table,
            "run_id": predictor.run_id,
            "protein_family": protein_family,
            "embedding_dataset": embedding_dataset_name,
            "embedding_model": embedding_model,
            "embedding_type": embedding_type,
            "work_dir": str(ctx.work_dir),
            "outputs_dir": str(ctx.outputs_dir),
        }

        return RuntimeResult(
            outputs={
                "predictions": predictions_df,
                "attention": attention_data,
                "property_table": property_table,
            },
            artifacts=[property_bundle],
            metadata=runtime_metadata,
        )

    def _resolve_resource(
        self,
        explicit: Optional[str],
        default_rel: Path,
        *,
        normalizer: Optional[Callable[[Path], Path]] = None,
    ) -> Path:
        import logging
        logger = logging.getLogger("LambdaAdapter")
        if explicit:
            candidate = Path(explicit)
            if not candidate.is_absolute():
                candidate = Path(self.paths.data_root) / candidate
        else:
            candidate = Path(self.paths.data_root) / default_rel

        if not candidate.exists():
            fallback = self._fallback_path(default_rel)
            logger.info("[lambda] Resolving resource '%s' -> fallback '%s' -> target '%s'", default_rel, fallback, candidate)
            candidate = copy_if_missing(fallback, candidate)

        if normalizer is not None:
            candidate = normalizer(candidate)
        return candidate

    def _resolve_positional_map(
        self,
        explicit: Optional[str],
        default_rel: Path,
        *,
        protein_family: str,
        binding_config: Optional[Path] = None,
    ) -> Path:
        import logging
        logger = logging.getLogger("LambdaAdapter")
        if explicit:
            candidate = Path(explicit)
            if not candidate.is_absolute():
                candidate = Path(self.paths.data_root) / candidate
        else:
            candidate = Path(self.paths.data_root) / default_rel

        if candidate.exists():
            return candidate

        fallback = self._fallback_path(default_rel)
        logger.info("[lambda] Resolving positional map '%s' -> fallback '%s' -> target '%s'", default_rel, fallback, candidate)
        if fallback.exists():
            return copy_if_missing(fallback, candidate)

        if binding_config is not None and binding_config.exists():
            return build_positional_map_from_binding(
                binding_config,
                candidate,
                protein_family=protein_family,
            )

        return ensure_positional_map(candidate, protein_family=protein_family)

    def _fallback_path(self, relative: Path) -> Path:
        trimmed = relative
        search_paths = [
            self.lambda_root / relative,
            self.lambda_root / trimmed,
            self._protos_data_root / trimmed,
            Path(__file__).resolve().parents[3] / trimmed,
        ]
        lambda_config_dir = self.lambda_root / "lmda" / "configs"
        search_paths.insert(0, lambda_config_dir / relative.name)
        search_paths.insert(1, lambda_config_dir / trimmed.name)
        for path in search_paths:
            if path.exists():
                return path
        return search_paths[0]


class ModelManager:
    """Coordinate model invocations based on ModelCards and adapters."""

    def __init__(self, data_root: Optional[Path] = None) -> None:
        self.paths = ProtosPaths(data_root=str(data_root) if data_root else None)
        self.cards: Dict[str, ModelCard] = {}
        self.adapters: Dict[str, ModelAdapterBase] = {}
        self._artifact_providers: Dict[str, Any] = {
            "sequence_dataset": self._provide_sequence_dataset,
            "grn_table": self._provide_grn_table,
            "embedding_dataset": self._provide_embedding_dataset,
            # Optional providers for new cards
            "graph_entity": self._provide_graph_entity,
            "ligand_file": self._provide_ligand_file,
            "ligand_molecule": self._provide_ligand_file,
            "json_payload": self._provide_json_payload,
            "structure_entity": self._provide_structure_entity,
            "file_path": self._provide_file_path,
        }
        self._register_defaults()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def prepare(
        self,
        model_name: str,
        *,
        inputs: Optional[MutableMapping[str, Any]] = None,
        config: Optional[MutableMapping[str, Any]] = None,
        metadata: Optional[MutableMapping[str, Any]] = None,
    ) -> ModelInvocation:
        import logging
        logger = logging.getLogger("ModelManager")
        card = self.cards.get(model_name)
        if card is None:
            raise ValueError(f"No model card registered for '{model_name}'")
        adapter = self.adapters.get(model_name)
        if adapter is None:
            raise ValueError(f"No adapter registered for '{model_name}'")

        request = ModelRequest(inputs=inputs, config=config, metadata=metadata)
        logger.info("[manager] Preparing model '%s' mode=%s", model_name, card.execution.mode)
        invocation = adapter.prepare(card, request)
        logger.info("[manager] Prepared invocation for '%s' external=%s runtime=%s", model_name, bool(invocation.job), bool(invocation.runtime))
        invocation.metadata.update(request.metadata)
        return invocation

    def prepare_boltz_mutations(
        self,
        dataset_name: str,
        mutations: List[Dict[str, Any]],
        *,
        base_config: Optional[MutableMapping[str, Any]] = None,
    ) -> List[ModelInvocation]:
        """Prepare Boltz jobs for a list of mutation configurations.

        Args:
            dataset_name: Name of the sequence dataset containing base entities.
            mutations: List of dictionaries with at minimum an "entity" key and
                optional "mutations" / "config" overrides.
            base_config: Shared configuration options applied to every job before
                per-entry overrides.
        """

        invocations: List[ModelInvocation] = []
        for entry in mutations:
            entity = entry.get("entity")
            if not entity:
                raise ValueError("Each mutation entry must include an 'entity'")

            config = deepcopy(base_config) if base_config else {}
            entry_config = entry.get("config", {})
            if entry_config:
                config.update(entry_config)
            if "mutations" in entry:
                config["mutations"] = entry["mutations"]

            invocation = self.prepare(
                "boltz2",
                inputs={"sequence_dataset": dataset_name, "entity": entity},
                config=config,
                metadata={"mutation_entry": entry},
            )
            invocations.append(invocation)

        return invocations

    # ------------------------------------------------------------------
    # Card/adapter registration
    # ------------------------------------------------------------------

    def _register_defaults(self) -> None:
        boltz_card = ModelCard(
            name="boltz2",
            version="2",
            description="Boltz structure prediction",
            execution=ExecutionSpec(mode="external_config", entrypoint="boltz predict"),
            input_spec=[
                ArtifactSpec(
                    name="sequence_dataset",
                    kind="sequence",
                    provider="sequence_dataset",
                    format="fasta",
                )
            ],
            output_spec=[],
        )
        self.register_model(boltz_card, BoltzAdapter(self))

        # Embedding runtime cards (delegated to embedding_runtime.run_embedding)
        for emb_model in ("esm2_t12_35m", "ankh_large", "esm2_t33_650m"):
            emb_card = ModelCard(
                name=f"embedding_{emb_model}",
                version="1.0",
                description=f"Embeddings with {emb_model}",
                execution=ExecutionSpec(
                    mode="runtime",
                    entrypoint="protos.models.embedding_runtime:run_embedding",
                    environment={"gpu": True},
                ),
                input_spec=[
                    ArtifactSpec(
                        name="sequence_dataset",
                        kind="sequence",
                        provider="sequence_dataset",
                        format="fasta",
                    )
                ],
                output_spec=[
                    ArtifactSpec(
                        name="embedding_artifact",
                        kind="embedding",
                        provider="embedding_runtime",
                        format="npz",
                    )
                ],
                metadata={"model_name": emb_model},
            )
            # Use the ConfigurableRuntimeAdapter by passing None
            self.register_model(emb_card, None)

        lambda_card = ModelCard(
            name="lambda",
            version="1.0",
            description="Lambda graph-based property predictor",
            execution=ExecutionSpec(mode="runtime", entrypoint="lambda_predict"),
            input_spec=[
                ArtifactSpec(
                    name="sequence_dataset",
                    kind="sequence",
                    provider="sequence_dataset",
                    format="fasta",
                ),
                ArtifactSpec(
                    name="grn_table",
                    kind="grn",
                    provider="grn_table",
                    format="csv",
                ),
                ArtifactSpec(
                    name="embedding_dataset",
                    kind="embedding",
                    provider="embedding_dataset",
                    format="npz",
                    optional=True,
                ),
            ],
            output_spec=[
                ArtifactSpec(
                    name="property_table",
                    kind="property",
                    provider="property_processor",
                    format="csv",
                    optional=True,
                )
            ],
            metadata={"default_run_id": "007061"},
        )
        self.register_model(lambda_card, LambdaAdapter(self))

        # ------------------------------------------------------------------
        # Placeholder model cards for bundled, yet-unwired models
        # ------------------------------------------------------------------

        ligand_mpnn_card = ModelCard(
            name="ligand_mpnn",
            version="0.1",
            description=(
                "LigandMPNN: design sequences conditioned on a protein PDB"
            ),
            execution=ExecutionSpec(
                mode="external_config",
                entrypoint="python",
                environment={"gpu": True},
            ),
            input_spec=[
                ArtifactSpec(
                    name="structure_pdb",
                    kind="structure",
                    provider="file_path",
                    format="pdb",
                ),
                ArtifactSpec(
                    name="ligand_molecule",
                    kind="ligand",
                    provider="ligand_file",
                    format="sdf",
                    optional=True,
                ),
                ArtifactSpec(
                    name="constraints",
                    kind="json",
                    provider="json_payload",
                    optional=True,
                ),
            ],
            output_spec=[],
        )
        self.register_model(ligand_mpnn_card, LigandMpnnAdapter(self))

        # Uni-Dock: GPU-accelerated docking (via unidock CLI)
        unidock_card = ModelCard(
            name="unidock",
            version="1.1",
            description="Uni-Dock: GPU-accelerated molecular docking",
            execution=ExecutionSpec(
                mode="external_config",
                entrypoint="unidock",
                environment={"gpu": True},
            ),
            input_spec=[
                ArtifactSpec(
                    name="receptor_pdb",
                    kind="structure",
                    provider="file_path",
                    format="pdb",
                ),
                ArtifactSpec(
                    name="ligand_file",
                    kind="ligand",
                    provider="ligand_file",
                    format="sdf",
                ),
                ArtifactSpec(
                    name="bbox",
                    kind="json",
                    provider="json_payload",
                    optional=True,
                ),
            ],
            output_spec=[],
        )
        self.register_model(unidock_card, UniDockAdapter(self))

        equibind_card = ModelCard(
            name="equibind",
            version="1.0",
            description=(
                "EquiBind: fast docking to predict ligand poses for a receptor"
            ),
            execution=ExecutionSpec(mode="external_config", entrypoint="python"),
            input_spec=[
                ArtifactSpec(
                    name="receptor_pdb",
                    kind="structure",
                    provider="file_path",
                    format="pdb",
                ),
                ArtifactSpec(
                    name="ligand_file",
                    kind="ligand",
                    provider="ligand_file",
                    format="sdf",
                ),
            ],
            output_spec=[
                ArtifactSpec(
                    name="docked_ligands",
                    kind="ligand",
                    provider="structure_exporter",
                    format="sdf",
                ),
                ArtifactSpec(
                    name="docking_scores",
                    kind="property",
                    provider="property_processor",
                    format="csv",
                    optional=True,
                ),
            ],
        )
        self.register_model(equibind_card, EquiBindAdapter(self))

        pocketdta_card = ModelCard(
            name="pocketdta",
            version="1.0",
            description=(
                "PocketDTA: multimodal DTA prediction using pocket and ligand"
            ),
            execution=ExecutionSpec(
                mode="external_config",
                entrypoint="python",
                environment={"gpu": True},
            ),
            input_spec=[
                ArtifactSpec(
                    name="dataset_dir",
                    kind="file",
                    provider="file_path",
                    format="dir",
                ),
                ArtifactSpec(
                    name="config_yaml",
                    kind="config",
                    provider="file_path",
                    format="yaml",
                ),
                ArtifactSpec(
                    name="checkpoint",
                    kind="file",
                    provider="file_path",
                    format="pth",
                    optional=True,
                ),
            ],
            output_spec=[],
        )
        self.register_model(pocketdta_card, PocketDtaAdapter(self))

        graphscore_dta_card = ModelCard(
            name="graphscoredta",
            version="1.0",
            description=(
                "GraphscoreDTA: GNN-based affinity prediction for pocket+ligand"
            ),
            execution=ExecutionSpec(
                mode="external_config",
                entrypoint="python",
                environment={"gpu": True},
            ),
            input_spec=[
                ArtifactSpec(
                    name="graphs_bin",
                    kind="file",
                    provider="file_path",
                    format="bin",
                ),
                ArtifactSpec(
                    name="labels_csv",
                    kind="file",
                    provider="file_path",
                    format="csv",
                ),
                ArtifactSpec(
                    name="vina_pkl",
                    kind="file",
                    provider="file_path",
                    format="pkl",
                ),
                ArtifactSpec(
                    name="model_path",
                    kind="file",
                    provider="file_path",
                    format="pth",
                ),
            ],
            output_spec=[],
        )
        self.register_model(graphscore_dta_card, GraphscoreDtaAdapter(self))

        pocket2mol_card = ModelCard(
            name="pocket2mol",
            version="1.0",
            description="Pocket2Mol: generate molecules for a protein pocket",
            execution=ExecutionSpec(mode="external_config", entrypoint="python", environment={"gpu": True}),
            input_spec=[
                ArtifactSpec(
                    name="structure_pdb",
                    kind="structure",
                    provider="file_path",
                    format="pdb",
                ),
                ArtifactSpec(
                    name="ligand_molecule",
                    kind="ligand",
                    provider="ligand_file",
                    format="sdf",
                    optional=True,
                ),
                ArtifactSpec(
                    name="bbox",
                    kind="json",
                    provider="json_payload",
                    optional=True,
                ),
            ],
            output_spec=[],
        )
        self.register_model(pocket2mol_card, Pocket2MolAdapter(self))

    def register_model(
        self, card: ModelCard, adapter: Optional[ModelAdapterBase]
    ) -> None:
        self.cards[card.name] = card
        if adapter is None:
            mode = (card.execution.mode or "").lower()
            if mode.startswith("runtime"):
                adapter = ConfigurableRuntimeAdapter(self)
            else:
                adapter = ConfigurableExternalAdapter(self)
        self.adapters[card.name] = adapter

    # ------------------------------------------------------------------
    # Artifact assembly
    # ------------------------------------------------------------------

    def assemble_inputs(
        self,
        card: ModelCard,
        request: ModelRequest,
    ) -> List[ArtifactBundle]:
        bundles: List[ArtifactBundle] = []
        for spec in card.input_spec:
            provider = self._artifact_providers.get(spec.provider)
            if provider is None:
                if spec.optional:
                    continue
                raise ValueError(f"Unsupported artifact provider: {spec.provider}")
            bundle = provider(spec, request)
            if bundle is None:
                if spec.optional:
                    continue
                raise ValueError(f"Required artifact '{spec.name}' could not be resolved")
            bundles.append(bundle)
        return bundles

    def _provide_sequence_dataset(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        dataset_name = request.get_input(spec.name)
        if not dataset_name:
            dataset_name = request.get_input("dataset")
        if not dataset_name:
            raise ValueError("Sequence dataset name not provided")

        seq_proc = SequenceProcessor()
        sequences = seq_proc.load_dataset(dataset_name)
        dataset_info = seq_proc.get_dataset_info(dataset_name)

        metadata = {
            "dataset": dataset_name,
            "sequences": sequences,
            "dataset_info": dataset_info,
        }

        artifact_rel = None
        if dataset_info:
            metadata_section = dataset_info.get("metadata", {})
            artifact_rel = metadata_section.get("artifact_path")
        if artifact_rel:
            path = Path(self.paths.data_root) / artifact_rel
        else:
            path = Path(seq_proc.path_fasta_dir) / f"{dataset_name}.fasta"

        return ArtifactBundle(spec=spec, path=path, metadata=metadata)

    def _provide_graph_entity(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        from protos.processing.graph.graph_processor import GraphProcessor

        entity = request.get_input(spec.name) or request.get_input("graph_entity")
        if not entity:
            raise ValueError("Graph entity name not provided")

        gp = GraphProcessor()
        info = gp.entity_registry.find_entity(entity, gp.processor_type)
        if not info:
            raise ValueError(f"Graph entity '{entity}' not found")
        path = Path(self.paths.data_root) / info.file_path
        payload = gp.load_graph(entity)
        metadata = {"entity": entity, "graph_metadata": payload.get("graph_metadata", {})}
        return ArtifactBundle(spec=spec, path=path, metadata=metadata)

    def _provide_ligand_file(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        value = request.get_input(spec.name) or request.get_input("ligand") or request.get_input("ligand_file")
        if not value:
            if spec.optional:
                return None  # type: ignore[return-value]
            raise ValueError("Ligand file path not provided")
        p = Path(value)
        if not p.is_absolute():
            if p.exists():
                path = p.resolve()
            else:
                path = Path(self.paths.data_root) / p
        else:
            path = p
        if not path.exists():
            raise FileNotFoundError(f"Ligand file not found: {path}")
        metadata = {"source": "user", "format": spec.format or path.suffix.lstrip(".")}
        return ArtifactBundle(spec=spec, path=path, metadata=metadata)

    def _provide_json_payload(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        value = request.get_input(spec.name) or request.get_input("payload") or request.get_input("metadata")
        payload_dir = Path(self.paths.data_root) / "models" / "_payloads"
        payload_dir.mkdir(parents=True, exist_ok=True)

        if isinstance(value, (dict, list)):
            path = payload_dir / f"{spec.name}.json"
            with open(path, "w", encoding="utf-8") as fh:
                json.dump(value, fh, indent=2)
            return ArtifactBundle(spec=spec, path=path, metadata={"source": "inline"})

        if isinstance(value, str):
            # JSON string → write to file
            trimmed = value.strip()
            if trimmed.startswith("{") or trimmed.startswith("["):
                path = payload_dir / f"{spec.name}.json"
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(value)
                return ArtifactBundle(spec=spec, path=path, metadata={"source": "inline"})

            # Treat as filesystem path
            p = Path(value)
            if not p.is_absolute():
                if p.exists():
                    path = p.resolve()
                else:
                    path = Path(self.paths.data_root) / p
            else:
                path = p
            if not path.exists():
                raise FileNotFoundError(f"JSON payload not found: {path}")
            return ArtifactBundle(spec=spec, path=path, metadata={"source": "file"})

        if spec.optional:
            return None  # type: ignore[return-value]
        raise ValueError(f"Unsupported json_payload for '{spec.name}'")

    def _provide_structure_entity(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        from protos.processing.structure import StructureProcessor

        name = request.get_input(spec.name) or request.get_input("structure") or request.get_input("receptor_structure")
        if not name:
            raise ValueError("Structure entity name not provided")

        sp = StructureProcessor()
        info = sp.entity_registry.find_entity(name, sp.processor_type)
        if not info:
            # Try auto-discovered PKL in cache
            pkl_path = sp.path_pkl_dir / f"{name}.pkl"
            if pkl_path.exists():
                sp._register_entity(name, pkl_path, {"auto_discovered": True})  # type: ignore[attr-defined]
                info = sp.entity_registry.find_entity(name, sp.processor_type)
        if not info:
            raise ValueError(f"Structure entity '{name}' not found")
        path = Path(self.paths.data_root) / info.file_path
        if not path.exists():
            raise FileNotFoundError(f"Structure file not found: {path}")
        return ArtifactBundle(spec=spec, path=path, metadata={"entity": name, "format": "pkl"})

    def _provide_file_path(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        value = request.get_input(spec.name)
        if not value:
            if spec.optional:
                return None  # type: ignore[return-value]
            raise ValueError(f"File path for '{spec.name}' not provided")
        p = Path(value)
        if not p.is_absolute():
            if p.exists():
                path = p.resolve()
            else:
                path = Path(self.paths.data_root) / p
        else:
            path = p
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        return ArtifactBundle(spec=spec, path=path, metadata={"source": "file"})

    def _provide_grn_table(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        table_name = request.get_input(spec.name) or request.get_input("grn_table")
        if not table_name:
            raise ValueError("GRN table name not provided")

        grn_proc = GRNProcessor()
        table = grn_proc.load_table(table_name)
        table_path = grn_proc.tables_dir / f"{table_name}.csv"

        metadata = {
            "table_name": table_name,
            "table": table,
        }

        return ArtifactBundle(spec=spec, path=table_path, metadata=metadata)

    def _provide_embedding_dataset(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        dataset_name = request.get_input(spec.name) or request.get_input("embedding_dataset")
        if not dataset_name:
            if spec.optional:
                return None  # type: ignore[return-value]
            raise ValueError("Embedding dataset name not provided")

        # Heuristic: dataset names often include model between double-underscores,
        # e.g., "<base>__ankh_large__per_residue". Use it to avoid a redundant
        # default-model initialization.
        initial_model: Optional[str] = None
        parts = dataset_name.split("__")
        if len(parts) >= 3:
            initial_model = parts[1] or None

        emb_proc = EmbeddingProcessor(model_name=initial_model) if initial_model else EmbeddingProcessor()
        dataset_info = emb_proc.get_dataset_info(dataset_name)
        if not dataset_info:
            raise ValueError(f"Embedding dataset '{dataset_name}' not found")

        meta = dataset_info.get("metadata", {})
        model_name = meta.get("models") or meta.get("model_name") or emb_proc.model_name
        embedding_type = meta.get("embedding_type", "per_residue")

        if model_name != emb_proc.model_name:
            emb_proc = EmbeddingProcessor(model_name=model_name)

        embeddings = emb_proc.load_embeddings(dataset_name)
        numpy_embeddings: Dict[str, np.ndarray] = {}
        for key, tensor in embeddings.items():
            if hasattr(tensor, "detach"):
                numpy_embeddings[key] = tensor.detach().cpu().numpy()
            else:
                numpy_embeddings[key] = np.asarray(tensor)

        artifact_rel = meta.get("artifact_path")
        if artifact_rel:
            path = Path(self.paths.data_root) / artifact_rel
        else:
            path = Path(emb_proc.datasets_dir) / dataset_name

        metadata = {
            "dataset": dataset_name,
            "model_name": model_name,
            "embedding_type": embedding_type,
            "embeddings": numpy_embeddings,
        }

        return ArtifactBundle(spec=spec, path=path, metadata=metadata)

    # ------------------------------------------------------------------
    # Output ingestion (registration)
    # ------------------------------------------------------------------
    def ingest_outputs(self, invocation: ModelInvocation) -> Dict[str, Any]:
        """Register model outputs into Protos datasets/entities.

        - Property tables: ensure datasets are recorded via PropertyProcessor
        - No-ops for other kinds until dedicated ingesters are added
        """
        summary: Dict[str, Any] = {"model": invocation.model, "ingested": []}

        # Ingest explicit artifacts bundled by adapters
        for bundle in invocation.outputs or []:
            if bundle.spec.kind == "property" and bundle.spec.provider == "property_processor":
                table_path = Path(bundle.path)
                table_name = table_path.stem
                try:
                    prop = PropertyProcessor()
                    # Load and save to update dataset metadata/index
                    df = prop.load_property_table(table_name)
                    prop.save_property_table(table_name)
                    summary["ingested"].append({
                        "type": "property_table",
                        "name": table_name,
                        "rows": int(len(df)) if df is not None else None,
                    })
                except Exception:
                    continue
            elif bundle.spec.kind == "ligand":
                # Register SDF artifact as a molecule record for discovery later
                try:
                    from protos.processing.molecule import MoleculeProcessor
                    mp = MoleculeProcessor()
                    rel = Path(bundle.path).relative_to(self.paths.data_root)
                    ent_name = Path(bundle.path).stem
                    mp.save_entity(ent_name, {"kind": "sdf", "file_path": str(rel)}, metadata={"source_model": invocation.model})
                    summary["ingested"].append({
                        "type": "ligand_sdf",
                        "name": ent_name,
                        "path": str(rel),
                    })
                except Exception:
                    pass

        # Also check adapter-provided metadata for common shortcuts
        if invocation.runtime and invocation.runtime.metadata:
            table_name = invocation.runtime.metadata.get("property_table")
            if table_name:
                try:
                    prop = PropertyProcessor()
                    df = prop.load_property_table(table_name)
                    prop.save_property_table(table_name)
                    summary["ingested"].append({
                        "type": "property_table",
                        "name": table_name,
                        "rows": int(len(df)) if df is not None else None,
                    })
                except Exception:
                    pass

        return summary


__all__ = [
    "ModelManager",
    "ModelAdapterBase",
    "ExternalJobAdapter",
    "RuntimeAdapter",
    "ConfigurableRuntimeAdapter",
    "ConfigurableExternalAdapter",
    "ModelRequest",
    "BoltzAdapter",
    "LambdaAdapter",
]

class UniDockAdapter(ExternalJobAdapter):
    """Prepare a Uni-Dock docking job.

    Inputs:
    - receptor_pdb: PDB/PDBQT/CIF path (if CIF, converted to PDB; if PDB, will attempt PDBQT)
    - ligand_file: SDF/MOL2/PDB file for docking (passed via --gpu_batch)

    Optional config keys:
    - center: [cx, cy, cz] (float list). If absent and ligand provided, centroid is used when possible.
    - size or bbox_size: single float for cubic box (default 22.5)
    - size_xyz: [sx, sy, sz] explicit sizes
    - scoring: vina|vinardo|ad4 (default vina)
    - search_mode: fast|balanced|detailed (if absent, uses exhaustiveness/max_step)
    - exhaustiveness: int (default 256)
    - max_step: int (default 10)
    - num_modes: int (default 10)
    - energy_range: float (default 3.0)
    - refine_step: int (default 5)
    - seed: int (default 181129)
    """

    def _add_unidock_tools_to_path(self) -> None:
        try:
            import sys  # noqa: F401
            base = Path(__file__).resolve().parent / "Uni-Dock" / "unidock_tools" / "src"
            if str(base) not in sys.path:
                sys.path.insert(0, str(base))
        except Exception:
            pass

    def _maybe_to_pdbqt(self, pdb_src: Path, ctx: "ModelRunContext") -> Path:
        # If already PDBQT, just copy into inputs dir for provenance
        if pdb_src.suffix.lower() == ".pdbqt":
            target = ctx.inputs_dir / pdb_src.name
            try:
                shutil.copy2(pdb_src, target)
                return target
            except Exception:
                return pdb_src

        # If PDB, attempt conversion via bundled unidock_tools (RDKit-based)
        if pdb_src.suffix.lower() == ".pdb":
            try:
                self._add_unidock_tools_to_path()
                from unidock_tools.modules.protein_prep.pdb2pdbqt import pdb2pdbqt  # type: ignore
                target = ctx.inputs_dir / f"{pdb_src.stem}.pdbqt"
                # First attempt: direct conversion
                try:
                    pdb2pdbqt(str(pdb_src), str(target))
                except Exception:
                    pass
                if target.exists():
                    return target

                # Fallback: sanitize PDB to protein-only ATOM records, then convert
                sanitized = ctx.inputs_dir / f"{pdb_src.stem}_protein_only.pdb"
                try:
                    with open(pdb_src, "r", encoding="utf-8", errors="ignore") as fh, open(
                        sanitized, "w", encoding="utf-8"
                    ) as out:
                        for line in fh:
                            if line.startswith("ATOM"):
                                out.write(line)
                        out.write("END\n")
                    pdb2pdbqt(str(sanitized), str(target))
                except Exception:
                    pass
                if target.exists():
                    return target
            except Exception:
                pass
            return pdb_src

        # Any other extension: return as-is
        return pdb_src

    def _center_from_ligand(self, ligand_path: Path) -> Optional[List[float]]:
        try:
            from rdkit import Chem  # type: ignore
            mol = None
            if ligand_path.suffix.lower() == ".sdf":
                suppl = Chem.SDMolSupplier(str(ligand_path), removeHs=False)
                mol = next((m for m in suppl if m is not None), None)
            elif ligand_path.suffix.lower() == ".mol2":
                mol = Chem.MolFromMol2File(str(ligand_path), removeHs=False)
            elif ligand_path.suffix.lower() == ".pdb":
                mol = Chem.MolFromPDBFile(str(ligand_path), removeHs=False)
            if mol is None or mol.GetNumAtoms() == 0:
                return None
            conf = mol.GetConformer()
            xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
            ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
            zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
            return [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
        except Exception:
            return None

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        receptor = self._require_bundle(inputs, "receptor_pdb")
        ligand = self._require_bundle(inputs, "ligand_file")
        bbox_bundle = next((b for b in inputs if b.spec.name == "bbox"), None)

        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        # Ensure receptor PDB (handles CIF→PDB) then convert to PDBQT when feasible
        pdb_src = Path(receptor.path)
        if pdb_src.suffix.lower() in {".cif", ".mmcif"}:
            pdb_src = self._ensure_pdb_file(pdb_src, ctx, purpose="unidock")
        else:
            # Copy file into inputs for provenance
            target = ctx.inputs_dir / pdb_src.name
            try:
                shutil.copy2(pdb_src, target)
                pdb_src = target
            except Exception:
                pass
        receptor_pdbqt = self._maybe_to_pdbqt(pdb_src, ctx)

        lig_src = Path(ligand.path)
        lig_target = ctx.inputs_dir / lig_src.name
        try:
            shutil.copy2(lig_src, lig_target)
        except Exception:
            lig_target = lig_src

        # Box center determination
        cfg = request.config or {}
        # If bbox artifact provided, allow it to override center/size
        if bbox_bundle is not None:
            try:
                with open(bbox_bundle.path, "r", encoding="utf-8") as fh:
                    bbox_cfg = json.load(fh)
                if isinstance(bbox_cfg, dict):
                    cfg = {**cfg, **bbox_cfg}
            except Exception:
                pass
        center = cfg.get("center")
        if not center:
            center = self._center_from_ligand(lig_target)
        if not center:
            # Fallback: CA centroid via gemmi if receptor is PDB/PDBQT with matching PDB available
            try:
                import gemmi  # type: ignore
                pdb_for_center = pdb_src
                if pdb_for_center.suffix.lower() not in {".pdb", ".cif", ".mmcif"}:
                    pdb_for_center = pdb_src.with_suffix(".pdb") if pdb_src.with_suffix(".pdb").exists() else pdb_src
                st = gemmi.read_structure(str(pdb_for_center))
                coords = []
                for model in st:
                    for chain in model:
                        for res in chain:
                            a = res.find_atom("CA")
                            if a is not None:
                                c = a.pos
                                coords.append((c.x, c.y, c.z))
                if coords:
                    cx = sum(c[0] for c in coords) / len(coords)
                    cy = sum(c[1] for c in coords) / len(coords)
                    cz = sum(c[2] for c in coords) / len(coords)
                    center = [cx, cy, cz]
            except Exception:
                center = None
        if not center:
            center = [0.0, 0.0, 0.0]

        # Box size
        size_xyz = cfg.get("size_xyz")
        if isinstance(size_xyz, (list, tuple)) and len(size_xyz) == 3:
            sx, sy, sz = [float(x) for x in size_xyz]
        else:
            # Support bbox.size as well
            size_val = (
                (cfg.get("size") if not isinstance(cfg.get("size"), (list, tuple)) else None)
                or cfg.get("bbox_size")
                or cfg.get("box_size")
                or 22.5
            )
            sx = sy = sz = float(size_val)

        scoring = str(cfg.get("scoring", "vina"))
        num_modes = int(cfg.get("num_modes", 10))
        energy_range = float(cfg.get("energy_range", 3.0))
        refine_step = int(cfg.get("refine_step", 5))
        seed = int(cfg.get("seed", 181129))

        command: List[str] = ["unidock"]
        if scoring.lower() == "ad4":
            try:
                self._add_unidock_tools_to_path()
                from unidock_tools.modules.docking.gen_grid import generate_ad4_grid  # type: ignore
                map_dir = ctx.inputs_dir / "mapdir"
                map_dir.mkdir(parents=True, exist_ok=True)
                map_prefix = generate_ad4_grid(
                    str(receptor_pdbqt), str(map_dir),
                    (float(center[0]), float(center[1]), float(center[2])),
                    (float(sx), float(sy), float(sz)),
                )
                command += ["--maps", str(map_prefix)]
            except Exception:
                # Fallback: still try receptor path
                command += ["--receptor", str(receptor_pdbqt)]
        else:
            # Require PDBQT for vina/vinardo
            if Path(receptor_pdbqt).suffix.lower() != ".pdbqt":
                raise ValueError(
                    "Uni-Dock (vina/vinardo) requires PDBQT receptor. Install RDKit and ensure unidock_tools available for PDB→PDBQT conversion."
                )
            command += ["--receptor", str(receptor_pdbqt)]

        command += [
            "--gpu_batch", str(lig_target),
            "--dir", str(ctx.outputs_dir),
            "--center_x", str(center[0]),
            "--center_y", str(center[1]),
            "--center_z", str(center[2]),
            "--size_x", str(sx),
            "--size_y", str(sy),
            "--size_z", str(sz),
            "--scoring", scoring,
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range),
            "--refine_step", str(refine_step),
            "--seed", str(seed),
            "--verbosity", "2",
            "--keep_nonpolar_H",
        ]

        search_mode = cfg.get("search_mode")
        if search_mode:
            command += ["--search_mode", str(search_mode)]
        else:
            exhaustiveness = int(cfg.get("exhaustiveness", 256))
            max_step = int(cfg.get("max_step", 10))
            command += ["--exhaustiveness", str(exhaustiveness), "--max_step", str(max_step)]

        # Package minimal manifest
        manifest = {
            "receptor": str(receptor_pdbqt),
            "ligand": str(lig_target),
            "center": center,
            "size": [sx, sy, sz],
            "scoring": scoring,
        }
        manifest_path = ctx.work_dir / "inputs.json"
        try:
            with open(manifest_path, "w", encoding="utf-8") as fh:
                json.dump(manifest, fh, indent=2)
            ctx.config_path = manifest_path
        except Exception:
            pass

        job = PreparedJob(
            command=command,
            working_dir=ctx.work_dir,
            artifacts=[receptor, ligand],
            metadata={
                "context": ctx.as_metadata(),
                "args": command[1:],
                "inputs": manifest,
            },
        )
        return job
class EquiBindAdapter(ExternalJobAdapter):
    """Prepare EquiBind input layout and config."""

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        receptor = self._require_bundle(inputs, "receptor_pdb")
        ligand = self._require_bundle(inputs, "ligand_file")

        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        complex_name = request.get_input("complex_name") or Path(receptor.path).stem
        complex_dir = ctx.inputs_dir / complex_name
        complex_dir.mkdir(parents=True, exist_ok=True)

        receptor_target = complex_dir / f"{complex_name}_protein.pdb"
        ligand_ext = Path(ligand.path).suffix.lower()
        ligand_target = complex_dir / f"{complex_name}_ligand{ligand_ext}"
        shutil.copy2(receptor.path, receptor_target)
        shutil.copy2(ligand.path, ligand_target)

        cfg = {
            "inference_path": str(ctx.inputs_dir),
            "output_directory": str(ctx.outputs_dir),
        }
        config_path = ctx.work_dir / "config.yaml"
        with open(config_path, "w", encoding="utf-8") as fh:
            yaml.safe_dump(cfg, fh, sort_keys=False)
        ctx.config_path = config_path

        inference_script = Path(__file__).resolve().parent / "EquiBind" / "inference.py"
        command = ["python", str(inference_script), "--config", str(config_path)]

        job = PreparedJob(
            command=command,
            working_dir=ctx.work_dir,
            artifacts=[
                receptor,
                ligand,
                ArtifactBundle(
                    spec=ArtifactSpec(
                        name="equibind_config",
                        kind="config",
                        provider="equibind_adapter",
                        format="yaml",
                    ),
                    path=config_path,
                    metadata={"inference_path": str(ctx.inputs_dir)},
                ),
            ],
            metadata={
                "context": ctx.as_metadata(),
                "inputs": {
                    "receptor_pdb": str(receptor_target),
                    "ligand_file": str(ligand_target),
                },
            },
        )
        return job

class PocketDtaAdapter(ExternalJobAdapter):
    """Prepare PocketDTA training/inference run using repo-expected inputs.

    Required inputs:
    - dataset_dir: directory containing process.csv and required assets (file_path)
    - config_yaml: path to PocketDTA OmegaConf YAML (file_path)
    Optional inputs:
    - checkpoint: model checkpoint path (file_path)
    """

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        dataset = self._require_bundle(inputs, "dataset_dir")
        config = self._require_bundle(inputs, "config_yaml")
        checkpoint = next((b for b in inputs if b.spec.name == "checkpoint"), None)

        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        # Package dataset and config into inputs for portability
        ds_src = Path(dataset.path)
        ds_dst = ctx.inputs_dir / "dataset"
        if ds_src.is_dir():
            if ds_dst.exists():
                shutil.rmtree(ds_dst)
            shutil.copytree(ds_src, ds_dst)
        else:
            raise ValueError("dataset_dir must be a directory including process.csv")

        cfg_src = Path(config.path)
        cfg_dst = ctx.inputs_dir / cfg_src.name
        shutil.copy2(cfg_src, cfg_dst)
        ctx.config_path = cfg_dst

        main_script = Path(__file__).resolve().parent / "PocketDTA" / "main.py"
        command = ["python", str(main_script), "--config", str(cfg_dst)]
        task = request.config.get("task")
        if task:
            command.extend(["--task", str(task)])

        job_artifacts = [dataset, config]
        if checkpoint is not None:
            # Copy checkpoint near inputs for portability
            c_src = Path(checkpoint.path)
            c_dst = ctx.inputs_dir / c_src.name
            shutil.copy2(c_src, c_dst)
            job_artifacts.append(checkpoint)

        job = PreparedJob(
            command=command,
            working_dir=ctx.work_dir,
            artifacts=job_artifacts,
            metadata={
                "context": ctx.as_metadata(),
                "dataset": str(ds_dst),
            },
        )
        return job

class GraphscoreDtaAdapter(ExternalJobAdapter):
    """Stage GraphscoreDTA test_set layout and build a predict job.

    Expects file paths for:
    - graphs_bin (all_in_one_graph_test2016.bin)
    - labels_csv (labels_test2016.csv)
    - vina_pkl (Vina_terms2016.pkl)
    - model_path (modelp.pth)
    """

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        graphs = self._require_bundle(inputs, "graphs_bin")
        labels = self._require_bundle(inputs, "labels_csv")
        vina = self._require_bundle(inputs, "vina_pkl")
        model = self._require_bundle(inputs, "model_path")

        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        # Predict.py expects ../test_set and ../model relative to its src dir
        repo_root = Path(__file__).resolve().parent / "GraphscoreDTA"
        test_set = repo_root / "test_set"
        model_dir = repo_root / "model"
        test_set.mkdir(parents=True, exist_ok=True)
        model_dir.mkdir(parents=True, exist_ok=True)
        graphs_target = test_set / "all_in_one_graph_test2016.bin"
        labels_target = test_set / "labels_test2016.csv"
        vina_target = test_set / "Vina_terms2016.pkl"
        shutil.copy2(graphs.path, graphs_target)
        shutil.copy2(labels.path, labels_target)
        shutil.copy2(vina.path, vina_target)
        model_target = model_dir / "modelp.pth"
        shutil.copy2(model.path, model_target)

        # Write a small manifest for tooling/debug
        manifest = {
            "test_set": {
                "graphs_bin": str(graphs_target),
                "labels_csv": str(labels_target),
                "vina_pkl": str(vina_target),
            },
            "model": str(model_target),
        }
        manifest_path = ctx.work_dir / "inputs.json"
        with open(manifest_path, "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)
        ctx.config_path = manifest_path

        predict_script_dir = repo_root / "src"
        predict_script = predict_script_dir / "predict.py"
        command = ["python", str(predict_script)]

        job = PreparedJob(
            command=command,
            working_dir=predict_script_dir,
            artifacts=[graphs, labels, vina, model],
            metadata={
                "context": ctx.as_metadata(),
                "staged": manifest,
            },
        )
        return job

class Pocket2MolAdapter(ExternalJobAdapter):
    """Prepare Pocket2Mol sampling run from protein structure and optional ligand.

    Inputs:
    - structure_pdb: file path (can be CIF; converted to PDB via gemmi)
    - ligand_molecule (optional): SDF/MOL2/PDB; used to compute pocket center
    - bbox/center (optional, via request.config): overrides auto center/size
    """

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        protein = self._require_bundle(inputs, "structure_pdb")
        ligand = next((b for b in inputs if b.spec.name in ("ligand_molecule", "ligand_file")), None)

        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        # Ensure PDB
        pdb_src = self._ensure_pdb_file(Path(protein.path), ctx, purpose="pocket2mol")

        # Compute pocket center
        center = request.config.get("center")
        bbox_size = float(request.config.get("bbox_size", 23.0))
        if center is None and ligand is not None:
            try:
                from rdkit import Chem  # type: ignore
                lig_path = Path(ligand.path)
                mol = None
                if lig_path.suffix.lower() == ".sdf":
                    suppl = Chem.SDMolSupplier(str(lig_path), removeHs=False)
                    mol = next((m for m in suppl if m is not None), None)
                elif lig_path.suffix.lower() == ".mol2":
                    mol = Chem.MolFromMol2File(str(lig_path), removeHs=False)
                elif lig_path.suffix.lower() == ".pdb":
                    mol = Chem.MolFromPDBFile(str(lig_path), removeHs=False)
                if mol is not None and mol.GetNumAtoms() > 0:
                    conf = mol.GetConformer()
                    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
                    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
                    zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
                    center = [sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)]
            except Exception:
                center = None
        if center is None:
            # Fallback: CA centroid via gemmi
            try:
                import gemmi  # type: ignore
                st = gemmi.read_structure(str(pdb_src))
                coords = []
                for model in st:
                    for chain in model:
                        for res in chain:
                            a = res.find_atom("CA")
                            if a is not None:
                                c = a.pos
                                coords.append((c.x, c.y, c.z))
                if coords:
                    cx = sum(c[0] for c in coords) / len(coords)
                    cy = sum(c[1] for c in coords) / len(coords)
                    cz = sum(c[2] for c in coords) / len(coords)
                    center = [cx, cy, cz]
                else:
                    center = [0.0, 0.0, 0.0]
            except Exception:
                center = [0.0, 0.0, 0.0]

        # Build command
        sample_script = Path(__file__).resolve().parent / "Pocket2Mol" / "sample_for_pdb.py"
        command = [
            "python",
            str(sample_script),
            "--pdb_path",
            str(pdb_src),
            "--center",
            f"{center[0]},{center[1]},{center[2]}",
            "--bbox_size",
            str(bbox_size),
            "--outdir",
            str(ctx.outputs_dir),
        ]
        # Config path in repo
        config_path = Path(__file__).resolve().parent / "Pocket2Mol" / "configs" / "sample_for_pdb.yml"
        if config_path.exists():
            command.extend(["--config", str(config_path)])
            ctx.config_path = config_path

        job = PreparedJob(
            command=command,
            working_dir=ctx.work_dir,
            artifacts=[protein] + ([ligand] if ligand else []),
            metadata={
                "context": ctx.as_metadata(),
                "center": center,
                "bbox_size": bbox_size,
            },
        )
        return job

class LigandMpnnAdapter(ExternalJobAdapter):
    """Prepare CLI invocation for LigandMPNN run.py with staged outputs.

    Required input spec:
    - structure_pdb (file_path)

    Optional config keys mapped to CLI:
    - seed, batch_size, number_of_batches, temperature, model_type
    - fixed_residues (list[str]), redesigned_residues (list[str])
    - bias_AA (str), omit_AA (str)
    """

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        pdb_bundle = self._require_bundle(inputs, "structure_pdb")
        ligand_bundle = next((b for b in inputs if b.spec.name in ("ligand_molecule", "ligand_file")), None)
        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        # Ensure PDB input exists; convert CIF→PDB if necessary
        pdb_src = self._ensure_pdb_file(Path(pdb_bundle.path), ctx, purpose="ligand_mpnn")

        # If ligand provided, merge into combined PDB
        if ligand_bundle is not None:
            try:
                from rdkit import Chem  # type: ignore
                lig_path = Path(ligand_bundle.path)
                mol = None
                if lig_path.suffix.lower() == ".sdf":
                    suppl = Chem.SDMolSupplier(str(lig_path), removeHs=False)
                    mol = next((m for m in suppl if m is not None), None)
                elif lig_path.suffix.lower() == ".mol2":
                    mol = Chem.MolFromMol2File(str(lig_path), removeHs=False)
                elif lig_path.suffix.lower() == ".pdb":
                    mol = Chem.MolFromPDBFile(str(lig_path), removeHs=False)
                if mol is None:
                    raise RuntimeError(f"Unable to read ligand file: {lig_path}")
                lig_block = Chem.MolToPDBBlock(mol)

                combined = ctx.inputs_dir / f"{pdb_src.stem}_with_ligand.pdb"
                with open(pdb_src, "r", encoding="utf-8", errors="ignore") as ph, open(
                    combined, "w", encoding="utf-8"
                ) as out:
                    for line in ph:
                        if not line.startswith("END"):
                            out.write(line)
                    # Append ligand atoms
                    for line in lig_block.splitlines(True):
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            out.write(line)
                    out.write("END\n")
                pdb_src = combined
            except Exception as exc:
                raise RuntimeError(f"Failed merging ligand into PDB: {exc}")

        # Compose CLI
        run_script = Path(__file__).resolve().parent / "LigandMPNN" / "run.py"
        cmd: List[str] = [
            "python",
            str(run_script),
            "--pdb_path",
            str(pdb_src),
            "--out_folder",
            str(ctx.outputs_dir),
        ]

        cfg = request.config or {}
        def add_flag(name: str, flag: str):
            val = cfg.get(name)
            if val is not None:
                cmd.extend([flag, str(val)])

        add_flag("seed", "--seed")
        add_flag("batch_size", "--batch_size")
        add_flag("number_of_batches", "--number_of_batches")
        add_flag("temperature", "--temperature")
        add_flag("model_type", "--model_type")

        # Residue selections as space-separated strings
        for key, flag in (("fixed_residues", "--fixed_residues"), ("redesigned_residues", "--redesigned_residues")):
            vals = cfg.get(key)
            if isinstance(vals, list):
                cmd.extend([flag, " ".join(str(v) for v in vals)])
            elif isinstance(vals, str):
                cmd.extend([flag, vals])

        add_flag("bias_AA", "--bias_AA")
        add_flag("omit_AA", "--omit_AA")

        # Build job
        job = PreparedJob(
            command=cmd,
            working_dir=ctx.work_dir,
            artifacts=[pdb_bundle],
            metadata={
                "context": ctx.as_metadata(),
                "args": cmd[2:],
            },
        )
        return job
