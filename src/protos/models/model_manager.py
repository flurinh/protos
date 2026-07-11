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
    IngestionSpec,
    JobResult,
    JobState,
    JobStatus,
    ModelBatch,
    ModelCard,
    ModelInvocation,
    PreparedJob,
    RuntimeResult,
)
from .job_client import (
    JobClient,
    ServerConfig,
    JobState as APIJobState,
    JobResult as APIJobResult,
    JobStatus as APIJobStatus,
)

# Lambda runtime imports are deprecated - models now run via remote/container execution
try:
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
except ImportError:
    _lambda_runtime = None
    InMemoryEmbeddingAdapter = None
    InMemoryGRNAdapter = None
    align_embeddings_to_grn = None
    build_grn_assignments = None
    build_property_rows = None
    copy_if_missing = None
    ensure_positional_map = None
    prepare_predictor = None
    build_positional_map_from_binding = None
    normalize_binding_config = None


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

    def _ensure_pdb_file(
        self, source: Path, ctx: "ModelRunContext", *, purpose: str
    ) -> Path:
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


# ---------------------------------------------------------------------------
# Job Execution Framework
# ---------------------------------------------------------------------------


class JobExecutor(ABC):
    """Abstract interface for executing prepared jobs.

    Implementations handle the actual execution of jobs on different backends:
    - DockerJobExecutor: Local Docker containers
    - Future: RemoteClusterExecutor for HPC/cloud clusters
    """

    @abstractmethod
    def submit(self, job: PreparedJob, model: str) -> JobState:
        """Submit a job for execution.

        Args:
            job: The prepared job to execute
            model: Name of the model (for tracking)

        Returns:
            JobState with job_id and initial status
        """

    @abstractmethod
    def status(self, job_id: str) -> JobState:
        """Get the current status of a job."""

    @abstractmethod
    def result(self, job_id: str) -> Optional[JobResult]:
        """Get the result of a completed job, or None if not ready."""

    @abstractmethod
    def cancel(self, job_id: str) -> bool:
        """Cancel a running job. Returns True if successfully cancelled."""

    @abstractmethod
    def list_jobs(self, status: Optional[JobStatus] = None) -> List[JobState]:
        """List jobs, optionally filtered by status."""


class DockerJobExecutor(JobExecutor):
    """Execute jobs in local Docker containers.

    This executor runs PreparedJob commands inside Docker containers with:
    - Volume mounts for input/output directories
    - GPU support when requested
    - Configurable image per model

    Job state is stored in job.json within each run directory:
    data/models/<model>/<run_id>/job.json

    After ingestion, job directories are cleaned up unless persistent=True.
    """

    # Default Docker images per model
    DEFAULT_IMAGES: Dict[str, str] = {
        "boltzgen": "protos/boltzgen:latest",
        "boltz": "protos/boltz:latest",
        "boltz2": "protos/boltz:latest",
        "lambda": "protos/lambda:latest",
    }

    def __init__(
        self,
        models_dir: Optional[Path] = None,
        default_image: str = "protos/base:latest",
        use_gpu: Optional[bool] = None,
    ) -> None:
        # All jobs stored under data/models/
        self.models_dir = models_dir or Path.home() / ".protos" / "models"
        self.models_dir.mkdir(parents=True, exist_ok=True)
        self.default_image = default_image
        # Auto-detect GPU if not specified
        self.use_gpu = use_gpu if use_gpu is not None else self._detect_gpu()
        self._jobs: Dict[str, JobState] = {}
        self._load_jobs()

    @property
    def jobs_dir(self) -> Path:
        """Backwards compatibility - returns models_dir."""
        return self.models_dir

    def _detect_gpu(self) -> bool:
        """Detect if NVIDIA GPU is available for Docker."""
        import subprocess
        try:
            # Check if nvidia-smi exists and works
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                capture_output=True,
                timeout=5,
            )
            if result.returncode == 0 and result.stdout.strip():
                # Also verify Docker can access GPU
                docker_result = subprocess.run(
                    ["docker", "run", "--rm", "--gpus", "all", "nvidia/cuda:12.2.2-base-ubuntu22.04", "nvidia-smi", "-L"],
                    capture_output=True,
                    timeout=30,
                )
                return docker_result.returncode == 0
        except Exception:
            pass
        return False

    def _load_jobs(self) -> None:
        """Load job states by scanning model directories for job.json files.

        Supports both new format (job.json in run directories) and legacy
        format (global job_states.json) for backwards compatibility.
        """
        # First, load from legacy job_states.json if it exists
        legacy_state_file = self.models_dir / "job_states.json"
        if legacy_state_file.exists():
            try:
                with open(legacy_state_file, "r") as f:
                    data = json.load(f)
                for job_data in data.get("jobs", []):
                    job_id = job_data["job_id"]
                    self._jobs[job_id] = self._deserialize_job_state(job_data)
            except Exception:
                pass

        # Then scan for job.json files in run directories
        if not self.models_dir.exists():
            return

        for model_dir in self.models_dir.iterdir():
            if not model_dir.is_dir():
                continue
            for run_dir in model_dir.iterdir():
                if not run_dir.is_dir():
                    continue
                job_json = run_dir / "job.json"
                if job_json.exists():
                    try:
                        state = self._read_job_state(job_json)
                        # Prefer job.json over legacy state
                        self._jobs[state.job_id] = state
                    except Exception:
                        pass

    def _write_job_state(self, path: Path, state: JobState) -> None:
        """Write job state to job.json in run directory."""
        data = self._serialize_job_state(state)
        with open(path, "w") as f:
            json.dump(data, f, indent=2, default=str)

    def _read_job_state(self, path: Path) -> JobState:
        """Read job state from job.json."""
        with open(path) as f:
            data = json.load(f)
        return self._deserialize_job_state(data)

    def _save_jobs(self) -> None:
        """Persist job states to their respective job.json files."""
        for state in self._jobs.values():
            job_json = state.prepared_job.working_dir / "job.json"
            if job_json.parent.exists():
                self._write_job_state(job_json, state)

    def _serialize_job_state(self, state: JobState) -> Dict[str, Any]:
        """Serialize JobState for persistence."""
        return {
            "job_id": state.job_id,
            "model": state.model,
            "status": state.status.value,
            "created_at": state.created_at.isoformat(),
            "submitted_at": state.submitted_at.isoformat() if state.submitted_at else None,
            "completed_at": state.completed_at.isoformat() if state.completed_at else None,
            "executor": state.executor,
            "error": state.error,
            "metadata": state.metadata,
            "working_dir": str(state.prepared_job.working_dir),
            "command": state.prepared_job.command,
            "run_id": state.prepared_job.run_id,
        }

    def _deserialize_job_state(self, data: Dict[str, Any]) -> JobState:
        """Deserialize JobState from persistence."""
        # For backwards compatibility, use job_id as run_id if run_id is missing
        run_id = data.get("run_id") or data["job_id"]
        return JobState(
            job_id=data["job_id"],
            model=data["model"],
            status=JobStatus(data["status"]),
            prepared_job=PreparedJob(
                run_id=run_id,
                command=data.get("command", []),
                working_dir=Path(data["working_dir"]),
                artifacts=[],
                metadata={},
            ),
            created_at=datetime.fromisoformat(data["created_at"]),
            submitted_at=datetime.fromisoformat(data["submitted_at"]) if data.get("submitted_at") else None,
            completed_at=datetime.fromisoformat(data["completed_at"]) if data.get("completed_at") else None,
            executor=data.get("executor", "docker"),
            error=data.get("error"),
            metadata=data.get("metadata", {}),
        )

    def _get_docker_image(self, model: str, config: Dict[str, Any]) -> str:
        """Get the Docker image for a model."""
        # Check config override first
        if "docker_image" in config:
            return config["docker_image"]
        return self.DEFAULT_IMAGES.get(model, self.default_image)

    def submit(self, job: PreparedJob, model: str, persistent: bool = False) -> JobState:
        """Submit a job to run in a Docker container.

        Args:
            job: The prepared job to execute (must have run_id set)
            model: Name of the model (for tracking)
            persistent: If False (default), job directory will be cleaned up
                       after successful ingestion. If True, keep all files.

        Returns:
            JobState with job_id (= run_id) and initial status
        """
        import subprocess
        import threading

        # Use run_id from PreparedJob as the canonical job identifier
        job_id = job.run_id
        docker_image = self._get_docker_image(model, job.metadata)

        # Create job state - working_dir is already under data/models/<model>/<run_id>/
        state = JobState(
            job_id=job_id,
            model=model,
            status=JobStatus.PENDING,
            prepared_job=job,
            executor="docker",
            metadata={
                "docker_image": docker_image,
                "use_gpu": self.use_gpu,
                "persistent": persistent,
                "job_name": job.metadata.get("job_name", job_id),
            },
        )
        self._jobs[job_id] = state

        # Write job state to job.json in the run directory
        self._write_job_state(job.working_dir / "job.json", state)

        # Build Docker command
        docker_cmd = self._build_docker_command(job, docker_image)

        # Store logs in the job's working directory (already under data/models/)
        stdout_file = job.working_dir / "stdout.log"
        stderr_file = job.working_dir / "stderr.log"

        def run_job():
            state.status = JobStatus.RUNNING
            state.submitted_at = datetime.now()
            self._save_jobs()

            start_time = datetime.now()
            try:
                with open(stdout_file, "w") as out, open(stderr_file, "w") as err:
                    proc = subprocess.run(
                        docker_cmd,
                        cwd=str(job.working_dir),
                        stdout=out,
                        stderr=err,
                        timeout=3600 * 24,  # 24 hour timeout
                    )

                duration = (datetime.now() - start_time).total_seconds()

                # Collect output files
                output_dir = Path(job.metadata.get("output_dir", job.working_dir / "outputs"))
                output_files = list(output_dir.glob("**/*")) if output_dir.exists() else []

                state.result = JobResult(
                    exit_code=proc.returncode,
                    stdout=stdout_file.read_text() if stdout_file.exists() else "",
                    stderr=stderr_file.read_text() if stderr_file.exists() else "",
                    output_dir=output_dir,
                    output_files=[f for f in output_files if f.is_file()],
                    duration_seconds=duration,
                )

                if proc.returncode == 0:
                    state.status = JobStatus.COMPLETED
                else:
                    state.status = JobStatus.FAILED
                    state.error = f"Exit code: {proc.returncode}"

            except subprocess.TimeoutExpired:
                state.status = JobStatus.FAILED
                state.error = "Job timed out after 24 hours"
            except Exception as e:
                state.status = JobStatus.FAILED
                state.error = str(e)

            state.completed_at = datetime.now()
            self._save_jobs()

        # Run in background thread
        thread = threading.Thread(target=run_job, daemon=True)
        thread.start()
        state.metadata["thread_id"] = thread.ident

        return state

    def _build_docker_command(self, job: PreparedJob, image: str) -> List[str]:
        """Build the Docker run command."""
        working_dir = job.working_dir
        output_dir = Path(job.metadata.get("output_dir", working_dir / "outputs"))

        cmd = ["docker", "run", "--rm"]

        # GPU support
        if self.use_gpu:
            cmd.extend(["--gpus", "all"])

        # Shared memory for PyTorch multiprocessing (default 64MB is too small)
        cmd.extend(["--shm-size", "8g"])

        # Volume mounts - working directory
        cmd.extend(["-v", f"{working_dir}:/workspace"])

        # Mount HuggingFace cache for model weights (shared across jobs)
        hf_cache = Path.home() / ".cache" / "huggingface"
        hf_cache.mkdir(parents=True, exist_ok=True)
        cmd.extend(["-v", f"{hf_cache}:/cache"])
        cmd.extend(["-e", "HF_HOME=/cache"])

        # Working directory inside container
        cmd.extend(["-w", "/workspace"])

        # User mapping to avoid permission issues
        cmd.extend(["--user", f"{os.getuid()}:{os.getgid()}"])

        # Image
        cmd.append(image)

        # The original command (already uses relative paths for container)
        cmd.extend(job.command)
        return cmd

    def status(self, job_id: str) -> JobState:
        """Get the current status of a job."""
        if job_id not in self._jobs:
            raise KeyError(f"Job '{job_id}' not found")
        return self._jobs[job_id]

    def result(self, job_id: str) -> Optional[JobResult]:
        """Get the result of a completed job."""
        state = self.status(job_id)
        return state.result

    def cancel(self, job_id: str) -> bool:
        """Cancel a running job."""
        if job_id not in self._jobs:
            return False

        state = self._jobs[job_id]
        if state.status not in (JobStatus.PENDING, JobStatus.RUNNING):
            return False

        # For Docker, we'd need to track container ID to stop it
        # For now, just mark as cancelled
        state.status = JobStatus.CANCELLED
        state.completed_at = datetime.now()
        self._save_jobs()
        return True

    def list_jobs(self, status: Optional[JobStatus] = None) -> List[JobState]:
        """List jobs, optionally filtered by status."""
        jobs = list(self._jobs.values())
        if status is not None:
            jobs = [j for j in jobs if j.status == status]
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)

    def cleanup_job(self, job_id: str, force: bool = False) -> bool:
        """Clean up a job's working directory after ingestion.

        Args:
            job_id: The job ID to clean up
            force: If True, clean up even if persistent=True or not ingested

        Returns:
            True if cleanup was performed
        """
        if job_id not in self._jobs:
            return False

        state = self._jobs[job_id]

        # Check if job should be kept
        if not force:
            if state.metadata.get("persistent", False):
                return False
            if not state.metadata.get("ingested", False):
                return False

        # Only clean completed/failed jobs
        if state.status not in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED):
            return False

        # Remove working directory
        working_dir = state.prepared_job.working_dir
        if working_dir and working_dir.exists():
            try:
                shutil.rmtree(working_dir)
            except Exception:
                return False

        # Mark as cleaned up
        state.metadata["cleaned_up"] = True
        state.metadata["cleaned_up_at"] = datetime.now().isoformat()
        self._save_jobs()

        return True

    def cleanup_completed_jobs(self, model: Optional[str] = None) -> int:
        """Clean up all completed and ingested job directories.

        Args:
            model: Optional model name to filter by

        Returns:
            Number of jobs cleaned up
        """
        cleaned = 0
        for job_id, state in list(self._jobs.items()):
            if model and state.model != model:
                continue
            if state.metadata.get("cleaned_up"):
                continue
            if self.cleanup_job(job_id):
                cleaned += 1
        return cleaned

    def remove_job(self, job_id: str) -> bool:
        """Remove a job from tracking (and optionally clean up files).

        This removes the job from the state file. Use cleanup_job first
        if you want to remove the files.
        """
        if job_id not in self._jobs:
            return False

        del self._jobs[job_id]
        self._save_jobs()
        return True


class ApptainerJobExecutor(JobExecutor):
    """Execute jobs in Apptainer/Singularity containers.

    This executor runs PreparedJob commands inside Apptainer containers with:
    - Bind mounts for input/output directories
    - GPU support via --nv flag
    - SIF images per model

    Job state is stored in job.json within each run directory:
    data/models/<model>/<run_id>/job.json

    Similar to DockerJobExecutor but uses Apptainer for HPC compatibility.
    """

    # Default SIF images per model (paths relative to protos/models)
    DEFAULT_SIFS: Dict[str, str] = {
        "rfdiffusion2": "RFdiffusion2/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif",
    }

    def __init__(
        self,
        models_dir: Optional[Path] = None,
        protos_models_dir: Optional[Path] = None,
        use_gpu: Optional[bool] = None,
    ) -> None:
        """Initialize Apptainer executor.

        Args:
            models_dir: Directory for job outputs (data/models/)
            protos_models_dir: Directory containing model installations (src/protos/models/)
            use_gpu: Whether to use GPU (--nv flag). Auto-detects if None.
        """
        self.models_dir = models_dir or Path.home() / ".protos" / "models"
        self.models_dir.mkdir(parents=True, exist_ok=True)
        self.protos_models_dir = protos_models_dir or Path(__file__).parent
        self.use_gpu = use_gpu if use_gpu is not None else self._detect_gpu()
        self._jobs: Dict[str, JobState] = {}
        self._load_jobs()

    def _detect_gpu(self) -> bool:
        """Detect if NVIDIA GPU is available for Apptainer."""
        import subprocess
        try:
            result = subprocess.run(
                ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
                capture_output=True,
                timeout=5,
            )
            return result.returncode == 0 and bool(result.stdout.strip())
        except Exception:
            return False

    def _load_jobs(self) -> None:
        """Load job states by scanning model directories for job.json files.

        Supports both new format (job.json in run directories) and legacy
        format (global apptainer_job_states.json) for backwards compatibility.
        """
        # First, load from legacy apptainer_job_states.json if it exists
        legacy_state_file = self.models_dir / "apptainer_job_states.json"
        if legacy_state_file.exists():
            try:
                with open(legacy_state_file, "r") as f:
                    data = json.load(f)
                for job_data in data.get("jobs", []):
                    job_id = job_data["job_id"]
                    self._jobs[job_id] = self._deserialize_job_state(job_data)
            except Exception:
                pass

        # Then scan for job.json files in run directories
        if not self.models_dir.exists():
            return

        for model_dir in self.models_dir.iterdir():
            if not model_dir.is_dir():
                continue
            for run_dir in model_dir.iterdir():
                if not run_dir.is_dir():
                    continue
                job_json = run_dir / "job.json"
                if job_json.exists():
                    try:
                        state = self._read_job_state(job_json)
                        # Only load apptainer jobs
                        if state.executor == "apptainer":
                            self._jobs[state.job_id] = state
                    except Exception:
                        pass

    def _write_job_state(self, path: Path, state: JobState) -> None:
        """Write job state to job.json in run directory."""
        data = self._serialize_job_state(state)
        with open(path, "w") as f:
            json.dump(data, f, indent=2, default=str)

    def _read_job_state(self, path: Path) -> JobState:
        """Read job state from job.json."""
        with open(path) as f:
            data = json.load(f)
        return self._deserialize_job_state(data)

    def _save_jobs(self) -> None:
        """Persist job states to their respective job.json files."""
        for state in self._jobs.values():
            job_json = state.prepared_job.working_dir / "job.json"
            if job_json.parent.exists():
                self._write_job_state(job_json, state)

    def _serialize_job_state(self, state: JobState) -> Dict[str, Any]:
        """Serialize JobState for persistence."""
        return {
            "job_id": state.job_id,
            "model": state.model,
            "status": state.status.value,
            "created_at": state.created_at.isoformat(),
            "submitted_at": state.submitted_at.isoformat() if state.submitted_at else None,
            "completed_at": state.completed_at.isoformat() if state.completed_at else None,
            "executor": state.executor,
            "error": state.error,
            "metadata": state.metadata,
            "working_dir": str(state.prepared_job.working_dir),
            "command": state.prepared_job.command,
            "run_id": state.prepared_job.run_id,
        }

    def _deserialize_job_state(self, data: Dict[str, Any]) -> JobState:
        """Deserialize JobState from persistence."""
        # For backwards compatibility, use job_id as run_id if run_id is missing
        run_id = data.get("run_id") or data["job_id"]
        return JobState(
            job_id=data["job_id"],
            model=data["model"],
            status=JobStatus(data["status"]),
            prepared_job=PreparedJob(
                run_id=run_id,
                command=data.get("command", []),
                working_dir=Path(data["working_dir"]),
                artifacts=[],
                metadata={},
            ),
            created_at=datetime.fromisoformat(data["created_at"]),
            submitted_at=datetime.fromisoformat(data["submitted_at"]) if data.get("submitted_at") else None,
            completed_at=datetime.fromisoformat(data["completed_at"]) if data.get("completed_at") else None,
            executor=data.get("executor", "apptainer"),
            error=data.get("error"),
            metadata=data.get("metadata", {}),
        )

    def _get_sif_path(self, model: str, config: Dict[str, Any]) -> Path:
        """Get the SIF container path for a model."""
        # Check config override first
        if "sif_path" in config:
            return Path(config["sif_path"])
        # Use default SIF path
        if model in self.DEFAULT_SIFS:
            return self.protos_models_dir / self.DEFAULT_SIFS[model]
        raise ValueError(f"No SIF container configured for model '{model}'")

    def submit(self, job: PreparedJob, model: str, persistent: bool = False) -> JobState:
        """Submit a job to run in an Apptainer container.

        Args:
            job: The prepared job to execute (must have run_id set)
            model: Name of the model (for tracking)
            persistent: If False (default), job directory will be cleaned up
                       after successful ingestion

        Returns:
            JobState with job_id (= run_id) and initial status
        """
        import subprocess
        import threading

        # Use run_id from PreparedJob as the canonical job identifier
        job_id = job.run_id
        sif_path = self._get_sif_path(model, job.metadata)

        if not sif_path.exists():
            raise FileNotFoundError(f"SIF container not found: {sif_path}")

        state = JobState(
            job_id=job_id,
            model=model,
            status=JobStatus.PENDING,
            prepared_job=job,
            executor="apptainer",
            metadata={
                "sif_path": str(sif_path),
                "use_gpu": self.use_gpu,
                "persistent": persistent,
                "job_name": job.metadata.get("job_name", job_id),
            },
        )
        self._jobs[job_id] = state

        # Write job state to job.json in the run directory
        self._write_job_state(job.working_dir / "job.json", state)

        # Build Apptainer command
        apptainer_cmd = self._build_apptainer_command(job, sif_path)

        # Store logs
        stdout_file = job.working_dir / "stdout.log"
        stderr_file = job.working_dir / "stderr.log"

        def run_job():
            state.status = JobStatus.RUNNING
            state.submitted_at = datetime.now()
            self._save_jobs()

            start_time = datetime.now()
            try:
                with open(stdout_file, "w") as out, open(stderr_file, "w") as err:
                    proc = subprocess.run(
                        apptainer_cmd,
                        cwd=str(job.working_dir),
                        stdout=out,
                        stderr=err,
                        timeout=3600 * 24,  # 24 hour timeout
                    )

                duration = (datetime.now() - start_time).total_seconds()

                output_dir = Path(job.metadata.get("output_dir", job.working_dir / "outputs"))
                output_files = list(output_dir.glob("**/*")) if output_dir.exists() else []

                state.result = JobResult(
                    exit_code=proc.returncode,
                    stdout=stdout_file.read_text() if stdout_file.exists() else "",
                    stderr=stderr_file.read_text() if stderr_file.exists() else "",
                    output_dir=output_dir,
                    output_files=[f for f in output_files if f.is_file()],
                    duration_seconds=duration,
                )

                if proc.returncode == 0:
                    state.status = JobStatus.COMPLETED
                else:
                    state.status = JobStatus.FAILED
                    state.error = f"Exit code: {proc.returncode}"

            except subprocess.TimeoutExpired:
                state.status = JobStatus.FAILED
                state.error = "Job timed out after 24 hours"
            except Exception as e:
                state.status = JobStatus.FAILED
                state.error = str(e)

            state.completed_at = datetime.now()
            self._save_jobs()

        thread = threading.Thread(target=run_job, daemon=True)
        thread.start()
        state.metadata["thread_id"] = thread.ident

        return state

    def _build_apptainer_command(self, job: PreparedJob, sif_path: Path) -> List[str]:
        """Build the Apptainer exec command."""
        working_dir = job.working_dir.absolute()
        output_dir = Path(job.metadata.get("output_dir", working_dir / "outputs")).absolute()

        cmd = ["apptainer", "exec"]

        # GPU support
        if self.use_gpu:
            cmd.append("--nv")

        # Working directory inside container (for models that need to run from their install dir)
        container_workdir = job.metadata.get("container_workdir")
        if container_workdir:
            cmd.extend(["--pwd", container_workdir])

        # Environment variables
        for key, value in job.metadata.get("env", {}).items():
            cmd.extend(["--env", f"{key}={value}"])

        # Bind mounts - working directory and output
        cmd.extend(["-B", f"{working_dir}:{working_dir}"])
        if output_dir.parent != working_dir:
            cmd.extend(["-B", f"{output_dir}:{output_dir}"])

        # Additional bind mounts from metadata
        for bind in job.metadata.get("bind_mounts", []):
            cmd.extend(["-B", bind])

        # SIF container
        cmd.append(str(sif_path))

        # The command to run inside container
        cmd.extend(job.command)

        return cmd

    def status(self, job_id: str) -> JobState:
        """Get the current status of a job."""
        if job_id not in self._jobs:
            raise KeyError(f"Job '{job_id}' not found")
        return self._jobs[job_id]

    def result(self, job_id: str) -> Optional[JobResult]:
        """Get the result of a completed job."""
        state = self.status(job_id)
        return state.result

    def cancel(self, job_id: str) -> bool:
        """Cancel a running job."""
        if job_id not in self._jobs:
            return False
        state = self._jobs[job_id]
        if state.status not in (JobStatus.PENDING, JobStatus.RUNNING):
            return False
        state.status = JobStatus.CANCELLED
        state.completed_at = datetime.now()
        self._save_jobs()
        return True

    def list_jobs(self, status: Optional[JobStatus] = None) -> List[JobState]:
        """List jobs, optionally filtered by status."""
        jobs = list(self._jobs.values())
        if status is not None:
            jobs = [j for j in jobs if j.status == status]
        return sorted(jobs, key=lambda j: j.created_at, reverse=True)


class APIJobExecutor(JobExecutor):
    """Execute jobs via the Protos Job Server API.

    This executor sends jobs to a FastAPI server (running locally in Docker
    or remotely) instead of executing Docker containers directly. This allows:
    - Consistent interface for local development and production
    - Future support for remote execution APIs
    - Job persistence and status tracking via the server

    Configuration is read from data/models/server_config.yaml.
    """

    def __init__(
        self,
        data_root: Optional[Path] = None,
        config: Optional[ServerConfig] = None,
    ) -> None:
        self.data_root = data_root or Path.home() / ".protos"
        self._client = JobClient(data_root=self.data_root, config=config)

    def submit(self, job: PreparedJob, model: str, persistent: bool = False) -> JobState:
        """Submit a job to the API server.

        Uses the new submit_run_dir method to package and upload the entire
        run directory, ensuring the run_id is preserved as the job_id.

        Args:
            job: The prepared job to execute (must have run_id set)
            model: Name of the model (for tracking)
            persistent: If True, keep job files after completion

        Returns:
            JobState with job_id (= run_id) and initial status
        """
        # Write job.json to the run directory before uploading
        job_json_path = job.working_dir / "job.json"
        job_data = {
            "run_id": job.run_id,
            "command": job.command,
            "metadata": {
                **job.metadata,
                "persistent": persistent,
            },
        }
        with open(job_json_path, "w") as f:
            json.dump(job_data, f, indent=2, default=str)

        # Use submit_run_dir to package and upload the entire run directory
        api_state = self._client.submit_run_dir(job.working_dir, model)

        # Convert API response to internal JobState
        return self._convert_api_state(api_state, job)

    def _convert_api_state(
        self,
        api_state: APIJobState,
        job: Optional[PreparedJob] = None,
    ) -> JobState:
        """Convert API JobState to internal JobState."""
        # Map API status to internal status
        status_map = {
            APIJobStatus.PENDING: JobStatus.PENDING,
            APIJobStatus.RUNNING: JobStatus.RUNNING,
            APIJobStatus.COMPLETED: JobStatus.COMPLETED,
            APIJobStatus.FAILED: JobStatus.FAILED,
            APIJobStatus.CANCELLED: JobStatus.CANCELLED,
        }

        # Create a PreparedJob placeholder if not provided
        if job is None:
            working_dir = Path(api_state.working_dir) if api_state.working_dir else Path("/tmp")
            # Use job_id as run_id for API jobs
            run_id = api_state.metadata.get("run_id", api_state.job_id)
            job = PreparedJob(
                run_id=run_id,
                command=[],
                working_dir=working_dir,
                artifacts=[],
                metadata={},
            )

        return JobState(
            job_id=api_state.job_id,
            model=api_state.model,
            status=status_map.get(api_state.status, JobStatus.PENDING),
            prepared_job=job,
            created_at=datetime.fromisoformat(api_state.created_at),
            submitted_at=datetime.fromisoformat(api_state.submitted_at) if api_state.submitted_at else None,
            completed_at=datetime.fromisoformat(api_state.completed_at) if api_state.completed_at else None,
            executor="api",
            error=api_state.error,
            metadata=api_state.metadata,
        )

    def status(self, job_id: str) -> JobState:
        """Get the current status of a job from the API."""
        api_state = self._client.status(job_id)
        return self._convert_api_state(api_state)

    def result(self, job_id: str) -> Optional[JobResult]:
        """Get the result of a completed job from the API."""
        state = self._client.status(job_id)
        if state.status not in (APIJobStatus.COMPLETED, APIJobStatus.FAILED):
            return None

        api_result = self._client.result(job_id)
        return JobResult(
            exit_code=api_result.exit_code,
            stdout=api_result.stdout,
            stderr=api_result.stderr,
            output_dir=Path(state.working_dir) / "predictions" if state.working_dir else Path("."),
            output_files=[Path(f) for f in api_result.output_files],
            duration_seconds=api_result.duration_seconds,
        )

    def cancel(self, job_id: str) -> bool:
        """Cancel a running job via the API."""
        return self._client.cancel(job_id)

    def list_jobs(self, status: Optional[JobStatus] = None) -> List[JobState]:
        """List jobs from the API server."""
        # Map internal status to API status
        api_status = None
        if status is not None:
            status_map = {
                JobStatus.PENDING: APIJobStatus.PENDING,
                JobStatus.RUNNING: APIJobStatus.RUNNING,
                JobStatus.COMPLETED: APIJobStatus.COMPLETED,
                JobStatus.FAILED: APIJobStatus.FAILED,
                JobStatus.CANCELLED: APIJobStatus.CANCELLED,
            }
            api_status = status_map.get(status)

        api_jobs = self._client.list_jobs(status=api_status)
        return [self._convert_api_state(j) for j in api_jobs]

    def wait_for_completion(
        self,
        job_id: str,
        timeout: float = 3600,
        poll_interval: float = 5.0,
    ) -> JobState:
        """Wait for a job to complete via the API."""
        api_state = self._client.wait_for_completion(
            job_id, timeout=timeout, poll_interval=poll_interval
        )
        return self._convert_api_state(api_state)

    def download_outputs(self, job_id: str, dest_dir: Path) -> List[Path]:
        """Download all output files from a completed job.

        Args:
            job_id: The job ID
            dest_dir: Local directory to download files to

        Returns:
            List of downloaded file paths
        """
        api_result = self._client.result(job_id)
        downloaded = []
        for file_path in api_result.output_files:
            dest = dest_dir / file_path
            self._client.download_file(job_id, file_path, dest)
            downloaded.append(dest)
        return downloaded


class ModelRunContext:
    """Filesystem context for a single model run under data/models/<model>/.

    Creates a stable directory layout expected by our workflow:
    - work_dir:    data/models/<model>/<run_id>/
    - inputs_dir:  .../inputs/
    - outputs_dir: .../outputs/
    - config_path: optional, when a config file is created

    The run_id is the canonical identifier - either a timestamp (default) or
    a user-provided name. When run_prefix is provided, it's incorporated into
    the run_id as "{prefix}_{timestamp}" for better organization.
    """

    def __init__(
        self,
        paths: ProtosPaths,
        card: ModelCard,
        run_prefix: Optional[str] = None,
        run_id: Optional[str] = None,
    ) -> None:
        self.paths = paths
        self.card = card
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if run_id:
            self.run_id = run_id
        elif run_prefix:
            self.run_id = f"{run_prefix}_{timestamp}"
        else:
            self.run_id = timestamp
        self.work_dir = (
            Path(self.paths.data_root)
            / "models"
            / card.name
            / self.run_id
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

        kwargs = {
            "card": card,
            "request": request,
            "inputs": inputs,
            "paths": self.paths,
        }
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
        meta.update(
            {
                "work_dir": str(ctx.work_dir),
                "outputs_dir": str(ctx.outputs_dir),
            }
        )
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
                payload = self._build_config_payload(
                    card, request, inputs, packaged_inputs
                )
                with open(config_path, "w", encoding="utf-8") as fh:
                    json.dump(payload, fh, indent=2)
            else:
                config_path = working_dir / "config.yaml"
                payload = self._build_config_payload(
                    card, request, inputs, packaged_inputs
                )
                with open(config_path, "w", encoding="utf-8") as fh:
                    yaml.safe_dump(payload, fh, sort_keys=False)
            cmd.append(str(config_path))
        ctx.config_path = config_path

        job = PreparedJob(
            run_id=ctx.run_id,
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

    @staticmethod
    def _build_config_payload(
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
        packaged: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        artifacts = {
            b.spec.name: {
                "path": str(b.path),
                "kind": b.spec.kind,
                "format": b.spec.format,
            }
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

        # Create run context with descriptive prefix
        run_prefix = f"{entity_name}_{config_id}"
        ctx = ModelRunContext(self.paths, card, run_prefix=run_prefix)
        ctx.create()

        yaml_data = self._generate_yaml(entity_sequences, request.config)
        yaml_path = ctx.work_dir / "config.yaml"
        fasta_path = ctx.inputs_dir / "sequences.fasta"

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

        metadata_path = ctx.work_dir / "metadata.json"
        metadata_content = {
            "entity": entity_name,
            "config_id": config_id,
            "mutations": mutations,
            "dataset": sequence_bundle.metadata.get("dataset"),
        }
        with open(metadata_path, "w", encoding="utf-8") as fh:
            json.dump(metadata_content, fh, indent=2)

        # Build command with appropriate flags - use relative paths for container
        command = ["boltz", "predict", "config.yaml", "--out_dir", "outputs"]

        # Add --use_msa_server unless MSA paths are provided
        has_msa = self._config_has_msa(request.config)
        if not has_msa and request.config.get("use_msa_server", True):
            command.append("--use_msa_server")

        # Optional: use potentials for better physical quality
        if request.config.get("use_potentials", False):
            command.append("--use_potentials")

        # Optional: recycling and sampling steps
        if "recycling_steps" in request.config:
            command.extend(["--recycling_steps", str(request.config["recycling_steps"])])
        if "sampling_steps" in request.config:
            command.extend(["--sampling_steps", str(request.config["sampling_steps"])])
        if "diffusion_samples" in request.config:
            command.extend(["--diffusion_samples", str(request.config["diffusion_samples"])])

        # Optional: device configuration
        if "devices" in request.config:
            command.extend(["--devices", str(request.config["devices"])])
        if "accelerator" in request.config:
            command.extend(["--accelerator", request.config["accelerator"]])

        job = PreparedJob(
            run_id=ctx.run_id,
            command=command,
            working_dir=ctx.work_dir,
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
                "output_dir": str(ctx.outputs_dir),
            },
        )
        return job

    def _config_has_msa(self, config: MutableMapping[str, Any]) -> bool:
        """Check if config provides MSA paths."""
        # Check if any sequence override has msa field
        for override in config.get("sequence_overrides", {}).values():
            if "fields" in override and "msa" in override["fields"]:
                return True
        # Check if sequences section directly includes msa
        for seq in config.get("sequences", []):
            if isinstance(seq, dict):
                for seq_type, seq_data in seq.items():
                    if isinstance(seq_data, dict) and "msa" in seq_data:
                        return True
        return False

    def _config_identifier(self, mutations: Iterable[Dict[str, Any]]) -> str:
        labels = [
            mut.get("name") or f"{mut['position']}{mut['mutant']}" for mut in mutations
        ]
        return "_".join(labels)

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

        # Ligand handling - supports multiple ligands via 'ligands' list
        # or single ligand via 'ligand'/'ligand_smiles' for backwards compatibility
        ligands_config = config.get("ligands", [])
        affinity_binder = None

        # Backwards compatibility: single ligand config
        if not ligands_config:
            ligand_cfg = config.get("ligand")
            if isinstance(ligand_cfg, dict):
                ligand_id = ligand_cfg.get("id") or ligand_cfg.get("name") or "LIG"
                ligand_smiles = ligand_cfg.get("smiles")
                if ligand_smiles:
                    ligands_config.append({"id": ligand_id, "smiles": ligand_smiles})
                    affinity_binder = ligand_id
            else:
                ligand_id = config.get("ligand_id")
                ligand_smiles = config.get("ligand_smiles")
                if ligand_smiles:
                    ligands_config.append({"id": ligand_id or "LIG", "smiles": ligand_smiles})
                    affinity_binder = ligand_id or "LIG"

        # Add all ligands to sequences
        for lig in ligands_config:
            lig_id = lig.get("id") or lig.get("name") or "LIG"
            lig_smiles = lig.get("smiles")
            if lig_smiles:
                seq_entries.append(
                    {
                        "ligand": {
                            "id": lig_id,
                            "smiles": lig_smiles,
                        }
                    }
                )
                # First non-cofactor ligand becomes affinity binder
                if affinity_binder is None and not lig.get("is_cofactor", False):
                    affinity_binder = lig_id

        yaml_sections["sequences"] = seq_entries

        for section in ("constraints", "properties"):
            if section in config:
                yaml_sections[section] = deepcopy(config[section])

        additional = config.get("additional_sections", {})
        for key, value in additional.items():
            yaml_sections[key] = deepcopy(value)

        # Default properties: affinity when ligand present and no explicit properties provided
        if affinity_binder and "properties" not in yaml_sections:
            yaml_sections["properties"] = [
                {"affinity": {"binder": affinity_binder}}
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
            entry_body["id"] = (
                identifier if isinstance(identifier, list) else [identifier]
            )

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


class BoltzGenAdapter(ExternalJobAdapter):
    """Adapter for preparing BoltzGen protein design configuration jobs.

    BoltzGen uses a different YAML schema than Boltz, focused on protein design
    with support for:
    - Designed protein sequences with length ranges (e.g., "80..140")
    - Structure files (.cif/.pdb) as design targets
    - Binding site specifications
    - Secondary structure constraints
    - Bond constraints for stapled peptides and disulfides
    """

    MODEL_DIR = "boltzgen"

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        config = request.config
        job_name = config.get("job_name", "design_job")

        # Create run context with job_name as prefix
        ctx = ModelRunContext(self.paths, card, run_prefix=job_name)
        ctx.create()

        yaml_data = self._generate_yaml(config)
        yaml_path = ctx.work_dir / "config.yaml"

        # Copy any referenced structure files to the inputs directory
        self._copy_structure_files(yaml_data, ctx.inputs_dir, config)

        with open(yaml_path, "w", encoding="utf-8") as fh:
            yaml.dump(
                yaml_data,
                fh,
                Dumper=BoltzYamlDumper,
                default_flow_style=False,
                sort_keys=False,
                indent=2,
            )

        metadata_path = ctx.work_dir / "metadata.json"
        metadata_content = {
            "job_name": job_name,
            "config": {k: v for k, v in config.items() if k != "entities"},
        }
        with open(metadata_path, "w", encoding="utf-8") as fh:
            json.dump(metadata_content, fh, indent=2)

        # Use relative paths for container execution
        # config.yaml is in working_dir, output goes to outputs/
        command = [
            "boltzgen", "run",
            "config.yaml",
            "--output", "outputs",
        ]

        # Optional: number of designs (direct CLI arg)
        if "num_designs" in config:
            command.extend(["--num_designs", str(config["num_designs"])])

        # Optional: diffusion batch size (direct CLI arg)
        if "diffusion_batch_size" in config:
            command.extend(["--diffusion_batch_size", str(config["diffusion_batch_size"])])

        # Optional: design step config overrides via --config design key=value
        # These are internal model parameters, not direct CLI args
        design_overrides = []
        if "recycling_steps" in config:
            design_overrides.append(f"recycling_steps={config['recycling_steps']}")
        if "sampling_steps" in config:
            design_overrides.append(f"sampling_steps={config['sampling_steps']}")
        if "diffusion_samples" in config:
            design_overrides.append(f"diffusion_samples={config['diffusion_samples']}")
        if design_overrides:
            command.extend(["--config", "design"] + design_overrides)

        # Optional: device configuration
        if "devices" in config:
            command.extend(["--devices", str(config["devices"])])
        if "accelerator" in config:
            command.extend(["--accelerator", config["accelerator"]])

        # Protocol selection
        if "protocol" in config:
            command.extend(["--protocol", config["protocol"]])

        job = PreparedJob(
            run_id=ctx.run_id,
            command=command,
            working_dir=ctx.work_dir,
            artifacts=[
                ArtifactBundle(
                    spec=ArtifactSpec(
                        name="boltzgen_config",
                        kind="config",
                        provider="boltzgen_adapter",
                        format="yaml",
                    ),
                    path=yaml_path,
                    metadata={"job_name": job_name},
                ),
            ],
            metadata={
                "job_name": job_name,
                "output_dir": str(ctx.outputs_dir),
                "docker_image": config.get("docker_image", "protos/boltzgen:latest"),
            },
        )
        return job

    def _copy_structure_files(
        self,
        yaml_data: OrderedDict,
        input_dir: Path,
        config: MutableMapping[str, Any],
    ) -> None:
        """Copy referenced structure files to the input directory."""
        entities = yaml_data.get("entities", [])
        structure_dir = config.get("structure_dir")

        for entity in entities:
            if "file" in entity:
                file_spec = entity["file"]
                orig_path = file_spec.get("path", "")
                if orig_path:
                    # Try to resolve the structure file
                    src_path = None
                    if structure_dir:
                        candidate = Path(structure_dir) / orig_path
                        if candidate.exists():
                            src_path = candidate
                    if src_path is None:
                        candidate = Path(self.paths.data_root) / "structure" / "mmcif" / orig_path
                        if candidate.exists():
                            src_path = candidate
                    if src_path is None:
                        candidate = Path(orig_path)
                        if candidate.exists():
                            src_path = candidate

                    if src_path and src_path.exists():
                        dst_path = input_dir / Path(orig_path).name
                        if not dst_path.exists():
                            shutil.copy2(src_path, dst_path)
                        # Update the path in yaml to be relative
                        file_spec["path"] = Path(orig_path).name

    def _generate_yaml(
        self,
        config: MutableMapping[str, Any],
    ) -> OrderedDict:
        """Generate BoltzGen YAML configuration.

        Args:
            config: Configuration dict that may contain:
                - entities: List of entity specifications
                - designed_proteins: List of designed protein specs
                - target_structure: Path to target structure file
                - constraints: Bond constraints
                - And other BoltzGen-specific options

        Returns:
            OrderedDict with the BoltzGen YAML structure.
        """
        yaml_sections: OrderedDict[str, Any] = OrderedDict()

        # Build entities section
        entities = self._prepare_entities_section(config)
        if entities:
            yaml_sections["entities"] = entities

        # Add constraints if present
        constraints = config.get("constraints")
        if constraints:
            yaml_sections["constraints"] = deepcopy(constraints)

        return yaml_sections

    def _prepare_entities_section(
        self,
        config: MutableMapping[str, Any],
    ) -> List[Dict[str, Any]]:
        """Prepare the entities section of the BoltzGen YAML.

        Handles three types of entities:
        1. Designed proteins with sequence patterns (e.g., "80..140", "3..5C6C3")
        2. Fixed proteins with explicit sequences
        3. File references for target structures
        4. Ligands (CCD codes or SMILES)
        """
        entities: List[Dict[str, Any]] = []

        # If explicit entities are provided, use them directly
        if "entities" in config:
            return deepcopy(config["entities"])

        # Handle designed proteins
        designed = config.get("designed_proteins", [])
        for dp in designed:
            entity: Dict[str, Any] = OrderedDict()
            entity["id"] = dp.get("id", "A")

            # Handle sequence specification
            seq = dp.get("sequence")
            if seq:
                entity["sequence"] = seq
            else:
                # Generate length range pattern
                min_len = dp.get("min_length", 80)
                max_len = dp.get("max_length", 140)
                entity["sequence"] = f"{min_len}..{max_len}"

            entities.append({"protein": entity})

        # Handle fixed proteins (non-designed)
        fixed = config.get("fixed_proteins", [])
        for fp in fixed:
            entity = OrderedDict()
            entity["id"] = fp.get("id", "X")
            entity["sequence"] = fp.get("sequence", "")
            entities.append({"protein": entity})

        # Handle target structure files
        target = config.get("target_structure")
        if target:
            file_entity: Dict[str, Any] = OrderedDict()
            file_entity["path"] = target if isinstance(target, str) else target.get("path")

            # Include specific chains
            include = target.get("include") if isinstance(target, dict) else None
            if include:
                file_entity["include"] = include

            # Binding types specification
            binding = target.get("binding_types") if isinstance(target, dict) else None
            if binding:
                file_entity["binding_types"] = binding

            # Structure groups
            groups = target.get("structure_groups") if isinstance(target, dict) else None
            if groups:
                file_entity["structure_groups"] = groups

            # Design regions in target
            design = target.get("design") if isinstance(target, dict) else None
            if design:
                file_entity["design"] = design

            # Secondary structure constraints
            ss = target.get("secondary_structure") if isinstance(target, dict) else None
            if ss:
                file_entity["secondary_structure"] = ss

            entities.append({"file": file_entity})

        # Handle ligands
        ligands = config.get("ligands", [])
        for lig in ligands:
            lig_entity: Dict[str, Any] = OrderedDict()
            lig_entity["id"] = lig.get("id", "LIG")
            if "ccd" in lig:
                lig_entity["ccd"] = lig["ccd"]
            elif "smiles" in lig:
                lig_entity["smiles"] = lig["smiles"]
            entities.append({"ligand": lig_entity})

        return entities

    def generate_design_config(
        self,
        designed_protein_id: str = "B",
        min_length: int = 80,
        max_length: int = 140,
        target_structure: Optional[str] = None,
        target_chain: Optional[str] = None,
        ligand_ccd: Optional[str] = None,
        constraints: Optional[List[Dict[str, Any]]] = None,
    ) -> OrderedDict:
        """Convenience method to generate a simple design configuration.

        Args:
            designed_protein_id: Chain ID for the designed protein
            min_length: Minimum length for designed protein
            max_length: Maximum length for designed protein
            target_structure: Path to target structure file (.cif or .pdb)
            target_chain: Chain ID to use from target structure
            ligand_ccd: CCD code for ligand (if any)
            constraints: List of bond constraints

        Returns:
            OrderedDict with BoltzGen YAML configuration
        """
        config: Dict[str, Any] = {
            "designed_proteins": [{
                "id": designed_protein_id,
                "min_length": min_length,
                "max_length": max_length,
            }],
        }

        if target_structure:
            target_config: Dict[str, Any] = {"path": target_structure}
            if target_chain:
                target_config["include"] = [{"chain": {"id": target_chain}}]
            config["target_structure"] = target_config

        if ligand_ccd:
            config["ligands"] = [{"id": "LIG", "ccd": ligand_ccd}]

        if constraints:
            config["constraints"] = constraints

        return self._generate_yaml(config)


class LambdaAdapter(RuntimeAdapter):
    """Adapter executing Lambda predictions within the current process or via Docker.

    Supports two execution modes:
    - Runtime mode (default): Runs in-process, requires torch/torch_geometric locally
    - Docker mode: Prepares inputs and runs via Docker container

    Set config["use_docker"] = True to use Docker execution.
    """

    DEFAULT_RUN_ID = "007061"
    DEFAULT_BINDING_CONFIG = Path("grn/configs/binding_domain2.json")
    DEFAULT_POSITIONAL_MAP = Path("grn/configs/final_mapping7.csv")

    def __init__(self, manager: "ModelManager") -> None:
        super().__init__(manager)
        self.lambda_root = Path(__file__).resolve().parent / "lambda"
        self._protos_data_root = Path(__file__).resolve().parents[3] / "data"

    def _prepare(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> ModelInvocation:
        """Prepare Lambda invocation - runtime or Docker based on config."""
        use_docker = request.config.get("use_docker", False)

        if use_docker:
            # Use Docker-based execution
            job = self.build_job(card, request, inputs)
            return ModelInvocation(
                model=card.name,
                card=card,
                job=job,
                inputs=inputs,
                metadata=dict(job.metadata) if job.metadata else {},
            )
        else:
            # Use in-process runtime (original behavior)
            return super()._prepare(card, request, inputs)

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
        embedding_bundle = next(
            (b for b in inputs if b.spec.name == "embedding_dataset"), None
        )

        sequences: Dict[str, str] = sequence_bundle.metadata.get("sequences", {})
        grn_table: pd.DataFrame = grn_bundle.metadata.get("table")
        dataset_name: str = sequence_bundle.metadata.get("dataset") or "lambda_dataset"

        logger.info(
            "[lambda] Sequences: %d | GRN rows: %d",
            len(sequences),
            int(len(grn_table) if grn_table is not None else 0),
        )

        embeddings_map: Dict[str, np.ndarray] = {}
        embedding_dataset_name: Optional[str] = None
        embedding_model = str(request.config.get("embedding_model", "ankh_large"))
        embedding_type = str(request.config.get("embedding_type", "per_residue"))
        ingest_embeddings = bool(request.config.get("ingest_embeddings", False))

        if not sequences:
            raise ValueError("Lambda adapter requires sequence data")
        if grn_table is None or grn_table.empty:
            raise ValueError("Lambda adapter requires a populated GRN table")

        protein_family = request.get_input("protein_family") or request.config.get(
            "protein_family"
        )
        if not protein_family:
            raise ValueError("Lambda adapter requires 'protein_family' input")

        # Resolve embeddings: either from provided dataset or compute via embedding card
        if embedding_bundle is not None:
            logger.info(
                "[lambda] Using provided embedding dataset '%s'",
                embedding_bundle.metadata.get("dataset"),
            )
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

            status = (
                emb_inv.runtime.outputs.get("status")
                if emb_inv.runtime.outputs
                else None
            )
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
                logger.info(
                    "[lambda] Loaded embeddings: %d entities from NPZ",
                    len(embeddings_map),
                )
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
                    logger.warning(
                        "[lambda] Failed to ingest embeddings dataset: %s", exc
                    )

        if not embeddings_map:
            raise ValueError("Lambda adapter requires embedding tensors")

        assignments = build_grn_assignments(grn_table)
        grn_dict, aligned_embeddings = align_embeddings_to_grn(
            assignments, embeddings_map
        )
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
            logger.info(
                "[lambda] Resolving resource '%s' -> fallback '%s' -> target '%s'",
                default_rel,
                fallback,
                candidate,
            )
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
        logger.info(
            "[lambda] Resolving positional map '%s' -> fallback '%s' -> target '%s'",
            default_rel,
            fallback,
            candidate,
        )
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

    def build_job(
        self,
        card: ModelCard,
        request: "ModelRequest",
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        """Prepare Lambda inputs for Docker-based execution.

        Creates an input directory with:
        - aligned_embeddings.npz: Per-residue embeddings keyed by sequence ID
        - grn_table.csv: GRN annotations with entity_name as index
        - config.json: Runtime configuration

        Returns a PreparedJob with Docker run command.
        """
        import logging

        logger = logging.getLogger("LambdaAdapter")

        # Create run context
        job_name = request.config.get("job_name", "lambda_job")
        ctx = ModelRunContext(self.paths, card, run_prefix=job_name)
        ctx.create()
        job_dir = ctx.work_dir
        input_dir = ctx.inputs_dir
        output_dir = ctx.outputs_dir

        # Resolve inputs
        sequence_bundle = self._require_bundle(inputs, "sequence_dataset")
        grn_bundle = self._require_bundle(inputs, "grn_table")
        embedding_bundle = next(
            (b for b in inputs if b.spec.name == "embedding_dataset"), None
        )

        sequences: Dict[str, str] = sequence_bundle.metadata.get("sequences", {})
        grn_table: pd.DataFrame = grn_bundle.metadata.get("table")
        dataset_name: str = sequence_bundle.metadata.get("dataset") or "lambda_dataset"

        if not sequences:
            raise ValueError("Lambda adapter requires sequence data")
        if grn_table is None or grn_table.empty:
            raise ValueError("Lambda adapter requires a populated GRN table")

        protein_family = request.get_input("protein_family") or request.config.get(
            "protein_family", "gpcr_a"
        )

        # Get or compute embeddings
        embeddings_map: Dict[str, np.ndarray] = {}
        embedding_model = str(request.config.get("embedding_model", "ankh_large"))
        embedding_type = str(request.config.get("embedding_type", "per_residue"))

        if embedding_bundle is not None:
            embeddings_map = embedding_bundle.metadata.get("embeddings", {})
            logger.info("[lambda] Using provided embeddings: %d entities", len(embeddings_map))
        else:
            # Try to compute embeddings via embedding card
            logger.info("[lambda] Computing embeddings with %s", embedding_model)
            try:
                emb_inv = self.manager.prepare(
                    f"embedding_{embedding_model}",
                    inputs={"sequence_dataset": dataset_name},
                    config={"embedding_type": embedding_type},
                )
                if emb_inv.runtime:
                    for b in emb_inv.runtime.artifacts:
                        if b.spec.kind == "embedding" and b.path.exists():
                            npz = np.load(b.path, allow_pickle=False)
                            for key in npz.files:
                                embeddings_map[key] = np.asarray(npz[key])
                            break
            except Exception as exc:
                logger.warning("[lambda] Failed to compute embeddings: %s", exc)

        if not embeddings_map:
            raise ValueError(
                "Lambda adapter requires embeddings. Provide embedding_dataset "
                "or ensure embedding model dependencies are available."
            )

        # Align embeddings to GRN positions
        assignments = build_grn_assignments(grn_table)
        grn_dict, aligned_embeddings = align_embeddings_to_grn(assignments, embeddings_map)

        if not aligned_embeddings:
            raise RuntimeError("No embeddings aligned to GRN positions")

        # Write aligned embeddings
        embeddings_path = input_dir / "aligned_embeddings.npz"
        np.savez(embeddings_path, **aligned_embeddings)
        logger.info("[lambda] Wrote embeddings: %s", embeddings_path)

        # Write GRN table (filter to available entities)
        available_ids = sorted(grn_dict.keys())
        filtered_table = grn_table.loc[grn_table.index.intersection(available_ids)]
        grn_path = input_dir / "grn_table.csv"
        filtered_table.to_csv(grn_path)
        logger.info("[lambda] Wrote GRN table: %s (%d rows)", grn_path, len(filtered_table))

        # Write config
        run_config = {
            "run_id": str(request.config.get("run_id", self.DEFAULT_RUN_ID)),
            "protein_family": protein_family,
            "batch_size": int(request.config.get("batch_size", 8)),
            "collect_attention": bool(request.config.get("collect_attention", False)),
            "verbose": bool(request.config.get("verbose", False)),
        }
        config_path = input_dir / "config.json"
        with open(config_path, "w", encoding="utf-8") as fh:
            json.dump(run_config, fh, indent=2)

        # Build Docker command
        docker_image = request.config.get("docker_image", "protos/lambda:latest")
        command = [
            "docker", "run", "--rm",
            "-v", f"{input_dir}:/srv/run/input:ro",
            "-v", f"{output_dir}:/srv/run/output",
            docker_image,
            "--input", "/srv/run/input",
            "--outputs", "/srv/run/output",
            "--config", "/srv/run/input/config.json",
        ]

        # Add GPU support if requested
        if request.config.get("use_gpu", False):
            command.insert(3, "--gpus")
            command.insert(4, "all")

        job = PreparedJob(
            run_id=ctx.run_id,
            command=command,
            working_dir=job_dir,
            artifacts=[
                ArtifactBundle(
                    spec=ArtifactSpec(
                        name="lambda_input",
                        kind="config",
                        provider="lambda_adapter",
                        format="directory",
                    ),
                    path=input_dir,
                    metadata={
                        "job_name": job_name,
                        "embeddings_file": str(embeddings_path),
                        "grn_table_file": str(grn_path),
                        "config_file": str(config_path),
                    },
                ),
            ],
            metadata={
                "job_name": job_name,
                "docker_image": docker_image,
                "input_dir": str(input_dir),
                "output_dir": str(output_dir),
                "run_config": run_config,
            },
        )
        return job


class ModelManager:
    """Coordinate model invocations based on ModelCards and adapters."""

    _DATASET_FIELD_MAP = {
        "sequence": "sequence_dataset",
        "structure": "structure_dataset",
        "graph": "graph_dataset",
        "property": "property_table",
    }

    def __init__(
        self,
        data_root: Optional[Path] = None,
        executor: Optional[JobExecutor] = None,
        use_api: Optional[bool] = None,
    ) -> None:
        """Initialize the ModelManager.

        Args:
            data_root: Root directory for protos data (default: ~/.protos)
            executor: Custom JobExecutor instance (overrides use_api)
            use_api: If True, use APIJobExecutor (FastAPI server).
                    If False, use DockerJobExecutor (direct Docker).
                    If None (default), read from server_config.yaml.
        """
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

        # Job executor for running prepared jobs
        # All jobs stored under data/models/<model>/jobs/
        models_dir = Path(self.paths.data_root) / "models"
        if executor is not None:
            self._executor = executor
        else:
            # Determine executor type from config or parameter
            config_path = models_dir / "server_config.yaml"
            server_config = ServerConfig.from_file(config_path)

            # use_api parameter overrides config file
            should_use_api = use_api if use_api is not None else (server_config.mode == "remote")

            if should_use_api:
                self._executor = APIJobExecutor(
                    data_root=Path(self.paths.data_root),
                    config=server_config,
                )
            else:
                self._executor = DockerJobExecutor(models_dir=models_dir)

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
        logger.info(
            "[manager] Preparing model '%s' mode=%s", model_name, card.execution.mode
        )
        invocation = adapter.prepare(card, request)
        logger.info(
            "[manager] Prepared invocation for '%s' external=%s runtime=%s",
            model_name,
            bool(invocation.job),
            bool(invocation.runtime),
        )
        invocation.metadata.update(request.metadata)
        return invocation

    # ------------------------------------------------------------------
    # Job Execution API
    # ------------------------------------------------------------------

    def submit_job(
        self,
        invocation: ModelInvocation,
        persistent: bool = False,
    ) -> JobState:
        """Submit a prepared job for execution.

        Args:
            invocation: A ModelInvocation with a prepared job (from prepare())
            persistent: If True, keep job directory after ingestion.
                       If False (default), clean up after successful ingestion.

        Returns:
            JobState with job_id for tracking

        Raises:
            ValueError: If the invocation doesn't have a prepared job
        """
        if not invocation.is_external():
            raise ValueError(
                "Cannot submit job: invocation is not an external job. "
                "Use runtime execution instead."
            )

        # Check if job requires Apptainer execution
        executor_type = invocation.job.metadata.get("executor", "docker")
        if executor_type == "apptainer":
            # Use ApptainerJobExecutor for this job
            if not hasattr(self, "_apptainer_executor"):
                models_dir = Path(self.paths.data_root) / "models"
                self._apptainer_executor = ApptainerJobExecutor(
                    models_dir=models_dir,
                    protos_models_dir=Path(__file__).parent,
                )
            return self._apptainer_executor.submit(invocation.job, invocation.model, persistent=persistent)

        return self._executor.submit(invocation.job, invocation.model, persistent=persistent)

    def prepare_and_submit(
        self,
        model_name: str,
        *,
        inputs: Optional[MutableMapping[str, Any]] = None,
        config: Optional[MutableMapping[str, Any]] = None,
        metadata: Optional[MutableMapping[str, Any]] = None,
        persistent: bool = False,
    ) -> JobState:
        """Prepare and immediately submit a job for execution.

        Convenience method combining prepare() and submit_job().

        Args:
            persistent: If True, keep job directory after ingestion.

        Returns:
            JobState with job_id for tracking
        """
        invocation = self.prepare(model_name, inputs=inputs, config=config, metadata=metadata)
        return self.submit_job(invocation, persistent=persistent)

    def job_status(self, job_id: str) -> JobState:
        """Get the current status of a submitted job.

        Args:
            job_id: The job ID returned from submit_job()

        Returns:
            JobState with current status and metadata
        """
        return self._executor.status(job_id)

    def job_result(self, job_id: str) -> Optional[JobResult]:
        """Get the result of a completed job.

        Args:
            job_id: The job ID returned from submit_job()

        Returns:
            JobResult if job is complete, None otherwise
        """
        return self._executor.result(job_id)

    def cancel_job(self, job_id: str) -> bool:
        """Cancel a running job.

        Args:
            job_id: The job ID to cancel

        Returns:
            True if successfully cancelled
        """
        return self._executor.cancel(job_id)

    def list_jobs(
        self,
        status: Optional[JobStatus] = None,
        model: Optional[str] = None,
    ) -> List[JobState]:
        """List submitted jobs.

        Args:
            status: Filter by job status (optional)
            model: Filter by model name (optional)

        Returns:
            List of JobState objects, sorted by creation time (newest first)
        """
        jobs = self._executor.list_jobs(status)
        if model is not None:
            jobs = [j for j in jobs if j.model == model]
        return jobs

    def wait_for_job(
        self,
        job_id: str,
        timeout_seconds: float = 3600,
        poll_interval: float = 5.0,
    ) -> JobState:
        """Wait for a job to complete.

        Args:
            job_id: The job ID to wait for
            timeout_seconds: Maximum time to wait (default 1 hour)
            poll_interval: Time between status checks (default 5 seconds)

        Returns:
            Final JobState when job completes

        Raises:
            TimeoutError: If job doesn't complete within timeout
        """
        import time

        start = time.time()
        while time.time() - start < timeout_seconds:
            state = self._executor.status(job_id)
            if state.status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED):
                return state
            time.sleep(poll_interval)

        raise TimeoutError(f"Job {job_id} did not complete within {timeout_seconds} seconds")

    def run_and_ingest(
        self,
        invocation: ModelInvocation,
        timeout_seconds: float = 3600,
    ) -> Dict[str, Any]:
        """Submit a job, wait for completion, and ingest results.

        Convenience method for synchronous job execution with result registration.

        Args:
            invocation: A prepared ModelInvocation
            timeout_seconds: Maximum time to wait for completion

        Returns:
            Ingestion summary from ingest_outputs()

        Raises:
            RuntimeError: If job fails or times out
        """
        state = self.submit_job(invocation)
        final_state = self.wait_for_job(state.job_id, timeout_seconds=timeout_seconds)

        if final_state.status == JobStatus.FAILED:
            raise RuntimeError(f"Job failed: {final_state.error}")
        if final_state.status == JobStatus.CANCELLED:
            raise RuntimeError("Job was cancelled")

        # Update invocation with output artifacts
        if final_state.result and final_state.result.output_files:
            invocation.outputs = self._parse_job_outputs(invocation, final_state.result)

        return self.ingest_outputs(invocation)

    def _parse_job_outputs(
        self,
        invocation: ModelInvocation,
        result: JobResult,
    ) -> List[ArtifactBundle]:
        """Parse job output files into ArtifactBundles.

        Override in model-specific adapters for custom parsing.
        """
        bundles: List[ArtifactBundle] = []

        for output_file in result.output_files:
            suffix = output_file.suffix.lower()

            # Determine artifact kind from file extension
            if suffix in (".cif", ".pdb", ".mmcif"):
                kind = "structure"
            elif suffix == ".csv":
                kind = "property"
            elif suffix == ".sdf":
                kind = "ligand"
            elif suffix in (".npz", ".npy"):
                kind = "embedding"
            elif suffix == ".json":
                kind = "metadata"
            else:
                kind = "file"

            bundle = ArtifactBundle(
                spec=ArtifactSpec(
                    name=output_file.stem,
                    kind=kind,
                    provider=f"{invocation.model}_output",
                    format=suffix.lstrip("."),
                ),
                path=output_file,
                metadata={"source_job": invocation.model},
            )
            bundles.append(bundle)

        return bundles

    @property
    def executor(self) -> JobExecutor:
        """Access the job executor."""
        return self._executor

    @executor.setter
    def executor(self, executor: JobExecutor) -> None:
        """Set a custom job executor."""
        self._executor = executor

    def prepare_input(
        self,
        model_name: str,
        *,
        entity_name: Optional[str] = None,
        entity_format: str = "sequence",
        dataset_name: Optional[str] = None,
        dataset_input_key: Optional[str] = None,
        inputs: Optional[MutableMapping[str, Any]] = None,
        config: Optional[MutableMapping[str, Any]] = None,
        metadata: Optional[MutableMapping[str, Any]] = None,
    ) -> ModelInvocation:
        """Convenience wrapper for preparing a model for a single entity."""

        resolved_inputs: Dict[str, Any] = dict(inputs or {})
        if entity_name:
            resolved_inputs.setdefault("entity", entity_name)

        dataset_field = dataset_input_key or self._dataset_field_for_format(
            entity_format
        )
        if dataset_name:
            if not dataset_field:
                raise ValueError(
                    f"No dataset input mapping available for format '{entity_format}'"
                )
            resolved_inputs.setdefault(dataset_field, dataset_name)

        return self.prepare(
            model_name,
            inputs=resolved_inputs,
            config=config,
            metadata=metadata,
        )

    def prepare_batch(
        self,
        model_name: str,
        entity_configs: Iterable[MutableMapping[str, Any]],
        *,
        batch_name: Optional[str] = None,
        default_entity_format: str = "sequence",
        base_config: Optional[MutableMapping[str, Any]] = None,
        batch_metadata: Optional[MutableMapping[str, Any]] = None,
    ) -> ModelBatch:
        """Normalize a batch of entity configurations for downstream execution."""

        normalized_inputs: List[Dict[str, Any]] = []
        for index, entry in enumerate(entity_configs, start=1):
            entity = entry.get("entity")
            if not entity:
                raise ValueError(f"Batch entry {index} is missing 'entity'")

            entry_format = str(entry.get("format") or default_entity_format).lower()
            entry_inputs = dict(entry.get("inputs", {}))
            dataset_name = entry.get("dataset_name")
            dataset_field = entry.get(
                "dataset_input_key"
            ) or self._dataset_field_for_format(entry_format)
            if dataset_name and dataset_field:
                entry_inputs.setdefault(dataset_field, dataset_name)

            final_config = self._merge_dicts(base_config, entry.get("config"))
            entry_metadata = dict(entry.get("metadata", {}))

            normalized_inputs.append(
                {
                    "entity": entity,
                    "format": entry_format,
                    "inputs": entry_inputs,
                    "config": final_config,
                    "metadata": entry_metadata,
                }
            )

        batch_label = (
            batch_name
            or f"{model_name}_batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )
        payload_metadata = dict(batch_metadata or {})
        payload_metadata.setdefault("entry_count", len(normalized_inputs))

        return ModelBatch(
            name=batch_label,
            model=model_name,
            inputs=normalized_inputs,
            metadata=payload_metadata,
        )

    def _dataset_field_for_format(self, entity_format: Optional[str]) -> Optional[str]:
        if not entity_format:
            return None
        return self._DATASET_FIELD_MAP.get(str(entity_format).lower())

    @staticmethod
    def _merge_dicts(
        base: Optional[MutableMapping[str, Any]],
        override: Optional[MutableMapping[str, Any]],
    ) -> Dict[str, Any]:
        result: Dict[str, Any] = dict(base or {})
        if override:
            result.update(override)
        return result

    def prepare_boltz_mutations(
        self,
        dataset_name: str,
        mutations: List[Dict[str, Any]],
        *,
        base_config: Optional[MutableMapping[str, Any]] = None,
        model_name: str = "boltz2",
    ) -> List[ModelInvocation]:
        """Prepare Boltz jobs for a list of mutation configurations.

        Args:
            dataset_name: Name of the sequence dataset containing base entities.
            mutations: List of dictionaries with at minimum an "entity" key and
                optional "mutations" / "config" overrides.
            base_config: Shared configuration options applied to every job before
                per-entry overrides.
            model_name: Registered model name to delegate to (defaults to "boltz2").
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

            entry_model = entry.get("model_name", model_name) or "boltz2"

            entry_metadata = dict(entry.get("metadata") or {})
            entry_metadata.setdefault("entity", entity)
            entry_metadata.setdefault("dataset", dataset_name)
            entry_metadata["mutation_entry"] = entry

            invocation = self.prepare(
                entry_model,
                inputs={"sequence_dataset": dataset_name, "entity": entity},
                config=config,
                metadata=entry_metadata,
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

        boltzgen_card = ModelCard(
            name="boltzgen",
            version="1",
            description="BoltzGen protein design configuration",
            execution=ExecutionSpec(mode="external_config", entrypoint="boltz design"),
            input_spec=[
                ArtifactSpec(
                    name="structure_file",
                    kind="structure",
                    provider="structure_loader",
                    format="cif",
                    optional=True,
                )
            ],
            output_spec=[
                ArtifactSpec(
                    name="designed_structures",
                    kind="structure",
                    provider="boltzgen_output",
                    format="cif",
                ),
                ArtifactSpec(
                    name="design_metrics",
                    kind="property",
                    provider="boltzgen_output",
                    format="csv",
                ),
            ],
            ingestion_spec=[
                # Register final ranked designs (in subdirectories)
                IngestionSpec(
                    output_type="structure",
                    file_pattern="predictions/final_ranked_designs/final_*_designs/rank*.cif",
                    processor="structure",
                    name_template="{job_name}_final_{stem}",
                    params={"register_entity": True, "copy_to_mmcif": True},
                ),
                # Register intermediate ranked designs
                IngestionSpec(
                    output_type="structure",
                    file_pattern="predictions/final_ranked_designs/intermediate_*_designs/rank*.cif",
                    processor="structure",
                    name_template="{job_name}_ranked_{stem}",
                    params={"register_entity": True, "copy_to_mmcif": True},
                    required=False,
                ),
                # Register design metrics as property table
                IngestionSpec(
                    output_type="property",
                    file_pattern="predictions/final_ranked_designs/*_metrics*.csv",
                    processor="property",
                    name_template="{job_name}_metrics",
                    params={"merge_existing": False},
                ),
                # Register all design metrics
                IngestionSpec(
                    output_type="property",
                    file_pattern="predictions/final_ranked_designs/all_designs_metrics.csv",
                    processor="property",
                    name_template="{job_name}_all_metrics",
                    params={"merge_existing": False},
                ),
                # Capture inverse-folded designs (refolded structures)
                IngestionSpec(
                    output_type="structure",
                    file_pattern="predictions/intermediate_designs_inverse_folded/refold_cif/*.cif",
                    processor="structure",
                    name_template="{job_name}_refolded_{stem}",
                    params={"register_entity": True, "intermediate": True},
                    required=False,
                ),
            ],
        )
        self.register_model(boltzgen_card, BoltzGenAdapter(self))

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
            description=("LigandMPNN: design sequences conditioned on a protein PDB"),
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
            execution=ExecutionSpec(
                mode="external_config", entrypoint="python", environment={"gpu": True}
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
                    name="bbox",
                    kind="json",
                    provider="json_payload",
                    optional=True,
                ),
            ],
            output_spec=[],
        )
        self.register_model(pocket2mol_card, Pocket2MolAdapter(self))

        # RFdiffusion2 - all-atom protein generation via diffusion
        rfdiffusion2_card = ModelCard(
            name="rfdiffusion2",
            version="2.0",
            description="All-atom protein backbone generation via diffusion with motif scaffolding",
            execution=ExecutionSpec(
                mode="external_config",
                entrypoint="apptainer exec",
                environment={
                    "container": "singularity",
                    "sif_path": "RFdiffusion2/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif",
                },
            ),
            input_spec=[
                ArtifactSpec(
                    name="structure_pdb",
                    kind="structure",
                    provider="file_path",
                    format="pdb",
                    optional=True,
                ),
            ],
            output_spec=[
                ArtifactSpec(
                    name="designed_structures",
                    kind="structure",
                    provider="rfdiffusion2_adapter",
                    format="pdb",
                ),
                ArtifactSpec(
                    name="sequences",
                    kind="sequence",
                    provider="rfdiffusion2_adapter",
                    format="fasta",
                    optional=True,
                ),
            ],
            ingestion_spec=[
                IngestionSpec(
                    output_type="structure",
                    file_pattern="*.pdb",
                    processor="structure",
                    name_template="{job_name}_{stem}",
                ),
            ],
            metadata={
                "supports_ligands": True,
                "supports_partial_diffusion": True,
                "supports_motif_scaffolding": True,
                "container_type": "apptainer",
            },
        )
        self.register_model(rfdiffusion2_card, RFdiffusion2Adapter(self))

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
                raise ValueError(
                    f"Required artifact '{spec.name}' could not be resolved"
                )
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
        metadata = {
            "entity": entity,
            "graph_metadata": payload.get("graph_metadata", {}),
        }
        return ArtifactBundle(spec=spec, path=path, metadata=metadata)

    def _provide_ligand_file(
        self,
        spec: ArtifactSpec,
        request: ModelRequest,
    ) -> ArtifactBundle:
        value = (
            request.get_input(spec.name)
            or request.get_input("ligand")
            or request.get_input("ligand_file")
        )
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
        value = (
            request.get_input(spec.name)
            or request.get_input("payload")
            or request.get_input("metadata")
        )
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
                return ArtifactBundle(
                    spec=spec, path=path, metadata={"source": "inline"}
                )

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

        name = (
            request.get_input(spec.name)
            or request.get_input("structure")
            or request.get_input("receptor_structure")
        )
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
        return ArtifactBundle(
            spec=spec, path=path, metadata={"entity": name, "format": "pkl"}
        )

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
        dataset_name = request.get_input(spec.name) or request.get_input(
            "embedding_dataset"
        )
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

        emb_proc = (
            EmbeddingProcessor(model_name=initial_model)
            if initial_model
            else EmbeddingProcessor()
        )
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

        Uses the model card's ingestion_spec to determine how to process outputs:
        - Structures (CIF/PDB) → StructureProcessor (register as entity)
        - Properties (CSV) → PropertyProcessor (record in property table)
        - Embeddings (NPZ) → EmbeddingProcessor
        - Ligands (SDF) → MoleculeProcessor
        """
        import logging
        import glob

        logger = logging.getLogger("ModelManager.ingest")
        summary: Dict[str, Any] = {
            "model": invocation.model,
            "job_id": invocation.metadata.get("job_id"),
            "ingested": [],
            "errors": [],
        }

        # Get working directory from job or metadata
        working_dir = None
        if invocation.job:
            working_dir = invocation.job.working_dir
        elif invocation.metadata.get("context"):
            working_dir = Path(invocation.metadata["context"].get("work_dir", ""))

        if not working_dir or not Path(working_dir).exists():
            logger.warning("[ingest] No valid working directory found")
            return summary

        working_dir = Path(working_dir)
        job_name = invocation.metadata.get("job_name", invocation.model)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        # Process ingestion specs from model card
        if invocation.card and invocation.card.ingestion_spec:
            for spec in invocation.card.ingestion_spec:
                try:
                    results = self._process_ingestion_spec(
                        spec, working_dir, job_name, timestamp, invocation.model
                    )
                    summary["ingested"].extend(results)
                except Exception as e:
                    error_msg = f"Failed to process {spec.output_type} spec: {e}"
                    logger.warning("[ingest] %s", error_msg)
                    if spec.required:
                        summary["errors"].append(error_msg)

        # Fallback: process explicit artifacts bundled by adapters
        for bundle in invocation.outputs or []:
            if (
                bundle.spec.kind == "property"
                and bundle.spec.provider == "property_processor"
            ):
                table_path = Path(bundle.path)
                table_name = table_path.stem
                try:
                    prop = PropertyProcessor()
                    df = prop.load_property_table(table_name)
                    prop.save_property_table(table_name)
                    summary["ingested"].append(
                        {
                            "type": "property_table",
                            "name": table_name,
                            "rows": int(len(df)) if df is not None else None,
                        }
                    )
                except Exception:
                    continue
            elif bundle.spec.kind == "ligand":
                try:
                    from protos.processing.molecule import MoleculeProcessor

                    mp = MoleculeProcessor()
                    rel = Path(bundle.path).relative_to(self.paths.data_root)
                    ent_name = Path(bundle.path).stem
                    mp.save_entity(
                        ent_name,
                        {"kind": "sdf", "file_path": str(rel)},
                        metadata={"source_model": invocation.model},
                    )
                    summary["ingested"].append(
                        {
                            "type": "ligand_sdf",
                            "name": ent_name,
                            "path": str(rel),
                        }
                    )
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
                    summary["ingested"].append(
                        {
                            "type": "property_table",
                            "name": table_name,
                            "rows": int(len(df)) if df is not None else None,
                        }
                    )
                except Exception:
                    pass

        logger.info(
            "[ingest] Completed for %s: %d items ingested, %d errors",
            invocation.model,
            len(summary["ingested"]),
            len(summary["errors"]),
        )
        return summary

    def _process_ingestion_spec(
        self,
        spec: IngestionSpec,
        working_dir: Path,
        job_name: str,
        timestamp: str,
        model: str,
    ) -> List[Dict[str, Any]]:
        """Process a single ingestion spec and register outputs.

        Returns list of ingestion result dicts.
        """
        import glob
        import logging

        logger = logging.getLogger("ModelManager.ingest")
        results: List[Dict[str, Any]] = []

        # Find matching files
        pattern = str(working_dir / spec.file_pattern)
        matching_files = glob.glob(pattern, recursive=True)

        if not matching_files:
            logger.debug("[ingest] No files match pattern: %s", spec.file_pattern)
            return results

        logger.info(
            "[ingest] Found %d files for %s pattern",
            len(matching_files),
            spec.output_type,
        )

        for file_path in matching_files:
            file_path = Path(file_path)
            stem = file_path.stem

            # Generate entity name from template
            entity_name = spec.name_template.format(
                job_name=job_name,
                stem=stem,
                model=model,
                timestamp=timestamp,
            )

            try:
                if spec.processor == "structure":
                    result = self._ingest_structure(
                        file_path, entity_name, spec.params, model
                    )
                elif spec.processor == "property":
                    result = self._ingest_property(
                        file_path, entity_name, spec.params, model
                    )
                elif spec.processor == "embedding":
                    result = self._ingest_embedding(
                        file_path, entity_name, spec.params, model
                    )
                elif spec.processor == "molecule":
                    result = self._ingest_molecule(
                        file_path, entity_name, spec.params, model
                    )
                else:
                    logger.warning("[ingest] Unknown processor: %s", spec.processor)
                    continue

                if result:
                    results.append(result)

            except Exception as e:
                logger.warning("[ingest] Failed to ingest %s: %s", file_path, e)

        return results

    def _ingest_structure(
        self,
        file_path: Path,
        entity_name: str,
        params: Dict[str, Any],
        source_model: str,
    ) -> Optional[Dict[str, Any]]:
        """Ingest a structure file (CIF/PDB) into Protos.

        Args:
            file_path: Path to the structure file
            entity_name: Name to register the entity as
            params: Additional parameters (register_entity, copy_to_mmcif, etc.)
            source_model: Name of the model that produced this output

        Returns:
            Ingestion result dict or None
        """
        import logging
        from protos.processing.structure import StructureProcessor
        from protos.io.ingest.structure_loader import StructureLoader

        logger = logging.getLogger("ModelManager.ingest")

        try:
            sp = StructureProcessor()
            loader = StructureLoader(processor=sp)

            # Copy structure to mmcif directory if requested
            if params.get("copy_to_mmcif", True):
                mmcif_dir = Path(self.paths.data_root) / "structure" / "mmcif"
                mmcif_dir.mkdir(parents=True, exist_ok=True)
                dest_path = mmcif_dir / f"{entity_name}.cif"

                # Copy file
                shutil.copy2(file_path, dest_path)
                logger.info("[ingest] Copied structure to %s", dest_path)

                # Also copy to models directory for organized storage
                models_dir = Path(self.paths.data_root) / "models" / source_model / "structures"
                models_dir.mkdir(parents=True, exist_ok=True)
                model_dest = models_dir / f"{entity_name}.cif"
                shutil.copy2(file_path, model_dest)

            # Register as entity if requested
            if params.get("register_entity", True):
                # Register via loader (will parse and index)
                try:
                    # Import the structure file to register it
                    target_path = dest_path if params.get("copy_to_mmcif") else file_path
                    # The structure is already in mmcif dir, just load to verify
                    df = sp.load_entity(entity_name)
                    if df is None:
                        # Try loading directly from file
                        from protos.io.formats.cif_utils import read_cif_file
                        df = read_cif_file(str(target_path))
                        if df is not None and not df.empty:
                            sp.save_entity(entity_name, df)
                    logger.info("[ingest] Registered entity: %s", entity_name)
                except Exception as e:
                    logger.debug("[ingest] Entity registration: %s", e)

            return {
                "type": "structure",
                "name": entity_name,
                "source_file": str(file_path),
                "registered": params.get("register_entity", True),
                "intermediate": params.get("intermediate", False),
            }

        except Exception as e:
            logger.warning("[ingest] Structure ingestion failed: %s", e)
            return None

    def _ingest_property(
        self,
        file_path: Path,
        table_name: str,
        params: Dict[str, Any],
        source_model: str,
    ) -> Optional[Dict[str, Any]]:
        """Ingest a property table (CSV) into Protos.

        Args:
            file_path: Path to the CSV file
            table_name: Name for the property table
            params: Additional parameters (merge_existing, etc.)
            source_model: Name of the model that produced this output

        Returns:
            Ingestion result dict or None
        """
        import logging

        logger = logging.getLogger("ModelManager.ingest")

        try:
            # Load the CSV
            df = pd.read_csv(file_path)

            # Copy to property directory
            property_dir = Path(self.paths.data_root) / "property"
            property_dir.mkdir(parents=True, exist_ok=True)
            dest_path = property_dir / f"{table_name}.csv"

            if params.get("merge_existing", False) and dest_path.exists():
                # Merge with existing table
                existing = pd.read_csv(dest_path)
                df = pd.concat([existing, df], ignore_index=True)
                df = df.drop_duplicates()

            df.to_csv(dest_path, index=False)

            # Also keep a copy in models directory
            models_dir = Path(self.paths.data_root) / "models" / source_model / "properties"
            models_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(dest_path, models_dir / f"{table_name}.csv")

            logger.info("[ingest] Saved property table: %s (%d rows)", table_name, len(df))

            return {
                "type": "property_table",
                "name": table_name,
                "rows": len(df),
                "source_file": str(file_path),
            }

        except Exception as e:
            logger.warning("[ingest] Property ingestion failed: %s", e)
            return None

    def _ingest_embedding(
        self,
        file_path: Path,
        dataset_name: str,
        params: Dict[str, Any],
        source_model: str,
    ) -> Optional[Dict[str, Any]]:
        """Ingest embeddings (NPZ) into Protos.

        Args:
            file_path: Path to the NPZ file
            dataset_name: Name for the embedding dataset
            params: Additional parameters
            source_model: Name of the model that produced this output

        Returns:
            Ingestion result dict or None
        """
        import logging

        logger = logging.getLogger("ModelManager.ingest")

        try:
            # Copy to embeddings directory
            emb_dir = Path(self.paths.data_root) / "embeddings"
            emb_dir.mkdir(parents=True, exist_ok=True)
            dest_path = emb_dir / f"{dataset_name}.npz"

            shutil.copy2(file_path, dest_path)

            # Count entities in NPZ
            npz = np.load(file_path, allow_pickle=False)
            num_entities = len(npz.files)

            logger.info("[ingest] Saved embedding dataset: %s (%d entities)", dataset_name, num_entities)

            return {
                "type": "embedding",
                "name": dataset_name,
                "entities": num_entities,
                "source_file": str(file_path),
            }

        except Exception as e:
            logger.warning("[ingest] Embedding ingestion failed: %s", e)
            return None

    def _ingest_molecule(
        self,
        file_path: Path,
        entity_name: str,
        params: Dict[str, Any],
        source_model: str,
    ) -> Optional[Dict[str, Any]]:
        """Ingest a molecule file (SDF) into Protos.

        Args:
            file_path: Path to the SDF file
            entity_name: Name for the molecule entity
            params: Additional parameters
            source_model: Name of the model that produced this output

        Returns:
            Ingestion result dict or None
        """
        import logging

        logger = logging.getLogger("ModelManager.ingest")

        try:
            from protos.processing.molecule import MoleculeProcessor

            mp = MoleculeProcessor()

            # Copy to molecules directory
            mol_dir = Path(self.paths.data_root) / "molecules"
            mol_dir.mkdir(parents=True, exist_ok=True)
            dest_path = mol_dir / f"{entity_name}.sdf"

            shutil.copy2(file_path, dest_path)

            # Register as entity
            rel_path = dest_path.relative_to(self.paths.data_root)
            mp.save_entity(
                entity_name,
                {"kind": "sdf", "file_path": str(rel_path)},
                metadata={"source_model": source_model},
            )

            logger.info("[ingest] Registered molecule: %s", entity_name)

            return {
                "type": "molecule",
                "name": entity_name,
                "source_file": str(file_path),
            }

        except Exception as e:
            logger.warning("[ingest] Molecule ingestion failed: %s", e)
            return None

    def ingest_job_outputs(
        self,
        job_id: str,
        *,
        force: bool = False,
        cleanup: bool = True,
    ) -> Dict[str, Any]:
        """Ingest outputs from a completed job by job ID.

        This is the main entry point for ingesting outputs from Docker/external jobs.
        It packages the outputs and registers them in Protos, then optionally
        cleans up the job directory.

        Args:
            job_id: The job ID to ingest outputs from
            force: Force re-ingestion even if already processed
            cleanup: If True (default), clean up job directory after successful
                    ingestion (unless job was marked as persistent)

        Returns:
            Ingestion summary dict
        """
        import logging

        logger = logging.getLogger("ModelManager.ingest")

        # Get job state
        state = self._executor.status(job_id)

        if state.status != JobStatus.COMPLETED:
            raise ValueError(f"Job {job_id} is not completed (status: {state.status})")

        # Check if already ingested
        if state.metadata.get("ingested") and not force:
            logger.info("[ingest] Job %s already ingested", job_id)
            return state.metadata.get("ingestion_summary", {})

        # Get model card
        card = self.cards.get(state.model)
        if not card:
            raise ValueError(f"No model card found for '{state.model}'")

        # Create a synthetic invocation for ingestion
        invocation = ModelInvocation(
            model=state.model,
            card=card,
            job=state.prepared_job,
            metadata={
                "job_id": job_id,
                "job_name": state.metadata.get("job_name", state.model),
            },
        )

        # Run ingestion
        summary = self.ingest_outputs(invocation)
        summary["job_id"] = job_id

        # Update job state to mark as ingested
        state.metadata["ingested"] = True
        state.metadata["ingestion_summary"] = summary
        state.metadata["ingested_at"] = datetime.now().isoformat()
        self._executor._save_jobs()

        logger.info("[ingest] Job %s ingestion complete: %d items", job_id, len(summary.get("ingested", [])))

        # Cleanup job directory if requested (and not persistent)
        if cleanup and not state.metadata.get("persistent", False):
            if self._executor.cleanup_job(job_id):
                summary["cleaned_up"] = True
                logger.info("[ingest] Cleaned up job directory for %s", job_id)
            else:
                summary["cleaned_up"] = False

        return summary

    def cleanup_job(self, job_id: str, force: bool = False) -> bool:
        """Clean up a job's working directory.

        Args:
            job_id: The job ID to clean up
            force: If True, clean up even if persistent=True or not ingested

        Returns:
            True if cleanup was performed
        """
        return self._executor.cleanup_job(job_id, force=force)

    def cleanup_completed_jobs(self, model: Optional[str] = None) -> int:
        """Clean up all completed and ingested job directories.

        Args:
            model: Optional model name to filter by

        Returns:
            Number of jobs cleaned up
        """
        return self._executor.cleanup_completed_jobs(model=model)


def prepare_mutation_screen(
    *,
    seq_proc: Optional[SequenceProcessor] = None,
    dataset_name: str,
    grn_positions: Iterable[str],
    mutations: Iterable[str],
    grn_table_name: Optional[str] = None,
    protein_family: str = "gpcr_a",
    reference_table: str = "gpcrdb_ref",
    manager: Optional[ModelManager] = None,
    model_name: str = "boltz2",
    base_config: Optional[MutableMapping[str, Any]] = None,
    metadata: Optional[MutableMapping[str, Any]] = None,
) -> List[ModelInvocation]:
    """Prepare Boltz-style mutation predictions across GRN positions."""

    if not dataset_name:
        raise ValueError("dataset_name is required for mutation screening")

    seq_proc = seq_proc or SequenceProcessor()
    manager = manager or ModelManager(data_root=Path(seq_proc.paths.data_root))

    sequences = seq_proc.load_dataset(dataset_name)
    if not sequences:
        raise ValueError(f"Sequence dataset '{dataset_name}' is empty or missing")

    if grn_table_name:
        grn_proc = GRNProcessor()
        grn_table = grn_proc.load_table(grn_table_name)
    else:
        grn_table = seq_proc.annotate_with_grn(
            dataset_name=dataset_name,
            reference_table=reference_table,
            protein_family=protein_family,
            output_table=None,
            allow_create=True,
            return_summary=False,
        )

    if grn_table is None or grn_table.empty:
        raise ValueError("GRN annotations are required to prepare mutations")

    invocations: List[ModelInvocation] = []
    normalized_labels = [str(label) for label in grn_positions]

    for entity_name, sequence in sequences.items():
        mapping = GRNProcessor.build_grn_to_seq_index(
            grn_table, sequence_id=entity_name
        )
        if not mapping:
            continue

        for grn_label in normalized_labels:
            seq_idx = mapping.get(grn_label)
            if not seq_idx:
                continue

            position = int(seq_idx)
            if position < 1 or position > len(sequence):
                continue
            original = sequence[position - 1]

            for mutant in mutations:
                mutant_residue = str(mutant).strip()
                if len(mutant_residue) != 1 or mutant_residue == original:
                    continue

                mutation_payload = {
                    "position": position,
                    "original": original,
                    "mutant": mutant_residue,
                    "name": f"{original}{position}{mutant_residue}",
                    "grn": grn_label,
                }

                entry_config = deepcopy(base_config) if base_config else {}
                entry_config["mutations"] = [mutation_payload]
                entry_config.setdefault(
                    "output_name", f"{entity_name}_{mutation_payload['name']}"
                )

                entry_metadata = dict(metadata or {})
                entry_metadata.update(
                    {
                        "entity": entity_name,
                        "grn_label": grn_label,
                        "mutations": [mutation_payload],
                        "dataset": dataset_name,
                    }
                )

                invocation = manager.prepare_input(
                    model_name,
                    entity_name=entity_name,
                    entity_format="sequence",
                    dataset_name=dataset_name,
                    config=entry_config,
                    metadata=entry_metadata,
                )
                invocations.append(invocation)

    return invocations


__all__ = [
    # Core manager
    "ModelManager",
    # Adapters
    "ModelAdapterBase",
    "ExternalJobAdapter",
    "RuntimeAdapter",
    "ConfigurableRuntimeAdapter",
    "ConfigurableExternalAdapter",
    "ModelRequest",
    "BoltzAdapter",
    "LambdaAdapter",
    # Job execution
    "JobExecutor",
    "DockerJobExecutor",
    # Utilities
    "prepare_mutation_screen",
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

            base = (
                Path(__file__).resolve().parent / "Uni-Dock" / "unidock_tools" / "src"
            )
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
                    with (
                        open(pdb_src, "r", encoding="utf-8", errors="ignore") as fh,
                        open(sanitized, "w", encoding="utf-8") as out,
                    ):
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
                    pdb_for_center = (
                        pdb_src.with_suffix(".pdb")
                        if pdb_src.with_suffix(".pdb").exists()
                        else pdb_src
                    )
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
                (
                    cfg.get("size")
                    if not isinstance(cfg.get("size"), (list, tuple))
                    else None
                )
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
                    str(receptor_pdbqt),
                    str(map_dir),
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
            "--gpu_batch",
            str(lig_target),
            "--dir",
            str(ctx.outputs_dir),
            "--center_x",
            str(center[0]),
            "--center_y",
            str(center[1]),
            "--center_z",
            str(center[2]),
            "--size_x",
            str(sx),
            "--size_y",
            str(sy),
            "--size_z",
            str(sz),
            "--scoring",
            scoring,
            "--num_modes",
            str(num_modes),
            "--energy_range",
            str(energy_range),
            "--refine_step",
            str(refine_step),
            "--seed",
            str(seed),
            "--verbosity",
            "2",
            "--keep_nonpolar_H",
        ]

        search_mode = cfg.get("search_mode")
        if search_mode:
            command += ["--search_mode", str(search_mode)]
        else:
            exhaustiveness = int(cfg.get("exhaustiveness", 256))
            max_step = int(cfg.get("max_step", 10))
            command += [
                "--exhaustiveness",
                str(exhaustiveness),
                "--max_step",
                str(max_step),
            ]

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
            run_id=ctx.run_id,
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
            run_id=ctx.run_id,
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
            run_id=ctx.run_id,
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
            run_id=ctx.run_id,
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
        ligand = next(
            (b for b in inputs if b.spec.name in ("ligand_molecule", "ligand_file")),
            None,
        )

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
        sample_script = (
            Path(__file__).resolve().parent / "Pocket2Mol" / "sample_for_pdb.py"
        )
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
        config_path = (
            Path(__file__).resolve().parent
            / "Pocket2Mol"
            / "configs"
            / "sample_for_pdb.yml"
        )
        if config_path.exists():
            command.extend(["--config", str(config_path)])
            ctx.config_path = config_path

        job = PreparedJob(
            run_id=ctx.run_id,
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
        ligand_bundle = next(
            (b for b in inputs if b.spec.name in ("ligand_molecule", "ligand_file")),
            None,
        )
        ctx = ModelRunContext(self.manager.paths, card)
        ctx.create()

        # Ensure PDB input exists; convert CIF→PDB if necessary
        pdb_src = self._ensure_pdb_file(
            Path(pdb_bundle.path), ctx, purpose="ligand_mpnn"
        )

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
                with (
                    open(pdb_src, "r", encoding="utf-8", errors="ignore") as ph,
                    open(combined, "w", encoding="utf-8") as out,
                ):
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
        for key, flag in (
            ("fixed_residues", "--fixed_residues"),
            ("redesigned_residues", "--redesigned_residues"),
        ):
            vals = cfg.get(key)
            if isinstance(vals, list):
                cmd.extend([flag, " ".join(str(v) for v in vals)])
            elif isinstance(vals, str):
                cmd.extend([flag, vals])

        add_flag("bias_AA", "--bias_AA")
        add_flag("omit_AA", "--omit_AA")

        # Build job
        job = PreparedJob(
            run_id=ctx.run_id,
            command=cmd,
            working_dir=ctx.work_dir,
            artifacts=[pdb_bundle],
            metadata={
                "context": ctx.as_metadata(),
                "args": cmd[2:],
            },
        )
        return job


class RFdiffusion2Adapter(ExternalJobAdapter):
    """Adapter for RFdiffusion2 protein backbone generation.

    RFdiffusion2 generates protein backbones via diffusion with support for:
    - Unconditional generation (de novo design)
    - Partial diffusion (inpainting)
    - Motif scaffolding with atomic constraints
    - Ligand-aware design
    - Multi-chain complexes

    Uses Apptainer/Singularity container for execution.
    """

    MODEL_DIR = "rfdiffusion2"
    SIF_RELATIVE = "RFdiffusion2/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif"
    PIPELINE_SCRIPT = "rf_diffusion/benchmark/pipeline.py"

    def build_job(
        self,
        card: ModelCard,
        request: ModelRequest,
        inputs: List[ArtifactBundle],
    ) -> PreparedJob:
        """Build RFdiffusion2 job from request config.

        Config keys:
            input_pdb: Path to input PDB/CIF structure
            contigs: Contig specification (e.g., "A1-100/0 B1-50")
            num_designs: Number of designs to generate
            ligands: Comma-separated ligand names (e.g., "RET,ATP")
            contig_atoms: Dict mapping residue IDs to atoms (e.g., {"A230": ["OG"]})
            guidepost: Use unindexed scaffolding (default True)
            partial_T: Partial diffusion timesteps (for inpainting)
            stop_step: Pipeline stop point ("design", "ligandmpnn", "end")
        """
        config = request.config or {}
        job_name = config.get("job_name", "rfd2_job")

        # Create run context with job_name as prefix
        ctx = ModelRunContext(self.paths, card, run_prefix=job_name)
        ctx.create()

        # Copy/convert input structure
        input_pdb = config.get("input_pdb")
        if not input_pdb:
            # Try to get from input artifacts
            for bundle in inputs:
                if bundle.spec.kind == "structure":
                    input_pdb = str(bundle.path)
                    break
        if not input_pdb:
            raise ValueError("RFdiffusion2 requires 'input_pdb' in config or structure input")

        # Ensure PDB exists (convert CIF if needed)
        pdb_path = self._ensure_pdb_file(Path(input_pdb), ctx, purpose="rfdiffusion2")

        # Build inference arguments
        args = self._build_inference_args(config, pdb_path, ctx.outputs_dir)

        # Build command for Apptainer execution
        # The pipeline.py script is run inside the container
        command = [
            "python",
            self.PIPELINE_SCRIPT,
        ] + args

        # Save job metadata
        metadata_path = ctx.work_dir / "metadata.json"
        metadata_content = {
            "job_name": job_name,
            "input_pdb": str(pdb_path),
            "contigs": config.get("contigs", ""),
            "num_designs": config.get("num_designs", 10),
            "config": {k: v for k, v in config.items() if k != "input_pdb"},
        }
        with open(metadata_path, "w", encoding="utf-8") as fh:
            json.dump(metadata_content, fh, indent=2)

        # Determine SIF path and RFdiffusion2 install directory
        rfd2_install_dir = Path(__file__).parent / "RFdiffusion2"
        sif_path = rfd2_install_dir / "rf_diffusion" / "exec" / "bakerlab_rf_diffusion_aa.sif"

        job = PreparedJob(
            run_id=ctx.run_id,
            command=command,
            working_dir=ctx.work_dir,
            artifacts=[
                ArtifactBundle(
                    spec=ArtifactSpec(
                        name="rfd2_input",
                        kind="structure",
                        provider="rfdiffusion2_adapter",
                        format="pdb",
                    ),
                    path=pdb_path,
                    metadata={"job_name": job_name},
                ),
            ],
            metadata={
                "job_name": job_name,
                "output_dir": str(ctx.outputs_dir),
                "sif_path": str(sif_path),
                "container_workdir": str(rfd2_install_dir),  # Run from RFdiffusion2 install dir
                "bind_mounts": [
                    f"{pdb_path.parent}:{pdb_path.parent}",
                    f"{ctx.outputs_dir}:{ctx.outputs_dir}",
                    f"{rfd2_install_dir}:{rfd2_install_dir}",  # Bind RFdiffusion2 install
                ],
                "env": {
                    "PYTHONPATH": str(rfd2_install_dir),  # Add RFdiffusion2 to Python path
                },
                "executor": "apptainer",
            },
        )
        return job

    def _build_inference_args(
        self,
        config: MutableMapping[str, Any],
        pdb_path: Path,
        output_dir: Path,
    ) -> List[str]:
        """Build RFdiffusion2 inference arguments."""
        args = []

        # Input PDB
        args.append(f"inference.input_pdb={pdb_path}")

        # Contigs specification
        contigs = config.get("contigs", "")
        if contigs:
            args.append(f"contigmap.contigs=['{contigs}']")

        # Number of designs
        num_designs = config.get("num_designs", 10)
        args.append(f"inference.num_designs={num_designs}")

        # Ligands
        ligands = config.get("ligands")
        if ligands:
            args.append(f"inference.ligand='{ligands}'")

        # Guidepost mode (unindexed scaffolding)
        guidepost = config.get("guidepost", True)
        args.append(f"inference.contig_as_guidepost={guidepost}")

        # Contig atoms (atomic constraints for motif scaffolding)
        contig_atoms = config.get("contig_atoms")
        if contig_atoms:
            atoms_str = json.dumps(contig_atoms).replace('"', "'")
            args.append(f'contigmap.contig_atoms="{atoms_str}"')

        # Partial diffusion (for inpainting)
        partial_T = config.get("partial_T")
        if partial_T is not None:
            args.append(f"diffuser.partial_T={partial_T}")

        # Stop step (design, ligandmpnn, or end)
        stop_step = config.get("stop_step", "end")
        args.append(f"stop_step='{stop_step}'")

        # Output directory
        args.append(f"outdir={output_dir}")

        # Additional inference parameters
        extra_params = config.get("extra_params", {})
        for key, value in extra_params.items():
            if isinstance(value, str):
                args.append(f"{key}='{value}'")
            else:
                args.append(f"{key}={value}")

        return args
