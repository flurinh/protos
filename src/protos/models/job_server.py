"""FastAPI Job Server for Protos model execution.

This server provides a REST API for submitting and managing model jobs.
It can run locally in Docker (for development) or be deployed remotely.

Endpoints:
    POST   /jobs           - Submit a new job
    GET    /jobs           - List all jobs
    GET    /jobs/{id}      - Get job status
    GET    /jobs/{id}/result - Get job result with output files
    DELETE /jobs/{id}      - Cancel a job
    POST   /jobs/{id}/ingest - Trigger output ingestion
"""

from __future__ import annotations

import asyncio
import io
import json
import os
import shutil
import subprocess
import tarfile
import tempfile
import uuid
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, File, Form, HTTPException, UploadFile, BackgroundTasks
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from pydantic import BaseModel, Field

app = FastAPI(
    title="Protos Job Server",
    description="REST API for submitting and managing model execution jobs",
    version="1.0.0",
)


# =============================================================================
# Models
# =============================================================================


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobSubmission(BaseModel):
    """Request body for job submission."""
    model: str = Field(..., description="Model name (e.g., 'boltzgen', 'boltz')")
    command: List[str] = Field(..., description="Command to execute")
    config: Optional[Dict[str, Any]] = Field(default=None, description="Job configuration")
    persistent: bool = Field(default=False, description="Keep job files after completion")


class JobState(BaseModel):
    """Job state response."""
    job_id: str
    model: str
    status: JobStatus
    created_at: str
    submitted_at: Optional[str] = None
    completed_at: Optional[str] = None
    error: Optional[str] = None
    working_dir: Optional[str] = None
    metadata: Dict[str, Any] = Field(default_factory=dict)


class JobResult(BaseModel):
    """Job result response."""
    job_id: str
    exit_code: int
    stdout: str
    stderr: str
    duration_seconds: float
    output_files: List[str] = Field(default_factory=list)


class JobListResponse(BaseModel):
    """Response for job listing."""
    jobs: List[JobState]
    total: int


# =============================================================================
# Job Storage
# =============================================================================


class JobStore:
    """In-memory job storage with disk persistence."""

    def __init__(self, data_dir: Optional[Path] = None):
        self.data_dir = data_dir or Path("/data/models")
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self._jobs: Dict[str, Dict[str, Any]] = {}
        self._load_jobs()

    def _state_file(self) -> Path:
        return self.data_dir / "job_server_state.json"

    def _load_jobs(self) -> None:
        state_file = self._state_file()
        if state_file.exists():
            try:
                with open(state_file) as f:
                    data = json.load(f)
                self._jobs = data.get("jobs", {})
            except Exception:
                self._jobs = {}

    def _save_jobs(self) -> None:
        state_file = self._state_file()
        with open(state_file, "w") as f:
            json.dump({"jobs": self._jobs}, f, indent=2, default=str)

    def create_job(
        self,
        model: str,
        command: List[str],
        config: Optional[Dict[str, Any]] = None,
        persistent: bool = False,
    ) -> str:
        job_id = f"{model}_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"

        # Create working directory under data/models/<model>/jobs/<job_id>/
        working_dir = self.data_dir / model / "jobs" / job_id
        working_dir.mkdir(parents=True, exist_ok=True)

        self._jobs[job_id] = {
            "job_id": job_id,
            "model": model,
            "status": JobStatus.PENDING.value,
            "command": command,
            "config": config or {},
            "persistent": persistent,
            "created_at": datetime.now().isoformat(),
            "submitted_at": None,
            "completed_at": None,
            "error": None,
            "working_dir": str(working_dir),
            "exit_code": None,
            "stdout": "",
            "stderr": "",
            "duration_seconds": 0.0,
            "output_files": [],
        }
        self._save_jobs()
        return job_id

    def get_job(self, job_id: str) -> Optional[Dict[str, Any]]:
        return self._jobs.get(job_id)

    def update_job(self, job_id: str, **updates) -> bool:
        if job_id not in self._jobs:
            return False
        self._jobs[job_id].update(updates)
        self._save_jobs()
        return True

    def list_jobs(
        self,
        model: Optional[str] = None,
        status: Optional[JobStatus] = None,
    ) -> List[Dict[str, Any]]:
        jobs = list(self._jobs.values())
        if model:
            jobs = [j for j in jobs if j["model"] == model]
        if status:
            jobs = [j for j in jobs if j["status"] == status.value]
        return sorted(jobs, key=lambda j: j["created_at"], reverse=True)

    def delete_job(self, job_id: str) -> bool:
        if job_id not in self._jobs:
            return False

        job = self._jobs[job_id]
        working_dir = Path(job["working_dir"])

        # Remove working directory if not persistent
        if not job.get("persistent") and working_dir.exists():
            try:
                shutil.rmtree(working_dir)
            except Exception:
                pass

        del self._jobs[job_id]
        self._save_jobs()
        return True


# Global job store - initialized on startup
job_store: Optional[JobStore] = None


@app.on_event("startup")
async def startup():
    global job_store
    data_dir = Path(os.environ.get("PROTOS_DATA_DIR", "/data/models"))
    job_store = JobStore(data_dir)


# =============================================================================
# Job Execution
# =============================================================================


# Docker images per model
DOCKER_IMAGES = {
    "boltzgen": "protos/boltzgen:latest",
    "boltz": "protos/boltz:latest",
    "lambda": "protos/lambda:latest",
}


async def run_job_in_docker(job_id: str) -> None:
    """Execute a job in a Docker container."""
    job = job_store.get_job(job_id)
    if not job:
        return

    job_store.update_job(
        job_id,
        status=JobStatus.RUNNING.value,
        submitted_at=datetime.now().isoformat(),
    )

    model = job["model"]
    working_dir = Path(job["working_dir"])
    command = job["command"]

    # Get Docker image
    docker_image = job.get("config", {}).get("docker_image") or DOCKER_IMAGES.get(model, "protos/base:latest")

    # Build Docker command
    docker_cmd = [
        "docker", "run", "--rm",
        "--gpus", "all",
        "--shm-size", "8g",
        "-v", f"{working_dir}:/workspace",
        "-v", f"{Path.home()}/.cache/huggingface:/cache",
        "-e", "HF_HOME=/cache",
        "-w", "/workspace",
        "--user", f"{os.getuid()}:{os.getgid()}",
        docker_image,
    ] + command

    stdout_file = working_dir / "stdout.log"
    stderr_file = working_dir / "stderr.log"

    start_time = datetime.now()
    try:
        with open(stdout_file, "w") as out, open(stderr_file, "w") as err:
            proc = await asyncio.create_subprocess_exec(
                *docker_cmd,
                stdout=out,
                stderr=err,
                cwd=str(working_dir),
            )
            await proc.wait()

        duration = (datetime.now() - start_time).total_seconds()

        # Collect output files
        output_dir = working_dir / "predictions"
        output_files = []
        if output_dir.exists():
            output_files = [
                str(f.relative_to(working_dir))
                for f in output_dir.rglob("*")
                if f.is_file()
            ]

        job_store.update_job(
            job_id,
            status=JobStatus.COMPLETED.value if proc.returncode == 0 else JobStatus.FAILED.value,
            completed_at=datetime.now().isoformat(),
            exit_code=proc.returncode,
            stdout=stdout_file.read_text() if stdout_file.exists() else "",
            stderr=stderr_file.read_text() if stderr_file.exists() else "",
            duration_seconds=duration,
            output_files=output_files,
            error=f"Exit code: {proc.returncode}" if proc.returncode != 0 else None,
        )

    except Exception as e:
        job_store.update_job(
            job_id,
            status=JobStatus.FAILED.value,
            completed_at=datetime.now().isoformat(),
            error=str(e),
        )


# =============================================================================
# API Endpoints
# =============================================================================


@app.post("/jobs", response_model=JobState)
async def submit_job(
    submission: JobSubmission,
    background_tasks: BackgroundTasks,
):
    """Submit a new job for execution."""
    job_id = job_store.create_job(
        model=submission.model,
        command=submission.command,
        config=submission.config,
        persistent=submission.persistent,
    )

    # Run job in background
    background_tasks.add_task(run_job_in_docker, job_id)

    job = job_store.get_job(job_id)
    return JobState(
        job_id=job["job_id"],
        model=job["model"],
        status=JobStatus(job["status"]),
        created_at=job["created_at"],
        working_dir=job["working_dir"],
        metadata={"persistent": job["persistent"]},
    )


@app.post("/jobs/upload", response_model=JobState)
async def submit_job_with_files(
    model: str = Form(...),
    command: str = Form(...),  # JSON-encoded list
    persistent: bool = Form(False),
    files: List[UploadFile] = File(default=[]),
    background_tasks: BackgroundTasks = None,
):
    """Submit a job with file uploads."""
    command_list = json.loads(command)

    job_id = job_store.create_job(
        model=model,
        command=command_list,
        persistent=persistent,
    )

    job = job_store.get_job(job_id)
    working_dir = Path(job["working_dir"])

    # Save uploaded files
    for file in files:
        file_path = working_dir / file.filename
        file_path.parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, "wb") as f:
            content = await file.read()
            f.write(content)

    # Run job in background
    background_tasks.add_task(run_job_in_docker, job_id)

    return JobState(
        job_id=job["job_id"],
        model=job["model"],
        status=JobStatus(job["status"]),
        created_at=job["created_at"],
        working_dir=job["working_dir"],
        metadata={"persistent": job["persistent"], "uploaded_files": [f.filename for f in files]},
    )


@app.post("/jobs/run/{run_id}", response_model=JobState)
async def submit_run_package(
    run_id: str,
    model: str = Form(...),
    run_package: UploadFile = File(...),
    background_tasks: BackgroundTasks = None,
):
    """Receive and execute a packaged run directory.

    The run_id becomes the job_id - the canonical identifier for this job.
    The uploaded tar.gz contains inputs/, config.yaml, and job.json.
    """
    # Create working directory using run_id
    run_dir = job_store.data_dir / model / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    # Unpack tar
    content = await run_package.read()
    tar_buffer = io.BytesIO(content)
    with tarfile.open(fileobj=tar_buffer, mode='r:gz') as tar:
        tar.extractall(run_dir)

    # Create outputs directory
    (run_dir / "outputs").mkdir(exist_ok=True)

    # Read job.json to get command and other settings
    job_json_path = run_dir / "job.json"
    if job_json_path.exists():
        with open(job_json_path) as f:
            job_data = json.load(f)
        command = job_data.get("command", [])
        persistent = job_data.get("metadata", {}).get("persistent", False)
    else:
        command = []
        persistent = False

    # Register job in job store with the run_id as job_id
    job_store._jobs[run_id] = {
        "job_id": run_id,
        "model": model,
        "status": JobStatus.PENDING.value,
        "command": command,
        "config": {},
        "persistent": persistent,
        "created_at": datetime.now().isoformat(),
        "submitted_at": None,
        "completed_at": None,
        "error": None,
        "working_dir": str(run_dir),
        "exit_code": None,
        "stdout": "",
        "stderr": "",
        "duration_seconds": 0.0,
        "output_files": [],
        "run_id": run_id,
    }
    job_store._save_jobs()

    # Run job in background
    background_tasks.add_task(run_job_in_docker, run_id)

    return JobState(
        job_id=run_id,
        model=model,
        status=JobStatus.PENDING,
        created_at=job_store._jobs[run_id]["created_at"],
        working_dir=str(run_dir),
        metadata={"persistent": persistent, "run_id": run_id},
    )


@app.get("/jobs/{job_id}/outputs")
async def get_outputs(job_id: str):
    """Return outputs directory as tar archive.

    Used to download results back to the local run directory.
    """
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    working_dir = Path(job["working_dir"])
    outputs_dir = working_dir / "outputs"

    if not outputs_dir.exists():
        raise HTTPException(status_code=404, detail=f"Outputs directory not found for job {job_id}")

    # Create tar of outputs
    tar_buffer = io.BytesIO()
    with tarfile.open(fileobj=tar_buffer, mode='w:gz') as tar:
        for item in outputs_dir.rglob("*"):
            if item.is_file():
                tar.add(item, arcname=str(item.relative_to(outputs_dir)))
    tar_buffer.seek(0)

    return StreamingResponse(
        tar_buffer,
        media_type="application/gzip",
        headers={"Content-Disposition": f"attachment; filename={job_id}_outputs.tar.gz"}
    )


@app.get("/jobs", response_model=JobListResponse)
async def list_jobs(
    model: Optional[str] = None,
    status: Optional[JobStatus] = None,
):
    """List all jobs, optionally filtered by model or status."""
    jobs = job_store.list_jobs(model=model, status=status)
    return JobListResponse(
        jobs=[
            JobState(
                job_id=j["job_id"],
                model=j["model"],
                status=JobStatus(j["status"]),
                created_at=j["created_at"],
                submitted_at=j["submitted_at"],
                completed_at=j["completed_at"],
                error=j["error"],
                working_dir=j["working_dir"],
                metadata={"persistent": j.get("persistent", False)},
            )
            for j in jobs
        ],
        total=len(jobs),
    )


@app.get("/jobs/{job_id}", response_model=JobState)
async def get_job(job_id: str):
    """Get the status of a specific job."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    return JobState(
        job_id=job["job_id"],
        model=job["model"],
        status=JobStatus(job["status"]),
        created_at=job["created_at"],
        submitted_at=job["submitted_at"],
        completed_at=job["completed_at"],
        error=job["error"],
        working_dir=job["working_dir"],
        metadata={"persistent": job.get("persistent", False)},
    )


@app.get("/jobs/{job_id}/result", response_model=JobResult)
async def get_job_result(job_id: str):
    """Get the result of a completed job."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job["status"] not in (JobStatus.COMPLETED.value, JobStatus.FAILED.value):
        raise HTTPException(status_code=400, detail=f"Job {job_id} is not complete")

    return JobResult(
        job_id=job["job_id"],
        exit_code=job.get("exit_code", -1),
        stdout=job.get("stdout", ""),
        stderr=job.get("stderr", ""),
        duration_seconds=job.get("duration_seconds", 0.0),
        output_files=job.get("output_files", []),
    )


@app.get("/jobs/{job_id}/files/{file_path:path}")
async def download_file(job_id: str, file_path: str):
    """Download an output file from a job."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    working_dir = Path(job["working_dir"])
    full_path = working_dir / file_path

    if not full_path.exists():
        raise HTTPException(status_code=404, detail=f"File {file_path} not found")

    if not full_path.is_file():
        raise HTTPException(status_code=400, detail=f"{file_path} is not a file")

    return FileResponse(full_path)


@app.delete("/jobs/{job_id}")
async def cancel_job(job_id: str):
    """Cancel a running job or delete a completed job."""
    job = job_store.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job {job_id} not found")

    if job["status"] == JobStatus.RUNNING.value:
        # TODO: Actually stop the Docker container
        job_store.update_job(
            job_id,
            status=JobStatus.CANCELLED.value,
            completed_at=datetime.now().isoformat(),
        )
        return {"message": f"Job {job_id} cancelled"}
    else:
        job_store.delete_job(job_id)
        return {"message": f"Job {job_id} deleted"}


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}


# =============================================================================
# Main
# =============================================================================


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
