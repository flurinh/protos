"""Job Client for communicating with Protos Job Server.

This client provides a unified interface for job submission, whether
running locally (via Docker) or remotely (via API).

Configuration is read from data/models/server_config.yaml:

```yaml
mode: local  # or 'remote'

local:
  # Local mode runs a FastAPI server in Docker
  docker_socket: /var/run/docker.sock
  server_image: protos/job-server:latest
  auto_start: true
  port: 8765

remote:
  # Remote mode connects to an external API
  url: https://api.protos.example.com
  api_key: ${PROTOS_API_KEY}
  timeout: 3600
```
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import requests
import yaml


class JobStatus(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class JobState:
    """Represents the state of a job."""
    job_id: str
    model: str
    status: JobStatus
    created_at: str
    submitted_at: Optional[str] = None
    completed_at: Optional[str] = None
    error: Optional[str] = None
    working_dir: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "JobState":
        return cls(
            job_id=data["job_id"],
            model=data["model"],
            status=JobStatus(data["status"]),
            created_at=data["created_at"],
            submitted_at=data.get("submitted_at"),
            completed_at=data.get("completed_at"),
            error=data.get("error"),
            working_dir=data.get("working_dir"),
            metadata=data.get("metadata", {}),
        )


@dataclass
class JobResult:
    """Result from a completed job."""
    job_id: str
    exit_code: int
    stdout: str
    stderr: str
    duration_seconds: float
    output_files: List[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "JobResult":
        return cls(
            job_id=data["job_id"],
            exit_code=data["exit_code"],
            stdout=data["stdout"],
            stderr=data["stderr"],
            duration_seconds=data["duration_seconds"],
            output_files=data.get("output_files", []),
        )

    @property
    def success(self) -> bool:
        return self.exit_code == 0


@dataclass
class ServerConfig:
    """Configuration for the job server."""
    mode: str = "local"  # 'local' or 'remote'

    # Local mode settings
    local_port: int = 8765
    local_host: str = "localhost"
    auto_start: bool = True
    server_image: str = "protos/job-server:latest"

    # Remote mode settings
    remote_url: str = ""
    api_key: str = ""
    timeout: int = 3600

    @classmethod
    def from_file(cls, path: Path) -> "ServerConfig":
        """Load config from YAML file."""
        if not path.exists():
            return cls()

        with open(path) as f:
            data = yaml.safe_load(f) or {}

        config = cls(mode=data.get("mode", "local"))

        local = data.get("local", {})
        config.local_port = local.get("port", 8765)
        config.local_host = local.get("host", "localhost")
        config.auto_start = local.get("auto_start", True)
        config.server_image = local.get("server_image", "protos/job-server:latest")

        remote = data.get("remote", {})
        config.remote_url = remote.get("url", "")
        config.api_key = os.environ.get("PROTOS_API_KEY", remote.get("api_key", ""))
        config.timeout = remote.get("timeout", 3600)

        return config

    def to_file(self, path: Path) -> None:
        """Save config to YAML file."""
        data = {
            "mode": self.mode,
            "local": {
                "port": self.local_port,
                "host": self.local_host,
                "auto_start": self.auto_start,
                "server_image": self.server_image,
            },
            "remote": {
                "url": self.remote_url,
                "api_key": "${PROTOS_API_KEY}",
                "timeout": self.timeout,
            },
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            yaml.safe_dump(data, f, default_flow_style=False)


class JobClient:
    """Client for interacting with Protos Job Server.

    Supports both local (Docker) and remote (API) modes.
    """

    def __init__(
        self,
        data_root: Optional[Path] = None,
        config: Optional[ServerConfig] = None,
    ):
        self.data_root = data_root or Path.home() / ".protos"
        self.models_dir = self.data_root / "models"
        self.models_dir.mkdir(parents=True, exist_ok=True)

        # Load or use provided config
        config_path = self.models_dir / "server_config.yaml"
        self.config = config or ServerConfig.from_file(config_path)

        # Save default config if it doesn't exist
        if not config_path.exists():
            self.config.to_file(config_path)

        self._session = requests.Session()
        self._local_server_running = False

    @property
    def base_url(self) -> str:
        """Get the base URL for the job server."""
        if self.config.mode == "remote":
            return self.config.remote_url.rstrip("/")
        return f"http://{self.config.local_host}:{self.config.local_port}"

    def _ensure_local_server(self) -> None:
        """Ensure the local job server is running."""
        if self.config.mode != "local":
            return

        if self._local_server_running:
            return

        # Check if server is already running
        try:
            resp = self._session.get(f"{self.base_url}/health", timeout=2)
            if resp.status_code == 200:
                self._local_server_running = True
                return
        except Exception:
            pass

        if not self.config.auto_start:
            raise RuntimeError(
                f"Local job server not running at {self.base_url}. "
                "Set auto_start=true or start it manually."
            )

        # Start the server
        self._start_local_server()

    def _start_local_server(self) -> None:
        """Start the local job server in Docker."""
        import subprocess

        # Check if container already exists
        result = subprocess.run(
            ["docker", "ps", "-a", "--filter", "name=protos-job-server", "--format", "{{.Names}}"],
            capture_output=True,
            text=True,
        )

        if "protos-job-server" in result.stdout:
            # Container exists, try to start it
            subprocess.run(["docker", "start", "protos-job-server"], capture_output=True)
        else:
            # Create and start new container
            cmd = [
                "docker", "run", "-d",
                "--name", "protos-job-server",
                "--restart", "unless-stopped",
                "-p", f"{self.config.local_port}:8000",
                "-v", f"{self.models_dir}:/data/models",
                "-v", "/var/run/docker.sock:/var/run/docker.sock",
                "-v", f"{Path.home()}/.cache/huggingface:/cache",
                "-e", "PROTOS_DATA_DIR=/data/models",
                self.config.server_image,
            ]
            subprocess.run(cmd, capture_output=True)

        # Wait for server to be ready
        for _ in range(30):
            try:
                resp = self._session.get(f"{self.base_url}/health", timeout=2)
                if resp.status_code == 200:
                    self._local_server_running = True
                    return
            except Exception:
                pass
            time.sleep(1)

        raise RuntimeError("Failed to start local job server")

    def _headers(self) -> Dict[str, str]:
        """Get request headers."""
        headers = {"Content-Type": "application/json"}
        if self.config.mode == "remote" and self.config.api_key:
            headers["Authorization"] = f"Bearer {self.config.api_key}"
        return headers

    def submit(
        self,
        model: str,
        command: List[str],
        config: Optional[Dict[str, Any]] = None,
        persistent: bool = False,
        files: Optional[Dict[str, Path]] = None,
    ) -> JobState:
        """Submit a job to the server.

        Args:
            model: Model name (e.g., 'boltzgen')
            command: Command to execute
            config: Optional job configuration
            persistent: Keep job files after completion
            files: Optional dict of {filename: local_path} to upload

        Returns:
            JobState with job_id and initial status
        """
        self._ensure_local_server()

        if files:
            # Use multipart form for file uploads
            form_data = {
                "model": model,
                "command": json.dumps(command),
                "persistent": str(persistent).lower(),
            }
            file_uploads = [
                ("files", (name, open(path, "rb")))
                for name, path in files.items()
            ]

            try:
                resp = self._session.post(
                    f"{self.base_url}/jobs/upload",
                    data=form_data,
                    files=file_uploads,
                    timeout=60,
                )
            finally:
                for _, (_, f) in file_uploads:
                    f.close()
        else:
            resp = self._session.post(
                f"{self.base_url}/jobs",
                headers=self._headers(),
                json={
                    "model": model,
                    "command": command,
                    "config": config or {},
                    "persistent": persistent,
                },
                timeout=60,
            )

        resp.raise_for_status()
        return JobState.from_dict(resp.json())

    def status(self, job_id: str) -> JobState:
        """Get the current status of a job."""
        self._ensure_local_server()

        resp = self._session.get(
            f"{self.base_url}/jobs/{job_id}",
            headers=self._headers(),
            timeout=30,
        )
        resp.raise_for_status()
        return JobState.from_dict(resp.json())

    def result(self, job_id: str) -> JobResult:
        """Get the result of a completed job."""
        self._ensure_local_server()

        resp = self._session.get(
            f"{self.base_url}/jobs/{job_id}/result",
            headers=self._headers(),
            timeout=30,
        )
        resp.raise_for_status()
        return JobResult.from_dict(resp.json())

    def cancel(self, job_id: str) -> bool:
        """Cancel a running job."""
        self._ensure_local_server()

        resp = self._session.delete(
            f"{self.base_url}/jobs/{job_id}",
            headers=self._headers(),
            timeout=30,
        )
        return resp.status_code == 200

    def list_jobs(
        self,
        model: Optional[str] = None,
        status: Optional[JobStatus] = None,
    ) -> List[JobState]:
        """List jobs, optionally filtered."""
        self._ensure_local_server()

        params = {}
        if model:
            params["model"] = model
        if status:
            params["status"] = status.value

        resp = self._session.get(
            f"{self.base_url}/jobs",
            headers=self._headers(),
            params=params,
            timeout=30,
        )
        resp.raise_for_status()
        data = resp.json()
        return [JobState.from_dict(j) for j in data.get("jobs", [])]

    def download_file(self, job_id: str, file_path: str, dest: Path) -> Path:
        """Download an output file from a job."""
        self._ensure_local_server()

        resp = self._session.get(
            f"{self.base_url}/jobs/{job_id}/files/{file_path}",
            headers=self._headers(),
            stream=True,
            timeout=300,
        )
        resp.raise_for_status()

        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
        return dest

    def submit_run_dir(self, run_dir: Path, model: str) -> JobState:
        """Submit a prepared run directory to remote server.

        Packages inputs/, config.yaml, job.json into tar and uploads.
        The run_id is derived from the run directory name.

        Args:
            run_dir: Path to the prepared run directory (data/models/<model>/<run_id>/)
            model: Model name for tracking

        Returns:
            JobState with job_id (= run_id) and initial status
        """
        import io
        import tarfile

        self._ensure_local_server()

        # Create tar of the run directory (excluding outputs/)
        tar_buffer = io.BytesIO()
        with tarfile.open(fileobj=tar_buffer, mode='w:gz') as tar:
            for item in run_dir.iterdir():
                if item.name != 'outputs':
                    tar.add(item, arcname=item.name)
        tar_buffer.seek(0)

        # Upload
        run_id = run_dir.name
        resp = self._session.post(
            f"{self.base_url}/jobs/run/{run_id}",
            files={"run_package": (f"{run_id}.tar.gz", tar_buffer, "application/gzip")},
            data={"model": model},
            timeout=120,
        )
        resp.raise_for_status()
        return JobState.from_dict(resp.json())

    def download_outputs(self, job_id: str, run_dir: Path) -> None:
        """Download outputs/ directory from remote job back to local run_dir.

        Args:
            job_id: The job ID (= run_id)
            run_dir: Local run directory to download outputs to
        """
        import io
        import tarfile

        self._ensure_local_server()

        resp = self._session.get(
            f"{self.base_url}/jobs/{job_id}/outputs",
            stream=True,
            timeout=600,
        )
        resp.raise_for_status()

        # Extract tar to run_dir/outputs/
        outputs_dir = run_dir / "outputs"
        outputs_dir.mkdir(parents=True, exist_ok=True)

        tar_buffer = io.BytesIO(resp.content)
        with tarfile.open(fileobj=tar_buffer, mode='r:gz') as tar:
            tar.extractall(outputs_dir)

    def wait_for_completion(
        self,
        job_id: str,
        timeout: float = 3600,
        poll_interval: float = 5.0,
    ) -> JobState:
        """Wait for a job to complete.

        Args:
            job_id: Job ID to wait for
            timeout: Maximum time to wait (seconds)
            poll_interval: Time between status checks

        Returns:
            Final JobState

        Raises:
            TimeoutError: If job doesn't complete within timeout
        """
        start = time.time()
        while time.time() - start < timeout:
            state = self.status(job_id)
            if state.status in (
                JobStatus.COMPLETED,
                JobStatus.FAILED,
                JobStatus.CANCELLED,
            ):
                return state
            time.sleep(poll_interval)

        raise TimeoutError(f"Job {job_id} did not complete within {timeout} seconds")


# Convenience function for quick access
def get_client(data_root: Optional[Path] = None) -> JobClient:
    """Get a JobClient instance."""
    return JobClient(data_root=data_root)
