# Experimental job server and client

`job_server.py` and `job_client.py` implement an experimental Docker-backed job
transport. They are not a production deployment architecture.

## Important limitations

- The server accepts commands and executes them in model containers.
- The current server has no authentication or authorization enforcement.
- Uploaded archives and file paths have not been hardened for untrusted input.
- Job state is an in-memory dictionary persisted to one JSON file.
- A DELETE request marks a running job cancelled but does not stop its Docker
  process.
- The client expects a `protos/job-server:latest` image, but this repository
  does not currently contain a maintained Dockerfile that builds that image.

Do not expose this service to an untrusted network. Treat it as development
code until those constraints are addressed.

## Python dependencies

The server imports FastAPI, Uvicorn, Pydantic, and multipart upload support; the
client imports Requests and PyYAML. FastAPI/Uvicorn/python-multipart are not in
ProtOS's base dependencies.

## Configuration

`JobClient(data_root=...)` reads `<data_root>/models/server_config.yaml`. Its
default `data_root` is `~/.protos`, independently of `PROTOS_DATA_ROOT`. If the
file is absent, the client writes a default local configuration.
`ServerConfig` supports:

```yaml
mode: remote
local:
  host: localhost
  port: 8765
  auto_start: false
  server_image: protos/job-server:latest
remote:
  url: https://jobs.example.invalid
  api_key: ${PROTOS_API_KEY}
  timeout: 3600
```

The `${PROTOS_API_KEY}` text is only a serialization placeholder. On load,
`ServerConfig.from_file()` uses the `PROTOS_API_KEY` environment variable when
set. Although the client sends a Bearer header in remote mode, the current
server does not validate it.

## Implemented HTTP routes

| Method and path | Behavior |
| --- | --- |
| `GET /health` | return server health/timestamp |
| `POST /jobs` | create and background a Docker job |
| `POST /jobs/upload` | create a job with uploaded files |
| `POST /jobs/run/{run_id}` | upload a prepared run archive |
| `GET /jobs` | list jobs, optionally by model/status |
| `GET /jobs/{job_id}` | return state |
| `GET /jobs/{job_id}/result` | return logs/result after completion/failure |
| `GET /jobs/{job_id}/outputs` | download the outputs directory as tar.gz |
| `GET /jobs/{job_id}/files/{path}` | download one file |
| `DELETE /jobs/{job_id}` | mark running state cancelled or delete stored state |

The server maps `boltzgen`, `boltz`, and `lambda` to Docker image names and uses
`protos/base:latest` for other model strings unless `config.docker_image` is
provided.

## Client surface

`JobClient` provides `submit`, `status`, `result`, `list_jobs`, `cancel`,
`download_file`, `submit_run_dir`, `download_outputs`, and
`wait_for_completion`. These methods require a reachable configured server.
Constructing a client alone does not start Docker; the first request calls the
local auto-start logic when `mode` is `local`.
