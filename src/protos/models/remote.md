# Remote Model Deployment Architecture

## 🔑 Core Idea

Keep lightweight services always running:
- **API server** (FastAPI) - CPU-only, minimal resources
- **Job queue** (Redis/RabbitMQ) - Message passing and coordination
- **Controller** (K8s/Slurm/compose) - Orchestration logic

Spin up heavy GPU containers only when jobs arrive, then tear them down automatically when idle.

## 1. Docker + Orchestration Layer

You already have adapters (`Boltz2Adapter`, `RFdiffusionAdapter`) in `model_manager.py`. Wrap each model in its own Docker image:

- `boltz2:latest`
- `rfdiffusion:latest`
- `proteinmpnn:latest`
- `alphafold:latest`

### Launch Strategy

**Local Development**
```bash
docker run --gpus all --rm boltz2:latest
```

**Production Options**

1. **Kubernetes + GPU nodes**
   ```bash
   kubectl scale deployment boltz2=1  # Scale up on job arrival
   kubectl scale deployment boltz2=0  # Auto-scale back when idle
   ```

2. **Slurm/HPC**
   - Submit containerized jobs to GPU queue
   - Use Singularity/Apptainer for HPC compatibility

3. **Cloud Providers**
   - AWS Batch / GCP Vertex / Azure Batch
   - Spin up GPU VM only for job duration
   - Automatic shutdown after completion

## 2. Idle Cost Avoidance

GPU costs are only incurred when containers run on GPU nodes.

### Cost Management Options

**On-Premise GPU (Lab Server)**
- Idle GPUs = sunk cost
- Focus on maximizing utilization

**Cloud GPU (Pay-Per-Use)**
- Configure autoscaling with 0 GPU nodes by default
- Workflow: Job arrives → Scale to 1 GPU node → Run container → Scale back to 0
- Result: GPU hours = actual model runtime only

Example AWS setup:
```yaml
# EKS autoscaling group configuration
minSize: 0
maxSize: 4
desiredCapacity: 0
instanceTypes: ["g4dn.xlarge"]  # GPU instances
```

## 3. Student Workflow

From the student's perspective, nothing changes:

```python
result = predict("boltz2",
                 entity_name="MyProtein",
                 backend=ModelBackend.API)
```

### Behind the Scenes

1. **API server** receives request
2. **Job queue** stores job details (model: boltz2, inputs: data)
3. **Controller** logic:
   - Checks if boltz2 GPU service is running
   - If not, spins up the container
   - Container runs model, stores results in `/predictions`
4. **API server** returns job result or async job ID
5. **Auto-shutdown** after X minutes of inactivity

## 4. Implementation Options

### Simple (Fastest Setup)
- Keep one GPU server available
- Run containers on-demand: `docker run --rm ...`
- GPU always "reserved" but simple to manage

### Efficient (Cloud Cost-Aware)
- Deploy Kubernetes with GPU autoscaler
- Configure each model service with `replicas=0` default
- Scale up to 1 on job submission
- Automatic scale-down on idle timeout

### Academic/HPC Cluster
- Use SLURM jobs for each model execution
- Containerized with Singularity/Apptainer
- Protos API submits to SLURM, retrieves outputs
- Zero cost until GPU job is actually running

## 5. Model Characteristics

Different models require different handling:

| Model | Characteristics | Best Approach |
|-------|----------------|---------------|
| **Boltz2, RFdiffusion** | Heavy, long-running jobs | Async batch jobs with job ID system |
| **ProteinMPNN** | Lighter, faster | Interactive mode, possibly keep warm on CPU |
| **AlphaFold** | Extremely heavy | Definitely async, not suitable for interactive use |

## ✅ Recommended Architecture

For academic/student use cases:

1. **Deploy FastAPI inference server** (CPU-only)
2. **Add job queue and worker containers**
3. **Workers launch on-demand** via Kubernetes autoscaler or SLURM
4. **Students always use the same API** - GPU management is transparent

### Benefits
- Students never need local GPUs
- Pay only for actual prediction runtime
- Support both fast/light (ProteinMPNN) and heavy jobs (RFdiffusion, AlphaFold)
- Unified interface for all models

## 🔑 Why Separate Docker Images Per Model

### 1. **No Dependency Conflicts**
Boltz, RFdiffusion, ProteinMPNN, and AlphaFold have very different requirements:
- Different CUDA versions
- Incompatible PyTorch/JAX versions
- Conflicting Python package versions

### 2. **Smaller Images & Faster Pulls**
- Students needing only ProteinMPNN don't download 20GB of Boltz+AlphaFold+JAX
- Faster startup times for specific models

### 3. **Easier Scaling**
- Run only needed model containers
- Boltz container for design jobs
- RFdiffusion container for structure generation
- Scale independently based on demand

### 4. **Independent Updates**
- Update Boltz to new version without rebuilding AlphaFold
- Test new model versions in isolation

### 5. **Cost Control**
- Kubernetes/Slurm can schedule per-container GPU jobs
- Idle models = no running container = no GPU costs

## ⚙️ Project Structure

```
docker/
├── boltz/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── entrypoint.py
├── rfdiffusion/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── entrypoint.py
├── proteinmpnn/
│   ├── Dockerfile
│   ├── requirements.txt
│   └── entrypoint.py
└── alphafold/
    ├── Dockerfile
    ├── requirements.txt
    └── entrypoint.py
```

### Dockerfile Template

```dockerfile
FROM nvidia/cuda:<version>-cudnn-runtime

# Install Python and models-specific dependencies
RUN apt-get update && apt-get install -y python3-pip
COPY requirements.txt .
RUN pip3 install -r requirements.txt

# Model weights strategy:
# Option 1: Copy into image (larger image, stable)
COPY weights/ /models/weights/

# Option 2: Mount at runtime (smaller image, flexible)
# docker run -v /path/to/weights:/models/weights

# Copy application code
COPY entrypoint.py /app/
WORKDIR /app

# Run models in server mode
CMD ["python3", "entrypoint.py"]
```

### Entrypoint Template

Each `entrypoint.py` runs the model as a service:
- FastAPI app exposing `/predict` endpoint
- Reads from `BaseModel` class (reuse existing `base_model.py`)
- Handles model loading and inference

## 🌀 Deployment Example (Kubernetes)

### Architecture
- **1 FastAPI gateway** (always running, CPU-only)
- **N model services** (scale 0→1 on demand):
  - `boltz-service`
  - `rfdiffusion-service`
  - `proteinmpnn-service`
  - `alphafold-service`

### Autoscaling Flow
```
Student Request → UnifiedModelClient → API Gateway
                                           ↓
                                    Job Queue (Redis)
                                           ↓
                                    Controller checks:
                                    Is boltz-service running?
                                           ↓
                                    No: Scale up to 1
                                           ↓
                                    Run prediction
                                           ↓
                                    Return results
                                           ↓
                                    After idle timeout:
                                    Scale down to 0
```

This architecture ensures efficient resource utilization while maintaining a simple interface for end users.