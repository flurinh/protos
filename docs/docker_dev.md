# Docker Development Notes - BoltzGen Integration

## Current Status (2026-01-31)

### Working Components

1. **NVIDIA Container Toolkit** - Installed and configured
   - RTX 5090 (32GB VRAM) accessible via Docker
   - Test: `docker run --rm --gpus all nvidia/cuda:12.2.2-base-ubuntu22.04 nvidia-smi -L`

2. **BoltzGen Docker Image** - Built successfully
   - Image: `protos/boltzgen:latest`
   - Base: `nvidia/cuda:12.2.2-cudnn8-devel-ubuntu22.04`
   - Location: `src/protos/models/boltzgen/Dockerfile`

3. **BoltzGenAdapter CLI Fix** - Completed
   - File: `src/protos/models/model_manager.py:1154-1168`
   - Changed from invalid `--recycling_steps` to `--config design recycling_steps=X`

4. **MDM2 Binder Design** - Successfully completed
   - Script: `scripts/boltzgen_real_workflow.py`
   - Designed 4 peptide binders against MDM2 (PDB: 1YCR)

### Known Issue: Docker Shared Memory

**Error:**
```
RuntimeError: unable to allocate shared memory(shm) for file </torch_XXX>: Resource temporarily unavailable (11)
```

**Cause:**
Docker's default `/dev/shm` size is 64MB, which is insufficient for PyTorch's multiprocessing data loaders.

**Fix Required:**
Add `--shm-size` flag to Docker run command in `DockerJobExecutor._run_job()`:

```python
# In model_manager.py, around line 467
cmd = ["docker", "run", "--rm"]

# GPU support
if self.use_gpu:
    cmd.extend(["--gpus", "all"])

# FIX: Add shared memory size for PyTorch
cmd.extend(["--shm-size", "8g"])  # Add this line
```

**Location:** `src/protos/models/model_manager.py:467`

### Files Modified

| File | Changes |
|------|---------|
| `src/protos/models/model_manager.py` | Fixed BoltzGen CLI args (lines 1154-1168) |
| `scripts/test_boltzgen_docker.py` | Fixed `devices: 0` config |
| `scripts/boltzgen_minimal_example.py` | Created minimal workflow example |
| `scripts/boltzgen_real_workflow.py` | Created real MDM2 binder workflow |
| `thesis/workflows/boltzgen_rhodozyme_workflow.py` | Created rhodozyme design workflow |

### Pending Work

1. **Fix Docker shared memory** - Add `--shm-size=8g` to DockerJobExecutor
2. **Complete rhodozyme workflow** - Test with all 3 substrates:
   - Benzamidine (trypsin)
   - Phenylalanine (papain)
   - Leucine (subtilisin)
3. **Add structure ingestion to model_manager** - Auto-register CIF outputs

### Test Commands

```bash
# Verify GPU access
docker run --rm --gpus all nvidia/cuda:12.2.2-base-ubuntu22.04 nvidia-smi

# Build BoltzGen image
cd src/protos/models/boltzgen && docker build -t protos/boltzgen:latest .

# Run MDM2 binder design (working)
python scripts/boltzgen_real_workflow.py

# Run rhodozyme design (needs shm fix)
python thesis/workflows/boltzgen_rhodozyme_workflow.py
```

### Environment

- Host CUDA: 12.8.1
- Host Driver: 570.124.06
- Docker Image CUDA: 12.2.2
- GPU: NVIDIA GeForce RTX 5090 (32GB)
- Python: 3.10 (host), 3.11 (container)

### Architecture Notes

```
ModelManager
    └── DockerJobExecutor
            └── docker run --gpus all -v {working_dir}:/workspace protos/boltzgen:latest
                    └── boltzgen run config.yaml --output predictions --config design ...
```

**Volume Mounts:**
- `{working_dir}` → `/workspace` (job config + inputs)
- `~/.cache/huggingface` → `/cache` (model weights)

**BoltzGen Output Structure:**
```
predictions/
├── config.cif                    # Input visualization
├── steps.yaml                    # Pipeline steps
├── intermediate_designs/         # Raw designs
├── intermediate_designs_inverse_folded/
│   ├── refold_cif/              # Refolded structures ← Register these
│   └── refold_design_cif/       # Binder-only structures
└── final_ranked_designs/         # Filtered designs
```

### Key Code Locations

| Component | File | Line |
|-----------|------|------|
| DockerJobExecutor | `model_manager.py` | 216-528 |
| BoltzGenAdapter | `model_manager.py` | 1091-1405 |
| GPU detection | `model_manager.py` | 283-300 |
| Job submission | `model_manager.py` | 370-400 |
| Docker run cmd | `model_manager.py` | 465-510 |
