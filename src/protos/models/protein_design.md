# Protein-design integrations in this repository

This page records standalone protein-design utilities implemented in ProtOS.
It is not a general protein-design protocol and does not claim that external
models are installed or experimentally validated.

## RFdiffusion2 helper

The tracked RFdiffusion2 helper modules expose:

- `RFD2Config` and `run_rfdiffusion2()` in `rfdiffusion2_runner.py`;
- CIF-to-PDB/ORI preparation helpers in `rfdiffusion2_utils.py`; and
- `RFD2ConfigBuilder` plus convenience builders in
  `rfdiffusion2_configs.py`.

`run_rfdiffusion2()` executes Apptainer against this exact expected image:

`src/protos/models/RFdiffusion2/rf_diffusion/exec/bakerlab_rf_diffusion_aa.sif`

The image and model weights are workstation artifacts and are not tracked. The
function raises `FileNotFoundError` if the SIF or input PDB is missing and
raises `RuntimeError` if the Apptainer command fails.

The builder can be inspected without running the model:

```python
from pathlib import Path
from protos.models.rfdiffusion2_runner import RFD2Config

config = RFD2Config(
    input_pdb=Path("input.pdb"),
    output_dir=Path("designs"),
    num_designs=4,
    ligand="RET",
)

args = config.to_command_args()
assert "inference.num_designs=4" in args
assert "inference.ligand=RET" in args
```

Building arguments does not validate RFdiffusion2 semantics or execute a model.
Use the pinned upstream documentation for supported Hydra options.
