# Protein-design integrations in this repository

This page records implemented ProtOS integration points. It is not a general
protein-design protocol and does not claim that a registered model card is
installed or experimentally validated.

## ModelManager integrations

`ModelManager` currently registers these design-related external cards:

| Card | Implemented result of `prepare()` | External requirement |
| --- | --- | --- |
| `boltz2` | Boltz prediction YAML/FASTA and `boltz predict` command | Boltz CLI or maintained container |
| `boltzgen` | BoltzGen design configuration and `boltz design` command | pinned BoltzGen checkout/environment |
| `ligand_mpnn` | command payload from a protein PDB, optional ligand, and constraints | LigandMPNN checkout/environment |
| `pocket2mol` | command payload from a structure, optional ligand, and bounding box | Pocket2Mol checkout/environment |

Other registered cards concern docking or property prediction; see
`../../../docs/05_model_manager.md`.

There are no `RFdiffusionAdapter`, `ProteinDesignWorkflow`, binder preset, or
enzyme preset classes in ProtOS. Older examples using those names were design
sketches, not executable API, and have been removed.

## RFdiffusion2 helper

RFdiffusion2 is not a `ModelManager` card. The tracked helper module exposes:

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
