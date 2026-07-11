# ProtOS

ProtOS is a Python framework for managing structural-biology data with named
entities, datasets, processors, and a shared registry.

This is the archival `thesis` branch. The submitted-thesis material and its
reproduction notes are in [`thesis/README.md`](thesis/README.md). Ongoing ProtOS
development belongs on `master`.

## Installation

```bash
git clone --branch thesis https://github.com/flurinh/protos.git
cd protos
python -m pip install -e .
```

ProtOS stores managed data in `~/protos_data` by default. Set
`PROTOS_DATA_ROOT` before starting Python if you want another location.

## Verified structure example

Acquisition and processing are separate operations: `StructureLoader` fetches
and registers source files, while `StructureProcessor` loads their canonical
representations and manages datasets.

```python
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor

sp = StructureProcessor()
loader = StructureLoader(processor=sp)

# Fetch an RCSB mmCIF file, register it, and store its canonical representation.
name = loader.download_and_register("1ubq")
if name is None:
    raise RuntimeError("RCSB download failed")

structure = sp.load_entity(name)

# Structure datasets are stacked DataFrames by default; request a mapping here.
sp.create_dataset("example_structures", [name])
structures = sp.load_dataset("example_structures", return_format="dict")
```

This example was executed against the archived branch. It requires network
access to RCSB for the initial download.

## Reference data

The curated GRN reference tables are packaged in
`src/protos/reference_data/grn/` and copied into the configured data root during
project-data setup. See
[`src/protos/reference_data/grn/README.md`](src/protos/reference_data/grn/README.md)
and its checksum manifest for the exact archived bundle.

The longer files under [`docs/`](docs/) describe the historical API surface;
for reproducible thesis figures, use the branch-specific instructions in
[`thesis/README.md`](thesis/README.md).
