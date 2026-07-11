# Acquisition and processing

ProtOS separates acquisition from processing.

- A loader obtains a source artifact and registers it.
- A processor reads and writes ProtOS's managed representation.
- `EntityRegistry` records names, formats, file paths, metadata, and
  relationships.
- `DatasetManager` stores named collections for one processor type.

`load_entity()` does not imply a network request. For example,
`StructureProcessor.load_entity()` reads a registered or cached canonical PKL;
use `StructureLoader` first when the structure is not managed yet.

## Structure acquisition

`StructureLoader` supports these sources:

| Source | Accepted input | Network required |
| --- | --- | --- |
| RCSB | a four-character PDB ID such as `1ubq` | yes |
| AlphaFold DB | a UniProt accession or `AF-...-F1[-model_vN]` ID | yes |
| local | an existing CIF/mmCIF file path (optionally gzipped) | no |

Downloaded or imported CIF/mmCIF data is canonicalized through the associated
`StructureProcessor` and stored as PKL.

```python
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor

sp = StructureProcessor()
loader = StructureLoader(processor=sp)

name = loader.download_and_register("1ubq")
if name is None:
    raise RuntimeError("RCSB download failed")

structure = sp.load_entity(name)
assert structure is not None
```

The example performs an RCSB network request the first time it runs.

## Sequence acquisition

`SequenceLoader` accepts an existing FASTA path, a bare UniProt accession, or
an identifier such as `uniprot:P00533,P07550`. A UniProt identifier requires a
network request. For deterministic local data, call
`SequenceProcessor.save_entity()` or use an existing FASTA file.

## BaseLoader behavior

Concrete loaders inherit the following implemented operations from
`protos.io.core.base_loader.BaseLoader`:

| Method | Behavior |
| --- | --- |
| `download_and_register(...)` | fetch one artifact and register it |
| `download_batch(...)` | process several identifiers and optionally create a dataset |
| `get_download_path(...)` | resolve a managed loader path |
| `check_entity_exists(...)` | query the registry for the loader's format |

`fetch_entity()` and `parse_identifier()` are abstract source-specific
operations. Import concrete loaders from their modules under
`protos.io.ingest`; they are not exported from the top-level `protos` package.
