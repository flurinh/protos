# StructureProcessor

Import: `from protos.processing.structure import StructureProcessor`

`StructureProcessor` stores one canonical pandas DataFrame per structure as a
PKL. Its index is `(structure_id, atom_id)`. Source CIF/PDB files are acquisition
artifacts; `load_entity()` itself only reads an already registered or cached
canonical PKL and returns `None` when no managed structure exists.

## Acquire and load a structure

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
assert structure.index.names == ["structure_id", "atom_id"]
```

This example requires network access for the initial RCSB download.

## Storage and datasets

| Method | Current behavior |
| --- | --- |
| `save_entity(id, df, format="pkl", metadata=None)` | canonicalize, write PKL, and register; other save formats are rejected |
| `load_entity(id)` | return one canonical DataFrame or `None` |
| `export_entity(id, out_path, format=...)` | materialize CIF, PDB, or SDF through the structure exporter |
| `create_dataset(name, ids)` | store a logical dataset definition |
| `load_dataset(name)` | return one stacked DataFrame |
| `load_dataset(name, return_format="dict")` | return `{structure_id: DataFrame}` |
| `save_dataset(name, ids, metadata=None)` | ensure loaded frames are persisted, then create the dataset |

## Implemented analysis and editing groups

- ligand, water, and ion summaries/contacts;
- chain, residue-range, query, and HETATM filtering;
- atom addition/deletion, annotation, transforms, merge, and reindexing;
- chain-sequence extraction and registration;
- pairwise and batch structural alignment;
- sequence-based GRN annotation; and
- CIF/PDB export.

These operations work on managed DataFrames. Methods that need an entity call
`load_entity()` and therefore do not download a missing structure. Refer to the
source signature in
`src/protos/processing/structure/structure_processor.py` for the operation's
filters and output options.

Deprecated compatibility methods `load_structure`, `save_structure`, and
`load_structures` remain present, but new code should use the entity and dataset
methods above.
