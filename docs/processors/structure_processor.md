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

## Structure-level GRN annotation and Plotly output

`annotate_with_grn()` extracts one C-alpha record per polymer residue, annotates
that sequence, and maps the resulting GRNs back onto every atom of each matched
residue. Polymeric modified residues stored as `HETATM` (for example MSE) are
included. `label_seq_id` defines intrinsic polymer order; author numbering plus
insertion code is the fallback lookup key.

Experimental constructs often include T4 lysozyme, endolysin, antibodies, or
other fusion partners. Consequently, generated insertion coordinates are
disabled by default for structure annotation. Use `residue_ranges={"A":
(start, end)}` to isolate a biological component when multiple components share
one author chain, then pass `assign_insertions=True` to assign genuine
intra-segment insertions and directional flexible-region labels inside that
component. `return_summary=True` returns the selected reference, coverage,
indel counts, and status for each chain.

`plot_grn_ca_structure()` produces a Plotly figure with one marker per C-alpha,
GRN hover text, and edges between sequential C-alpha atoms within each chain.
Edges longer than 4.5 Å are omitted by default so unresolved structure gaps are
not bridged visually.

The following exercised script downloads and annotates a standalone GPCR,
arrestin, G-protein heterotrimer, GPCR–Gs complex, and GPCR–arrestin complex,
then writes individual self-contained HTML files and a combined gallery:

```bash
python scripts/visualize_grn_structures.py \
  --data-root data/grn_structure_demo \
  --output-dir data/visualizations/grn_structures
```

See [`scripts/visualize_grn_structures.py`](../../scripts/visualize_grn_structures.py)
for the exact PDB chains, component ranges, and reference tables used.
