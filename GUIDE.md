# Protos Workflow Guide

This guide explains how Protos organises data through *entities*,
the global *registry*, *datasets*, and zero‑configuration *ProtosPaths*. Each
section links to the runnable `test_*.py` scripts that demonstrate the
concepts in practice.

---

## 1. Zero‑Configuration Data Root

Every script calls `protos.set_data_path(<directory>)` (see for example
`test_sequence_workflow.py`, `test_structure_alignment.py`, `test_ligand_workflow.py`).
Once set, the singleton **ProtosPaths** object creates a deterministic
directory tree (`structure/`, `sequence/`, `ligand/`, `graph/`, `property/`, …)
inside that root. Processors should never hardcode file paths; instead they call
`self.get_subdirectory_path(...)` inherited from `BaseProcessor`.

ProtosPaths responsibilities:
- lazily create processor directories only when accessed;
- expose helper methods such as `get_processor_path('input')` used by scripts to
  stage FASTA/SDF files before ingestion (see `test_sequence_workflow.py`).

---

## 2. Entities and the Global Registry

An **entity** is the canonical representation of a data object (structure, SMILES,
graph, GRN table row, property row, …). Each entity has:
- a human‑readable name (`original_id`);
- a `format_type` identifying the owning processor (e.g. `sequence`, `ligand`,
  `graph`);
- optional metadata and an optional relative file path.

When processors persist data they must call
`self.entity_registry.register_entity(...)`. The registry automatically tracks
relationships between entities. For example:

- `test_structure_alignment.py` registers aligned structures and records
  `derived_from` links to the original PDB codes.
- `test_graph_workflow.py` generates graphs and records `graph -> structure`
  relationships.
- `test_ligand_workflow.py` registers both SMILES ligands and structure‑derived
  ligands, again linking them back to the source structure with
  `derived_from`.

You can inspect relationships via helper methods. For instance
`StructureProcessor.list_dataset_related_sequences()` (see
`test_cross_processor_annotation.py`) surfaces sequence entities linked to each
structure.

---

## 3. Datasets and `DatasetManager`

A **dataset** is a named collection of entities plus metadata. `BaseProcessor`
exposes high‑level helpers (`create_dataset`, `add_to_dataset`, `load_dataset`,
`export_dataset`) that wrap the underlying `DatasetManager`.

- `test_sequence_workflow.py` shows FASTA ingestion that creates two datasets:
  `single_sequence` (materialised entity) and `demo_sequences` (dataset pointing
  to multiple sequences). Later it calls `export_dataset` to write FASTA files.
- `test_graph_workflow.py` generates an atom‑level graph dataset named
  `graph_demo_structures__atom_graphs`.
- `test_ligand_workflow.py` registers a SMILES dataset and, when hetero atoms are
  present, creates a structure‑derived ligand dataset `<pdb>__ligands`.
- `test_property_workflow.py` records two property tables that are themselves
  datasets stored under the PropertyProcessor.

Always access dataset contents via `processor.load_dataset(<name>)` so that the
registry can resolve current entity names even if files move.

---

## 4. Processor Overview & Test Scripts

### StructureProcessor
- Loads canonical PKL frames, aligns structures, and exports CIFs.
- Demonstrated in `test_structure_alignment.py` (CEalign alignment pipeline) and
  `test_cross_processor_annotation.py` (extract chains → register sequences).

### SequenceProcessor
- Handles single and multi‑sequence FASTA ingestion, alignment utilities, mutant
  libraries, and dataset export.
- See `test_sequence_workflow.py` for ingestion/export, `test_sequence_alignment.py`
  for reference alignment, and `test_sequence_real_workflow.py` for UniProt
  downloads.

### GRNProcessor & Sequence integration
- `test_grn_workflow.py` aligns GPCR sequences against the `gpcrdb_ref` table and
  records a GRN dataset. The TODO entry explains current limitations (MMseqs
  dependency, alignment fallbacks).

### GraphProcessor *(new)*
- Converts structures to PyTorch Geometric–style graphs; stores metadata about
  the featurisation and geometry.
- Demonstrated in `test_graph_workflow.py`. PyG is optional—when missing the
  script still persists NumPy payloads and reports node/edge counts.

### LigandProcessor *(new)*
- Registers ligands from SMILES or structure extractions, stores optional SDF
  artefacts, and computes primitive ligand–protein contacts.
- Demonstrated in `test_ligand_workflow.py`.
- Interactions can be logged to `PropertyProcessor` for downstream analysis and
  QSAR modelling (see TODO section).

### EmbeddingProcessor
- Generates per-residue or global embeddings (ESM/Ankh). Optional torch +
  transformers dependencies.
- `test_embedding_workflow.py` iterates over all registered models and gracefully
  skips when dependencies are absent.

### PropertyProcessor
- Stores tabular annotations keyed by entities. Each row becomes a property
  entity linked to the annotated scope (structure, sequence, ligand…).
- `test_property_workflow.py` demonstrates recording sequence alignment scores
  and structure classifications in property tables.

---

## 5. Workflow Index by Script

| Script | What it demonstrates | Key processors |
| --- | --- | --- |
| `test_sequence_workflow.py` | FASTA ingestion, dataset export, mutant library generation | SequenceProcessor |
| `test_sequence_alignment.py` | Pairwise alignment utilities | SequenceProcessor |
| `test_sequence_real_workflow.py` | Remote UniProt downloads and dataset persistence | SequenceProcessor |
| `test_structure_alignment.py` | GPCR structure downloads, CEalign alignment, CIF export | StructureProcessor |
| `test_cross_processor_annotation.py` | Structure → sequence chain extraction and annotation | StructureProcessor, SequenceProcessor |
| `test_grn_workflow.py` | Align sequences to GRN reference tables and record annotations | SequenceProcessor, GRNProcessor |
| `test_graph_workflow.py` | Build PyG-compatible graphs from structures, register datasets | GraphProcessor |
| `test_ligand_workflow.py` | Register SMILES ligands, extract structure ligands, compute contacts | LigandProcessor, StructureProcessor, PropertyProcessor |
| `test_property_workflow.py` | Record property tables for sequences and structures | PropertyProcessor |
| `test_embedding_workflow.py` | Iterate across embedding models, skip when dependencies missing | EmbeddingProcessor |

Use these scripts as reference templates when developing new workflows: set the
data root, instantiate the required processor(s), load or ingest data through
their zero-config APIs, and rely on the registry to persist entities, datasets,
and relationships.

