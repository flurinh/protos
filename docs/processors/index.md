# Processors

Processors are the core data manipulation components of Protos. Each processor manages a specific data type, providing operations for data registration, loading, saving, and analysis.

## Design Principles

1. **Zero Configuration** - All processors work out of the box with no path parameters
2. **Human-Readable Names** - All public APIs use names, never UUIDs
3. **Automatic Registration** - Data is automatically registered with EntityRegistry
4. **ProtosPaths Integration** - All path operations managed through the singleton

## Available Processors

| Processor | Type | Description |
|-----------|------|-------------|
| [StructureProcessor](structure_processor.md) | `structure` | Protein structures as DataFrames |
| [SequenceProcessor](sequence_processor.md) | `sequence` | Sequences with alignment capabilities |
| [GRNProcessor](grn_processor.md) | `grn` | Generic Residue Numbering |
| [EmbeddingProcessor](embedding_processor.md) | `embedding` | Protein embeddings (ESM-2, Ankh) |
| [GraphProcessor](graph_processor.md) | `graph` | PyTorch Geometric graphs |
| [PropertyProcessor](property_processor.md) | `property` | Tabular properties linked to entities |
| [MoleculeProcessor](molecule_processor.md) | `molecule` | Small molecule descriptors |

## BaseProcessor Architecture

All processors inherit from `BaseProcessor` (`protos.io.core.base_processor`).

### Core Attributes

| Attribute | Description |
|-----------|-------------|
| `processor_type` | Type identifier (e.g., "structure", "sequence") |
| `paths` | ProtosPaths singleton instance |
| `entity_registry` | EntityRegistry singleton instance |
| `dataset_manager` | DatasetManager for this processor type |
| `data_path` | Root path for processor data |

### Abstract Methods (must implement)

```python
@abstractmethod
def load_entity(self, name: str) -> Any:
    """Load entity by human-readable name."""
    pass

@abstractmethod
def save_entity(self, name: str, data: Any, metadata: dict = None):
    """Save entity with human-readable name."""
    pass
```

### Common Methods

| Method | Description |
|--------|-------------|
| `list_entities()` | List all entities of this processor type |
| `entity_exists(name)` | Check if entity exists |
| `delete_entity(name)` | Delete entity from registry and storage |
| `export_entity(name, path, format)` | Export entity to file |
| `create_dataset(name, entities)` | Create dataset from entities |
| `load_dataset(name)` | Load all entities in a dataset |
| `list_datasets()` | List available datasets |

---

## Data Registration

All processors automatically register entities with the EntityRegistry when saving:

```python
proc = StructureProcessor()

# This registers the entity automatically
proc.save_entity("my_protein", df, metadata={"source": "custom"})

# Entity is now accessible
from protos.io.core import get_registry
registry = get_registry()
info = registry.find_entity("my_protein", "structure")
```

### Bulk Registration

Register multiple files from a directory:

```python
report = proc.register_directory(
    Path("/path/to/files"),
    extensions=[".cif", ".pdb"],
    recursive=False,
    dry_run=False
)

print(f"Registered: {len(report['registered'])}")
print(f"Skipped: {len(report['skipped'])}")
print(f"Errors: {len(report['errors'])}")
```

---

## Loading Entities and Datasets

### Single Entity

```python
# Load by name
data = proc.load_entity("my_protein")
```

### Dataset Operations

```python
# Create dataset
proc.create_dataset("study_set", ["protein_a", "protein_b", "protein_c"])

# Load all entities in dataset
data_dict = proc.load_dataset("study_set")
for name, data in data_dict.items():
    print(f"Processing {name}")

# Modify dataset
proc.add_to_dataset("study_set", ["protein_d"])
proc.remove_from_dataset("study_set", ["protein_a"])

# Get dataset info
info = proc.get_dataset_info("study_set")
print(f"Contains {info['entity_count']} entities")
```

---

## Processor Dependencies

Some processors depend on others for cross-format analysis:

| Processor | Dependencies |
|-----------|--------------|
| GRNProcessor | SequenceProcessor (for alignment) |
| EmbeddingProcessor | SequenceProcessor (sequence data) |
| GraphProcessor | StructureProcessor (structure data) |
| PropertyProcessor | Any (annotates existing entities) |

### Example: Cross-Processor Pipeline

```python
from protos import (
    StructureProcessor, SequenceProcessor,
    GRNProcessor, EmbeddingProcessor
)

# Load structure
struct_proc = StructureProcessor()
df = struct_proc.load_entity("3sn6")

# Extract sequence from chain
chain_seq = "..."  # Extract from structure

# Get sequence processor
seq_proc = SequenceProcessor()
seq_proc.save_entity("3sn6_chain_R", chain_seq)

# Annotate with GRN
grn_proc = GRNProcessor()
annotations, _ = grn_proc.annotate_sequences(
    {"3sn6_chain_R": chain_seq},
    reference_table="class_a_gpcr",
    protein_family="gpcr_a"
)

# Generate embeddings
emb_proc = EmbeddingProcessor(model_name="esm2_t30_150m")
embedding = emb_proc.generate_embedding(chain_seq)
emb_proc.save_entity("3sn6_chain_R", embedding)
```

---

## Relationship Resolution

Processors can resolve related entities:

```python
# Get related structures for a sequence
related = seq_proc.list_source_structures(["seq_1", "seq_2"])

# Get related entities for a dataset
related = struct_proc.resolve_related_entities_for_dataset(
    "my_dataset",
    rel_type="derived_from",
    direction="outgoing"
)
```
