# Protos Documentation

Protos is a unified data access and processing framework for computational biology. It provides zero-configuration data management for protein structures, sequences, molecules, and related biological data with automatic entity registration and relationship tracking.

## Documentation Structure

### [1. Unified Data Access](01_unified_data_access.md)
Loaders for acquiring data from external sources:
- **StructureLoader** - Fetch structures from RCSB PDB, AlphaFold, or local files
- **SequenceLoader** - Fetch sequences from UniProt or local FASTA files
- **LigandLoader** - Import molecules from SDF files or SMILES strings

### [2. Entity, Registry and Datasets](02_entity_registry_datasets.md)
Centralized entity management and dataset operations:
- **EntityRegistry** - Track entities across formats with human-readable names
- **Relationships** - Define connections between entities (derived_from, subset_of, etc.)
- **DatasetManager** - Create and manage collections of entities

### [3. Zero-Configuration / ProtosPaths](03_zero_configuration.md)
Automatic path management with no configuration required:
- **ProtosPaths** - Singleton path manager for all data directories
- **Environment Variables** - Optional customization via `PROTOS_DATA_ROOT`
- **Directory Structure** - Standard layout for each processor type

### [4. Processors](processors/index.md)
Data registration, loading, operations, and analytical capabilities:
- **[StructureProcessor](processors/structure_processor.md)** - Manage protein structures with DataFrame representation
- **[SequenceProcessor](processors/sequence_processor.md)** - Manage sequences with alignment capabilities
- **[GRNProcessor](processors/grn_processor.md)** - Generic Residue Numbering assignment and mapping
- **[EmbeddingProcessor](processors/embedding_processor.md)** - Generate and store protein embeddings
- **[GraphProcessor](processors/graph_processor.md)** - Create PyTorch Geometric graphs from structures
- **[PropertyProcessor](processors/property_processor.md)** - Store tabular properties linked to entities
- **[MoleculeProcessor](processors/molecule_processor.md)** - Manage ligand descriptors (SMILES, InChI)

### [5. Model Manager](05_model_manager.md)
Model orchestration and job submission:
- **ModelCard** - Declarative model specifications
- **ArtifactSpec** - Input/output artifact definitions
- **Data Packaging** - Assemble artifacts from processors
- **Job Submission** - Create execution packages for external models

## Quick Start

```python
from protos import (
    StructureLoader, SequenceLoader,
    StructureProcessor, SequenceProcessor,
    set_data_path
)

# Optional: Set custom data path (must call BEFORE creating loaders/processors)
# set_data_path("/my/custom/data/root")

# Load structures from PDB
loader = StructureLoader()
loader.download_and_register("1ubq")  # RCSB PDB
loader.download_and_register("P00533", source="alphafold")  # AlphaFold

# Work with structures
proc = StructureProcessor()
df = proc.load_entity("1ubq")
print(df.columns)  # structure_id, atom_id, chain_id, residue_number, x, y, z, ...

# Create datasets
proc.create_dataset("my_structures", ["1ubq", "P00533"])
```

## Design Principles

1. **Zero Configuration** - Works out of the box with sensible defaults
2. **Human-Readable Names** - All public APIs use names, not UUIDs
3. **Automatic Registration** - Entities tracked automatically across formats
4. **Separation of Concerns** - Loaders acquire data; Processors manipulate data
5. **Lazy Initialization** - Resources created only when needed
