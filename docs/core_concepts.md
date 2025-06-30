# Core Concepts

Protos is built on several fundamental abstractions that enable seamless management of heterogeneous biological data. Understanding these concepts is essential for effective use of the framework.

## The Protos Philosophy

**Biological names, not file paths.** Researchers work with familiar identifiers (PDB IDs, UniProt accessions, custom names) while Protos handles all file system operations internally. This abstraction provides:

- Platform independence (Windows/Linux/Mac)
- Centralized data management
- Automatic format handling
- Simplified workflows

## Key Abstractions

### 1. Entities

An **entity** represents a single biological data item that can exist in multiple formats. Each entity receives a unique, deterministic hash ID that remains consistent across all representations.

```python
# Example: Same protein in different formats
"1ubq" (PDB ID) → entity_id "3f8a9c2d1e"
  ├── structure: /structure/mmcif/3f8a9c2d1e.cif
  ├── sequence: /sequence/fasta/3f8a9c2d1e.fasta
  └── embedding: /embedding/3f8a9c2d1e.pkl
```

Key properties:
- Deterministic ID generation (same input → same hash)
- Multi-format support
- Metadata tracking
- Cross-format queries

### 2. Datasets

A **dataset** is a named collection of related entities. Datasets enable:
- Logical grouping of data
- Batch operations
- Reproducible analyses
- Metadata management

```python
# Creating a dataset
processor.create_dataset(
    dataset_id="kinase_study",
    name="Human kinase structures",
    description="Crystal structures of human kinases",
    content=["1atp", "2src", "3erk"]  # Entity names
)
```

### 3. Processors

**Processors** are specialized classes that handle specific data types. All processors inherit from `BaseProcessor` and share a common interface:

```python
# Common operations across all processors
processor.list_entities()           # List available entities
processor.load_entity(name)         # Load by name
processor.save_entity(name, data)   # Save with name
processor.create_dataset(...)       # Create dataset
processor.load_dataset(name)        # Load dataset
```

### 4. Registry System

The **registry** provides centralized tracking of all entities and datasets:

- **Entity Registry**: Maps names to hash IDs, tracks formats
- **Dataset Registry**: Manages dataset definitions
- **Global Registry**: Unified access to all registries

```python
# Registry automatically handles:
- Name → ID resolution
- Format tracking
- Metadata storage
- Cross-format queries
```

### 5. Path Management

**ProtosPaths** handles all file system operations:

```python
# Users never construct paths
# ❌ WRONG
path = Path("data/structure/mmcif/1ubq.cif")

# ✅ CORRECT  
structure = processor.load_structure("1ubq")
```

## Data Organization

Protos organizes data in a hierarchical structure:

```
protos_data/
├── structure/          # 3D structures (CifBaseProcessor)
│   ├── mmcif/         # PDB/CIF files
│   └── alignments/    # Structure alignments
├── sequence/          # Sequences (SeqProcessor)
│   ├── fasta/         # FASTA files
│   └── alignments/    # Sequence alignments
├── grn/               # GRN annotations (GRNBaseProcessor)
│   ├── tables/        # GRN tables
│   └── configs/       # GRN configurations
├── property/          # Properties (PropertyProcessor)
└── embedding/         # Embeddings (EmbeddingProcessor)
```

## Common Patterns

### Loading Data
```python
# By biological identifier
structure = cp.load_structure("1ubq")      # PDB ID
sequence = sp.load_sequence("P12345")      # UniProt ID
grn_table = gp.load_grn_table("rhodopsins") # Table name

# By entity ID (automatic resolution)
data = processor.load_entity("1ubq")       # Resolves to hash ID
```

### Saving Data
```python
# Save with biological name
sp.save_sequences({"MY_PROTEIN": "MKTAY..."}, "my_protein.fasta")

# Automatic entity registration
# Entity ID generated from "MY_PROTEIN"
# Registered in entity registry
# Available for cross-format operations
```

### Cross-Format Operations
```python
# Extract sequences from structures
sequences = struct_processor.get_seq_dict()

# Use in another processor
seq_processor.save_sequences(sequences, "extracted.fasta")

# Generate embeddings
embeddings = emb_processor.embed_sequences(sequences)
```

### Dataset Operations
```python
# Create dataset
processor.create_dataset(
    dataset_id="my_study",
    name="My protein study",
    content=["1ubq", "2gb1", "1crn"]
)

# Load dataset
processor.load_dataset("my_study")

# Iterate over dataset
for entity in processor.iter_dataset("my_study"):
    data = processor.load_entity(entity)
    # Process data...
```

## Best Practices

### 1. Use Biological Names
Always work with meaningful biological identifiers:
```python
# Good
processor.load_structure("1ubq")
processor.load_sequence("EGFR_HUMAN")

# Avoid
processor.load_entity("3f8a9c2d1e")  # Hash ID
```

### 2. Leverage Datasets
Group related data for batch processing:
```python
# Create logical groupings
processor.create_dataset("kinases", content=kinase_list)
processor.create_dataset("phosphatases", content=phosphatase_list)
```

### 3. Enable Cross-Format Tracking
Register entities consistently across formats:
```python
# Same name across formats enables tracking
struct_proc.save_structure("MY_PROTEIN", structure_data)
seq_proc.save_sequence("MY_PROTEIN", sequence_data)
# Both automatically linked via entity system
```

### 4. Use Type Hints
Processors enforce data types for reliability:
```python
# Load with type enforcement
structure = cp.load_structure("1ubq", apply_dtypes=True)
```

### 5. Handle Missing Data
Check for empty results:
```python
structure = cp.load_structure("1ubq")
if structure is not None and not structure.empty:
    # Process structure
    ca_atoms = structure[structure['atom_name'] == 'CA']
```

## Advanced Concepts

### Entity Metadata
Each entity can store format-specific metadata:
```python
{
    "entity_id": "3f8a9c2d1e",
    "formats": {
        "structure": {
            "resolution": 1.8,
            "method": "X-RAY DIFFRACTION",
            "chains": ["A", "B"]
        },
        "sequence": {
            "length": 150,
            "organism": "Homo sapiens"
        }
    }
}
```

### Processor Inheritance
Create custom processors by inheriting from base classes:
```python
class CustomProcessor(BaseProcessor):
    def __init__(self, name="custom"):
        super().__init__(name, processor_type="custom")
    
    def process_entity(self, entity_name):
        # Custom processing logic
        pass
```

### Registry Queries
Advanced queries across the registry:
```python
# Find all entities with structures
struct_entities = registry.get_entities_by_format("structure")

# Find entities present in multiple formats
multi_format = registry.get_multi_format_entities()

# Query by metadata
high_res = registry.query_entities({
    "format": "structure",
    "resolution": {"lt": 2.0}
})
```

## Summary

Protos provides a powerful abstraction layer that:
1. Eliminates file path manipulation
2. Enables cross-format data integration
3. Maintains data relationships automatically
4. Provides a consistent interface across data types
5. Supports complex research workflows

The combination of entities, datasets, processors, and the registry system creates a flexible framework for managing the complexity of structural biology data.