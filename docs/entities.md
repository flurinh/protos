# Entity System

The entity system is the cornerstone of Protos' data management architecture. It provides universal tracking of biological data across multiple formats, enabling seamless integration and cross-format operations.

## Overview

An **entity** represents a single biological item (protein, sequence, structure) that can exist in multiple data formats. Each entity receives a unique, deterministic identifier that remains constant across all representations.

## Core Concepts

### Entity Identity

Every entity has:
- **Original Name**: The biological identifier (e.g., "1ubq", "P12345", "EGFR_HUMAN")
- **Entity ID**: A deterministic hash generated from the name
- **Format Entries**: References to the entity's data in different formats

```python
# Example entity structure
{
    "entity_id": "3f8a9c2d1e",
    "original_ids": ["1ubq", "UBIQ_HUMAN"],  # Can have aliases
    "formats": {
        "structure": {
            "path": "structure/mmcif/3f8a9c2d1e.cif",
            "metadata": {"resolution": 1.8, "method": "X-RAY"}
        },
        "sequence": {
            "path": "sequence/fasta/3f8a9c2d1e.fasta",
            "metadata": {"length": 76, "organism": "Homo sapiens"}
        },
        "embedding": {
            "path": "embedding/3f8a9c2d1e.pkl",
            "metadata": {"model": "esm2", "dimension": 1280}
        }
    }
}
```

### Deterministic ID Generation

Entity IDs are generated using SHA-256 hashing:

```python
from protos.io.data_access import generate_entity_id

# Same input always produces same ID
entity_id = generate_entity_id("1ubq")  # → "3f8a9c2d1e"
entity_id = generate_entity_id("1ubq")  # → "3f8a9c2d1e" (identical)

# Different inputs produce different IDs
id1 = generate_entity_id("1ubq")  # → "3f8a9c2d1e"
id2 = generate_entity_id("2gb1")  # → "7b5c2a9f4d"
```

Properties:
- Deterministic: Same name → same ID
- Unique: Different names → different IDs
- Short: First 10 characters of SHA-256 hash
- Stable: IDs never change

## Entity Registry

The EntityRegistry manages all entity tracking:

### Registry Operations

```python
from protos.io.data_access import GlobalRegistry

registry = GlobalRegistry()

# Register new entity
entity_id = registry.register_entity(
    original_id="MY_PROTEIN",
    format_type="sequence",
    file_path="sequence/fasta/abc123.fasta",
    metadata={"length": 250, "organism": "E. coli"}
)

# Resolve name to ID
entity_id = registry.resolve_identifier("MY_PROTEIN")

# Get entity information
entity_info = registry.get_entity(entity_id)

# List entities by format
structure_entities = registry.list_entities_by_format("structure")

# Query entities
high_res = registry.query_entities({
    "format": "structure",
    "resolution": {"lt": 2.0}
})
```

### Multi-Format Tracking

Entities can exist in multiple formats simultaneously:

```python
# Same entity, multiple formats
entity_name = "RHODOPSIN"

# Add structure format
struct_proc.save_structure(entity_name, structure_data)

# Add sequence format  
seq_proc.save_sequence(entity_name, "MNGTEGPNFYVPFS...")

# Add embedding format
emb_proc.save_embedding(entity_name, embedding_array)

# Registry tracks all formats
info = registry.get_entity_info(entity_name)
print(info['formats'])  # ['structure', 'sequence', 'embedding']
```

## Entity Operations

### Creating Entities

Entities are created automatically when saving data:

```python
# Saving automatically creates/updates entity
processor.save_entity("NEW_PROTEIN", data)

# Behind the scenes:
# 1. Generate entity ID from name
# 2. Save data to format-specific location
# 3. Register in entity registry
# 4. Update format metadata
```

### Loading Entities

Load entities by name, with automatic ID resolution:

```python
# Load by biological name
structure = cp.load_structure("1ubq")
sequence = sp.load_sequence("P12345")

# Load by entity ID (less common)
data = processor.load_entity("3f8a9c2d1e")

# Check if entity exists
if processor.has_entity("MY_PROTEIN"):
    data = processor.load_entity("MY_PROTEIN")
```

### Listing Entities

List operations return human-readable names:

```python
# List all entities (returns names, not IDs)
entities = processor.list_entities()
# ["1ubq", "2gb1", "P12345", "MY_PROTEIN"]

# List with metadata
detailed = processor.list_entities(include_metadata=True)
# [
#   {"name": "1ubq", "id": "3f8a9c2d1e", "formats": ["structure", "sequence"]},
#   {"name": "2gb1", "id": "7b5c2a9f4d", "formats": ["structure"]}
# ]
```

### Cross-Format Access

Access the same entity across different processors:

```python
# Structure processor
structure = struct_proc.load_entity("RHODOPSIN")

# Sequence processor - same entity
sequence = seq_proc.load_entity("RHODOPSIN")

# Property processor - same entity
properties = prop_proc.get_entity_properties("RHODOPSIN")

# All accessing the same biological entity
```

## Advanced Features

### Entity Aliases

Entities can have multiple names:

```python
# Register aliases
registry.add_alias("1ubq", "UBIQ_HUMAN")
registry.add_alias("1ubq", "Ubiquitin")

# All names resolve to same entity
id1 = registry.resolve_identifier("1ubq")        # → "3f8a9c2d1e"
id2 = registry.resolve_identifier("UBIQ_HUMAN")  # → "3f8a9c2d1e"
id3 = registry.resolve_identifier("Ubiquitin")   # → "3f8a9c2d1e"
```

### Entity Metadata

Rich metadata tracking per format:

```python
# Structure metadata
{
    "resolution": 1.8,
    "method": "X-RAY DIFFRACTION",
    "chains": ["A"],
    "ligands": ["ZN", "SO4"],
    "deposition_date": "1995-07-10"
}

# Sequence metadata  
{
    "length": 76,
    "organism": "Homo sapiens",
    "uniprot_id": "P0CG48",
    "sequence_version": 2
}

# Custom metadata
processor.update_entity_metadata("MY_PROTEIN", {
    "experiment_id": "EXP001",
    "conditions": {"pH": 7.4, "temp": 298},
    "notes": "Crystallized with new method"
})
```

### Entity Relationships

Track relationships between entities:

```python
# Parent-child relationships
registry.add_relationship(
    parent="1ubq",
    child="1ubq_A",
    relationship_type="chain_extraction"
)

# Mutation relationships
registry.add_relationship(
    parent="WT_PROTEIN",
    child="MUT_L234A",
    relationship_type="mutation",
    metadata={"mutation": "L234A"}
)

# Query relationships
children = registry.get_children("1ubq")
mutations = registry.get_related("WT_PROTEIN", type="mutation")
```

### Entity Versioning

Track entity versions over time:

```python
# Save new version
processor.save_entity(
    "MY_PROTEIN",
    updated_data,
    version="2.0",
    version_notes="Improved refinement"
)

# Access specific version
data_v1 = processor.load_entity("MY_PROTEIN", version="1.0")
data_v2 = processor.load_entity("MY_PROTEIN", version="2.0")

# List versions
versions = registry.get_entity_versions("MY_PROTEIN")
```

## Entity Lifecycle

### 1. Creation
```python
# Entity created on first save
sp.save_sequence("NEW_PROT", "MKTAYIA...")
# → Entity ID generated
# → File saved
# → Registry updated
```

### 2. Multi-Format Addition
```python
# Add more formats to existing entity
struct_proc.save_structure("NEW_PROT", struct_data)
emb_proc.save_embedding("NEW_PROT", embedding)
# → Same entity ID used
# → New formats registered
```

### 3. Discovery
```python
# Find entity across system
formats = registry.get_entity_formats("NEW_PROT")
# → ["sequence", "structure", "embedding"]

# Load from any processor
seq = sp.load_entity("NEW_PROT")
struct = cp.load_entity("NEW_PROT")
```

### 4. Deletion
```python
# Delete from specific format
sp.delete_entity("NEW_PROT")
# → Removes sequence format only

# Delete entirely
registry.delete_entity("NEW_PROT", all_formats=True)
# → Removes from all formats
```

## Best Practices

### 1. Use Meaningful Names

```python
# ✅ GOOD - Descriptive names
processor.save_entity("EGFR_HUMAN_WT", data)
processor.save_entity("kinase_inhibitor_complex", data)

# ❌ BAD - Cryptic names
processor.save_entity("test1", data)
processor.save_entity("x", data)
```

### 2. Consistent Naming Across Formats

```python
# Use same name for same entity
struct_proc.save_structure("RHO_HUMAN", structure)
seq_proc.save_sequence("RHO_HUMAN", sequence)
# Automatically linked through entity system
```

### 3. Add Rich Metadata

```python
# Include useful metadata
processor.save_entity(
    "KINASE_STUDY_001",
    data,
    metadata={
        "date": "2024-01-15",
        "experimenter": "J. Smith",
        "conditions": {
            "buffer": "PBS",
            "pH": 7.4,
            "additives": ["ATP", "Mg2+"]
        },
        "purpose": "ATP binding study"
    }
)
```

### 4. Leverage Cross-Format Capabilities

```python
# Complete analysis workflow
# 1. Load structure
structure = cp.load_entity("1ATP")

# 2. Extract sequence
sequence = cp.extract_sequence(structure)

# 3. Find similar sequences
similar = sp.find_similar(sequence)

# 4. Load their structures
for seq_id in similar:
    if cp.has_entity(seq_id):
        similar_struct = cp.load_entity(seq_id)
        # Analyze...
```

## Entity System Patterns

### Pattern 1: Batch Registration

```python
# Register multiple entities efficiently
entities_to_register = {
    "PROT1": "path/to/prot1.fasta",
    "PROT2": "path/to/prot2.fasta",
    "PROT3": "path/to/prot3.fasta"
}

for name, path in entities_to_register.items():
    with open(path) as f:
        sequence = f.read()
    sp.save_sequence(name, sequence)
```

### Pattern 2: Entity Migration

```python
# Migrate from old naming to entity system
old_files = Path("old_data").glob("*.fasta")

for old_file in old_files:
    # Extract meaningful name
    name = old_file.stem.upper()
    
    # Load data
    with open(old_file) as f:
        sequence = f.read()
    
    # Save with entity system
    sp.save_sequence(name, sequence)
    
    print(f"Migrated {old_file.name} → Entity: {name}")
```

### Pattern 3: Cross-Format Validation

```python
# Validate sequence matches structure
struct_sequence = cp.extract_sequence("1ubq")
fasta_sequence = sp.load_sequence("1ubq")

if struct_sequence != fasta_sequence:
    print("Warning: Sequence mismatch between formats")
    # Update to ensure consistency
    sp.save_sequence("1ubq", struct_sequence)
```

## Troubleshooting

### Common Issues

1. **Entity not found**
```python
# Check if entity exists
if not processor.has_entity("MY_PROT"):
    print("Entity does not exist in this format")
    
# Check all formats
formats = registry.get_entity_formats("MY_PROT")
print(f"Entity exists in: {formats}")
```

2. **Duplicate entities**
```python
# Entities are unique by name
# Saving with same name updates existing
sp.save_sequence("PROT1", new_sequence)  # Updates
```

3. **Name conflicts**
```python
# Use namespacing for projects
processor.save_entity("PROJECT1_KINASE1", data)
processor.save_entity("PROJECT2_KINASE1", data)
```

## Summary

The entity system provides:
- Universal tracking across formats
- Deterministic ID generation
- Automatic registration and discovery
- Rich metadata support
- Cross-format integration
- Relationship tracking

By abstracting data management through entities, Protos enables complex multi-format workflows while maintaining data integrity and traceability.