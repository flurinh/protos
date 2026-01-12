# Entity, Registry and Datasets

Protos provides centralized entity management through the EntityRegistry and DatasetManager. These components track all biological entities across multiple formats and organize them into datasets for batch operations.

## Design Principles

1. **Human-Readable Names** - All public APIs use names, never UUIDs (UUIDs are internal only)
2. **Multi-Format Support** - Same entity can exist in multiple formats (structure, sequence, etc.)
3. **Relationship Tracking** - Entities can have typed relationships to other entities
4. **Singleton Pattern** - One registry instance shared across all components

---

## EntityRegistry

The `EntityRegistry` is the central tracking system for all biological entities.

**Location:** `protos.io.core.entity_registry`

### Getting the Registry

```python
from protos.io.core import get_registry

registry = get_registry()  # Returns singleton instance
```

### Data Structures

#### EntityInfo
Information about an entity in a specific format:

```python
@dataclass
class EntityInfo:
    entity_id: str       # UUID (internal use only)
    original_id: str     # Human-readable name
    format_type: str     # e.g., 'structure', 'sequence'
    file_path: str       # Path to data file
    metadata: Dict       # Additional metadata
    created: str         # ISO timestamp
    modified: str        # ISO timestamp
```

#### Entity
Complete entity with all formats:

```python
@dataclass
class Entity:
    entity_id: str                    # UUID (internal)
    original_id: str                  # Human-readable name
    aliases: List[str]                # Alternative names
    formats: Dict[str, Dict]          # format_type -> format data
    relationships: List[Dict]         # Relationships to other entities
    created: str
    modified: str
```

### Basic Operations

#### Register an Entity

```python
# Register entity with format and file path
name = registry.register_entity(
    name="ubiquitin",
    format_type="structure",
    file_path="/path/to/ubiquitin.cif",
    metadata={"source": "pdb", "pdb_id": "1ubq"}
)
```

#### Find an Entity

```python
# Find entity by name
info = registry.find_entity("ubiquitin")
print(info.file_path)
print(info.format_type)
print(info.metadata)

# Find specific format
info = registry.find_entity("ubiquitin", format_type="sequence")
```

#### List Entities

```python
# List all entities
all_entities = registry.list_entities()

# List entities with specific format
structures = registry.list_entities(format_type="structure")
sequences = registry.list_entities(format_type="sequence")
```

#### Check Entity Existence

```python
# Check if entity exists
exists = registry.entity_exists("ubiquitin")

# Check specific format
has_structure = registry.entity_exists("ubiquitin", format_type="structure")
```

### Multi-Format Entities

A single entity can have multiple formats (e.g., both structure and sequence):

```python
# Register structure format
registry.register_entity("protein_a", "structure", "/path/to/structure.cif")

# Register sequence format for same entity
registry.register_entity("protein_a", "sequence", "/path/to/sequence.fasta")

# Get all formats
formats = registry.get_entity_formats("protein_a")
# ['structure', 'sequence']
```

### Aliases

Entities can have multiple names (aliases):

```python
# Add alias
registry.add_alias("ubiquitin", "UBI")
registry.add_alias("ubiquitin", "P0CG48")

# Find by any name
info = registry.find_entity("UBI")  # Same as find_entity("ubiquitin")
```

### Metadata Operations

```python
# Get metadata
metadata = registry.get_entity_metadata("ubiquitin", "structure")

# Update metadata (merges with existing)
registry.update_metadata("ubiquitin", "structure", {
    "resolution": 1.8,
    "method": "X-ray"
})
```

### Entity Renaming

```python
# Rename entity (preserves ID and relationships)
registry.rename_entity("old_name", "new_name")
```

### Remove Formats

```python
# Remove a format from entity
registry.remove_format("protein_a", "structure")

# If no formats remain, entity is deleted entirely
```

---

## Relationships

Entities can have directed relationships to other entities.

### Relationship Types

| Type | Inverse | Symmetric | Description |
|------|---------|-----------|-------------|
| `derived_from` | `derives_to` | No | Entity derived from another |
| `subset_of` | `contains` | No | Entity is subset of another |
| `merged_from` | `merged_into` | No | Entity merged from others |
| `version_of` | `has_version` | No | Entity is version of another |
| `aligned_to` | `aligned_to` | Yes | Entities are aligned |
| `annotated_by` | `annotates` | No | Entity annotated by another |

### Creating Relationships

```python
# Add relationship
registry.add_relationship(
    source_name="variant_1",
    target_name="wildtype",
    rel_type="derived_from",
    metadata={"mutation": "K48R"}
)
```

### Querying Relationships

```python
# Get all relationships
relationships = registry.get_relationships("wildtype")

# Filter by type
derived = registry.get_relationships("wildtype", rel_type="derived_from")

# Filter by direction
incoming = registry.get_relationships("wildtype", direction="incoming")
outgoing = registry.get_relationships("wildtype", direction="outgoing")

# Get related entity names
related = registry.get_related_entities("wildtype", rel_type="derived_from")
# ['variant_1', 'variant_2']
```

### Relationship Response Format

```python
{
    'type': 'derived_from',
    'source_name': 'variant_1',
    'target_name': 'wildtype',
    'direction': 'incoming',
    'metadata': {'mutation': 'K48R'},
    'created': '2024-01-15T10:30:00'
}
```

### Remove Relationships

```python
registry.remove_relationship(
    source_name="variant_1",
    target_name="wildtype",
    rel_type="derived_from"
)
```

---

## DatasetManager

The `DatasetManager` organizes entities into named collections for batch operations.

**Location:** `protos.io.core.dataset_manager`

### Getting a DatasetManager

DatasetManager is processor-specific and typically accessed through a processor:

```python
from protos import StructureProcessor

proc = StructureProcessor()
dm = proc.dataset_manager

# Or create directly
from protos.io.core.dataset_manager import DatasetManager
dm = DatasetManager("structure")
```

### Creating Datasets

```python
# Create dataset with entities
dm.create_dataset(
    name="gpcr_structures",
    entities=["3sn6", "4lde", "5c1m"],
    metadata={
        "description": "GPCR crystal structures",
        "family": "class_a"
    }
)
```

### Loading Datasets

```python
# Load dataset configuration
dataset = dm.load_dataset("gpcr_structures")
print(dataset['entities'])
print(dataset['metadata'])

# Get entity list only
entities = dm.get_dataset_entities("gpcr_structures")
# ['3sn6', '4lde', '5c1m']
```

### Listing Datasets

```python
# List all datasets for this processor
datasets = dm.list_datasets()

# Check if dataset exists
exists = dm.dataset_exists("gpcr_structures")
```

### Modifying Datasets

```python
# Add entities
dm.add_to_dataset("gpcr_structures", ["6cmo", "7dh8"])

# Remove entities
dm.remove_from_dataset("gpcr_structures", ["4lde"])

# Update metadata
dm.update_metadata("gpcr_structures", {"updated_by": "user"})
```

### Dataset Operations

```python
# Copy dataset
dm.copy_dataset("gpcr_structures", "gpcr_structures_backup")

# Merge datasets
dm.merge_datasets(
    ["dataset_a", "dataset_b", "dataset_c"],
    target_name="merged_dataset"
)

# Delete dataset
dm.delete_dataset("old_dataset")
```

### Dataset Information

```python
# Get detailed info
info = dm.get_dataset_info("gpcr_structures")

print(info['name'])
print(info['entity_count'])
print(info['created'])
print(info['modified'])

# Entity details include format information
for entity in info['entities']:
    print(f"{entity['name']}: {entity['formats']}")
    if entity.get('missing'):
        print("  (entity has been deleted)")
```

### Handling Entity Renames

When entities are renamed, datasets automatically track via stable IDs:

```python
# Entity names update automatically
registry.rename_entity("old_name", "new_name")

# Force refresh if needed
dm.refresh_dataset_entities("my_dataset")

# Refresh all datasets
dm.refresh_all_datasets()
```

---

## Dataset Storage Format

Datasets are stored as JSON files in `{processor_path}/datasets/{name}.json`:

```json
{
  "name": "gpcr_structures",
  "processor_type": "structure",
  "entities": ["3sn6", "4lde", "5c1m"],
  "entity_ids": ["uuid1", "uuid2", "uuid3"],
  "metadata": {
    "description": "GPCR crystal structures"
  },
  "created": "2024-01-15T10:00:00",
  "modified": "2024-01-15T11:30:00"
}
```

---

## Registry Storage Format

The registry is stored as JSON at `{data_root}/registry.json`:

```json
{
  "entities": {
    "uuid-1234": {
      "original_id": "ubiquitin",
      "aliases": ["UBI", "P0CG48"],
      "formats": {
        "structure": {
          "file_path": "structure/pkl/ubiquitin.pkl",
          "metadata": {"source": "pdb", "pdb_id": "1ubq"},
          "created": "2024-01-15T10:00:00"
        },
        "sequence": {
          "file_path": "sequence/fasta/ubiquitin.fasta",
          "metadata": {"source": "uniprot"},
          "created": "2024-01-15T10:05:00"
        }
      },
      "relationships": [
        {
          "type": "derived_from",
          "source": "uuid-5678",
          "target": "uuid-1234",
          "metadata": {},
          "created": "2024-01-15T10:10:00"
        }
      ],
      "created": "2024-01-15T10:00:00",
      "modified": "2024-01-15T10:10:00"
    }
  },
  "name_index": {
    "ubiquitin": "uuid-1234",
    "UBI": "uuid-1234",
    "P0CG48": "uuid-1234"
  }
}
```

---

## Complete Example

```python
from protos import StructureLoader, SequenceLoader, StructureProcessor
from protos.io.core import get_registry

# Initialize components
registry = get_registry()
struct_loader = StructureLoader()
seq_loader = SequenceLoader()
proc = StructureProcessor()

# Download and register multiple structures
struct_loader.download_batch(
    ["1ubq", "2w9s", "3sn6"],
    dataset_name="test_structures"
)

# Add sequence data for one entity
seq_loader.download_and_register("P0CG48", name="1ubq")

# Check formats
formats = registry.get_entity_formats("1ubq")
print(f"1ubq has formats: {formats}")
# ['structure', 'sequence']

# Add alias
registry.add_alias("1ubq", "ubiquitin")

# Create relationship
registry.add_relationship(
    source_name="2w9s",
    target_name="1ubq",
    rel_type="aligned_to"
)

# Work with datasets
dm = proc.dataset_manager
entities = dm.get_dataset_entities("test_structures")

# Iterate over dataset
for entity_name in entities:
    df = proc.load_entity(entity_name)
    print(f"{entity_name}: {len(df)} atoms")
```

---

## Thread Safety and Multi-Process Access

The registry uses atomic writes with retry logic for safe concurrent access:

```python
# Refresh to see updates from other processes
registry.refresh()

# Reset registry (creates backup)
backup_path = registry.reset(backup=True)
```

The registry file uses atomic rename operations to prevent corruption during concurrent writes.
