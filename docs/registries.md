# Registry System

The registry system provides centralized metadata management and discovery capabilities for Protos. It tracks all entities, datasets, and their relationships across the entire framework.

## Overview

Protos uses a hierarchical registry architecture:

```
GlobalRegistry
├── EntityRegistry      # Tracks all entities and formats
├── DatasetRegistry     # Manages dataset definitions
└── ProcessorRegistries # Format-specific registries
    ├── structure/registry.json
    ├── sequence/registry.json
    ├── grn/registry.json
    └── property/registry.json
```

## Registry Components

### 1. Entity Registry

The EntityRegistry is the core tracking system for all biological data:

```python
from protos.io.data_access import EntityRegistry

# Access through GlobalRegistry
registry = GlobalRegistry().entity_registry

# Direct operations
entity_id = registry.register_entity(
    original_id="1ubq",
    format_type="structure",
    file_path="structure/mmcif/1ubq.cif",
    metadata={"resolution": 1.8}
)
```

Key features:
- Maps biological names to entity IDs
- Tracks all formats per entity
- Stores format-specific metadata
- Enables cross-format queries

### 2. Dataset Registry

Manages collections of related entities:

```python
# Register dataset
registry.register_dataset(
    dataset_id="kinase_study",
    name="Human kinase structures",
    content=["1atp", "2src", "3erk"],
    metadata={
        "organism": "Homo sapiens",
        "purpose": "Drug discovery"
    }
)

# Query datasets
datasets = registry.list_datasets()
dataset_info = registry.get_dataset("kinase_study")
```

### 3. Processor Registries

Each processor maintains format-specific registries:

```python
# Structure registry tracks PDB-specific data
structure_registry = {
    "entities": {
        "3f8a9c2d1e": {
            "pdb_id": "1ubq",
            "chains": ["A"],
            "resolution": 1.8
        }
    },
    "datasets": {
        "test_structures": {
            "content": ["1ubq", "2gb1", "1crn"]
        }
    }
}
```

## Registry Operations

### Registration

Register new entities and datasets:

```python
# Register entity
entity_id = registry.register_entity(
    original_id="MY_PROTEIN",
    format_type="sequence",
    file_path="sequence/fasta/my_protein.fasta",
    metadata={
        "length": 350,
        "organism": "E. coli",
        "date_added": "2024-01-15"
    }
)

# Register dataset
registry.register_dataset(
    dataset_id="ecoli_proteins",
    name="E. coli protein collection",
    content=["PROT1", "PROT2", "PROT3"],
    description="Essential E. coli proteins"
)
```

### Discovery

Find and query registry contents:

```python
# List all entities
all_entities = registry.list_entities()

# Filter by format
structures = registry.list_entities_by_format("structure")
sequences = registry.list_entities_by_format("sequence")

# Get entity details
info = registry.get_entity("3f8a9c2d1e")
original_name = registry.get_original_id("3f8a9c2d1e")

# Query by metadata
high_res = registry.query_entities({
    "format": "structure",
    "resolution": {"lt": 2.0}
})
```

### Relationships

Track relationships between entities:

```python
# Register relationship
registry.add_relationship(
    source="WT_PROTEIN",
    target="MUT_L123A",
    relationship_type="mutation",
    metadata={"position": 123, "mutation": "L->A"}
)

# Query relationships
mutations = registry.get_related_entities(
    "WT_PROTEIN",
    relationship_type="mutation"
)

# Get lineage
lineage = registry.get_entity_lineage("DERIVED_PROTEIN")
```

## Registry Schema

### Entity Registry Schema

```json
{
    "entities": {
        "3f8a9c2d1e": {
            "original_ids": ["1ubq", "UBIQ_HUMAN"],
            "formats": {
                "structure": {
                    "path": "structure/mmcif/3f8a9c2d1e.cif",
                    "metadata": {
                        "resolution": 1.8,
                        "method": "X-RAY DIFFRACTION",
                        "chains": ["A"]
                    }
                },
                "sequence": {
                    "path": "sequence/fasta/3f8a9c2d1e.fasta",
                    "metadata": {
                        "length": 76,
                        "organism": "Homo sapiens"
                    }
                }
            },
            "relationships": {
                "parent": null,
                "children": ["1ubq_chainA"],
                "derived_from": null
            }
        }
    },
    "aliases": {
        "UBIQ_HUMAN": "3f8a9c2d1e",
        "Ubiquitin": "3f8a9c2d1e"
    }
}
```

### Dataset Registry Schema

```json
{
    "datasets": {
        "kinase_study": {
            "name": "Human kinase structures",
            "description": "Crystal structures for drug discovery",
            "content": ["1atp", "2src", "3erk"],
            "metadata": {
                "organism": "Homo sapiens",
                "created": "2024-01-15",
                "curator": "J. Smith",
                "tags": ["kinase", "drug-discovery", "human"]
            },
            "statistics": {
                "entity_count": 3,
                "format_distribution": {
                    "structure": 3,
                    "sequence": 3
                }
            }
        }
    }
}
```

## Advanced Registry Features

### 1. Metadata Queries

Complex queries using metadata:

```python
# Query with multiple conditions
results = registry.query_entities({
    "format": "structure",
    "resolution": {"lt": 2.5},
    "method": "X-RAY DIFFRACTION",
    "organism": {"in": ["Homo sapiens", "Mus musculus"]}
})

# Query with custom function
def has_ligand(metadata):
    return "ligands" in metadata and len(metadata["ligands"]) > 0

structures_with_ligands = registry.query_entities(
    {"format": "structure"},
    custom_filter=has_ligand
)
```

### 2. Bulk Operations

Efficient bulk registration:

```python
# Bulk register entities
entities_to_register = [
    {
        "original_id": "PROT1",
        "format_type": "sequence",
        "file_path": "seq1.fasta",
        "metadata": {"length": 100}
    },
    {
        "original_id": "PROT2",
        "format_type": "sequence",
        "file_path": "seq2.fasta",
        "metadata": {"length": 200}
    }
]

registered_ids = registry.bulk_register_entities(entities_to_register)
```

### 3. Registry Export/Import

Export and import registry data:

```python
# Export registry
registry.export_to_file("registry_backup.json")

# Import registry
registry.import_from_file("registry_backup.json")

# Merge registries
registry.merge_registry("external_registry.json")
```

### 4. Registry Validation

Validate registry integrity:

```python
# Check registry consistency
issues = registry.validate()
for issue in issues:
    print(f"Issue: {issue['type']} - {issue['description']}")

# Repair common issues
registry.repair_registry()

# Verify file references
missing_files = registry.verify_file_references()
```

## Registry Integration

### With Processors

Processors automatically interact with registries:

```python
class CustomProcessor(BaseProcessor):
    def save_entity(self, name: str, data: Any) -> str:
        # Save data
        file_path = self._save_data(data)
        
        # Auto-register in registry
        entity_id = self.registry.register_entity(
            original_id=name,
            format_type=self.processor_type,
            file_path=file_path,
            metadata=self._extract_metadata(data)
        )
        
        return entity_id
```

### With Entity System

Seamless entity tracking:

```python
# Entity operations update registry
sp.save_sequence("NEW_PROTEIN", "MKTAYIA...")
# Registry automatically updated

# Query finds entity
info = registry.get_entity_info("NEW_PROTEIN")
print(info["formats"])  # ["sequence"]
```

### With Datasets

Dataset management through registry:

```python
# Create dataset updates registry
cp.create_dataset(
    "pdb_2024",
    content=["7xyz", "8abc", "9def"]
)

# Registry tracks dataset
dataset = registry.get_dataset("pdb_2024")
print(f"Contains {len(dataset['content'])} structures")
```

## Registry Patterns

### Pattern 1: Cross-Format Discovery

```python
def find_entities_with_all_formats(formats):
    """Find entities that exist in all specified formats."""
    all_entities = registry.list_entities()
    matching = []
    
    for entity_id in all_entities:
        entity_info = registry.get_entity(entity_id)
        entity_formats = set(entity_info["formats"].keys())
        
        if set(formats).issubset(entity_formats):
            matching.append(entity_info["original_ids"][0])
    
    return matching

# Find proteins with structure, sequence, and embeddings
complete_proteins = find_entities_with_all_formats(
    ["structure", "sequence", "embedding"]
)
```

### Pattern 2: Metadata Aggregation

```python
def get_resolution_statistics():
    """Get statistics on structure resolutions."""
    structures = registry.query_entities({"format": "structure"})
    resolutions = []
    
    for entity_id in structures:
        info = registry.get_entity(entity_id)
        if "structure" in info["formats"]:
            res = info["formats"]["structure"]["metadata"].get("resolution")
            if res:
                resolutions.append(res)
    
    return {
        "count": len(resolutions),
        "mean": np.mean(resolutions),
        "median": np.median(resolutions),
        "best": min(resolutions),
        "worst": max(resolutions)
    }
```

### Pattern 3: Lineage Tracking

```python
def trace_mutation_lineage(protein_name):
    """Trace all mutations derived from a protein."""
    lineage_tree = {"name": protein_name, "mutations": []}
    
    # Get direct mutations
    mutations = registry.get_related_entities(
        protein_name,
        relationship_type="mutation"
    )
    
    # Recursively trace each mutation
    for mut in mutations:
        mut_info = registry.get_entity_metadata(mut)
        lineage_tree["mutations"].append({
            "name": mut,
            "mutation": mut_info.get("mutation"),
            "mutations": trace_mutation_lineage(mut)["mutations"]
        })
    
    return lineage_tree
```

## Registry Maintenance

### Cleanup Operations

```python
# Remove orphaned entries
orphans = registry.find_orphaned_entries()
for entity_id in orphans:
    registry.remove_entity(entity_id)

# Remove duplicate entries
duplicates = registry.find_duplicates()
registry.merge_duplicates(duplicates)

# Compact registry files
registry.compact()
```

### Backup Strategy

```python
# Regular backups
from datetime import datetime

backup_name = f"registry_backup_{datetime.now():%Y%m%d_%H%M%S}.json"
registry.export_to_file(backup_name)

# Incremental backups
changes = registry.get_changes_since(last_backup_time)
registry.export_changes(changes, "incremental_backup.json")
```

### Migration Support

```python
# Migrate from old format
def migrate_old_registry(old_registry_path):
    """Migrate from legacy registry format."""
    with open(old_registry_path) as f:
        old_data = json.load(f)
    
    for old_id, old_info in old_data.items():
        # Convert to new format
        registry.register_entity(
            original_id=old_info["name"],
            format_type=old_info["type"],
            file_path=old_info["path"],
            metadata=old_info.get("metadata", {})
        )
```

## Best Practices

### 1. Regular Validation

```python
# Schedule regular validation
def validate_registry_health():
    issues = registry.validate()
    if issues:
        logger.warning(f"Registry has {len(issues)} issues")
        registry.repair_registry()
```

### 2. Metadata Standards

```python
# Define metadata schema
STRUCTURE_METADATA_SCHEMA = {
    "resolution": float,
    "method": str,
    "chains": list,
    "ligands": list,
    "organism": str
}

# Validate before registration
def validate_metadata(metadata, schema):
    for key, expected_type in schema.items():
        if key in metadata:
            if not isinstance(metadata[key], expected_type):
                raise TypeError(f"{key} must be {expected_type}")
```

### 3. Efficient Queries

```python
# Cache frequent queries
@lru_cache(maxsize=128)
def get_high_quality_structures():
    return registry.query_entities({
        "format": "structure",
        "resolution": {"lt": 2.0},
        "method": "X-RAY DIFFRACTION"
    })
```

## Summary

The registry system provides:
- Centralized entity and dataset tracking
- Rich metadata management
- Relationship tracking
- Cross-format discovery
- Query capabilities
- Data integrity validation

By maintaining comprehensive registries, Protos enables sophisticated data management and discovery across all supported formats and processors.