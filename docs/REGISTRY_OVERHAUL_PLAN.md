# Registry System Overhaul Plan

## Overview

This document outlines the plan to overhaul Protos' registry system to support proper entity/dataset management with standardized hash-based entity IDs.

## Core Concepts

### Entity IDs as Hashes

All entities in Protos will use standardized hash-based IDs:
- **Hash length**: 10 characters (provides ~1 trillion unique IDs)
- **Hash algorithm**: SHA-256 truncated to 10 chars
- **Benefits**:
  - Consistent filenames across all data types
  - No special characters or path issues
  - Deterministic ID generation from content
  - Easy to reference and share

### Entity Types by Processor

1. **Structure Entities** (CifBaseProcessor)
   - Entity = Single PDB/CIF file
   - Hash input: PDB ID + chain (if specific) or file content
   - Example: `1ABC` → `a3f2d8c91b`
   - Stored as: `structure/entities/a3f2d8c91b.cif`

2. **Sequence Entities** (SeqProcessor)
   - Entity = Single sequence
   - Hash input: Sequence content + identifier
   - Example: `>P12345\nMLKLLI...` → `b7e4a2f3c8`
   - Stored as: `sequence/entities/b7e4a2f3c8.fasta`

3. **GRN Table Entities** (GRNBaseProcessor)
   - Entity = Collection of GRN-assigned sequences
   - Each row in table gets its own entity ID
   - Hash input: Sequence ID + GRN assignments
   - Table stored as: `grn/entities/c9d5e1a7f2.csv`
   - Table contains: `entity_id | sequence_id | 1.50 | 2.50 | ...`

4. **Embedding Entities** (EmbeddingProcessor)
   - Entity = Single protein embedding
   - Hash input: Protein ID + model name + sequence
   - Example: `P12345_esm2_t33` → `d8f3b2c1a9`
   - Stored as: `embedding/entities/d8f3b2c1a9.pkl`

## Data Structure

### Entity Registry Structure

```json
{
  "entities": {
    "a3f2d8c91b": {
      "type": "structure",
      "original_id": "1ABC",
      "metadata": {
        "pdb_id": "1ABC",
        "chains": ["A", "B"],
        "resolution": 2.5,
        "method": "X-RAY"
      },
      "datasets": ["experimental_structures", "test_set"],
      "created": "2024-01-15T10:30:00Z",
      "modified": "2024-01-15T10:30:00Z"
    },
    "b7e4a2f3c8": {
      "type": "sequence",
      "original_id": "P12345",
      "metadata": {
        "sequence_id": "P12345",
        "organism": "Homo sapiens",
        "length": 350
      },
      "datasets": ["human_proteins", "kinases"],
      "created": "2024-01-15T10:31:00Z",
      "modified": "2024-01-15T10:31:00Z"
    }
  },
  "datasets": {
    "experimental_structures": {
      "name": "Experimental Structures",
      "description": "X-ray and NMR structures",
      "type": "structure",
      "entities": ["a3f2d8c91b", "e5f6g7h8i9"],
      "metadata": {
        "source": "PDB",
        "filters": "resolution < 3.0"
      }
    }
  }
}
```

### GRN Table Format

For GRN tables, each sequence becomes an entity:

```csv
entity_id,sequence_id,1.50,2.50,3.50,4.50,5.50,6.50,7.50
f1a2b3c4d5,BR1_HUMAN,L45,V87,I123,W165,F203,W241,K279
g2b3c4d5e6,BR2_MOUSE,L46,V88,I124,W166,F204,W242,K280
h3c4d5e6f7,BR3_BOVIN,L47,V89,I125,W167,F205,W243,K281
```

## Implementation Details

### Hash Generation Function

```python
import hashlib

def generate_entity_id(content: str, prefix: str = "") -> str:
    """
    Generate a standardized 10-character entity ID from content.
    
    Args:
        content: String to hash (could be sequence, PDB ID, etc.)
        prefix: Optional prefix to make hash unique across types
        
    Returns:
        10-character hash ID
    """
    # Combine prefix and content for uniqueness
    to_hash = f"{prefix}:{content}" if prefix else content
    
    # Generate SHA-256 hash
    hash_obj = hashlib.sha256(to_hash.encode('utf-8'))
    
    # Take first 10 characters of hex digest
    entity_id = hash_obj.hexdigest()[:10]
    
    return entity_id
```

### Entity Storage Pattern

```
protos_data/
├── structure/
│   ├── entities/
│   │   ├── a3f2d8c91b.cif
│   │   ├── e5f6g7h8i9.cif
│   │   └── metadata/
│   │       ├── a3f2d8c91b.json
│   │       └── e5f6g7h8i9.json
│   └── registry.json
├── sequence/
│   ├── entities/
│   │   ├── b7e4a2f3c8.fasta
│   │   └── metadata/
│   │       └── b7e4a2f3c8.json
│   └── registry.json
├── grn/
│   ├── entities/
│   │   ├── c9d5e1a7f2.csv  # Full GRN table
│   │   └── metadata/
│   │       └── c9d5e1a7f2.json
│   └── registry.json
└── global_registry.json
```

## API Design

### BaseProcessor Entity Methods

```python
class BaseProcessor:
    def list_entities(self, dataset: Optional[str] = None) -> List[str]:
        """List all entity IDs, optionally filtered by dataset."""
        
    def load_entity(self, entity_id: str) -> Any:
        """Load a single entity by its hash ID."""
        
    def save_entity(self, data: Any, original_id: str = None, metadata: Dict = None) -> str:
        """Save entity and return its hash ID."""
        
    def entity_exists(self, entity_id: str) -> bool:
        """Check if an entity exists."""
        
    def get_entity_metadata(self, entity_id: str) -> Dict:
        """Get metadata for an entity."""
        
    def find_entity_by_original_id(self, original_id: str) -> Optional[str]:
        """Find entity hash ID by original ID (e.g., PDB ID)."""
```

### GRNBaseProcessor Specific Methods

```python
class GRNBaseProcessor:
    def load_grn_entity(self, entity_id: str) -> pd.Series:
        """Load a single GRN assignment by entity ID."""
        
    def save_grn_entities(self, grn_table: pd.DataFrame, table_name: str) -> List[str]:
        """Save GRN table and return list of entity IDs."""
        
    def get_grn_table_entities(self, table_name: str) -> List[str]:
        """Get all entity IDs in a GRN table."""
```

## Migration Strategy

### Phase 1: Parallel System
1. Implement new entity registry alongside existing system
2. Add entity ID generation to save operations
3. Maintain backward compatibility with original IDs

### Phase 2: Migration Tools
```python
def migrate_processor_data(processor_type: str):
    """Migrate existing data to entity-based system."""
    # 1. Scan existing files
    # 2. Generate entity IDs
    # 3. Create entity registry entries
    # 4. Copy/link files to new structure
    # 5. Update datasets to reference entity IDs
```

### Phase 3: Deprecation
1. Update all code to use entity IDs
2. Provide lookup by original ID
3. Deprecate direct file access
4. Remove old registry system

## Example Usage

### Loading Structures
```python
# Old way
processor.load_structure("1ABC")

# New way - by entity ID
processor.load_entity("a3f2d8c91b")

# Convenience method - lookup by original ID
entity_id = processor.find_entity_by_original_id("1ABC")
processor.load_entity(entity_id)
```

### Working with GRN Tables
```python
# Save GRN table with entity IDs
grn_processor = GRNBaseProcessor()
entity_ids = grn_processor.save_grn_entities(grn_df, "mo_ref")

# Load specific GRN assignment
grn_data = grn_processor.load_grn_entity("f1a2b3c4d5")
print(f"Sequence {grn_data['sequence_id']} has 3.50 = {grn_data['3.50']}")

# List all entities in a GRN table
entities = grn_processor.get_grn_table_entities("mo_ref")
```

### Dataset Operations
```python
# Create dataset from entities
dataset_manager = DatasetManager("structure")
dataset_manager.create_dataset(
    "experimental_structures",
    entity_ids=["a3f2d8c91b", "e5f6g7h8i9"],
    metadata={"resolution_cutoff": 3.0}
)

# List entities in dataset
entities = processor.list_entities(dataset="experimental_structures")
```

## Benefits

1. **Consistency**: All entities use same ID format
2. **Portability**: Hash IDs are platform-independent
3. **Deduplication**: Same content always gets same ID
4. **Traceability**: Can track entity across datasets
5. **Scalability**: No naming conflicts or special characters
6. **Performance**: Direct lookup by hash is fast

## Testing Strategy

1. **Unit Tests**:
   - Hash generation consistency
   - Entity CRUD operations
   - Migration tools

2. **Integration Tests**:
   - Cross-processor entity references
   - Dataset-entity relationships
   - Backward compatibility

3. **Performance Tests**:
   - Entity lookup speed
   - Large dataset handling
   - Migration performance

## Timeline

- **Week 1**: Implement EntityRegistry and hash generation
- **Week 2**: Update BaseProcessor with entity methods
- **Week 3**: Implement processor-specific entity handling
- **Week 4**: Create migration tools and tests
- **Week 5**: Documentation and examples
- **Week 6**: Deprecation planning and cleanup