# Protos Data Management V3: Enhanced Dual Identifier System

## Table of Contents
1. [Overview](#overview)
2. [Core Improvements](#core-improvements)
3. [Entity System](#entity-system)
4. [Dataset System](#dataset-system)
5. [Registry Architecture](#registry-architecture)
6. [Relationship System](#relationship-system)
7. [Metadata Schemas](#metadata-schemas)
8. [Implementation Details](#implementation-details)

---

## Overview

Protos V3 addresses key weaknesses in the dual identifier system:

- **Stable UUIDs**: Entity IDs are UUIDv4, not name-dependent
- **Canonical Aliases**: Normalized alias resolution with priority rules
- **Dataset IDs**: Datasets also have stable UUIDs
- **Multi-parent Relationships**: Support complex derivations
- **Scalable Registry**: SQLite backend with JSON export
- **Schema Validation**: Enforced metadata consistency
- **Configurable Paths**: Format-to-path mapping in config

---

## Core Improvements

### 1. UUID-Based Entity IDs

```python
import uuid

def generate_entity_id() -> str:
    """Generate stable entity ID independent of name."""
    # Use UUIDv4 for true stability
    return str(uuid.uuid4())

# Alternative: Content-based for deduplication
def generate_content_id(file_path: Path) -> str:
    """Generate ID from file content for deduplication."""
    hasher = hashlib.sha256()
    with open(file_path, 'rb') as f:
        # Read in chunks for large files
        while chunk := f.read(8192):
            hasher.update(chunk)
    return hasher.hexdigest()[:16]  # First 16 chars
```

### 2. Alias Resolution System

```python
class AliasResolver:
    """Canonical alias resolution with priority rules."""
    
    def __init__(self):
        self.rules = [
            # Rule 1: Case normalization
            lambda x: x.upper(),
            # Rule 2: Remove organism suffixes
            lambda x: re.sub(r'_HUMAN$|_MOUSE$|_RAT$', '', x),
            # Rule 3: UniProt to gene name
            self._uniprot_to_gene,
            # Rule 4: PDB chain notation
            lambda x: re.sub(r'_CHAIN_[A-Z]$', '', x)
        ]
        
    def normalize(self, alias: str) -> List[str]:
        """Generate normalized variations of alias."""
        normalized = [alias]  # Original always first
        
        for rule in self.rules:
            try:
                variant = rule(alias)
                if variant not in normalized:
                    normalized.append(variant)
            except:
                pass
                
        return normalized
    
    def resolve(self, query: str, aliases: List[str]) -> Optional[str]:
        """Find best matching alias."""
        # Direct match first
        if query in aliases:
            return query
            
        # Try normalized versions
        query_normalized = self.normalize(query)
        for q_norm in query_normalized:
            for alias in aliases:
                alias_normalized = self.normalize(alias)
                if q_norm in alias_normalized:
                    return alias
                    
        return None
```

### 3. Format Configuration

```yaml
# formats.yaml - Configurable format mappings
formats:
  structure:
    extensions: [".cif", ".pdb", ".mmcif"]
    subdirectory: "structures/pdb"
    processor: "CifBaseProcessor"
    metadata_schema: "structure_metadata_v1"
    
  sequence:
    extensions: [".fasta", ".fa", ".seq"]
    subdirectory: "sequences/fasta"
    processor: "SeqProcessor"
    metadata_schema: "sequence_metadata_v1"
    
  embedding:
    extensions: [".h5", ".npy", ".pkl"]
    subdirectory: "embeddings/{model_name}"
    processor: "EmbeddingProcessor"
    metadata_schema: "embedding_metadata_v1"
    
  custom_format:  # Easy to add new formats
    extensions: [".xyz"]
    subdirectory: "custom/data"
    processor: "CustomProcessor"
    metadata_schema: "custom_metadata_v1"
```

---

## Entity System

### Enhanced Entity Structure

```json
{
  "entity_id": "550e8400-e29b-41d4-a716-446655440000",  // UUIDv4
  "entity_name": "1UBQ",                                 // Current primary name
  "canonical_name": "UBIQ_HUMAN",                        // Canonical form
  "aliases": [
    {
      "name": "1UBQ",
      "type": "pdb_id",
      "priority": 1
    },
    {
      "name": "UBIQ_HUMAN", 
      "type": "gene_name",
      "priority": 2
    },
    {
      "name": "P62988",
      "type": "uniprot_id", 
      "priority": 3
    }
  ],
  "formats": {
    "structure": {
      "file_path": "structures/pdb/1UBQ.cif",
      "file_hash": "sha256:abcdef123456...",
      "file_size": 125840,
      "metadata": {
        // Validated against structure_metadata_v1 schema
        "resolution": 1.8,
        "method": "X-RAY",
        "chains": ["A"],
        "ligands": ["ZN", "SO4"]
      }
    }
  },
  "relationships": {
    "derived_from": ["entity_id_1", "entity_id_2"],  // Multi-parent
    "derived_entities": ["entity_id_3", "entity_id_4"],
    "version_chain": {
      "previous": "entity_id_0",
      "next": ["entity_id_5", "entity_id_6"]  // Branching versions
    }
  },
  "tags": ["reviewed", "high_quality", "drug_target"],
  "created": "2024-01-15T10:30:00Z",
  "modified": "2024-01-20T14:22:00Z"
}
```

### Entity Operations

```python
class EntityRegistry:
    def __init__(self, backend="sqlite"):
        if backend == "sqlite":
            self.backend = SQLiteBackend("entity_registry.db")
        else:
            self.backend = JSONBackend("entity_registry.json")
    
    def register_entity(self, name: str, format_type: str, 
                       file_path: str, metadata: Dict) -> str:
        """Register new entity with UUID."""
        # Generate stable UUID
        entity_id = str(uuid.uuid4())
        
        # Validate metadata against schema
        schema = self.get_metadata_schema(format_type)
        validate_metadata(metadata, schema)
        
        # Compute file hash for deduplication
        file_hash = compute_file_hash(file_path)
        
        # Check if content already exists
        existing = self.find_by_content_hash(file_hash)
        if existing:
            # Add as new format to existing entity
            entity_id = existing.entity_id
        
        # Create or update entity
        entity = Entity(
            entity_id=entity_id,
            entity_name=name,
            canonical_name=self.canonicalize(name),
            formats={format_type: {
                "file_path": file_path,
                "file_hash": file_hash,
                "metadata": metadata
            }}
        )
        
        self.backend.save_entity(entity)
        return entity_id
```

---

## Dataset System

### Dataset with UUID

```json
{
  "dataset_id": "660e8400-e29b-41d4-a716-446655440001",
  "dataset_name": "kinase_study",
  "format_type": "structure",  // Scoped by format
  "description": "Human kinase structures for drug design",
  "entities": [
    {
      "entity_id": "550e8400-e29b-41d4-a716-446655440000",
      "entity_name": "EGFR",  // Cached for performance
      "added": "2024-01-15T10:30:00Z"
    },
    {
      "entity_id": "770e8400-e29b-41d4-a716-446655440002",
      "entity_name": "ABL1",
      "added": "2024-01-16T14:20:00Z"
    }
  ],
  "metadata": {
    "organism": "homo_sapiens",
    "protein_family": "kinase",
    "curation_status": "reviewed"
  },
  "parent_datasets": [],  // Datasets can be derived too
  "child_datasets": ["880e8400-e29b-41d4-a716-446655440003"],
  "created": "2024-01-15T10:00:00Z",
  "modified": "2024-01-20T16:00:00Z"
}
```

### Dataset Operations

```python
class DatasetManager:
    def create_dataset(self, name: str, format_type: str,
                      entity_names: List[str], metadata: Dict) -> str:
        """Create dataset with UUID and format scoping."""
        # Generate dataset UUID
        dataset_id = str(uuid.uuid4())
        
        # Check for name collision within format
        existing = self.find_dataset(name, format_type)
        if existing:
            raise ValueError(
                f"Dataset '{format_type}:{name}' already exists"
            )
        
        # Resolve entity names to IDs
        entities = []
        for entity_name in entity_names:
            entity_id = self.registry.resolve_to_id(entity_name)
            if entity_id:
                entities.append({
                    "entity_id": entity_id,
                    "entity_name": entity_name,  # Cache current name
                    "added": datetime.now().isoformat()
                })
        
        dataset = Dataset(
            dataset_id=dataset_id,
            dataset_name=name,
            format_type=format_type,
            entities=entities,
            metadata=metadata
        )
        
        self.backend.save_dataset(dataset)
        return dataset_id
    
    def get_dataset_identifier(self, name: str, format_type: str) -> str:
        """Get fully qualified dataset identifier."""
        return f"{format_type}:{name}"
```

---

## Registry Architecture

### SQLite Backend for Scalability

```python
class SQLiteBackend:
    """Scalable backend for large registries."""
    
    def __init__(self, db_path: str):
        self.conn = sqlite3.connect(db_path)
        self._init_schema()
    
    def _init_schema(self):
        """Create tables with proper indices."""
        self.conn.executescript("""
            -- Entities table
            CREATE TABLE IF NOT EXISTS entities (
                entity_id TEXT PRIMARY KEY,
                entity_name TEXT NOT NULL,
                canonical_name TEXT,
                created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            CREATE INDEX IF NOT EXISTS idx_entity_name 
                ON entities(entity_name);
            CREATE INDEX IF NOT EXISTS idx_canonical 
                ON entities(canonical_name);
            
            -- Aliases table (normalized)
            CREATE TABLE IF NOT EXISTS aliases (
                entity_id TEXT,
                alias TEXT,
                alias_type TEXT,
                priority INTEGER,
                FOREIGN KEY (entity_id) REFERENCES entities(entity_id)
            );
            CREATE INDEX IF NOT EXISTS idx_alias ON aliases(alias);
            
            -- Formats table
            CREATE TABLE IF NOT EXISTS formats (
                entity_id TEXT,
                format_type TEXT,
                file_path TEXT,
                file_hash TEXT,
                metadata JSON,
                FOREIGN KEY (entity_id) REFERENCES entities(entity_id)
            );
            CREATE INDEX IF NOT EXISTS idx_format 
                ON formats(entity_id, format_type);
            CREATE INDEX IF NOT EXISTS idx_file_hash 
                ON formats(file_hash);
            
            -- Relationships table
            CREATE TABLE IF NOT EXISTS relationships (
                from_entity_id TEXT,
                to_entity_id TEXT,
                relationship_type TEXT,
                metadata JSON,
                FOREIGN KEY (from_entity_id) REFERENCES entities(entity_id),
                FOREIGN KEY (to_entity_id) REFERENCES entities(entity_id)
            );
            CREATE INDEX IF NOT EXISTS idx_relationships 
                ON relationships(from_entity_id, relationship_type);
            
            -- Datasets table
            CREATE TABLE IF NOT EXISTS datasets (
                dataset_id TEXT PRIMARY KEY,
                dataset_name TEXT NOT NULL,
                format_type TEXT NOT NULL,
                metadata JSON,
                created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                modified TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            );
            CREATE UNIQUE INDEX IF NOT EXISTS idx_dataset_name 
                ON datasets(dataset_name, format_type);
            
            -- Dataset members table
            CREATE TABLE IF NOT EXISTS dataset_members (
                dataset_id TEXT,
                entity_id TEXT,
                added TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (dataset_id) REFERENCES datasets(dataset_id),
                FOREIGN KEY (entity_id) REFERENCES entities(entity_id)
            );
            CREATE INDEX IF NOT EXISTS idx_dataset_members 
                ON dataset_members(dataset_id);
        """)
    
    def find_entity(self, name: str) -> Optional[Entity]:
        """Fast entity lookup with alias resolution."""
        # Try direct name match
        cursor = self.conn.execute("""
            SELECT entity_id FROM entities 
            WHERE entity_name = ? OR canonical_name = ?
        """, (name, name))
        
        if row := cursor.fetchone():
            return self.get_entity(row[0])
        
        # Try alias match
        cursor = self.conn.execute("""
            SELECT entity_id FROM aliases 
            WHERE alias = ?
            ORDER BY priority
            LIMIT 1
        """, (name,))
        
        if row := cursor.fetchone():
            return self.get_entity(row[0])
        
        return None
    
    def export_to_json(self, output_path: str):
        """Export database to JSON for portability."""
        # Implementation exports full registry as JSON
        pass
```

---

## Relationship System

### Multi-Parent Support

```json
{
  "relationships": {
    "derived_from": [
      {
        "entity_id": "550e8400-e29b-41d4-a716-446655440000",
        "relationship_type": "template",
        "metadata": {"coverage": "90%"}
      },
      {
        "entity_id": "660e8400-e29b-41d4-a716-446655440001",
        "relationship_type": "ligand_source",
        "metadata": {"ligands": ["ATP", "MG"]}
      }
    ],
    "derived_entities": [
      {
        "entity_id": "770e8400-e29b-41d4-a716-446655440002",
        "relationship_type": "filtered",
        "metadata": {"filter": "chain_A_only"}
      }
    ],
    "related_entities": [
      {
        "entity_id": "880e8400-e29b-41d4-a716-446655440003",
        "relationship_type": "homolog",
        "metadata": {"identity": 0.85}
      }
    ]
  }
}
```

### Relationship Queries

```python
class RelationshipManager:
    def add_relationship(self, from_id: str, to_id: str, 
                        rel_type: str, metadata: Dict = None):
        """Add typed relationship between entities."""
        self.backend.add_relationship(
            from_entity_id=from_id,
            to_entity_id=to_id,
            relationship_type=rel_type,
            metadata=metadata or {}
        )
    
    def get_lineage(self, entity_id: str) -> Dict:
        """Get complete lineage tree."""
        lineage = {
            "parents": self.get_parents(entity_id),
            "children": self.get_children(entity_id),
            "siblings": self.get_siblings(entity_id)
        }
        return lineage
    
    def find_related(self, entity_id: str, 
                    rel_type: str = None,
                    max_depth: int = 1) -> List[str]:
        """Find related entities up to max_depth."""
        # Graph traversal implementation
        pass
```

---

## Metadata Schemas

### JSON Schema Validation

```json
{
  "structure_metadata_v1": {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "required": ["method", "resolution", "chains"],
    "properties": {
      "method": {
        "type": "string",
        "enum": ["X-RAY", "NMR", "EM", "NEUTRON", "MODEL"]
      },
      "resolution": {
        "type": "number",
        "minimum": 0
      },
      "chains": {
        "type": "array",
        "items": {"type": "string", "pattern": "^[A-Za-z0-9]$"}
      },
      "ligands": {
        "type": "array",
        "items": {"type": "string"}
      },
      "organism": {
        "type": "string"
      },
      "expression_system": {
        "type": "string"
      }
    }
  },
  
  "sequence_metadata_v1": {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "required": ["length", "sequence_type"],
    "properties": {
      "length": {
        "type": "integer",
        "minimum": 1
      },
      "sequence_type": {
        "type": "string",
        "enum": ["protein", "dna", "rna"]
      },
      "organism": {
        "type": "string"
      },
      "source": {
        "type": "string",
        "enum": ["uniprot", "genbank", "custom", "extracted"]
      }
    }
  }
}
```

### Schema Validation

```python
import jsonschema

class MetadataValidator:
    def __init__(self, schema_dir: Path):
        self.schemas = self._load_schemas(schema_dir)
    
    def validate(self, metadata: Dict, format_type: str):
        """Validate metadata against format schema."""
        schema_name = f"{format_type}_metadata_v1"
        
        if schema_name not in self.schemas:
            raise ValueError(f"No schema for format '{format_type}'")
        
        try:
            jsonschema.validate(metadata, self.schemas[schema_name])
        except jsonschema.ValidationError as e:
            raise ValueError(f"Invalid metadata: {e.message}")
    
    def migrate_metadata(self, metadata: Dict, 
                        from_version: str, to_version: str):
        """Migrate metadata between schema versions."""
        # Version migration logic
        pass
```

---

## Implementation Details

### Configuration Management

```python
class ProtosConfig:
    """Central configuration for paths and formats."""
    
    def __init__(self, config_path: Path = None):
        self.config_path = config_path or Path.home() / ".protos" / "config.yaml"
        self.config = self._load_config()
    
    def get_format_config(self, format_type: str) -> Dict:
        """Get configuration for a format type."""
        return self.config["formats"].get(format_type, {})
    
    def get_format_path(self, format_type: str) -> Path:
        """Get subdirectory for format."""
        config = self.get_format_config(format_type)
        template = config.get("subdirectory", f"{format_type}/data")
        
        # Support dynamic paths like "embeddings/{model_name}"
        return self._resolve_path_template(template)
    
    def register_format(self, format_type: str, config: Dict):
        """Register new format type."""
        self.config["formats"][format_type] = config
        self._save_config()
```

### Migration Tools

```python
class RegistryMigration:
    """Tools for migrating between registry versions."""
    
    @staticmethod
    def migrate_v2_to_v3(old_registry: Dict) -> None:
        """Migrate from hash-based to UUID-based IDs."""
        backend = SQLiteBackend("entity_registry.db")
        
        for old_id, entity_data in old_registry["entities"].items():
            # Generate new UUID
            new_id = str(uuid.uuid4())
            
            # Migrate entity
            entity = Entity(
                entity_id=new_id,
                entity_name=entity_data["entity_name"],
                canonical_name=canonicalize(entity_data["entity_name"]),
                aliases=[{"name": a, "type": "legacy", "priority": 5} 
                        for a in entity_data.get("aliases", [])],
                formats=entity_data.get("formats", {}),
                created=entity_data.get("created"),
                modified=entity_data.get("modified")
            )
            
            backend.save_entity(entity)
            
            # Migrate relationships if present
            if "relationships" in entity_data:
                # Map old IDs to new UUIDs
                pass
```

---

## Summary

Protos V3 provides a robust, scalable data management system:

1. **Stable UUIDs**: True stability independent of names
2. **Smart Aliases**: Canonical resolution with priority rules
3. **Scoped Datasets**: Format-specific with UUID identification
4. **Multi-Parent**: Support complex derivation relationships
5. **Scalable Backend**: SQLite for performance, JSON for portability
6. **Schema Validation**: Consistent metadata with migrations
7. **Configurable Paths**: Easy format extension without code changes

This design scales from small projects to enterprise-level data management while maintaining simplicity for users.