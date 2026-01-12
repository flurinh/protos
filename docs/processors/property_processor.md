# PropertyProcessor

Manages tabular properties linked to existing entities.

**Location:** `protos.processing.property.property_processor`

**Processor Type:** `property`

## Overview

The PropertyProcessor provides:
- Registry-aware tabular property storage
- Properties that annotate existing entities
- Scope-based linking to structures, sequences, etc.
- CSV-based storage with optional row materialization
- Index files for efficient entity lookup

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(name)` | Load property entry by ID. Returns dict or None. |
| `save_entity(name, data, metadata)` | Save property entry with scope. |

### Table Operations

| Method | Description |
|--------|-------------|
| `record_properties(table_name, rows, metadata, allow_create, materialize_entries)` | Insert rows into property table with entity linking. |
| `load_table(table_name)` | Load property table as DataFrame. |
| `list_tables()` | List all property tables. |
| `list_datasets()` | List property datasets (alias). |

### Property Table CRUD

| Method | Description |
|--------|-------------|
| `create_property_table(table_name, columns, metadata)` | Create empty table with schema. |
| `load_property_table(table_name)` | Load table (alias for load_table). |
| `save_property_table(table_name, df, metadata)` | Save/overwrite entire table. |
| `add_property_column(table_name, column_name, default_value)` | Add new column to existing table. |

### Querying

| Method | Description |
|--------|-------------|
| `get_properties(entity_name, table_name, format_type)` | Get properties for specific entity. |
| `load_dataset_rows(table_name, entity_names, format_type)` | Load rows for multiple entities. |
| `filter_by_property(table_name, column, operator, value)` | Filter table by property value. |

### Path Properties

| Property | Description |
|----------|-------------|
| `tables_dir` | CSV property tables directory |
| `datasets_dir` | Dataset definitions directory |

---

## Key Concepts

### Scope Column

The `scope` column links property rows to existing entities:

```python
# Single scope - annotates one entity
scope = [{"format": "structure", "name": "1ubq"}]

# Multiple scopes - annotates multiple entities
scope = [
    {"format": "structure", "name": "1ubq"},
    {"format": "sequence", "name": "P0CG48"}
]
```

### Materialization Modes

| Mode | `materialize_entries` | Description |
|------|----------------------|-------------|
| Compact | `False` (default) | Single dataset entity, index file for lookup |
| Full | `True` | Each row is a separate entity |

### Special Columns

| Column | Description |
|--------|-------------|
| `property_entry_id` | Unique row ID (if materialized) |
| `entity_name` | Primary entity name (from last scope) |
| `scope` | List of entity references |

---

## Usage Examples

### Basic Property Recording

```python
from protos.processing.property import PropertyProcessor

proc = PropertyProcessor()

# Record properties for structures
rows = [
    {
        "scope": [{"format": "structure", "name": "1ubq"}],
        "binding_affinity": 5.2,
        "ic50": 100.0,
        "assay": "SPR"
    },
    {
        "scope": [{"format": "structure", "name": "3sn6"}],
        "binding_affinity": 3.8,
        "ic50": 50.0,
        "assay": "ITC"
    }
]

table = proc.record_properties(
    "binding_data",
    rows,
    metadata={"experiment": "binding_assay"}
)
```

### Loading Properties

```python
# Load entire table
df = proc.load_table("binding_data")
print(df)

# Get properties for specific entity
props = proc.get_properties("1ubq", table_name="binding_data")

# Filter by value
high_affinity = proc.filter_by_property(
    "binding_data",
    column="binding_affinity",
    operator=">",
    value=5.0
)
```

### Recording with DataFrame

```python
import pandas as pd

df = pd.DataFrame([
    {
        "scope": [{"format": "structure", "name": "protein_a"}],
        "property_1": 1.5,
        "property_2": "active"
    },
    {
        "scope": [{"format": "structure", "name": "protein_b"}],
        "property_1": 2.3,
        "property_2": "inactive"
    }
])

proc.record_properties("my_properties", df)
```

### Multiple Scope Annotations

```python
# Property linking structure and sequence
rows = [{
    "scope": [
        {"format": "structure", "name": "1ubq"},
        {"format": "sequence", "name": "P0CG48"}
    ],
    "cross_reference": "verified",
    "confidence": 0.95
}]

proc.record_properties("cross_refs", rows)
```

### Materialization Options

```python
# Compact mode (default) - efficient for large tables
proc.record_properties(
    "residue_properties",
    residue_rows,
    materialize_entries=False  # Single entity, index file
)

# Full materialization - each row is an entity
proc.record_properties(
    "important_findings",
    rows,
    materialize_entries=True  # Individual entity per row
)
```

### Entity Creation Control

```python
# Fail if entity doesn't exist
try:
    proc.record_properties("data", rows, allow_create=False)
except ValueError as e:
    print(f"Entity not found: {e}")

# Auto-create placeholder entities
proc.record_properties("data", rows, allow_create=True)
```

### Creating and Managing Tables

```python
# Create empty table with schema
proc.create_property_table(
    "experiment_results",
    columns=["entity_name", "scope", "value", "error", "method"],
    metadata={"created_by": "pipeline"}
)

# Add column to existing table
proc.add_property_column(
    "experiment_results",
    column_name="validated",
    default_value=False
)

# Overwrite entire table
proc.save_property_table("experiment_results", new_df)
```

### Querying Multiple Entities

```python
# Load rows for specific entities
rows = proc.load_dataset_rows(
    "binding_data",
    entity_names=["1ubq", "3sn6", "4lde"],
    format_type="structure"
)

# Filter with comparison
active = proc.filter_by_property(
    "screening_results",
    column="activity",
    operator="==",
    value="active"
)
```

### Checking Relationships

```python
from protos.io.core import get_registry

registry = get_registry()

# Find what properties annotate an entity
related = registry.get_relationships(
    "1ubq",
    rel_type="annotated_by",
    direction="incoming"
)

for rel in related:
    print(f"Property table: {rel['source_name']}")
    print(f"Table name: {rel['metadata'].get('table')}")
```

### Pipeline Example

```python
from protos import StructureProcessor
from protos.processing.property import PropertyProcessor

struct_proc = StructureProcessor()
prop_proc = PropertyProcessor()

# Calculate properties for dataset
rows = []
for entity in struct_proc.dataset_manager.get_dataset_entities("my_study"):
    df = struct_proc.load_entity(entity)

    rows.append({
        "scope": [{"format": "structure", "name": entity}],
        "atom_count": len(df),
        "chain_count": df.reset_index()["auth_chain_id"].nunique(),
        "has_ligand": (df["group"] == "HETATM").any()
    })

# Record computed properties
prop_proc.record_properties(
    "structure_stats",
    rows,
    metadata={"pipeline": "structure_analysis"}
)
```
