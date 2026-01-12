# GRNProcessor

Manages Generic Residue Numbering (GRN) tables and sequence annotations.

**Location:** `protos.processing.grn.grn_processor`

**Processor Type:** `grn`

## Overview

The GRNProcessor provides:
- Storage and management of GRN reference tables
- Sequence annotation with GRN positions
- Alignment-based GRN assignment
- Reference tables bundled with Protos

## GRN Concepts

Generic Residue Numbering (GRN) provides a universal numbering system for homologous residues across protein families. For example, GPCRs use the Ballesteros-Weinstein numbering (e.g., "3.50x50" for a conserved position in TM3).

**GRN Table Format:**
- Rows: Sequence identifiers
- Columns: GRN position labels (e.g., "1.50x50", "3.50x50")
- Cells: Values like "M123" (residue M at sequence position 123)

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(name)` | Load GRN data for entity from any registered table. Returns dict. |
| `save_entity(name, data, metadata)` | Save GRN row for entity. Requires table name in metadata. |

### Table Operations

| Method | Description |
|--------|-------------|
| `record_table(table_name, table, metadata, per_entity_metadata, link_entities)` | Persist GRN table and create entity relationships. |
| `load_table(table_name)` | Load a recorded GRN table as DataFrame. |
| `list_tables()` | List all registered GRN tables. |
| `load_reference_table(reference_name)` | Load bundled reference GRN table. |

### Sequence Annotation

| Method | Description |
|--------|-------------|
| `annotate_sequences(sequences, reference_table, protein_family, search)` | Annotate sequences with GRN using reference alignment. Returns (DataFrame, summary). |

### Static Utility Methods

| Method | Description |
|--------|-------------|
| `parse_grn_value(value)` | Parse cell value like "M123" into (residue, position) tuple. |
| `build_grn_to_seq_index(grn_table, sequence_id)` | Build mapping from GRN labels to sequence positions. |

### Path Properties

| Property | Description |
|----------|-------------|
| `tables_dir` | GRN mapping tables directory |
| `reference_dir` | Bundled reference data directory |
| `configs_dir` | Configuration files directory |

---

## Data Format

### GRN Table Structure

```
              1.50x50  2.50x50  3.50x50  7.50x50
entity_name
seq_a         N45      V98      R123     P312
seq_b         N47      V100     R125     -
seq_c         -        V95      R121     P308
```

- Index: `entity_name` (sequence identifiers)
- Columns: GRN position labels
- Values: Residue + sequence position (e.g., "R123") or "-" for missing

### Cell Value Format

| Value | Meaning |
|-------|---------|
| `M123` | Methionine at sequence position 123 |
| `R45` | Arginine at sequence position 45 |
| `-` | Position not present in this sequence |

---

## Usage Examples

### Basic Table Operations

```python
from protos.processing.grn import GRNProcessor

proc = GRNProcessor()

# Load bundled reference
ref = proc.load_reference_table("class_a_gpcr")
print(f"Reference has {len(ref)} sequences")
print(f"GRN columns: {list(ref.columns)[:5]}")

# List available tables
tables = proc.list_tables()

# Load recorded table
table = proc.load_table("my_grn_table")
```

### Recording GRN Tables

```python
import pandas as pd

# Create GRN table
data = {
    "1.50x50": ["N45", "N47", "-"],
    "3.50x50": ["R123", "R125", "R121"],
    "7.50x50": ["P312", "-", "P308"]
}
df = pd.DataFrame(data, index=["seq_a", "seq_b", "seq_c"])
df.index.name = "entity_name"

# Record with entity linking
table = proc.record_table(
    "my_receptor_grn",
    df,
    metadata={"family": "gpcr", "class": "A"},
    per_entity_metadata={
        "seq_a": {"organism": "human"},
        "seq_b": {"organism": "mouse"}
    },
    link_entities=True,
    allow_create=False
)
```

### Sequence Annotation

```python
# Define sequences to annotate
sequences = {
    "my_gpcr_1": "MTDKYRVFFVNVITNTMVM...",
    "my_gpcr_2": "MGSRVRDYFFVNVITNTMV..."
}

# Annotate using reference
annotations, summary = proc.annotate_sequences(
    sequences,
    reference_table="class_a_gpcr",
    protein_family="gpcr_a",
    search="pairwise"
)

# Check results
print(f"Annotated: {summary['global']['annotated']}/{summary['global']['total']}")

for seq_name, info in summary['per_sequence'].items():
    print(f"\n{seq_name}:")
    print(f"  Best reference: {info['reference']}")
    print(f"  Coverage: {info['coverage']:.2%}")
    print(f"  Status: {info['status']}")

# Save annotations
proc.record_table("my_annotations", annotations)
```

### GRN Parsing Utilities

```python
# Parse single cell value
parsed = GRNProcessor.parse_grn_value("M123")
print(parsed)  # ('M', 123)

parsed = GRNProcessor.parse_grn_value("-")
print(parsed)  # None

# Build position mapping for sequence
mapping = GRNProcessor.build_grn_to_seq_index(
    annotations,
    sequence_id="my_gpcr_1"
)

# Use mapping
print(f"3.50x50 is at position: {mapping.get('3.50x50')}")
print(f"7.50x50 is at position: {mapping.get('7.50x50')}")
```

### Entity Operations

```python
# Load GRN data for specific entity
grn_data = proc.load_entity("my_gpcr_1")
if grn_data:
    print(f"3.50x50: {grn_data.get('3.50x50')}")

# Save GRN data for entity
proc.save_entity(
    "new_sequence",
    {"1.50x50": "N45", "3.50x50": "R123"},
    metadata={"table": "my_grn_table"}
)
```

### Integration with StructureProcessor

```python
from protos import StructureProcessor
from protos.processing.grn import GRNProcessor

struct_proc = StructureProcessor()
grn_proc = GRNProcessor()

# Get GRN mapping
grn_table = grn_proc.load_table("my_grn_table")
mapping = GRNProcessor.build_grn_to_seq_index(grn_table, sequence_id="my_seq")

# Annotate structure with GRN
struct_proc.annotate_with_grn(
    "my_structure",
    grn_table=grn_table,
    sequence_id="my_seq"
)
```

### Checking Entity Relationships

```python
from protos.io.core import get_registry

registry = get_registry()

# Find entities annotated by GRN table
related = registry.get_relationships(
    "my_grn_table",
    rel_type="annotated_by",
    direction="outgoing"
)

for rel in related:
    print(f"Annotates: {rel['target_name']}")
```
