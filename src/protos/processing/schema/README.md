# Standard Schema Definitions for Protos

This module provides standardized schema definitions, interfaces, and conversion utilities for the Protos package. These standards ensure consistent data exchange between different processors and utilities.

## Components

### 1. Schema Definitions (`schema_definitions.py`)

Defines standard DataFrame schemas and data formats for:

- **Structure Data**: Atom-level structure information
  - Core columns: `pdb_id`, `auth_chain_id`, `atom_name`, `x`, `y`, `z`, etc.
  - Extended columns: `phi`, `psi`, `omega`, `b_factor`, etc.
  - Standard multi-index: `[structure_idx, atom_idx]`

- **GRN Data**: GRN assignments and tables
  - Core columns: `id`, `name`, `species`, `family`, etc.
  - GRN patterns: Standard (`1x50`), N-terminal (`n.10`), C-terminal (`c.5`), Loop (`2.45`)
  - Gap and unknown symbols: `-` and `X`

- **Sequence Data**: Sequence alignments and annotations
  - Alignment columns: `query_id`, `target_id`, `sequence_identity`, etc.
  - Feature columns: `seq_id`, `position`, `residue`, `feature_type`, etc.

Includes validation and creation functions for each schema type.

### 2. Interface Definitions (`interface_definitions.py`)

Defines abstract interfaces for common operations:

- **StructureInterface**: Standard methods for structure operations
  - `load_structure`, `save_structure`, `get_atoms`, `get_chains`, etc.

- **GRNInterface**: Standard methods for GRN operations
  - `load_grn_table`, `assign_grns`, `map_grns_to_structure`, etc.

- **SequenceInterface**: Standard methods for sequence operations
  - `align_sequences`, `load_fasta`, `calculate_identity`, etc.

Includes validation decorators and standard return value types.

### 3. Conversion Utilities (`conversion_utilities.py`)

Provides utility functions for converting between different data formats:

- **Structure Conversions**:
  - CIF ↔ Structure DataFrame
  - Structure DataFrame ↔ Dictionary
  - 3-letter ↔ 1-letter amino acid codes

- **GRN Conversions**:
  - GRN string ↔ float (for sorting)
  - GRN mapping ↔ DataFrame
  - Sorting GRNs in standard order

- **Sequence Conversions**:
  - FASTA ↔ Dictionary
  - Alignment results ↔ DataFrame

## Usage Examples

### Structure Schema

```python
from protos.processing.schema import (
    STRUCTURE_CORE_COLUMNS,
    create_empty_structure_df,
    validate_structure_df
)

# Create a new structure DataFrame
df = create_empty_structure_df()

# Add data to the DataFrame
df.loc[(0, 0), :] = {
    'pdb_id': '1abc',
    'auth_chain_id': 'A',
    # ... other required columns
}

# Validate the DataFrame
validate_structure_df(df)
```

### GRN Schema

```python
from protos.processing.schema import (
    sort_grns,
    grn_str_to_float,
    grn_float_to_str
)

# Sort GRNs in standard order
sorted_grns = sort_grns(['3x50', '1x50', '7x50', 'n.10', 'c.5', '2.45'])

# Convert between GRN representations
float_val = grn_str_to_float('1x50')  # 1.5
grn_str = grn_float_to_str(1.5)       # '1x50'
```

### Sequence Schema

```python
from protos.processing.schema import (
    fasta_to_dict,
    dict_to_fasta,
    alignment_result_to_df
)

# Convert between FASTA and dictionary
sequences = fasta_to_dict(fasta_content)
fasta_content = dict_to_fasta(sequences, width=80)

# Convert alignment results to DataFrame
df = alignment_result_to_df(
    query_id='seq1',
    target_id='seq2',
    aligned_query='ACDEFG-HI',
    aligned_target='ACDE-GYHI',
    score=0.8
)
```

## Implementation Guide

When implementing these standards in your code:

1. **Use Standard Schemas**: Always use the standard schemas for storing and exchanging data
2. **Implement Interfaces**: When creating processor classes, implement the appropriate interfaces
3. **Validate Data**: Use the validation functions to ensure data conforms to the schemas
4. **Use Conversion Utilities**: When transforming data between formats, use the provided utilities

## Development

This module is part of the standardization effort for Protos and serves as the foundation for more consistent and maintainable code. Future updates will include:

- Additional validation rules for specialized data types
- Performance optimizations for large datasets
- Extended schema definitions for new data types
- More conversion utilities for other data formats
