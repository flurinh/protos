# Protos Data Usage Guide

## Overview

Protos is a closed system that manages its own data structure. **DO NOT use input/output folders** - all data should be organized within the Protos data structure.

## Data Organization

### For Reference Data (included with Protos)
Location: `src/protos/reference_data/`

```
src/protos/reference_data/
├── grn/
│   ├── ref/              # Reference GRN tables (e.g., mo_ref.csv)
│   ├── tables/           # Generated GRN tables
│   └── configs/          # GRN configuration files
├── sequence/
│   ├── fasta/            # FASTA sequence files
│   └── alignments/       # Sequence alignments
└── structure/
    ├── mmcif/            # PDB/CIF structure files
    └── structure_dataset/ # Dataset definitions
```

### For User Data (during runtime)
Location: `~/protos_data/` or `$PROTOS_DATA_ROOT`

The same structure is automatically created when processors are used.

### For Test Data
Location: `tests/test-data/`

Follows the same structure for testing purposes.

## Proper Usage Examples

### Loading Sequences
```python
# WRONG - Don't use input folder
sequences = read_fasta("input/sequences.fasta")

# CORRECT - Use Protos data structure
sequences = read_fasta("src/protos/reference_data/sequence/fasta/sequences.fasta")
```

### Saving GRN Tables
```python
# WRONG - Don't use output folder
output_df.to_csv("output/results.csv")

# CORRECT - Save to GRN tables directory
output_file = Path("src/protos/reference_data/grn/tables/results.csv")
output_file.parent.mkdir(exist_ok=True, parents=True)
output_df.to_csv(output_file)
```

### Using Processors
```python
from protos.processing.grn.grn_base_processor import GRNBaseProcessor

# Processors automatically handle paths
processor = GRNBaseProcessor(name="my_analysis")
# This creates ~/protos_data/grn/ structure automatically
```

## Key Points

1. **No input/output folders** - Protos manages its own data structure
2. **Reference data** goes in `src/protos/reference_data/`
3. **User data** goes in `~/protos_data/` (or `$PROTOS_DATA_ROOT`)
4. **Test data** goes in `tests/test-data/`
5. **Always use the appropriate subdirectory** based on data type (grn/, sequence/, structure/)
6. **Let processors handle paths** - they know where to find and save data

## Migration from input/output

If you have data in old input/output folders:

1. FASTA files from `input/` → `src/protos/reference_data/sequence/fasta/`
2. GRN tables from `output/` → `src/protos/reference_data/grn/tables/`
3. Structure files → `src/protos/reference_data/structure/mmcif/`

Then remove the input/output folders - they are not part of the Protos system.