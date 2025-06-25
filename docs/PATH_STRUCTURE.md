# Protos Path Structure

## Overview

Protos organizes data in a hierarchical directory structure within a single data directory. This unified approach ensures consistent file organization across different data types and processor modules while simplifying data management.

## Directory Structure

The standard Protos data directory structure:

```
~/protos_data/              # Default location (or $PROTOS_DATA_ROOT)
├── structure/              # Structure-related data
│   ├── mmcif/              # PDB/mmCIF structure files
│   ├── alignments/         # Structure alignment results
│   ├── structure_dataset/  # Structure dataset definitions
│   ├── temp_cif/           # Temporary structure files
│   └── registry.json       # Structure dataset registry
│
├── grn/                    # GRN (Generic Residue Numbering) data
│   ├── ref/                # Reference GRN tables
│   ├── tables/             # Calculated GRN tables
│   ├── configs/            # GRN configuration files
│   ├── assignments/        # GRN assignment results
│   └── registry.json       # GRN dataset registry
│
├── sequence/               # Sequence-related data
│   ├── fasta/              # FASTA sequence files
│   ├── alignments/         # Sequence alignment results
│   ├── metadata/           # Sequence metadata
│   └── registry.json       # Sequence dataset registry
│
├── graph/                  # Graph-related data
│   └── registry.json
│
├── property/               # Property and metadata
│   └── registry.json
│
├── embedding/              # Protein embeddings
│   └── registry.json
│
└── global_registry.json    # Global dataset registry
```

Each processor type has its own subdirectory with a local `registry.json` file for backward compatibility. The global registry provides a unified view of all datasets.

## Default Location

The data directory is located at:

1. **Environment Variable**: If `PROTOS_DATA_ROOT` is set, that location is used
2. **Home Directory**: `~/protos_data/` (default)
3. **Custom Path**: When explicitly specified during processor initialization

The directory is automatically created when any processor or registry is first used.

## Registry System

The registry system provides a catalog of all datasets:

### Global Registry

Located at `data_root/global_registry.json`, maintains a unified view of all datasets across all processors.

### Processor Registries

Located at `data_root/processor_type/registry.json`, maintained for backward compatibility.

### Registry Format

```json
{
  "dataset_id": {
    "path": "/home/user/protos_data/structure/mmcif/1abc.cif",
    "metadata": {
      "processor_type": "structure",
      "dataset_type": "cif",
      "description": "Crystal structure of protein ABC",
      "resolution": 2.5,
      "chains": ["A", "B"]
    },
    "timestamp": "2023-01-01T00:00:00"
  }
}
```

## File Organization Examples

### Structure Files
```
~/protos_data/structure/
├── mmcif/
│   ├── 1abc.cif        # Individual structure files
│   ├── 2xyz.cif
│   └── 3def.cif
└── structure_dataset/
    ├── kinases.json    # Dataset definition listing multiple structures
    └── gpcrs.json
```

### GRN Data
```
~/protos_data/grn/
├── ref/
│   ├── gpcr_ref.csv    # Reference GRN table for GPCRs
│   └── opsin_ref.csv   # Reference table for opsins
├── tables/
│   ├── my_proteins.csv # User-generated GRN table
│   └── project_x.csv
└── configs/
    └── config.json     # GRN numbering configuration
```

### Sequence Data
```
~/protos_data/sequence/
├── fasta/
│   ├── human_proteome.fasta
│   └── project_sequences.fasta
└── alignments/
    ├── family_alignment.aln
    └── pairwise_results.csv
```

## Best Practices

1. **Use processor APIs** - They handle all path resolution automatically
2. **Organize by type** - Keep files in their appropriate subdirectories
3. **Use descriptive names** - Make dataset IDs and filenames self-documenting
4. **Let the system create directories** - They're created automatically as needed
5. **One data directory** - Everything lives in one place for easy backup and sharing
6. **Use the registry** - Track datasets through the registry system, not hardcoded paths

## Automatic Directory Creation

Directories are created automatically when needed:

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# This automatically creates:
# ~/protos_data/
# ~/protos_data/structure/
# ~/protos_data/structure/mmcif/
# ~/protos_data/structure/alignments/
# ~/protos_data/structure/structure_dataset/
# ~/protos_data/global_registry.json
processor = CifBaseProcessor(name="my_processor")
```

No manual directory creation is required - the system handles it all.