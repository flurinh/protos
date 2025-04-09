# Protos Path Structure

## Overview

Protos organizes data in a hierarchical directory structure to ensure consistent file organization across different data types and processor modules. The path structure is designed to separate reference data (distributed with the package) from user data (created at runtime).

## Data Sources

Protos recognizes two distinct data sources:

1. **Reference Data**: Read-only data distributed with the package, containing essential reference files and default datasets
2. **User Data**: Read-write data created during runtime, containing user-specific datasets and analysis results

## Directory Structure

Both reference and user data follow the same basic directory structure:

```
data_root/
├── structure/               # Structure-related data
│   ├── mmcif/               # mmCIF structure files
│   ├── alignments/          # Structure alignment results
│   ├── structure_dataset/   # Structure dataset files
│   └── temp_cif/            # Temporary structure files
├── grn/                     # GRN-related data
│   ├── tables/              # GRN tables
│   ├── grn/                 # GRN definitions
│   ├── configs/             # GRN configuration files
│   └── assignments/         # GRN assignment results
├── sequence/                # Sequence-related data
│   ├── fasta/               # FASTA sequence files
│   ├── alignments/          # Sequence alignment results
│   └── metadata/            # Sequence metadata
├── graph/                   # Graph-related data
├── property/                # Property-related data
├── embedding/               # Embedding-related data
└── global_registry.json     # Global dataset registry (user data only)
```

Each processor-specific directory also contains a `registry.json` file for tracking datasets associated with that processor.

## Path Resolution

### Reference Data Paths

Reference data is located within the package at:

```
protos/reference_data/
```

This location is determined during package installation and is normally not directly accessed by users.

### User Data Paths

User data is located at either:

1. The directory specified by the `PROTOS_DATA_ROOT` environment variable
2. A `data/` directory in the current working directory (default)

## Registry System

The registry system maps dataset IDs to file paths and metadata:

### Global Registry

Located at `user_data_root/global_registry.json`, provides a unified view of all datasets (both reference and user).

### Processor Registries

Located at `data_root/processor_type/registry.json`, provides backward compatibility with older code.

### Registry Format

```json
{
  "dataset_id": {
    "path": "/path/to/file",
    "metadata": {
      "processor_type": "structure",
      "dataset_type": "cif",
      "source": "user",
      "description": "Description of the dataset"
    },
    "timestamp": "2023-01-01T00:00:00"
  }
}
```

## Configuration

The path structure can be configured through:

### Environment Variables

- `PROTOS_DATA_ROOT`: Sets the root directory for user data
- `PROTOS_REF_DATA_ROOT`: Overrides the reference data location (rarely needed)

### Programmatic Configuration

```python
from protos.io.paths import ProtosPaths, DataSource

# Initialize with custom user data path
paths = ProtosPaths(
    user_data_root="/path/to/user/data",
    ref_data_root=None,  # Defaults to package reference data
    create_dirs=True,
    validate=True
)

# Access paths
structure_path = paths.get_structure_subdir_path("structure_dir", DataSource.USER)
```

## Best Practices

1. Use the high-level API when possible, which handles path resolution automatically
2. For writing data, always use user data paths
3. Use dataset IDs and the registry system instead of hardcoding paths
4. Use the `GlobalRegistry` to manage datasets across processor types
5. Prefer relative paths in registry entries to improve portability
6. Use `pathlib.Path` for cross-platform compatibility
7. Always check both reference and user data locations when loading files
8. Let processor classes handle path resolution automatically when possible