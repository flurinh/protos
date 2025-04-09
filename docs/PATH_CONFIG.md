# Path Configuration in Protos

## Overview

Protos implements a flexible path management system that handles two distinct types of data:

1. **Reference Data**: Read-only data distributed with the Protos package
2. **User Data**: Read-write data created at runtime by the user

This document explains how the path system works and how to configure it for your needs.

## Data Types

### Reference Data

Reference data is distributed with the Protos package and includes:

- Standard structure files
- Reference GRN tables
- Standard sequence datasets
- Default configurations

Reference data is located within the package at `protos/reference_data/` and is read-only. It should not be modified directly by users.

### User Data

User data is created at runtime and includes:

- User-generated datasets
- Custom structures
- Analysis results
- User configurations

User data is stored in the directory specified by the `PROTOS_DATA_ROOT` environment variable, or defaults to a `data/` directory in the current working directory.

## Configuration

### Environment Variables

Protos uses these environment variables for path configuration:

- `PROTOS_DATA_ROOT`: Location for user data (default: `data/` in current directory)
- `PROTOS_REF_DATA_ROOT`: Override for reference data location (usually not needed)

### Programmatic Configuration

Paths can be configured programmatically:

```python
from protos.io.paths import ProtosPaths, DataSource

# Initialize with custom paths
paths = ProtosPaths(
    user_data_root="/path/to/user/data",
    ref_data_root="/path/to/reference/data",  # Optional override
    create_dirs=True,  # Create directories if they don't exist
    validate=True      # Validate directory structure
)

# Use the paths object
structure_path = paths.get_processor_path("structure", DataSource.USER)
grn_path = paths.get_grn_subdir_path("table_dir", DataSource.REFERENCE)
```

## Directory Structure

Both reference and user data follow the same directory structure:

```
data_root/
├── structure/
│   ├── mmcif/
│   ├── alignments/
│   ├── structure_dataset/
│   └── temp_cif/
├── grn/
│   ├── tables/
│   ├── grn/
│   ├── configs/
│   └── assignments/
├── sequence/
│   ├── fasta/
│   ├── alignments/
│   └── metadata/
└── global_registry.json
```

Each processor type directory also has its own `registry.json` file for backward compatibility.

## Processor Path Resolution

Processor classes automatically resolve file paths without requiring user intervention:

1. **Initialization**: When a processor is instantiated, it automatically determines its type from the class name
   ```python
   # CifProcessor sets processor_type = "structure"
   # GRNProcessor sets processor_type = "grn"
   ```

2. **Path Resolver Creation**: A path resolver is created using `ProtosPaths`
   ```python
   self.path_resolver = ProtosPaths()
   ```

3. **Default Paths**: Data paths are constructed following conventions
   ```python
   self.data_root = self.path_resolver.get_data_root()
   self.data_path = os.path.join(data_root, self.processor_type)
   ```

4. **Cross-Platform Handling**: Paths use `pathlib.Path` for cross-platform compatibility
   ```python
   # Instead of: os.path.join("data", "structure", "mmcif")
   # The code uses: Path(data_root) / "structure" / "mmcif"
   ```

5. **Dual-Source Checking**: File operations check both reference and user data sources
   ```python
   # First checks user data (e.g., data/structure/mmcif/1abc.cif)
   # Then checks reference data (e.g., protos/reference_data/structure/mmcif/1abc.cif)
   ```

## Registry System

Protos uses a registry system to track datasets:

- **Global Registry**: Located at `data_root/global_registry.json`, maintains a unified view of all datasets
- **Processor Registries**: Located at `data_root/processor_type/registry.json`, for backward compatibility

The registry maps dataset IDs to file paths and metadata:

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

## Usage

### Using Global Helper Functions

```python
from protos.io.paths import (
    get_structure_path, 
    get_grn_path,
    get_dataset_path,
    DataSource
)

# Get paths (will check both reference and user data)
structure_path = get_structure_path("1abc")
grn_path = get_grn_path("grn_table")

# Explicitly specify source
user_path = get_dataset_path("user_dataset", source=DataSource.USER)
ref_path = get_dataset_path("reference_dataset", source=DataSource.REFERENCE)
```

### Using the Global Registry

```python
from protos.io.data_access import GlobalRegistry
from protos.io.paths import DataSource

# Initialize the registry
registry = GlobalRegistry()

# Register a dataset
registry.register_dataset(
    dataset_id="my_dataset",
    file_path="/path/to/file.dat",
    processor_type="structure",
    dataset_type="cif",
    source=DataSource.USER,
    metadata={"description": "My custom dataset"}
)

# Get dataset path
path = registry.get_dataset_path("my_dataset")

# List datasets by type
structure_datasets = registry.list_datasets("structure")
cif_datasets = registry.get_datasets_by_type("cif")
user_datasets = registry.get_datasets_by_source(DataSource.USER)
```

## Best Practices

1. Use the high-level API when possible, which will check both reference and user data sources
2. For writing data, always use user data paths
3. Use dataset IDs and the registry system instead of hardcoding paths
4. Use the `GlobalRegistry` for cross-processor data relationships
5. Use environment variables for deployment-specific configurations
6. Use `pathlib.Path` for cross-platform path handling instead of string concatenation
7. Always check both reference and user data locations when reading files