# Path Configuration in Protos

## Overview

Protos implements a simplified path management system with a single data directory that contains all project data. This unified approach ensures consistency and simplifies data management across all processors and components.

This document explains how the path system works and how to configure it for your needs.

## Data Directory

### Default Location

The Protos data directory is automatically initialized in one of these locations (in order of precedence):

1. **Environment Variable**: If `PROTOS_DATA_ROOT` is set, that location is used
2. **Home Directory**: `~/protos_data/` (default location)
3. **Custom Path**: When explicitly specified via `ProtosPaths` initialization

### Automatic Initialization

The data directory is automatically created when:
- Any processor is instantiated
- The registry system is accessed
- `ProtosPaths` is initialized

This ensures that the data directory always exists when needed, without requiring manual setup.

### Data Organization

All data is organized within the single data directory:

```
~/protos_data/                # Default location
├── structure/               # Structure-related data
│   ├── mmcif/              # Structure files
│   ├── alignments/         # Alignment results
│   └── structure_dataset/   # Structure datasets
├── grn/                    # GRN-related data
│   ├── ref/                # Reference GRN tables
│   ├── tables/             # Calculated GRN tables  
│   └── configs/            # GRN configurations
├── sequence/               # Sequence-related data
│   ├── fasta/              # FASTA files
│   └── alignments/         # Sequence alignments
└── global_registry.json    # Global dataset registry
```

## Configuration

### Environment Variables

Protos uses a single environment variable for path configuration:

- `PROTOS_DATA_ROOT`: Location for the data directory (default: `~/protos_data/`)

### Programmatic Configuration

Paths can be configured programmatically:

```python
from protos.io.paths import ProtosPaths

# Use default location (~protos_data/)
paths = ProtosPaths()

# Initialize with custom path
paths = ProtosPaths(
    data_root="/path/to/my/data",
    create_dirs=True,  # Create directories if they don't exist (default: True)
    validate=True      # Validate directory structure (default: True)
)

# Use the paths object
structure_path = paths.get_processor_path("structure")
grn_path = paths.get_grn_subdir_path("tables")
```

### Automatic Initialization by Processors

Processors automatically initialize paths if not already configured:

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# No need to manually configure paths - processor handles it
processor = CifBaseProcessor(name="my_processor")
# Data directory is automatically created at ~/protos_data/

# Or specify a custom location
processor = CifBaseProcessor(
    name="my_processor",
    data_root="/my/project/data"
)
```

## Processor Path Resolution

Processor classes automatically handle path resolution:

1. **Automatic Initialization**: When a processor is instantiated without a data_root, it uses the default location
   ```python
   # Automatically uses ~/protos_data/
   processor = CifBaseProcessor(name="my_processor")
   ```

2. **Type Detection**: Processors automatically determine their subdirectory from the class name
   ```python
   # CifBaseProcessor → "structure" subdirectory
   # GRNBaseProcessor → "grn" subdirectory
   # SeqProcessor → "sequence" subdirectory
   ```

3. **Directory Creation**: Required directories are automatically created
   ```python
   # These directories are created automatically:
   # ~/protos_data/structure/mmcif/
   # ~/protos_data/structure/alignments/
   # ~/protos_data/structure/structure_dataset/
   ```

4. **Cross-Platform Support**: All paths use `pathlib.Path` for compatibility
   ```python
   # Works on Windows, Linux, and macOS
   data_path = Path.home() / "protos_data" / "structure"
   ```

## Registry System

The registry system tracks all datasets in the data directory:

- **Global Registry**: Located at `data_root/global_registry.json`
- **Processor Registries**: Located at `data_root/processor_type/registry.json` (for backward compatibility)
- **Automatic Registration**: Datasets are automatically registered when created through processors

The registry maps dataset IDs to file paths and metadata:

```json
{
  "dataset_id": {
    "path": "/home/user/protos_data/structure/mmcif/1abc.cif",
    "metadata": {
      "processor_type": "structure",
      "dataset_type": "cif",
      "description": "Crystal structure of protein ABC"
    },
    "timestamp": "2023-01-01T00:00:00"
  }
}
```

## Usage Examples

### Basic Usage

```python
# Import a processor - paths are handled automatically
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Create processor - data directory is initialized automatically
processor = CifBaseProcessor(name="my_analysis")

# Load a structure
processor.load_structure("1abc")

# Work with data - all paths are resolved automatically
processor.save_dataset("my_results")
```

### Custom Data Location

```python
# Specify a custom data directory
processor = CifBaseProcessor(
    name="my_analysis",
    data_root="/my/project/data"
)

# Or set via environment variable before importing
import os
os.environ["PROTOS_DATA_ROOT"] = "/my/project/data"

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
grn_processor = GRNBaseProcessor(name="grn_analysis")
```

### Using the Registry

```python
from protos.io.data_access import GlobalRegistry

# Registry automatically uses the configured data directory
registry = GlobalRegistry()

# List all datasets
all_datasets = registry.list_datasets()

# Find specific types
structure_datasets = registry.list_datasets("structure")

# Get dataset information
info = registry.get_dataset_info("my_dataset")
print(f"Dataset location: {info['path']}")
print(f"Dataset type: {info['metadata']['dataset_type']}")
```

## Best Practices

1. **Let processors handle paths automatically** - They will use the correct default location
2. **Use dataset IDs** instead of hardcoding file paths
3. **Set `PROTOS_DATA_ROOT` once** at the start of your project if using a custom location
4. **Use the registry system** to track and discover datasets
5. **Avoid manual path construction** - Use the processor and registry APIs
6. **Single data directory** - All data (reference tables, calculated results, etc.) lives in one place

## Migration from Dual-Path System

If you're migrating from the old dual-path system (separate reference and user data):

1. **Copy reference data** to the new data directory structure:
   ```bash
   cp -r /old/reference/data/grn/ref/* ~/protos_data/grn/ref/
   ```

2. **Update environment variables**:
   ```bash
   # Old system
   export PROTOS_DATA_ROOT=/path/to/user/data
   export PROTOS_REF_DATA_ROOT=/path/to/reference/data
   
   # New system - only one variable
   export PROTOS_DATA_ROOT=~/protos_data
   ```

3. **Update code** - Remove DataSource specifications:
   ```python
   # Old code
   path = paths.get_processor_path("structure", DataSource.USER)
   
   # New code
   path = paths.get_processor_path("structure")
   ```