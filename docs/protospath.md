# ProtosPaths System

The ProtosPaths system is the foundation of Protos' file system abstraction. It provides centralized path management, ensuring that users never need to construct or manipulate file paths directly.

## Overview

ProtosPaths implements a singleton pattern that manages all path resolution for the framework. Once configured, all processors automatically use the same data root, ensuring consistency across the entire system.

## Architecture

```
ProtosPaths (Singleton)
    ├── Global Configuration
    │   └── set_data_root()
    ├── Path Resolution
    │   ├── get_data_root()
    │   ├── get_processor_path()
    │   └── resolve_path()
    └── Directory Structure
        ├── structure/
        ├── sequence/
        ├── grn/
        ├── property/
        └── embedding/
```

## Configuration

### Global Setup

ProtosPaths must be configured once at the start of a session:

```python
from protos.io.paths.path_config import ProtosPaths

# Set the global data root
ProtosPaths.set_data_root("/path/to/protos_data")

# All processors will now use this root automatically
```

### Environment Variables

ProtosPaths also supports environment variable configuration:

```bash
export PROTOS_DATA_ROOT=/path/to/protos_data
```

Priority order:
1. Explicit `set_data_root()` call
2. `PROTOS_DATA_ROOT` environment variable
3. Default location (`~/protos_data`)

## Directory Structure

ProtosPaths maintains a standardized directory hierarchy:

```
protos_data/
├── entity_registry.json      # Global entity registry
├── structure/               # Structure data (CifBaseProcessor)
│   ├── mmcif/              # PDB/CIF files
│   ├── pdb/                # Legacy PDB format
│   ├── alignments/         # Structure alignments
│   ├── datasets/           # Dataset definitions
│   └── registry.json       # Structure-specific registry
├── sequence/               # Sequence data (SeqProcessor)
│   ├── fasta/              # FASTA files
│   ├── processed/          # Processed sequences
│   ├── alignments/         # Sequence alignments
│   ├── datasets/           # Dataset definitions
│   └── registry.json       # Sequence-specific registry
├── grn/                    # GRN data (GRNBaseProcessor)
│   ├── tables/             # GRN annotation tables
│   ├── configs/            # GRN configurations
│   ├── ref/                # Reference GRN tables
│   └── registry.json       # GRN-specific registry
├── property/               # Property data (PropertyProcessor)
│   ├── datasets/           # Property datasets
│   └── registry.json       # Property-specific registry
└── embedding/              # Embeddings (EmbeddingProcessor)
    ├── models/             # Embedding models
    ├── embeddings/         # Generated embeddings
    └── registry.json       # Embedding-specific registry
```

## Usage Patterns

### Automatic Path Resolution

Processors automatically resolve paths without user intervention:

```python
# User provides name
processor.load_structure("1ubq")

# ProtosPaths resolves internally:
# 1. Get data root: /path/to/protos_data
# 2. Add processor directory: /path/to/protos_data/structure
# 3. Add file type: /path/to/protos_data/structure/mmcif
# 4. Resolve entity: /path/to/protos_data/structure/mmcif/3f8a9c2d1e.cif
```

### Processor Integration

All processors inherit path management from BaseProcessor:

```python
class CifBaseProcessor(BaseProcessor):
    def __init__(self, name="structures"):
        super().__init__(name, processor_type="structure")
        # self.data_path automatically set by ProtosPaths
        # No manual path configuration needed
```

### Path Properties

Processors expose convenient path properties:

```python
processor = CifBaseProcessor()

# Available path properties
processor.data_path          # /path/to/protos_data/structure
processor.path_mmcif         # /path/to/protos_data/structure/mmcif
processor.path_datasets      # /path/to/protos_data/structure/datasets
processor.path_registry      # /path/to/protos_data/structure/registry.json
```

## Advanced Features

### Custom Directory Creation

ProtosPaths can create custom subdirectories:

```python
# Create custom analysis directory
analysis_path = ProtosPaths.get_processor_path("structure", "analysis")
# Creates: /path/to/protos_data/structure/analysis/
```

### Path Validation

ProtosPaths validates paths and creates directories as needed:

```python
# Automatic directory creation
paths = ProtosPaths()
paths.ensure_directories()  # Creates all required directories
```

### Cross-Platform Support

ProtosPaths handles platform differences automatically:

```python
# Windows path
# C:\Users\name\protos_data\structure\mmcif\file.cif

# Linux/Mac path  
# /home/name/protos_data/structure/mmcif/file.cif

# User code remains the same
processor.load_structure("1ubq")
```

## Configuration Examples

### Research Project Setup

```python
# research_project.py
import os
from pathlib import Path
from protos.io.paths.path_config import ProtosPaths

# Configure project-specific data location
project_dir = Path(__file__).parent
data_dir = project_dir / "research_data"

# Set globally for all processors
ProtosPaths.set_data_root(str(data_dir))

# Now import and use processors
from protos.processing.structure.struct_base_processor import CifBaseProcessor
processor = CifBaseProcessor()  # Uses research_data automatically
```

### Testing Configuration

```python
# conftest.py (pytest configuration)
import pytest
from pathlib import Path
from protos.io.paths.path_config import ProtosPaths

@pytest.fixture(scope="session", autouse=True)
def configure_test_paths():
    """Configure ProtosPaths for testing."""
    test_data = Path(__file__).parent / "test-data"
    ProtosPaths.set_data_root(str(test_data))
    yield
    # Cleanup if needed
```

### Multi-User Environment

```python
# Respect user preferences
import os
from pathlib import Path
from protos.io.paths.path_config import ProtosPaths

# Check for user override
user_data = os.environ.get("PROTOS_USER_DATA")
if user_data:
    data_root = user_data
else:
    # Default to user home
    data_root = Path.home() / "protos_data"

ProtosPaths.set_data_root(str(data_root))
```

## Best Practices

### 1. Configure Once

Set the data root once at application startup:

```python
# main.py or __init__.py
from protos.io.paths.path_config import ProtosPaths
ProtosPaths.set_data_root("/path/to/data")

# All subsequent imports will use this configuration
```

### 2. Never Construct Paths

Let ProtosPaths handle all path operations:

```python
# ❌ WRONG - Manual path construction
file_path = Path(data_dir) / "structure" / "mmcif" / "1ubq.cif"
with open(file_path) as f:
    data = f.read()

# ✅ CORRECT - Use processor methods
structure = processor.load_structure("1ubq")
```

### 3. Use Path Properties

Access paths through processor properties:

```python
# Get paths when needed
mmcif_dir = processor.path_mmcif
dataset_dir = processor.path_datasets

# But prefer processor methods
datasets = processor.list_datasets()  # Better than os.listdir(dataset_dir)
```

### 4. Handle Missing Directories

ProtosPaths creates directories automatically, but verify in production:

```python
from protos.io.paths.path_config import ProtosPaths

# Ensure all directories exist
paths = ProtosPaths()
paths.ensure_directories()

# Or check specific directory
if not processor.path_mmcif.exists():
    processor.path_mmcif.mkdir(parents=True)
```

## Troubleshooting

### Common Issues

1. **"Data root not set" error**
   ```python
   # Solution: Configure before importing processors
   from protos.io.paths.path_config import ProtosPaths
   ProtosPaths.set_data_root("/path/to/data")
   ```

2. **Permission errors**
   ```python
   # Ensure write permissions
   data_root = ProtosPaths.get_data_root()
   os.access(data_root, os.W_OK)  # Check write permission
   ```

3. **Path not found**
   ```python
   # Verify path exists and create if needed
   if not processor.data_path.exists():
       processor.data_path.mkdir(parents=True, exist_ok=True)
   ```

### Debug Path Resolution

```python
# Debug path resolution
from protos.io.paths.path_config import ProtosPaths

paths = ProtosPaths()
print(f"Data root: {paths.get_data_root()}")
print(f"Structure path: {paths.get_processor_path('structure')}")

# Check processor paths
processor = CifBaseProcessor()
print(f"Processor data path: {processor.data_path}")
print(f"MMCIF path: {processor.path_mmcif}")
```

## Implementation Details

### Singleton Pattern

ProtosPaths uses a singleton to ensure consistent configuration:

```python
class ProtosPaths:
    _instance = None
    _data_root = None
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance
    
    @classmethod
    def set_data_root(cls, path: str):
        """Set global data root."""
        cls._data_root = Path(path).absolute()
```

### Path Resolution Logic

```python
def resolve_entity_path(self, entity_name: str, format_type: str) -> Path:
    """Resolve entity name to file path."""
    # 1. Generate entity ID
    entity_id = generate_entity_id(entity_name)
    
    # 2. Get processor path
    processor_path = self.get_processor_path(format_type)
    
    # 3. Determine file extension
    extension = self._get_extension(format_type)
    
    # 4. Construct full path
    return processor_path / f"{entity_id}{extension}"
```

## Summary

ProtosPaths provides:
- Centralized path configuration
- Automatic path resolution
- Cross-platform compatibility
- Standardized directory structure
- Seamless processor integration

By abstracting file system operations, ProtosPaths enables researchers to focus on their data and analysis rather than file management.