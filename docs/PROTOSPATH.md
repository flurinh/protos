# ProtosPaths: Complete Path Management Architecture

## Core Philosophy

**Protos completely abstracts the file system from users.** Users work exclusively with dataset/entity names while Protos handles all aspects of data management internally. This is not just a convenience - it's the fundamental design principle that enables Protos to provide a consistent, professional data analysis framework.

## Overview

ProtosPaths is the internal engine that makes this abstraction possible. It provides a unified, consistent way to handle all file paths, ensuring that users never need to construct paths manually and that all data access goes through processors. ProtosPaths is an implementation detail that users should never need to interact with directly.

## Why This Matters

By completely abstracting the file system, Protos provides:
- **Portability**: Code works identically across different systems and environments
- **Consistency**: All data is accessed the same way regardless of format or location
- **Safety**: No risk of path errors, overwrites, or missing directories
- **Simplicity**: Users focus on analysis, not file management
- **Professionalism**: Clean API that hides implementation complexity

## Core Principles

1. **Single Source of Truth**: ProtosPaths manages all path resolution internally
2. **Name-Based Access**: Users work exclusively with dataset/entity names
3. **Processor-Centric**: All data I/O goes through processor methods
4. **Format Agnostic**: Same interface regardless of underlying file format
5. **Completely Transparent**: File system is 100% hidden from users

## Architecture

### 1. ProtosPaths Class (`protos.io.paths.path_config.py`)

The central path configuration class that:
- Manages the global data root directory
- Provides path resolution for all data types
- Supports environment-based configuration
- Enables test isolation through data root switching

```python
from protos.io.paths.path_config import ProtosPaths

# Set global data root (typically done once at startup or in tests)
ProtosPaths.set_data_root("/path/to/data")

# Get paths instance
paths = ProtosPaths()
```

### 2. Data Directory Structure

```
$PROTOS_DATA_ROOT/
├── grn/                    # GRN data
│   ├── ref/               # Reference tables (mo_ref.csv, gpcrdb_ref.csv)
│   ├── tables/            # Processed GRN tables
│   ├── configs/           # Configuration files
│   └── registry.json      # Dataset registry
├── structure/             # Structure data
│   ├── mmcif/            # PDB/CIF files
│   ├── alignments/       # Structure alignments
│   ├── structure_dataset/# Dataset definitions
│   └── registry.json     # Dataset registry
├── sequence/              # Sequence data
│   ├── fasta/            # FASTA files
│   ├── alignments/       # Sequence alignments
│   └── registry.json     # Dataset registry
├── embedding/             # Embeddings
│   ├── models/           # Cached models
│   ├── embeddings/       # Computed embeddings
│   └── registry.json     # Dataset registry
└── property/              # Properties
    ├── datasets/         # Property datasets
    └── registry.json     # Dataset registry
```

### 3. BaseProcessor Integration

All processors inherit from BaseProcessor which:
- Automatically uses ProtosPaths for all path resolution
- Provides save_data() and load_data() methods
- Manages dataset registries
- Handles format inference

```python
class MyProcessor(BaseProcessor):
    def __init__(self, name="my_processor"):
        super().__init__(
            name=name,
            processor_data_dir="my_data"  # Subdirectory under data root
        )
```

## Usage Patterns

### For Users

Users NEVER construct paths. They only use dataset names:

```python
# ✅ CORRECT
processor = GRNBaseProcessor(name="analysis")
processor.load_grn_table("mo_ref")  # Just the dataset name
processor.save_data("results", df)   # Just the dataset name

# ❌ WRONG
processor.load_data("/path/to/mo_ref.csv")  # NO!
df.to_csv("/path/to/results.csv")           # NO!
```

### For Tests

Tests use the global configuration from conftest.py:

```python
# conftest.py handles this globally
@pytest.fixture(scope="session", autouse=True)
def configure_test_paths():
    """Configure ProtosPaths to use test-data directory for all tests."""
    test_data_root = Path(__file__).parent / "test-data"
    ProtosPaths.set_data_root(str(test_data_root))
    yield
    ProtosPaths.set_data_root(None)
```

Individual tests then just use processors normally:

```python
def test_something():
    # ✅ CORRECT - paths are already configured
    processor = GRNBaseProcessor(name="test")
    processor.load_grn_table("mo_ref")
    
    # ❌ WRONG - don't manipulate paths or environment
    os.environ["PROTOS_DATA_ROOT"] = "/tmp/test"  # NO!
    Path(__file__).parent / "data"                 # NO!
```

### For Processor Developers

When creating a new processor:

1. Inherit from BaseProcessor
2. Specify your processor_data_dir
3. Use self.save_data() and self.load_data()
4. Never construct paths manually

```python
class NewProcessor(BaseProcessor):
    def __init__(self, name="new_processor"):
        super().__init__(
            name=name,
            processor_data_dir="new_data"
        )
    
    def process_something(self, dataset_name):
        # Load data using dataset name
        data = self.load_data(dataset_name)
        
        # Process...
        result = do_processing(data)
        
        # Save using dataset name
        self.save_data(f"{dataset_name}_processed", result)
```

## Registry System

Each processor maintains a registry.json that tracks:
- Available datasets
- File formats
- Metadata
- Last modified times

The registry is automatically updated by save_data() and consulted by load_data().

## Environment Variables

- `PROTOS_DATA_ROOT`: Override default data directory
- Used primarily for testing and CI/CD
- ProtosPaths.set_data_root() takes precedence

## Migration Path

To fix existing code:

1. **Remove all Path() usage**: Replace with processor methods
2. **Remove all os.path usage**: Use dataset names instead
3. **Remove all environment manipulation**: Use ProtosPaths.set_data_root()
4. **Remove all .exists() checks**: Use processor.is_dataset_available()
5. **Remove all direct file I/O**: Use processor.save_data()/load_data()

## Best Practices

1. **Dataset Naming**: Use descriptive, consistent names
   - Good: "mo_ref", "7bmh_aligned", "esm2_embeddings"
   - Bad: "data1", "temp", "output"

2. **Format Inference**: Let the system infer formats from extensions
   - The system supports: csv, json, pkl, npy, npz, fasta

3. **Processor Isolation**: Each processor type has its own subdirectory
   - Prevents naming conflicts
   - Enables clear data organization

4. **Test Data**: Use setup_test_data_from_reference.py to populate test-data
   - Ensures consistent test environment
   - Copies minimal necessary data

## Common Errors and Solutions

### Error: "No such file or directory"
**Cause**: Manual path construction
**Solution**: Use processor.load_data(dataset_name)

### Error: "PROTOS_DATA_ROOT not set"  
**Cause**: Environment not configured
**Solution**: Use ProtosPaths.set_data_root() or rely on defaults

### Error: "Dataset not found"
**Cause**: Dataset not in registry or wrong name
**Solution**: Check processor.list_datasets() for available names

### Error: Tests work locally but fail in CI
**Cause**: Hardcoded paths or missing test data
**Solution**: Use only dataset names and ensure test data setup

## Summary

The ProtosPaths system ensures:
- **Consistency**: All paths handled the same way
- **Portability**: Code works across different environments  
- **Testability**: Easy to switch data roots for testing
- **Simplicity**: Users never deal with paths
- **Safety**: No accidental file overwrites or path errors

By following these patterns, Protos maintains a clean, professional API that handles all path complexity internally while presenting a simple interface to users.