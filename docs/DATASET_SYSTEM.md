# Dataset System in Protos

This document describes the standardized dataset system implemented in the Protos framework.

## Overview

Datasets in Protos are collections of related data items (e.g., structures, sequences, properties) with associated metadata. The dataset system provides a consistent interface for working with datasets across different processor types.

## Dataset Class

The core of the dataset system is the `Dataset` class, which provides:

- Standard metadata fields (id, name, description, etc.)
- Content storage for different types (lists, dictionaries, sets)
- Serialization and deserialization
- Methods for adding, removing, and updating items
- Support for custom metadata

## Integration with Processors

Each processor type (structure, GRN, sequence, etc.) can work with datasets using:

1. The legacy API in `BaseProcessor`:
   - `load_data()`, `save_data()`, `list_datasets()`, etc.

2. The new standardized API:
   - `create_standard_dataset()`, `load_standard_dataset()`, etc.
   - Access via `processor.dataset_manager`

## Dataset Manager

The `DatasetManager` class provides centralized dataset operations:

- Creating and loading datasets
- Saving and deleting datasets
- Listing available datasets
- Dataset conversion and merging
- Registry integration

## Directory Structure

Datasets are stored in a consistent directory structure:

```
<processor_root>/
  ├── datasets/
  │   ├── dataset1.json
  │   ├── dataset2.json
  │   └── ...
  └── registry.json
```

## Dataset Registry

The registry system tracks all available datasets:

- Each processor has its own registry
- The global registry provides cross-processor discoverability
- Registries store metadata and file locations

## Creating a Dataset

To create a new dataset:

```python
from protos.processing.structure import CifProcessor

# Initialize processor
processor = CifProcessor(name="my_processor")

# Create dataset from a list of PDB IDs
dataset = processor.create_standard_dataset(
    dataset_id="my_dataset",
    name="My Structure Dataset",
    description="A collection of important structures",
    content=["1abc", "2xyz", "3def"],
    metadata={"category": "important", "author": "John Doe"}
)
```

## Working with Datasets

```python
# Load an existing dataset
dataset = processor.load_standard_dataset("my_dataset")

# Access dataset content
for pdb_id in dataset.content:
    print(f"Processing {pdb_id}")

# Add a new item
dataset.add_item("4ghi")

# Save changes
processor.save_standard_dataset(dataset)

# List available datasets
available_datasets = processor.list_standard_datasets()
for info in available_datasets:
    print(f"{info['id']}: {info['name']}")
```

## Converting Legacy Datasets

To convert legacy datasets to the new format:

```python
# Convert a legacy dataset
converted_dataset = processor.convert_to_standard_dataset(
    legacy_dataset_id="old_dataset",
    new_dataset_id="new_dataset"
)
```

## Merging Datasets

```python
# Merge multiple datasets
merged_dataset = processor.dataset_manager.merge_datasets(
    dataset_ids=["dataset1", "dataset2", "dataset3"],
    new_dataset_id="merged_dataset",
    name="Merged Dataset",
    description="A combined dataset from multiple sources"
)
```

## Dataset Iteration

```python
# Iterate over a dataset
dataset = processor.load_standard_dataset("my_dataset")
for item in dataset:
    print(item)

# Check if an item is in the dataset
if "1abc" in dataset:
    print("1abc is in the dataset")
```

## Custom Metadata

Datasets can store arbitrary metadata:

```python
# Add custom metadata
dataset.update_metadata({
    "source": "PDB",
    "date_created": "2023-01-01",
    "tags": ["membrane", "protein", "structure"],
    "custom_field": {
        "nested": "value",
        "array": [1, 2, 3]
    }
})
```

## Best Practices

1. Use descriptive names and detailed descriptions for datasets
2. Include relevant metadata for searchability
3. Use the standardized API for new code
4. Convert legacy datasets when possible
5. Use consistent dataset IDs across the project