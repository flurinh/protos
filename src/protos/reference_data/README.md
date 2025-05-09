# Protos Reference Data

This directory contains reference data distributed with the Protos package. This data is essential for the framework to function properly and includes:

- Structure files (mmCIF)
- GRN tables (CSV)
- Reference FASTA sequences
- Standard datasets

## Directory Structure

The reference data follows the same directory structure as user data:

```
reference_data/
   structure/
      mmcif/
      alignments/
      structure_dataset/
   grn/
      tables/
      configs/
      assignments/
   sequence/
       fasta/
       alignments/
       metadata/
```

## Usage

Reference data is automatically loaded and made available through the `GlobalRegistry`. When requesting data through the Protos API, the system will first check for user-defined datasets in the user data directory, and if not found, will check the reference data.

## Extending

When distributing custom reference data, place it in the appropriate subdirectory and ensure that a registry.json file is present in each processor type directory to enable the framework to locate and use the data.

## Registry Format

Each processor directory should contain a registry.json file with the following format:

```json
{
  "dataset_id": {
    "path": "relative/path/to/file",
    "metadata": {
      "dataset_type": "type",
      "description": "description",
      "additional_metadata": "value"
    }
  }
}
```

The paths in the registry can be relative to the processor directory or absolute.