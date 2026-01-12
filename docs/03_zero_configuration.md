# Zero-Configuration / ProtosPaths

Protos provides automatic path management through the ProtosPaths system. All components work out of the box with no configuration required, while still allowing customization when needed.

## Design Principles

1. **Zero Configuration** - Works immediately with sensible defaults
2. **Lazy Initialization** - Directories created only when first accessed
3. **Singleton Pattern** - Single shared instance across all components
4. **Automatic Structure** - Standard directory layout for all data types

---

## Default Behavior

By default, Protos stores all data in `~/protos_data`:

```python
from protos import StructureProcessor

# No configuration needed!
proc = StructureProcessor()  # Uses ~/protos_data automatically
```

The directory structure is created automatically on first use.

---

## Customizing the Data Path

### Environment Variable

Set `PROTOS_DATA_ROOT` to override the default location:

```bash
export PROTOS_DATA_ROOT=/path/to/my/data
```

```python
# Now uses /path/to/my/data
proc = StructureProcessor()
```

### Programmatic Configuration

Use `set_data_path()` **before** creating any processors:

```python
from protos import set_data_path

# MUST be called before creating processors!
set_data_path("/path/to/my/data")

# Now processors will use the custom path
from protos import StructureProcessor
proc = StructureProcessor()
```

**Important:** Once any processor has been created, the path cannot be changed:

```python
from protos import StructureProcessor, set_data_path

proc = StructureProcessor()  # Initializes paths

set_data_path("/new/path")  # RuntimeError! Too late to change
```

---

## Directory Structure

Protos creates this directory layout automatically:

```
~/protos_data/                    # Data root
├── global_registry.json          # Entity registry
├── .protos_initialized           # Initialization marker
│
├── structure/                    # Structure processor
│   ├── mmcif/                   # Raw CIF files
│   ├── cache/                   # Processed structures (PKL)
│   ├── datasets/                # Dataset definitions
│   ├── sdf/                     # Ligand SDF files
│   ├── pdb/                     # Exported PDB files
│   ├── alignments/              # Structure alignments
│   └── temp_cif/                # Temporary files
│
├── sequence/                     # Sequence processor
│   ├── fasta/
│   │   ├── entities/            # Single-sequence FASTA
│   │   └── datasets/            # Multi-sequence FASTA
│   ├── alignments/
│   │   ├── pairwise/
│   │   ├── multiple/
│   │   └── mmseqs/
│   ├── metadata/
│   └── datasets/
│
├── grn/                          # GRN processor
│   ├── tables/                  # GRN mapping tables
│   ├── reference/               # Reference data (bundled)
│   ├── configs/                 # Configuration files
│   ├── assignments/             # GRN assignments
│   └── datasets/
│
├── embedding/                    # Embedding processor
│   ├── embeddings/              # Saved embeddings
│   └── datasets/
│
├── property/                     # Property processor
│   ├── tables/                  # Property tables (CSV)
│   └── datasets/
│
├── molecule/                     # Molecule processor
│   ├── records/                 # Molecule JSON records
│   └── datasets/
│
├── graph/                        # Graph processor
│   ├── graphs/                  # Graph representations
│   ├── analysis/                # Analysis results
│   └── datasets/
│
├── input/                        # Input staging area
│
└── temp/                         # Temporary files
```

---

## ProtosPaths API

### Getting the Instance

```python
from protos.io.paths import get_protos_paths

paths = get_protos_paths()  # Returns singleton

# Or with custom root (only works before initialization)
paths = get_protos_paths("/custom/path")
```

### Core Methods

| Method | Description |
|--------|-------------|
| `get_processor_path(type)` | Get processor directory |
| `get_subdir_path(type, subdir)` | Get subdirectory within processor |
| `get_dataset_path(type, name)` | Get dataset definition file path |
| `get_global_registry_path()` | Get registry file path |
| `resolve_path(path)` | Resolve relative path to absolute |
| `exists(path)` | Check if path exists |

### Examples

```python
from protos.io.paths import get_protos_paths

paths = get_protos_paths()

# Get processor directories
struct_dir = paths.get_processor_path("structure")
seq_dir = paths.get_processor_path("sequence")

# Get subdirectories
cache_dir = paths.get_subdir_path("structure", "cache_dir")
fasta_dir = paths.get_subdir_path("sequence", "entity_fasta_dir")

# Get dataset path
dataset_file = paths.get_dataset_path("structure", "my_dataset")
# ~/protos_data/structure/datasets/my_dataset.json

# Get registry path
registry_file = paths.get_global_registry_path()
# ~/protos_data/global_registry.json
```

---

## Subdirectory Keys

### Structure Subdirectories

| Key | Directory | Purpose |
|-----|-----------|---------|
| `structure_dir` | `mmcif/` | Raw CIF structure files |
| `cache_dir` | `cache/` | Processed structures (PKL) |
| `dataset_dir` | `structure_dataset/` | Dataset PKL files |
| `datasets_dir` | `datasets/` | Dataset JSON definitions |
| `sdf_dir` | `sdf/` | Ligand SDF files |
| `pdb_dir` | `pdb/` | Exported PDB files |
| `alignments_dir` | `alignments/` | Alignment results |
| `temp_dir` | `temp_cif/` | Temporary files |

### Sequence Subdirectories

| Key | Directory | Purpose |
|-----|-----------|---------|
| `entity_fasta_dir` | `fasta/entities/` | Single-sequence FASTA |
| `dataset_fasta_dir` | `fasta/datasets/` | Multi-sequence FASTA |
| `alignment_dir` | `alignments/` | General alignments |
| `pairwise_alignment_dir` | `alignments/pairwise/` | Pairwise alignments |
| `multiple_alignment_dir` | `alignments/multiple/` | MSA results |
| `mmseqs_alignment_dir` | `alignments/mmseqs/` | MMseqs2 results |
| `database_dir` | `databases/` | Sequence databases |
| `metadata_dir` | `metadata/` | Sequence metadata |
| `datasets_dir` | `datasets/` | Dataset definitions |

### GRN Subdirectories

| Key | Directory | Purpose |
|-----|-----------|---------|
| `table_dir` | `tables/` | GRN mapping tables |
| `reference_dir` | `reference/` | Reference GRN data |
| `configs_dir` | `configs/` | Configuration files |
| `assignment_dir` | `assignments/` | GRN assignments |
| `datasets_dir` | `datasets/` | Dataset definitions |

### Other Processors

| Processor | Key | Directory |
|-----------|-----|-----------|
| Embedding | `embeddings_dir` | `embeddings/` |
| Property | `tables_dir` | `tables/` |
| Molecule | `records_dir` | `records/` |
| Graph | `graphs_dir` | `graphs/` |

---

## Path Helper Functions

The paths module provides convenience functions:

```python
from protos.io.paths import (
    get_structure_path,
    get_sequence_path,
    get_grn_path,
    sanitize_storage_name,
)

paths = get_protos_paths()

# Get structure file path
cif_path = get_structure_path(paths, "1ubq")
# ~/protos_data/structure/mmcif/1ubq.cif

# Get sequence file path
fasta_path = get_sequence_path(paths, "protein_a")
# ~/protos_data/sequence/fasta/protein_a.fasta

# Get GRN table path
grn_path = get_grn_path(paths, "receptor_numbering")
# ~/protos_data/grn/tables/receptor_numbering.csv

# Sanitize names for filesystem
safe_name = sanitize_storage_name("My Protein (variant)")
# "My_Protein_variant"
```

---

## Resetting Data

To reset the data directory:

```python
from protos.io.paths import reset_protos_data

# Reset with backup (default)
reset_protos_data()

# Reset without reinstalling reference data
reset_protos_data(reinstall_reference=False)

# Complete wipe and recreate
reset_protos_data(wipe=True)

# Reset without registry backup
reset_protos_data(backup_registry=False)
```

### Reinitialize ProtosPaths

```python
paths = get_protos_paths()

# Reinitialize directory structure
paths.reinitialize()

# Complete wipe and reinstall
paths.reinitialize(wipe=True, reinstall_reference=True)
```

---

## Reference Data

Protos includes bundled reference data (e.g., GRN configurations) that is automatically installed on first use:

```
~/protos_data/grn/
├── configs/
│   └── config.json       # GRN configuration
└── reference/
    └── ...               # Reference GRN tables
```

The `.protos_initialized` marker file tracks whether reference data has been installed.

---

## Cross-Platform Support

ProtosPaths handles platform differences automatically:

- Uses `os.path.normpath()` for cross-platform paths
- Handles path separators correctly on Windows/Unix
- Expands `~` to user home directory
- Supports WSL path conversion

---

## Best Practices

### 1. Set Path Early

```python
# At the start of your script
from protos import set_data_path
set_data_path("/my/data/path")

# Then import and use processors
from protos import StructureProcessor
proc = StructureProcessor()
```

### 2. Use Environment Variable for Deployment

```bash
# In production/deployment
export PROTOS_DATA_ROOT=/data/protos

# Code works without modification
python my_script.py
```

### 3. Access Paths Through Processors

```python
proc = StructureProcessor()

# Use processor's path properties
cif_dir = proc.path_cif_dir
cache_dir = proc.path_pkl_dir
```

### 4. Let Protos Manage Paths

```python
# DON'T construct paths manually
# bad_path = "/home/user/protos_data/structure/mmcif/1ubq.cif"

# DO use ProtosPaths
paths = get_protos_paths()
good_path = get_structure_path(paths, "1ubq")
```

---

## Troubleshooting

### "Cannot change ProtosPaths after initialization"

This error occurs when trying to set the data path after processors have been created:

```python
# Wrong order
proc = StructureProcessor()
set_data_path("/new/path")  # Error!

# Correct order
set_data_path("/new/path")
proc = StructureProcessor()  # OK
```

### Missing Directories

Directories are created lazily. If a directory seems missing:

```python
paths = get_protos_paths()
paths.reinitialize()  # Recreates all directories
```

### Permission Issues

Ensure the data directory has appropriate permissions:

```bash
chmod -R 755 ~/protos_data
```
