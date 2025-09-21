# Protos Data Management: Zero-Configuration Design

## Overview

Protos provides a zero-configuration data management system that works immediately after installation. The system is built on three core components that work together seamlessly:

1. **ProtosPaths**: Automatic path management (no setup required)
2. **EntityRegistry**: Transparent entity tracking (users see only human names)  
3. **Processors**: Format-specific interfaces (the only component users interact with)

## Core Principle: Zero Configuration

```python
# This is ALL you need - no setup, no configuration, no paths
from protos.processing.structure import StructureProcessor

processor = StructureProcessor()
processor.load_structure("1ubq")  # Works immediately!
```

Everything else happens automatically behind the scenes.

## Architecture

### Component Hierarchy

```
User Code
    ↓
Processors (User Interface)
    ↓
ProtosPaths (Automatic)
    ↓
EntityRegistry (Invisible)
    ↓
File System (Human-Readable)
```

### Key Design Decisions

1. **Default Data Location**: `./data` in current working directory
2. **No Environment Variables Required**: Works out of the box
3. **No Configuration Files**: Everything has sensible defaults
4. **Human-Readable Everywhere**: No hash IDs in filenames or user-facing APIs
5. **UUID-Based Registry**: Stable entity tracking with human-friendly interface

## ProtosPaths: Zero-Configuration Path Management

ProtosPaths automatically manages all file paths without any user configuration.

### Automatic Initialization

```python
class ProtosPaths:
    """Zero-configuration path management."""
    
    def __init__(self, base_path: Optional[Path] = None):
        # Automatic default: ./data in current directory
        self.base_path = Path(base_path or "data").absolute()
        
        # Create base directory if needed
        self.base_path.mkdir(exist_ok=True)
        
        # Standard directory structure
        self.structure = {
            "structure": ["mmcif", "cache", "datasets", "alignments"],
            "sequence": ["fasta", "alignments", "databases", "datasets"],
            "grn": ["tables", "ref", "configs", "datasets", "assignments"],
            "property": ["tables", "datasets", "cache"],
            "embedding": ["embeddings", "models", "datasets"],
            "ligand": ["sdf", "cache", "datasets"],
            "graph": ["networks", "analysis", "datasets"],
            "input": ["pending", "processed", "rejected"],
            "temp": []
        }
        
        # Auto-create all directories
        self._ensure_directories()
```

### Directory Structure (Created Automatically)

```
./data/                          # Auto-created in working directory
├── entity_registry.json         # Central entity tracking (UUID-based)
├── structure/                   # Structure data
│   ├── mmcif/                  # CIF/PDB files (human names: 1ubq.cif)
│   ├── cache/                  # Processed structures (1ubq.pkl)
│   ├── datasets/               # Dataset definitions (kinases.json)
│   └── alignments/             # Structure alignments
├── sequence/                   # Sequence data
│   ├── fasta/                  # FASTA files (EGFR_HUMAN.fasta)
│   ├── alignments/             # MSAs and pairwise alignments
│   ├── databases/              # MMseqs2 databases
│   └── datasets/               # Dataset definitions
├── grn/                        # Generic Residue Numbering
│   ├── tables/                 # GRN alignment tables (CSV)
│   ├── ref/                    # Reference GRN tables
│   ├── configs/                # GRN configurations
│   ├── datasets/               # Dataset definitions
│   └── assignments/            # Assignment results
├── property/                   # Experimental properties
│   ├── tables/                 # Property data (CSV)
│   ├── datasets/               # Dataset definitions
│   └── cache/                  # Property cache
├── input/                      # Safe input workflow
│   ├── pending/                # Files awaiting processing
│   ├── processed/              # Successfully processed files
│   └── rejected/               # Failed validation
└── temp/                       # Temporary files
```

## Entity Registry: Invisible to Users

The EntityRegistry uses UUIDs internally but presents only human-readable names to users.

### Internal Structure (Users Never See This)

```json
{
  "entities": {
    "550e8400-e29b-41d4-a716-446655440000": {  // UUID (internal only)
      "entity_name": "1ubq",                    // Current primary name
      "aliases": ["UBIQ_HUMAN", "P62988"],      // Alternative names
      "formats": {
        "structure": {
          "file_path": "structure/mmcif/1ubq.cif",  // Human-readable path
          "metadata": {
            "resolution": 1.8,
            "method": "X-RAY",
            "chains": ["A"]
          },
          "created": "2024-01-15T10:30:00"
        }
      },
      "created": "2024-01-15T10:30:00",
      "modified": "2024-01-20T14:22:00"
    }
  },
  "name_index": {  // Fast name-to-UUID lookup
    "1ubq": "550e8400-e29b-41d4-a716-446655440000",
    "UBIQ_HUMAN": "550e8400-e29b-41d4-a716-446655440000",
    "P62988": "550e8400-e29b-41d4-a716-446655440000"
  }
}
```

### What Users See

```python
# Users ONLY work with human names
processor.list_entities()
# Returns: ['1ubq', '2gb1', 'EGFR_HUMAN', 'my_protein']

processor.load_structure("1ubq")  # Human name
processor.entity_exists("EGFR_HUMAN")  # Human name
processor.find_entity("P62988")  # Returns: "1ubq" (primary name)
```

## Processors: The Only User Interface

Users interact ONLY with processors. Everything else is automatic.

### Base Processor Pattern

```python
class BaseProcessor:
    """All processors inherit this zero-config base."""
    
    def __init__(self, name: str = None, paths: ProtosPaths = None):
        # Automatic setup cascade
        self.paths = paths or ProtosPaths()  # Auto-creates ./data
        self.entity_registry = EntityRegistry(self.paths)  # Auto-loads registry
        self.processor_type = self._get_processor_type()
        self.data_path = self.paths.get_processor_path(self.processor_type)
        
        # Everything ready to use!
```

### Processor Examples

#### Structure Processor

```python
from protos.processing.structure import StructureProcessor

# Zero configuration required
processor = StructureProcessor()

# Load structure (automatic path resolution)
structure = processor.load_structure("1ubq")

# Save structure (automatic registration)
processor.save_structure("my_protein", structure_data)

# Create dataset
processor.create_dataset("kinases", ["EGFR", "ABL1", "SRC"])

# Load dataset
structures = processor.load_dataset("kinases")
for name, structure in structures.items():
    print(f"{name}: {len(structure)} atoms")
```

#### Sequence Processor

```python
from protos.processing.sequence import SequenceProcessor

# Zero configuration
processor = SequenceProcessor()

# Work with sequences
sequence = processor.load_sequence("EGFR_HUMAN")
processor.save_sequence("my_sequence", "MTEYKLVVVG...")

# Align sequences
alignment = processor.align_sequences(["seq1", "seq2", "seq3"])
```

#### GRN Processor

```python
from protos.processing.grn import GRNProcessor

# Zero configuration
processor = GRNProcessor()

# Load GRN table (each row is an entity)
grn_table = processor.load_grn_table("gpcr_alignment")

# Get GRN for specific protein
grn_positions = processor.get_entity_grn("BACR_HALSA")
# Returns: {"1.50": "R", "2.50": "L", "3.50": "V", ...}

# Save new GRN table
processor.save_grn_table("my_family", grn_dataframe)
```

#### Property Processor

```python
from protos.processing.property import PropertyProcessor

# Zero configuration
processor = PropertyProcessor()

# Assign properties to any entity
processor.assign_property("BACR_HALSA", "lambda_max", 568)
processor.assign_property("1ubq", "experiment_date", "2024-01-15")

# Get entity properties
props = processor.get_entity_properties("BACR_HALSA")
# Returns: {"lambda_max": 568, "photocycle": "fast", ...}

# Save property dataset
processor.save_property_table("opsin_properties", properties_df)
```

## Common Workflows

### 1. Drop Files and Go

```bash
# User drops files into data directory
cp ~/downloads/*.cif ./data/structure/mmcif/
cp ~/sequences/*.fasta ./data/sequence/fasta/
```

```python
# Immediately available in code
processor = StructureProcessor()
processor.list_entities()  # Shows all CIF files
```

### 2. Cross-Format Entity Tracking

```python
# Load structure
struct_proc = StructureProcessor()
structure = struct_proc.load_structure("1ubq")

# Extract and save sequence
seq_proc = SequenceProcessor()
sequence = struct_proc.extract_sequence(structure, chain="A")
seq_proc.save_sequence("1ubq_A", sequence)

# Both tracked as same entity
# Registry knows 1ubq exists in both structure and sequence formats
```

### 3. Property Assignment

```python
# Assign properties to any entity type
prop_proc = PropertyProcessor()

# To structures
prop_proc.assign_property("1ubq", "quality", "high")
prop_proc.assign_property("2gb1", "quality", "medium")

# To sequences
prop_proc.assign_property("EGFR_HUMAN", "organism", "human")
prop_proc.assign_property("ABL1_MOUSE", "organism", "mouse")

# Query by property
high_quality = prop_proc.find_entities_with_property("quality", "high")
human_proteins = prop_proc.find_entities_with_property("organism", "human")
```

### 4. Dataset Creation and Sharing

```python
# Create dataset
processor = StructureProcessor()
processor.create_dataset("covid_proteins", 
    ["6M0J", "6LU7", "6W63"],
    metadata={
        "description": "SARS-CoV-2 protein structures",
        "created_by": "researcher_name",
        "purpose": "drug_design"
    }
)

# Dataset saved as human-readable JSON
# data/structure/datasets/covid_proteins.json

# Share dataset - just share the JSON file
# Colleague can load immediately
their_proc = StructureProcessor()
structures = their_proc.load_dataset("covid_proteins")
```

## Advanced Features (Still Zero Config)

### Safe Input Processing

```python
from protos.processing.input import InputProcessor

# Monitor input folder for new files
input_proc = InputProcessor()
input_proc.process_pending()  # Validates and registers new files

# Automatic workflow:
# 1. Files in data/input/pending/ are validated
# 2. Valid files are registered and moved to appropriate directories
# 3. Invalid files go to data/input/rejected/ with error logs
# 4. Successful files archived in data/input/processed/
```

### Custom Data Location (Optional)

```python
# Only if you need non-default location
from protos.io.paths import ProtosPaths

# Set custom location
custom_paths = ProtosPaths(base_path="/large/storage/protos_data")

# Pass to processors
processor = StructureProcessor(paths=custom_paths)
```

### Multi-Project Setup

```python
# Project 1
paths1 = ProtosPaths(base_path="./project1_data")
proc1 = StructureProcessor(paths=paths1)

# Project 2  
paths2 = ProtosPaths(base_path="./project2_data")
proc2 = StructureProcessor(paths=paths2)

# Completely isolated data
```

## Design Principles

### 1. Progressive Disclosure
- Simple tasks are simple
- Complex tasks are possible
- Advanced features don't complicate basic usage

### 2. Fail-Safe Defaults
- Data in current directory (always writeable)
- Auto-create directories
- Human-readable names everywhere
- No configuration required

### 3. Explicit Over Implicit
- No hidden environment variables
- No config files to find
- Clear data organization
- Predictable behavior

### 4. User-Centric Design
- Users never see UUIDs
- File operations are transparent
- Errors are helpful
- Common tasks are easy

## Implementation Guidelines

### For Processor Developers

```python
class NewProcessor(BaseProcessor):
    def __init__(self, name: str = "new_processor", paths: ProtosPaths = None):
        super().__init__(name, paths)
        # That's it! Everything is set up
        
    def save_entity(self, name: str, data: Any):
        # Always use human names
        file_path = self.data_path / f"{self._sanitize_filename(name)}.ext"
        
        # Save file
        self._write_file(file_path, data)
        
        # Auto-register (UUID created internally)
        self.entity_registry.register_entity(
            name=name,  # Human name
            format_type=self.processor_type,
            file_path=str(file_path),
            metadata={...}
        )
        
    def load_entity(self, name: str) -> Any:
        # Registry resolves name to path
        entity_info = self.entity_registry.find_entity(name)
        if entity_info:
            return self._read_file(entity_info.file_path)
        
        # Fallback for unregistered files
        file_path = self.data_path / f"{self._sanitize_filename(name)}.ext"
        if file_path.exists():
            return self._read_file(file_path)
```

### For Users

```python
# Just import and use - no setup required
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.property import PropertyProcessor

# Create processors (automatic initialization)
struct = StructureProcessor()
seq = SequenceProcessor()
prop = PropertyProcessor()

# Everything just works
structure = struct.load_structure("1ubq")
sequence = seq.load_sequence("EGFR_HUMAN")
prop.assign_property("1ubq", "validated", True)
```

## Summary

Protos achieves true zero-configuration through:

1. **Automatic Paths**: ProtosPaths creates `./data` and all subdirectories
2. **Hidden Complexity**: UUIDs used internally, humans see only friendly names
3. **Single Interface**: Users only interact with processors
4. **Smart Defaults**: Everything works out of the box
5. **Optional Overrides**: Advanced users can customize if needed

The result is a system that:
- Works immediately after installation
- Requires no setup or configuration
- Handles complex data relationships transparently
- Scales from simple scripts to large projects
- Maintains human-readable data at all times