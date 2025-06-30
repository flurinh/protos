# Protos Documentation

Protos is a comprehensive structural biology framework that provides a unified interface for managing, analyzing, and integrating diverse biological data types including protein structures, sequences, annotations, and experimental properties.

## Core Philosophy

Protos abstracts away file system complexity - researchers work with biological identifiers (names, IDs) while the framework handles all path resolution, file format conversions, and data management operations internally.

## Documentation Overview

### Getting Started
- [Installation Guide](cli/installation.md) - Install Protos and set up the environment
- [Quick Start Tutorial](cli/quickstart.md) - Basic usage examples to get started
- [Core Concepts](core_concepts.md) - Understand the fundamental abstractions

### Core Infrastructure
- [ProtosPaths System](protospath.md) - Centralized path management
- [BaseProcessor Architecture](baseprocessor.md) - Common interface for all processors
- [Entity System](entities.md) - Universal data tracking across formats
- [Registry System](registries.md) - Metadata and dataset management
- [File Formats](fileformats.md) - Supported formats and conventions

### Data Processors
- [CifBaseProcessor](processors/cifbase_processor.md) - 3D structure management
- [GRNBaseProcessor](processors/grn_processor.md) - Generic Residue Numbering
- [SeqProcessor](processors/seq_processor.md) - Sequence data handling
- [EmbeddingProcessor](processors/embedding_processor.md) - ML embeddings
- [PropertyProcessor](processors/property_processor.md) - Experimental properties

### Processing Operations
- [Alignments](processing/alignments.md) - Sequence and structure alignment
- [Annotations](processing/annotations.md) - GRN and functional annotations
- [Conversions](processing/conversions.md) - Cross-format data conversion

### Command Line Interface
- [CLI Reference](cli/commands.md) - Complete command reference
- [Usage Examples](cli/examples.md) - Real-world usage patterns

### API Reference
- [API Documentation](api/index.md) - Detailed API reference
- [Code Examples](api/examples.md) - Programming examples

## Key Features

### Universal Entity System
Every biological entity receives a deterministic hash ID that remains consistent across all data formats:
```python
# Same protein, multiple representations
"P12345" → entity_id "a3f2d8c91b"
  ├── sequence: /sequence/fasta/a3f2d8c91b.fasta
  ├── structure: /structure/mmcif/a3f2d8c91b.cif
  ├── embedding: /embedding/a3f2d8c91b.pkl
  └── properties: tracked in property datasets
```

### Name-Based Access
```python
# Work with biological names, not file paths
processor.load_structure("1ubq")          # PDB ID
processor.load_sequence("P12345")         # UniProt ID
processor.load_grn_table("opsin_alignment")  # Table name
```

### Cross-Format Integration
```python
# Extract sequences from structures
sequences = struct_processor.get_seq_dict()

# Generate embeddings from sequences
embeddings = emb_processor.embed_sequences(sequences)

# Associate properties with entities
prop_processor.create_property_dataset(properties, entity_column="protein_id")
```

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        PROTOS FRAMEWORK                          │
├─────────────────────────────────────────────────────────────────┤
│ 1. CORE INFRASTRUCTURE                                          │
│    ├─ ProtosPaths: Centralized path management                 │
│    ├─ BaseProcessor: Abstract base for all processors          │
│    └─ Entity Registry: Universal entity tracking               │
│                                                                 │
│ 2. SPECIALIZED PROCESSORS                                       │
│    ├─ CifBaseProcessor: 3D protein structures                  │
│    ├─ GRNBaseProcessor: Generic Residue Numbering              │
│    ├─ SeqProcessor: Sequence data                              │
│    ├─ PropertyProcessor: Metadata and properties               │
│    └─ EmbeddingProcessor: ML embeddings                        │
│                                                                 │
│ 3. DATA ORGANIZATION                                            │
│    ├─ Entities: Individual data items                          │
│    ├─ Datasets: Collections of related entities                │
│    └─ Cross-format tracking: Same entity, multiple formats     │
└─────────────────────────────────────────────────────────────────┘
```

## Quick Examples

### Loading and Analyzing Structures
```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Initialize processor
cp = CifBaseProcessor(name="analysis")

# Load structure by PDB ID
structure = cp.load_structure("1ubq")

# Create dataset from multiple structures
cp.create_dataset("test_set", name="Test proteins", 
                  content=["1ubq", "2gb1", "1crn"])

# Load entire dataset
cp.load_dataset("test_set")
```

### Working with Sequences
```python
from protos.processing.sequence.seq_processor import SeqProcessor

# Initialize processor
sp = SeqProcessor(name="sequences")

# Load sequences
sequences = sp.load_sequences(["P12345", "Q9Y6K1"])

# Save new sequences
sp.save_sequences({"MY_PROTEIN": "MKTAYIAKQRQ..."}, "my_proteins.fasta")
```

### Cross-Format Operations
```python
# Extract sequences from structures
sequences = struct_processor.get_seq_dict()

# Save to sequence processor
seq_processor.save_sequences(sequences, "structure_derived.fasta")

# Generate embeddings
embeddings = emb_processor.embed_sequences(sequences)
```

## Next Steps

1. Follow the [Installation Guide](cli/installation.md)
2. Work through the [Quick Start Tutorial](cli/quickstart.md)
3. Explore specific [processors](processors/cifbase_processor.md) for your data type
4. Check [examples](cli/examples.md) for common workflows

## Support

- GitHub Issues: [Report bugs or request features](https://github.com/protos/protos/issues)
- Documentation: [This documentation](https://protos.readthedocs.io)
- PyPI: [Python package](https://pypi.org/project/protos)