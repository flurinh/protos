# Protos: Structural Biology Framework

![Protos-MCP](resources/logo.png)A no-nonsense Python framework for handling structural biology data: structures, sequences, annotations, properties, and ML embeddings—all with built-in entity tracking and zero configuration.

## Why Protos?
Structural biology workflows suck: scattered data formats, path hell, inconsistent IDs, and endless boilerplate. Protos fixes that. It's a developer kit that lets protein engineers focus on science, not plumbing.

## Zero Setup: Clone, pip install, done. 
ProtosPaths auto-manages your data—no configs, no env vars required (but customizable if you want).
- **Unified Entities**: Track the same protein across PDB, FASTA, GRN tables, properties, and embeddings.
- **Use human-readable names**: Protos handles the rest.
- **Batteries Included**: Processors for every major task, CLI for automation, parallel support for scale.
- **Extensible**: Base classes make adding custom processors trivial.
- **Research-Ready**: Built for real workflows—GPCR analysis, opsin engineering, ML featurization.

If you're tired of wrestling with Biopython wrappers or custom scripts, Protos streamlines it all.

## Key Features

- **🏗️ Structure Management**: PDB/mmCIF parsing, coordinate analysis, structural alignments, and intelligent caching
- **🧬 Sequence Processing**: FASTA handling, database downloads, and sequence analysis
- **📊 GRN System**: Generic Residue Numbering for protein family comparisons
- **📋 Property Management**: Associate experimental data with entities using familiar identifiers
- **🤖 ML Integration**: Protein embeddings generation and management
- **🔗 Universal Tracking**: Same entity across multiple data formats
- **⚙️ CLI Tools**: Command-line utilities for automation and batch processing

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/protos.git
cd protos

# Install in development mode
pip install -e .
```

### Basic Usage

```python
from protos.processing.structure import CifBaseProcessor
from protos.processing.sequence import SeqProcessor

# Initialize processors (paths handled automatically)
struct_proc = CifBaseProcessor(name="my_structures")
seq_proc = SeqProcessor(name="my_sequences")

# Download and process a structure
struct_proc.download_structures(["1ubq"])
structure_data = struct_proc.load_structure("1ubq")

# Extract sequence from structure
sequences = struct_proc.extract_sequence("1ubq")
seq_proc.save_sequence("1ubq", sequences)

# Same entity, multiple formats - automatically tracked
entity_info = struct_proc.entity_registry.find_entity("1ubq")
# Entity ID works across all processors for the same biological entity
```

```bash
protos init
```
## Architecture

### Core Principle: ProtosPaths Manages Everything

**ProtosPaths is the ONLY path management system in Protos.** No component ever specifies or manages paths directly. This ensures:

- **Zero configuration** - works immediately with sensible defaults
- **Consistent organization** - all data follows the same structure
- **Easy portability** - change data location in one place
- **Automatic discovery** - components know where to find data

### Data Organization (Managed by ProtosPaths)

```
working_dir/
└── data/                    # Default base (or user-specified)
    ├── entity_registry.json # Universal entity tracking
    ├── structure/           # Structure data
    │   ├── mmcif/          # PDB/CIF files (human-readable names)
    │   ├── cache/          # Preprocessed structure PKL files
    │   ├── structure_dataset/ # Preprocessed dataset PKL files
    │   ├── alignments/     # Structural alignments
    │   ├── datasets/       # Dataset JSON definitions
    │   └── registry.json   # Structure processor registry
    ├── sequence/           # Sequence data
    │   ├── fasta/         # FASTA files (human-readable names)
    │   ├── alignments/    # All alignment results
    │   │   ├── pairwise/  # Pairwise alignments
    │   │   ├── multiple/  # Multiple sequence alignments
    │   │   └── mmseqs/    # MMseqs2 alignments
    │   ├── databases/     # MMseqs2 databases
    │   ├── metadata/      # Sequence metadata
    │   ├── datasets/      # Dataset JSON definitions
    │   └── registry.json  # Sequence processor registry
    ├── grn/               # GRN data
    │   ├── tables/        # GRN annotation tables
    │   ├── reference/     # Reference GRN tables
    │   ├── configs/       # GRN configurations
    │   ├── assignments/   # GRN assignment results
    │   ├── temp/          # Temporary files
    │   ├── datasets/      # Dataset JSON definitions
    │   └── registry.json  # GRN processor registry
    ├── property/          # Property data
    │   ├── tables/        # Property tables (CSV)
    │   ├── datasets/      # Dataset JSON definitions
    │   └── registry.json  # Property processor registry
    ├── embedding/         # ML embeddings
    │   ├── embeddings/    # Saved embeddings
    │   ├── datasets/      # Dataset JSON definitions
    │   └── registry.json  # Embedding processor registry
    ├── ligand/            # Ligand data
    │   ├── sdf/           # SDF/MOL files
    │   ├── cache/         # Cached ligand data
    │   ├── datasets/      # Dataset JSON definitions
    │   └── registry.json  # Ligand processor registry
    ├── graph/             # Graph/network data
    │   ├── networks/      # Graph/network files
    │   ├── analysis/      # Graph analysis results
    │   ├── datasets/      # Dataset JSON definitions
    │   └── registry.json  # Graph processor registry
    └── temp/              # Temporary processing files
```

## Processors

### CifBaseProcessor

Manages 3D protein structures from PDB/mmCIF files with intelligent caching:

- Automatic structure downloads from PDB
- Coordinate parsing and analysis
- Chain and residue filtering
- Cross-format sequence extraction
- Caches preprocessed structures (PKL) for performance
- Supports dataset-level caching for large-scale analysis

```python
from protos.processing.structure import StructureProcessor

struct_proc = StructureProcessor()

# Annotate GPCR chains with GRNs and persist the table
grn_df = struct_proc.annotate_structures_with_grn(
    ["3sn6", "5d5a", "6b73"],
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
    output_table="gpcr_structure_grn",
)

print(grn_df.head())
```

### SeqProcessor

Handles protein sequences and FASTA files:

- UniProt database integration
- Sequence analysis and comparison
- Multiple sequence alignment
- Batch sequence processing

### GRNBaseProcessor

Generic Residue Numbering for protein families:

- Standardized position numbering across homologs
- Support for GPCRs, opsins, and custom schemes
- Position conservation analysis
- Functional site identification
- Table-based storage with reference tables and configurations

```python
from protos.processing.sequence import SequenceProcessor

seq_proc = SequenceProcessor()

# Annotate a sequence dataset using the bundled GPCR reference
grn_df = seq_proc.annotate_with_grn(
    dataset_name="gpcr_sequences",
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
    output_table="gpcr_sequences_grn",
)

print(grn_df.head())  # residue-by-residue GRN mapping
```

### PropertyProcessor

Associates experimental data with entities:

- Flexible CSV import with human-readable identifiers
- Property filtering and analysis
- Cross-format property queries
- Experimental metadata tracking

### EmbeddingProcessor

Machine learning embeddings for sequences:

- ESM-2, Ankh, and other protein language models
- Mean, CLS, and per-residue embeddings
- Automatic caching and reuse
- ML-ready dataset preparation

## CLI Tools

```bash
# Data management
protos init --path ./data      # Initialize (or refresh) a data directory
protos clear --force           # Delete and optionally reinitialize the active data root

# Legacy helpers (maintained for compatibility)
protos-init                    # Equivalent to `protos init`
protos-cleanup --reinit        # Equivalent to `protos clear`

# Analysis tools
protos assign-grns             # Assign GRN numbers to sequences
protos expand-annotation       # Expand GRN annotations
protos process-structures      # Batch structure analysis
protos generate-embeddings     # Batch embedding generation
```

## Examples

### Drag-and-Drop Workflow

Add PDB files to the `data/structure/mmcif/` folder and use them immediately:

```python
from protos.processing.structure import CifBaseProcessor

# Files in data/structure/mmcif/6xyz.cif
processor = CifBaseProcessor()
structure = processor.load_structure("6xyz")  # Loads directly

# Create dataset
processor.create_dataset("my_structures", ["6xyz", "7abc"])
structures = processor.load_dataset("my_structures")
```

### Cross-Format Workflow

Work with the same entity across multiple formats:

```python
from protos.processing.structure import CifBaseProcessor
from protos.processing.sequence import SeqProcessor

# Load structure
struct_proc = CifBaseProcessor()
structure = struct_proc.load_structure("1ubq")

# Extract and save sequence
sequence = struct_proc.extract_sequence("1ubq")
seq_proc = SeqProcessor()
seq_proc.save_sequence("1ubq", sequence)

# Entity tracked across formats in entity_registry.json
```

### Property-Based Analysis

Import and query experimental data:

```python
from protos.processing.property import PropertyProcessor
import pandas as pd

# Import properties from DataFrame (dataset-level storage by default)
properties = pd.DataFrame({
    "entity_name": ["BACR_HALSA", "ChR2", "NpHR"],
    "lambda_max": [568, 470, 590],
    "photocycle": ["fast", "slow", "fast"],
    "scope": [[{"format": "sequence", "name": seq}] for seq in ["BACR_HALSA", "ChR2", "NpHR"]],
})
prop_proc = PropertyProcessor()
prop_proc.record_properties("opsin_properties", properties)

# Retrieve rows associated with a given sequence when needed
blue_opsins = prop_proc.load_dataset_rows(
    "opsin_properties", entity_name="ChR2", format_type="sequence"
)
```

## Comprehensive Review

Run the complete framework demonstration:

```bash
# Full comprehensive review
python protos_review.py

# Quick overview only
python protos_review.py --quick

# Interactive demo mode
python protos_review.py --demo
```

Or explore the Jupyter notebook:

```bash
jupyter notebook protos_review.ipynb
```

## Research Applications

- **Protein Family Analysis**: Compare homologous proteins using standardized numbering
- **Structure-Function Studies**: Correlate experimental properties with structural features
- **Machine Learning**: Prepare datasets with sequences, structures, and experimental labels
- **Drug Discovery**: Identify and analyze therapeutic targets
- **Protein Engineering**: Guide rational design using conservation and property data

## Contributing

We welcome contributions! Please see our contributing guidelines and:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Documentation

- **Framework Overview**: See `protos_review.py` for comprehensive examples
- **API Documentation**: Generated from docstrings
- **Tutorials**: Available in the `examples/` directory
- **CLI Reference**: Run `protos --help` for command-line usage

## Support

- **Issues**: Report bugs and request features on GitHub
- **Discussions**: Join our community discussions
- **Documentation**: Check the comprehensive review and examples

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use Protos in your research, please cite:

```bibtex
@software{protos2024,
  title={Protos: A Structural Biology Framework for Universal Data Management},
  author={Your Team},
  year={2024},
  url={https://github.com/your-org/protos}
}
```

---

**Ready to accelerate your structural biology research? Start with** `python protos_review.py` **to see everything Protos can do!** 🚀
