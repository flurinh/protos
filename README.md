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

- **Structure Management**: PDB/mmCIF parsing, coordinate analysis, structural alignments, and intelligent caching
- **Sequence Processing**: FASTA handling, database downloads, and sequence analysis
- **GRN System**: Generic Residue Numbering for protein family comparisons
- **Property Management**: Associate experimental data with entities using familiar identifiers
- **ML Integration**: Protein embeddings generation and management (ESM-2, Ankh, etc.)
- **Universal Tracking**: Same entity across multiple data formats
- **CLI Tools**: Command-line utilities for automation and batch processing

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/protos.git
cd protos

# Install in development mode
pip install -e .

# With embedding support (optional)
pip install -e ".[embedding]"

# With GPU support (optional)
pip install -e ".[gpu]"
```

### Basic Usage

```python
import protos
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor

# Optional: Set custom data path before creating processors
# protos.set_data_path("/path/to/your/data")

# Initialize processors (paths handled automatically)
struct_proc = StructureProcessor(name="my_structures")
seq_proc = SequenceProcessor(name="my_sequences")

# Load a structure entity
structure_df = struct_proc.load_entity("1ubq")
print(f"Loaded structure with {len(structure_df)} atoms")

# Extract sequences from structure
sequences = struct_proc.get_all_sequences("1ubq")

# Save sequence entity
for name, seq in sequences.items():
    seq_proc.save_entity(name, seq)

# Same entity, multiple formats - automatically tracked
```

### Initialize Data Directory

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

### StructureProcessor

Manages 3D protein structures from PDB/mmCIF files with intelligent caching:

```python
from protos.processing.structure import StructureProcessor

processor = StructureProcessor(name="structures")

# Load a structure
structure_df = processor.load_entity("1ubq")

# Extract sequences from all chains
sequences = processor.get_all_sequences("1ubq")

# Get single chain sequence
sequence = processor.get_sequence("1ubq", chain_id="A")

# Align structures
processor.align_and_record(
    structure_ids=["1ubq", "2ubq"],
    reference_id="1ubq",
    method="cealign"
)
```

### SequenceProcessor

Handles protein sequences and FASTA files:

```python
from protos.processing.sequence import SequenceProcessor

processor = SequenceProcessor(name="sequences")

# Load a sequence entity
sequence = processor.load_entity("my_protein")

# Save a sequence
processor.save_entity("new_protein", "MVLSPADKTN...")

# Perform pairwise alignment
seq1 = processor.load_entity("protein_1")
seq2 = processor.load_entity("protein_2")
score, alignment = processor.align_sequences(seq1, seq2, "protein_1", "protein_2")

# Apply mutations to a sequence
mutant = processor.mutate_sequence(sequence, ["V91A", "T219F"], "my_mutant")

# Generate all variants at positions
variants = processor.generate_variants(
    sequence,
    positions=[91, 219],
    possible_aas=[["A", "V", "L"], ["F", "W"]],
    base_id="variant"
)
```

### GRNProcessor

Generic Residue Numbering for protein families:

```python
from protos.processing.grn import GRNProcessor

processor = GRNProcessor(name="grn_annotations")

# Load a GRN reference table
processor.load_grn_table("mo_ref")

# Get sequences from GRN table
sequences = processor.get_seq_dict()

# Save GRN annotations
processor.save_grn_table("my_annotations")
```

### PropertyProcessor

Associates experimental data with entities:

```python
from protos.processing.property import PropertyProcessor

processor = PropertyProcessor(name="properties")

# Assign properties to entities
processor.assign_property(
    entity_identifier="opsin_001",
    property_name="lambda_max",
    property_value=568,
    dataset_name="spectral_properties"
)

# Query properties
props = processor.get_entity_properties("opsin_001")

# Filter entities by property
blue_opsins = processor.filter_entities_by_property(
    "spectral_properties",
    {"lambda_max": {"lt": 500}}
)
```

### EmbeddingProcessor

Machine learning embeddings for sequences:

```python
from protos.processing.embedding import EmbeddingProcessor

processor = EmbeddingProcessor(
    name="embeddings",
    model_name="esm2_t6_8m",  # Small model for demo
    batch_size=2
)

# Check if dependencies are available
deps = processor.check_dependencies()

# Generate embeddings
embeddings = processor.embed_sequences(
    {"protein_1": "MVLSPAD...", "protein_2": "MTEYKL..."},
    embedding_type="mean"  # or "cls", "per_residue"
)

# Save embeddings as dataset
processor.embed_sequences(
    sequences,
    save_dataset="my_embeddings"
)
```

## CLI Tools

```bash
# Initialize data directory
protos init                     # Initialize at default location
protos init --path /custom/path # Initialize at custom location

# Clear data directory
protos clear                    # Clear with confirmation
protos clear --force            # Clear without confirmation
protos clear --no-reinit        # Clear without reinitializing
```

## Examples

### Drag-and-Drop Workflow

Add PDB files to the `data/structure/mmcif/` folder and use them immediately:

```python
from protos.processing.structure import StructureProcessor

# Files in data/structure/mmcif/6xyz.cif
processor = StructureProcessor(name="my_processor")
structure = processor.load_entity("6xyz")  # Loads directly

# Create and load a dataset
processor.dataset_manager.create_dataset(
    dataset_id="my_structures",
    name="My Structure Dataset",
    content=["6xyz", "7abc"]
)
structures = processor.load_dataset("my_structures")
```

### Cross-Format Workflow

Work with the same entity across multiple formats:

```python
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor

# Load structure
struct_proc = StructureProcessor(name="structures")
structure = struct_proc.load_entity("1ubq")

# Extract and save sequences
sequences = struct_proc.get_all_sequences("1ubq")
seq_proc = SequenceProcessor(name="sequences")
for name, seq in sequences.items():
    seq_proc.save_entity(name, seq)

# Entity tracked across formats in entity_registry.json
```

### GRN Annotation Workflow

Annotate sequences with Generic Residue Numbers:

```python
from protos.processing.grn import GRNProcessor
from protos.processing.sequence import SequenceProcessor

# Initialize processors
grn_proc = GRNProcessor(name="annotation")
seq_proc = SequenceProcessor(name="sequences")

# Load reference GRN table
grn_proc.load_grn_table("mo_ref")

# Get reference sequences
ref_sequences = grn_proc.get_seq_dict()

# Align query sequences and transfer GRN annotations
# (See examples/grn_basic_annotation.py for full workflow)
```

### Embedding Generation

Generate protein language model embeddings:

```python
from protos.processing.embedding import EmbeddingProcessor

# Initialize with desired model
processor = EmbeddingProcessor(
    name="esm_embeddings",
    model_name="esm2_t36_3b",
    batch_size=4
)

# Check dependencies
if processor.check_dependencies()['ready']:
    # Generate mean embeddings
    embeddings = processor.embed_sequences(
        {"protein_1": sequence_1, "protein_2": sequence_2},
        embedding_type="mean",
        save_dataset="my_embeddings"
    )
```

### Property-Based Analysis

Import and query experimental data:

```python
from protos.processing.property import PropertyProcessor

# Initialize processor
prop_proc = PropertyProcessor(name="opsin_properties")

# Batch assign properties
assignments = [
    {'entity_identifier': 'br_1', 'property_name': 'lambda_max', 'property_value': 568},
    {'entity_identifier': 'br_1', 'property_name': 'photocycle', 'property_value': 'fast'},
    {'entity_identifier': 'chr_1', 'property_name': 'lambda_max', 'property_value': 470},
]
prop_proc.assign_properties_batch(assignments, "opsin_data")

# Query and filter
high_lambda = prop_proc.filter_entities_by_property(
    "opsin_data",
    {"lambda_max": {"gt": 550}}
)
```

### Model Input Preparation

Prepare inputs for structure prediction models:

```python
from protos.models.model_manager import ModelManager
from protos.processing.sequence import SequenceProcessor

# Initialize model manager
manager = ModelManager()

# List available models
print(manager.list_models())  # ['boltz2', 'rfdiffusion', ...]

# Prepare Boltz-2 input
config = {
    "recycling": 5,
    "num_samples": 3,
    "device": "cuda"
}

model_input = manager.prepare_input(
    model_name="boltz2",
    entity_name="my_protein",
    entity_format="sequence",
    config=config
)
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

- **Examples**: Available in the `examples/` directory
- **API Documentation**: Generated from docstrings
- **CLI Reference**: Run `protos --help` for command-line usage

## Support

- **Issues**: Report bugs and request features on GitHub
- **Discussions**: Join our community discussions
- **Documentation**: Check the examples directory

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

**Ready to accelerate your structural biology research? Start with** `protos init` **to set up your data directory!**
