# Protos: Structural Biology Framework

![Protos-MCP](resources/logo.png)

A comprehensive Python framework for managing, analyzing, and integrating structural biology data with universal entity tracking and cross-format compatibility.

## Overview

Protos provides a unified interface for handling diverse biological data types including protein structures, sequences, residue numbering schemes, experimental properties, and machine learning embeddings. The framework eliminates common pain points in structural biology research through centralized path management, automatic entity registration, and cross-format data integration.

## Key Features

- **🏗️ Structure Management**: PDB/mmCIF parsing, coordinate analysis, and structural alignments
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

# Initialize data directory
python -m protos.cli.init_data
```

### Basic Usage

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.sequence.seq_processor import SeqProcessor

# Initialize processors (paths handled automatically)
struct_proc = CifBaseProcessor(name="my_structures")
seq_proc = SeqProcessor(name="my_sequences")

# Download and process a structure
struct_proc.download_structures(["1ubq"])
structure_data = struct_proc.load_structure("1ubq")

# Extract sequence from structure
sequences = struct_proc.get_seq_dict()
seq_proc.save_sequences(sequences, "extracted_sequences.fasta")

# Same entity, multiple formats - automatically tracked
entity_id = struct_proc.generate_entity_id("1ubq")
# This ID works across all processors for the same biological entity
```

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                           PROTOS FRAMEWORK                     │
├─────────────────────────────────────────────────────────────────┤
│ 1. CORE INFRASTRUCTURE                                          │
│    ├─ ProtosPaths: Centralized path management                  │
│    ├─ BaseProcessor: Abstract base for all processors          │
│    └─ Entity Registry: Universal entity tracking system        │
│                                                                 │
│ 2. SPECIALIZED PROCESSORS                                       │
│    ├─ CifBaseProcessor: 3D protein structure management        │
│    ├─ GRNBaseProcessor: Generic Residue Numbering system       │
│    ├─ SeqProcessor: Sequence data management                   │
│    ├─ PropertyProcessor: Metadata and properties               │
│    └─ EmbeddingProcessor: ML embeddings generation             │
│                                                                 │
│ 3. DATA ABSTRACTION                                             │
│    ├─ Entities: Individual data items                          │
│    ├─ Datasets: Collections of related entities                │
│    └─ Cross-format tracking: Same entity across formats        │
│                                                                 │
│ 4. APPLICATIONS                                                 │
│    ├─ CLI tools: Command-line utilities                        │
│    ├─ Analysis workflows: Multi-processor pipelines            │
│    └─ Research integration: Real-world use cases               │
└─────────────────────────────────────────────────────────────────┘
```

## Processors

### CifBaseProcessor
Manages 3D protein structures from PDB/mmCIF files:
- Automatic structure downloads from PDB
- Coordinate parsing and analysis
- Chain and residue filtering
- Cross-format sequence extraction

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

### PropertyProcessor
Associates experimental data with entities:
- CSV import with user-friendly identifiers
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
protos init-data              # Initialize data directory
protos cleanup-data           # Clean and organize data
protos list-entities          # Browse available entities
protos list-datasets          # Browse datasets

# Analysis tools
protos assign-grns            # Assign GRN numbers to sequences
protos expand-annotation      # Expand GRN annotations
protos process-structures     # Batch structure processing
protos generate-embeddings    # Batch embedding generation
```

## Examples

### Multi-Format Entity Workflow
```python
# Start with a structure
struct_proc = CifBaseProcessor(name="demo")
struct_proc.download_structures(["1ubq"])

# Extract sequence
sequences = struct_proc.get_seq_dict()
seq_proc = SeqProcessor(name="demo")
seq_proc.save_sequences(sequences, "1ubq_sequence.fasta")

# Add to GRN table
grn_proc = GRNBaseProcessor(name="demo")
# ... create GRN alignment including 1ubq sequence ...

# Add experimental properties
prop_proc = PropertyProcessor(name="demo")
prop_proc.assign_property("1ubq", "resolution", 1.8, "structural_data")
prop_proc.assign_property("1ubq", "organism", "Human", "structural_data")

# Same entity "1ubq" now exists in 4 formats with consistent tracking
```

### Property-Based Analysis
```python
# Import experimental data from CSV
prop_proc = PropertyProcessor(name="analysis")
prop_proc.create_property_dataset_from_file(
    "experimental_data.csv",
    "opsin_properties",
    entity_column='protein_id'  # Uses familiar names like PDB IDs
)

# Filter by properties
blue_opsins = prop_proc.filter_entities_by_property(
    "opsin_properties", 
    {"lambda_max": {"lt": 500}}  # Wavelength < 500nm
)

thermostable = prop_proc.filter_entities_by_property(
    "opsin_properties",
    {"thermal_stability": {"gt": 60}}  # Tm > 60°C
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

**Ready to accelerate your structural biology research? Start with `python protos_review.py` to see everything Protos can do!** 🚀
