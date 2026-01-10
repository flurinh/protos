# Protos: Structural Biology Framework

A no-nonsense Python framework for handling structural biology data: structures, sequences, annotations, properties, and ML embeddings—all with built-in entity tracking and zero configuration.

## Quick Start

```bash
# Clone and install
git clone https://github.com/flurinh/protos.git
cd protos
pip install -e .

# Optional: embedding support
pip install -e ".[embedding]"
```

## Core Concepts

### 1. Data Path Setup

All protos data lives under a single root directory. Set it **before** creating any processors:

```python
import protos

# Set custom data path (optional - defaults to ./data)
protos.set_data_path("/path/to/your/data")

# Get current path
print(protos.get_data_path())
```

### 2. Entity Registry

Every piece of data (structure, sequence, property, embedding) is an **entity** tracked in a central registry. Entities are identified by human-readable names and can be linked across formats.

```python
from protos.processing.structure import StructureProcessor

sp = StructureProcessor()

# Load creates/retrieves an entity
structure = sp.load_entity("1ubq")

# Check if entity exists
exists = sp.entity_exists("1ubq")

# List all entities
entities = sp.list_entities()
```

### 3. Datasets

Datasets group entities for batch operations. Each processor manages its own datasets.

```python
# Create a dataset
sp.create_dataset("my_structures", ["1ubq", "2ubq", "3ubq"])

# List datasets
datasets = sp.list_datasets()

# Load all entities in a dataset
structures = sp.load_dataset("my_structures")  # Returns Dict[str, DataFrame]
```

---

## Processors

### StructureProcessor

Manages 3D protein structures from PDB/mmCIF files.

#### Loader: StructureLoader

Data sources:
- **RCSB PDB**: `"1ubq"`, `"3sn6"`
- **AlphaFold**: `"AF-P00533-F1"` or UniProt ID `"P00533"`
- **Local files**: Drag-and-drop into `data/structure/mmcif/`

```python
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor

loader = StructureLoader()

# Download single structure
loader.download_and_register("1ubq")

# Download batch and create dataset
loader.download_batch(
    ["3sn6", "4ldo", "2rh1"],
    dataset_name="gpcr_structures",
    create_dataset=True,
)
```

#### Loading Data

```python
sp = StructureProcessor()

# Load single structure (returns pandas DataFrame)
df = sp.load_entity("1ubq")
print(f"Atoms: {len(df)}, Columns: {list(df.columns)}")

# Load dataset
structures = sp.load_dataset("gpcr_structures")
for pdb_id, df in structures.items():
    print(f"{pdb_id}: {len(df)} atoms")
```

#### Operations

```python
# Extract sequences from structure
sequences = sp.get_all_sequences("1ubq")  # Dict[chain_id, sequence]
sequence_a = sp.get_sequence("1ubq", chain_id="A")

# Register chain sequences as entities (cross-processor linking)
sp.register_chain_sequences(
    ["3sn6", "4ldo"],
    dataset_prefix="gpcr_chains",
    create_dataset=True,
)

# Summarize ligands
ligands = sp.summarize_ligands("5d5a")

# Annotate with GRN (Generic Residue Numbering)
sp.annotate_with_grn("3sn6", chains=["R"])
```

#### Analysis

```python
from protos.analysis.structure_ligand_analysis import calculate_ligand_interactions

# Get ligand atoms
df = sp.load_entity("5d5a").reset_index()
ligand_atoms = df[(df["group"] == "HETATM") & (df["res_name3l"] == "CAU")]

# Calculate interactions
interactions = calculate_ligand_interactions(
    sp, "5d5a", ligand_atoms,
    detailed=True,
    cutoff=4.0
)
print(f"Binding residues: {len(interactions['binding_residues'])}")
print(f"H-bonds: {len(interactions['hydrogen_bonds'])}")
```

---

### SequenceProcessor

Handles protein sequences in FASTA format.

#### Loader: SequenceLoader

Data sources:
- **UniProt**: `"uniprot:P00533"` or `"uniprot:P00533,P12345,Q9Y5N6"`
- **Local FASTA files**: Any `.fasta` or `.fa` file path

```python
from protos.io.ingest.sequence_loader import SequenceLoader
from protos.processing.sequence import SequenceProcessor

seq_proc = SequenceProcessor()
loader = SequenceLoader(processor=seq_proc)

# Download from UniProt
loader.download_and_register(
    "uniprot:P00533,P07550",
    name="gpcr_sequences",
    materialize_entities=True,  # Save each sequence as entity
)

# Load local FASTA
loader.download_and_register(
    "/path/to/sequences.fasta",
    name="my_sequences",
)
```

#### Loading Data

```python
seq_proc = SequenceProcessor()

# Load single sequence
sequence = seq_proc.load_entity("P00533")

# Load dataset
sequences = seq_proc.load_dataset("gpcr_sequences")  # Dict[id, sequence]
```

#### Operations

```python
# Save a sequence
seq_proc.save_entity("my_protein", "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH")

# Align two sequences
score, alignment = seq_proc.align_sequences(
    seq1, seq2, "protein_1", "protein_2"
)

# Apply mutations
mutant = seq_proc.mutate_sequence(sequence, ["V91A", "T219F"], "my_mutant")

# Generate mutant library
library, metadata = seq_proc.create_mutant_library(
    base_sequence_id="wild_type",
    mutation_map={91: ["A", "V", "L"], 219: ["F", "W"]},
    limit=10,
    return_metadata=True,
)

# Conservation analysis
conservation = seq_proc.compute_conservation(sequences=library)

# Export to FASTA
seq_proc.export_dataset("gpcr_sequences", export_name="gpcr_export")
```

#### Analysis: GRN Annotation

```python
# Annotate sequences with Generic Residue Numbering
grn_table, summary = seq_proc.annotate_with_grn(
    dataset_name="gpcr_sequences",
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
    output_table="my_grn_annotations",
    return_summary=True,
    allow_create=True,
)

print(f"Annotated: {summary['global']['annotated']}/{summary['global']['total']}")
```

---

### PropertyProcessor

Associates experimental data with entities using scoped properties.

#### Loading/Storing Properties

```python
from protos.processing.property import PropertyProcessor

prop_proc = PropertyProcessor()

# Record properties (creates table if needed)
rows = [
    {
        "scope": [{"format": "sequence", "name": "opsin_1"}],
        "entity_name": "opsin_1",
        "lambda_max": 568,
        "photocycle": "fast",
    },
    {
        "scope": [
            {"format": "structure", "name": "5d5a"},
            {"format": "sequence", "name": "5d5a_chain_A"},
        ],
        "entity_name": "5d5a_chain_A",
        "classification": "gpcr_like",
        "similarity_score": 0.85,
    },
]

prop_proc.record_properties(
    "opsin_properties",
    rows,
    metadata={"source": "experimental"},
    allow_create=True,
)
```

#### Querying Properties

```python
# Get properties for an entity
props = prop_proc.get_properties("opsin_1")
print(props.to_string())

# Get properties for a structure (across all linked sequences)
struct_props = prop_proc.get_properties("5d5a")
```

---

### GRNProcessor

Manages Generic Residue Numbering tables for protein family comparisons.

```python
from protos.processing.grn import GRNProcessor

grn_proc = GRNProcessor()

# Load reference table (for alignment-based annotation)
ref_table = grn_proc.load_reference_table("gpcrdb_ref")
print(f"Reference table: {ref_table.shape}")

# List available annotation tables
tables = grn_proc.list_tables()

# Load a custom annotation table
if tables:
    grn_table = grn_proc.load_table(tables[0])

# Annotate sequences directly
sequences = {"protein_1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"}
annotations, summary = grn_proc.annotate_sequences(
    sequences,
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
)
```

---

### EmbeddingProcessor

Generates protein language model embeddings.

#### Direct Usage

```python
from protos.processing.embedding import EmbeddingProcessor

# Initialize with model name
emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")

# Check available models
print(EmbeddingProcessor.MODEL_REGISTRY.keys())
# ['esm2_t6_8m', 'esm2_t12_35m', 'esm2_t30_150m', ..., 'ankh_base', 'ankh_large']

# Generate embeddings directly
sequences = {"protein_1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH"}
embeddings = emb_proc.embed_sequences(
    sequences,
    embedding_type="mean",  # or "per_residue", "sum", "cls"
    save_dataset="my_embeddings"  # optional: save as dataset
)

# Load saved embeddings (returns dict with original sequence IDs)
embeddings = emb_proc.load_embeddings("my_embeddings")
```

#### Via ModelManager

```python
from protos.processing.embedding import EmbeddingProcessor
from protos.models.model_manager import ModelManager

manager = ModelManager()

# List embedding cards (use .cards dict)
embedding_cards = [name for name in manager.cards.keys() if 'embedding' in name]
# ['embedding_esm2_t12_35m', 'embedding_ankh_large', ...]

# Prepare and run embedding
invocation = manager.prepare(
    "embedding_esm2_t12_35m",
    inputs={"sequence_dataset": "gpcr_sequences"},
    config={"embedding_type": "mean"},  # or "per_residue", "sum"
)

# Ingest results into processor
emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")
emb_proc.ingest_from_invocation(invocation, dataset_name="gpcr_embeddings")

# Load embeddings
embeddings = emb_proc.load_embeddings("gpcr_embeddings")
```

---

### MoleculeProcessor

Handles small molecule/ligand data.

```python
from protos.processing.molecule import MoleculeProcessor

mol_proc = MoleculeProcessor()

# Save ligand with SMILES
mol_proc.save_entity(
    "5d5a_CAU_A",
    {"smiles": "CC(C)NC1=NC(=O)N(C=C1)C2=CN=C(N)N=C2N", "kind": "smiles_record"},
    metadata={"source_structure": "5d5a"},
)
```

---

## Model Integration

The `ModelManager` provides a unified interface for structure prediction and analysis models.

```python
from protos.models.model_manager import ModelManager

manager = ModelManager()

# List available models via .cards dict
print(list(manager.cards.keys()))
# ['boltz2', 'boltzgen', 'embedding_esm2_t12_35m', 'embedding_ankh_large', 'lambda', ...]

# Prepare model input
invocation = manager.prepare(
    "boltz2",
    inputs={"sequence_dataset": "my_sequences", "entity": "protein_1"},
    config={
        "output_name": "protein_1_prediction",
        "recycling": 5,
        "num_samples": 3,
    },
)

# Access job details (for external models)
if invocation.job:
    job = invocation.job
    print(f"Command: {' '.join(job.command)}")
    print(f"Working dir: {job.working_dir}")

# Access runtime results (for in-process models)
if invocation.runtime:
    print(f"Status: {invocation.runtime.outputs.get('status')}")
```

---

## Complete Workflow Example

```python
import protos
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.property import PropertyProcessor

# Setup - must be called BEFORE creating any processors
protos.set_data_path("/path/to/data")  # or use default ./data

# Download structures
loader = StructureLoader()
loader.download_batch(["3sn6", "4ldo", "2rh1"], dataset_name="gpcr_study")

# Initialize processors
struct_proc = StructureProcessor()
seq_proc = SequenceProcessor()
prop_proc = PropertyProcessor()

# Extract and register sequences
struct_proc.register_chain_sequences(
    ["3sn6", "4ldo", "2rh1"],
    dataset_prefix="gpcr_chains",
    create_dataset=True,
)

# Annotate with GRN
grn_table, summary = seq_proc.annotate_with_grn(
    dataset_name="gpcr_chains_3sn6",
    reference_table="gpcrdb_ref",
    protein_family="gpcr_a",
    output_table="gpcr_grn",
    return_summary=True,
    allow_create=True,
)

# Record analysis results (scope is required for each property row)
prop_proc.record_properties(
    "gpcr_analysis",
    [
        {"scope": [{"format": "structure", "name": "3sn6"}], "entity_name": "3sn6", "state": "active", "ligand_type": "full_agonist"},
        {"scope": [{"format": "structure", "name": "4ldo"}], "entity_name": "4ldo", "state": "active", "ligand_type": "full_agonist"},
        {"scope": [{"format": "structure", "name": "2rh1"}], "entity_name": "2rh1", "state": "inactive", "ligand_type": "inverse_agonist"},
    ],
    allow_create=True,
)

print("Workflow complete!")
```

---

## Data Organization

```
data/
├── entity_registry.json      # Central entity tracking
├── structure/
│   ├── mmcif/               # PDB/CIF files
│   ├── cache/               # Preprocessed PKL files
│   ├── datasets/            # Dataset definitions
│   └── registry.json
├── sequence/
│   ├── fasta/               # FASTA files
│   ├── alignments/          # Alignment results
│   ├── datasets/
│   └── registry.json
├── grn/
│   ├── tables/              # GRN annotation tables
│   ├── reference/           # Reference GRN tables
│   └── registry.json
├── property/
│   ├── tables/              # Property tables (CSV)
│   └── registry.json
├── embedding/
│   ├── embeddings/          # Saved embeddings
│   └── registry.json
└── ligand/
    ├── sdf/                 # SDF/MOL files
    └── registry.json
```

---

## CLI Tools

```bash
protos init                     # Initialize data directory
protos init --path /custom/path # Custom location
protos clear                    # Clear with confirmation
protos clear --force            # Clear without confirmation
```

---

## Examples

See `examples/` for complete workflows:
- `examples/workflow_tests/gpcr_ligand_mechanism.py` - GPCR binding analysis
- `examples/workflow_tests/structure/test_structure_grn_annotation.py` - GRN annotation
- `examples/workflow_tests/sequence/test_sequence_workflow.py` - Sequence operations
- `examples/workflow_tests/embedding/test_embedding_workflow.py` - ML embeddings

---

## License

MIT License - see LICENSE file for details.
