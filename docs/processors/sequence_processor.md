# SequenceProcessor

Manages protein sequences with alignment capabilities.

**Location:** `protos.processing.sequence.sequence_processor`

**Processor Type:** `sequence`

## Overview

The SequenceProcessor provides:
- Loading and saving sequences from FASTA files
- Pairwise and multiple sequence alignments
- MMseqs2 integration for fast searches
- Mutation analysis and variant generation
- Sequence metadata and conservation analysis

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(name)` | Load sequence by name. Returns string (single) or dict (multi-sequence). |
| `save_entity(name, data, metadata)` | Save single sequence to FASTA and register. |
| `load_sequence(identifier, use_cache)` | Load sequence with optional caching. |
| `get_sequence(identifier)` | Get sequence string by identifier. |
| `export_entity(name, export_name, format, sequence_ids)` | Export sequence to file. |
| `export_dataset(dataset_name, output_dir, format)` | Export all sequences in dataset. |

### Dataset Operations

| Method | Description |
|--------|-------------|
| `load_dataset(dataset_name)` | Load all sequences in dataset as dict. |
| `load_sequences(fasta_file, dataset_name, register_entities)` | Load from FASTA file into dataset. |
| `save_sequences(sequences, output_file, dataset_name, materialize_entities)` | Save multiple sequences. |
| `list_entities(entity_type)` | List registered sequence entities. |
| `list_datasets()` | List available sequence datasets. |

### Alignment Operations

| Method | Description |
|--------|-------------|
| `align_sequences(seq1, seq2, return_alignment)` | Pairwise alignment using BioPython. |
| `align_sequences_mmseqs(query, targets, min_seq_id)` | Fast alignment using MMseqs2. |
| `multiple_sequence_alignment(sequences, method)` | Multiple sequence alignment (MSA). |
| `align_and_record(query_id, reference_id, method, dataset_name)` | Align and store results with metadata. |
| `find_best_match(query_seq, target_sequences, top_k)` | Find best matching sequences. |

### MMseqs2 Database Operations

| Method | Description |
|--------|-------------|
| `create_mmseqs_database(sequences, db_name)` | Create MMseqs2 search database. |
| `list_mmseqs_databases()` | List available MMseqs2 databases. |

### Mutation Analysis

| Method | Description |
|--------|-------------|
| `mutate_sequence(sequence, mutations, format)` | Apply mutations to sequence. |
| `generate_variants(sequence, positions, amino_acids)` | Generate single-point variants. |
| `create_mutant_library(sequence, mutations, combinatorial)` | Create systematic mutant library. |
| `generate_mutants_for_sequence(sequence_id, mutations, dataset_name)` | Generate and register mutant sequences. |

### Conservation Analysis

| Method | Description |
|--------|-------------|
| `compute_conservation(sequences, method, reference)` | Compute position-wise conservation scores. |
| `compute_linkage(sequences, positions, method)` | Compute co-evolution / linkage scores. |

### GRN Annotation

| Method | Description |
|--------|-------------|
| `annotate_with_grn(sequence_ids, grn_table, protein_family)` | Add GRN annotations to sequences. |

### Sequence Metadata

| Method | Description |
|--------|-------------|
| `get_sequence_metadata(sequence_ids)` | Get computed metadata (MW, pI, composition). |

### Alignment I/O

| Method | Description |
|--------|-------------|
| `save_alignment(alignment_data, output_file, format)` | Save alignment results. |
| `load_alignment(alignment_file, alignment_type)` | Load saved alignment. |

### Related Entities

| Method | Description |
|--------|-------------|
| `list_source_structures(sequence_ids, rel_type, direction)` | Find related structure entities. |
| `list_dataset_source_structures(dataset_name, rel_type)` | Find structures for dataset sequences. |

### Path Properties

| Property | Description |
|----------|-------------|
| `path_entity_dir` | Single-sequence FASTA files |
| `path_fasta_dir` | Dataset FASTA files |
| `path_alignments_dir` | General alignments |
| `path_pairwise_alignments_dir` | Pairwise alignment results |
| `path_multiple_alignments_dir` | MSA results |
| `path_mmseqs_alignments_dir` | MMseqs2 results |
| `path_databases_dir` | MMseqs2 databases |
| `path_metadata_dir` | Sequence metadata |

### Utility Properties

| Property | Description |
|----------|-------------|
| `alignment_engine` | Lazy-loaded SequenceAlignmentEngine instance |
| `aligner` | BioPython PairwiseAligner instance |

---

## Data Format

Sequences are stored as FASTA files:

```
>sequence_id
MKWVTFISLLLLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIA...
```

- **Single-sequence files**: `load_entity()` returns string
- **Multi-sequence files**: `load_entity()` returns `Dict[str, str]`

---

## Usage Examples

### Basic Operations

```python
from protos import SequenceProcessor

proc = SequenceProcessor()

# Load sequence
seq = proc.load_entity("EGFR_HUMAN")

# Save sequence
proc.save_entity("my_protein", "MKWVTFISLLLLFSSAYS...", metadata={
    "organism": "human",
    "gene": "EGFR"
})

# Get sequence metadata
metadata = proc.get_sequence_metadata(["my_protein"])
print(f"MW: {metadata['molecular_weight'].iloc[0]:.1f} Da")
```

### Pairwise Alignment

```python
# Align two sequences
result = proc.align_sequences(seq1, seq2, return_alignment=True)
print(f"Score: {result['score']}")
print(f"Identity: {result['identity']:.2%}")
print(f"Alignment:\n{result['alignment']}")

# Find best match in database
matches = proc.find_best_match(query_seq, target_sequences, top_k=5)
for match in matches:
    print(f"{match['name']}: {match['identity']:.2%}")
```

### Multiple Sequence Alignment

```python
sequences = {
    "seq1": "MKWVTFIS...",
    "seq2": "MLKFTISA...",
    "seq3": "MKVLTDIS..."
}

# Perform MSA
alignment = proc.multiple_sequence_alignment(sequences, method="blosum62")
```

### MMseqs2 Fast Search

```python
# Create database
proc.create_mmseqs_database(target_sequences, "my_database")

# Search
results = proc.align_sequences_mmseqs(
    query_seq,
    targets="my_database",
    min_seq_id=0.3
)
```

### Mutation Analysis

```python
# Apply single mutation
mutant = proc.mutate_sequence("MKWVTFIS...", ["K48R"])

# Generate all single-point variants at positions
variants = proc.generate_variants(
    sequence,
    positions=[48, 51, 55],
    amino_acids=["A", "R", "K", "D", "E"]
)

# Create combinatorial mutant library
library = proc.create_mutant_library(
    sequence,
    mutations=["K48R", "E51A", "D55N"],
    combinatorial=True  # All combinations
)
```

### Conservation Analysis

```python
# Load aligned sequences
sequences = proc.load_dataset("my_alignment")

# Compute conservation
conservation = proc.compute_conservation(
    sequences,
    method="shannon",
    reference="seq1"
)

# Compute co-evolution
linkage = proc.compute_linkage(
    sequences,
    positions=[48, 51, 55, 100],
    method="mutual_information"
)
```

### Working with Datasets

```python
# Load sequences from FASTA
sequences = proc.load_sequences(
    "proteins.fasta",
    dataset_name="my_proteins",
    register_entities=True
)

# Save multiple sequences
proc.save_sequences(
    {"seq1": "MKW...", "seq2": "MLK..."},
    output_file="combined.fasta",
    dataset_name="combined"
)

# Find related structures
related = proc.list_source_structures(
    ["seq1", "seq2"],
    rel_type="derived_from"
)
```
