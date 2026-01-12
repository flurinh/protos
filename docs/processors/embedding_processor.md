# EmbeddingProcessor

Generates and stores protein sequence embeddings using transformer models.

**Location:** `protos.processing.embedding.embedding_processor`

**Processor Type:** `embedding`

**Dependencies:** torch, transformers (optional)

```bash
# Install with GPU support
pip install --no-build-isolation -e ".[gpu]"

# Install CPU-only
pip install -e ".[embedding]"
```

## Overview

The EmbeddingProcessor provides:
- Embedding generation using ESM-2, Ankh, and other PLMs
- Multiple pooling strategies (mean, sum, cls, per-residue)
- Batch processing for efficiency
- Automatic GPU detection and usage

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(entity_id, format_type)` | Load saved embedding for entity. |
| `save_entity(entity_id, data, format_type)` | Save embedding and register. |
| `load_embedding_entity(identifier)` | Load embedding as torch tensor. |
| `list_entities(dataset)` | List registered embedding entities. |
| `list_embedding_entities(dataset)` | Alias for list_entities. |

### Embedding Generation

| Method | Description |
|--------|-------------|
| `embed_sequences(sequences, embedding_type, batch_size, dataset_name)` | Generate embeddings for sequence dict. |
| `embed_from_fasta(fasta_path, embedding_type, dataset_name)` | Generate embeddings from FASTA file. |
| `get_residue_embeddings(embeddings, include_special_tokens)` | Extract per-residue embeddings. |
| `collapse_per_residue(embeddings, method)` | Collapse per-residue to single vector. |

### Dataset Operations

| Method | Description |
|--------|-------------|
| `load_embeddings(dataset_name)` | Load all embeddings in dataset. |
| `ingest_from_artifact(artifact_path, dataset_name)` | Import embeddings from NPZ artifact. |
| `ingest_from_invocation(invocation, dataset_name)` | Import from ModelInvocation result. |

### Model Management

| Method | Description |
|--------|-------------|
| `available_models()` | Class method: list available model configurations. |
| `list_available_models()` | Instance method: list models with details. |
| `get_embedding_dim()` | Get embedding dimension for current model. |
| `check_dependencies()` | Check if torch/transformers available. |
| `clear_cache()` | Clear loaded model from memory. |

### Properties

| Property | Description |
|----------|-------------|
| `model` | Lazy-loaded transformer model. |
| `tokenizer` | Lazy-loaded tokenizer. |
| `dependencies_available` | Whether torch/transformers are installed. |
| `embeddings_dir` | Directory for saved embeddings. |
| `datasets_dir` | Directory for dataset definitions. |

---

## Available Models

| Model Name | Parameters | Embedding Dim | Description |
|------------|------------|---------------|-------------|
| `esm2_t6_8m` | 8M | 320 | ESM-2 tiny |
| `esm2_t12_35m` | 35M | 480 | ESM-2 small |
| `esm2_t30_150m` | 150M | 640 | ESM-2 medium |
| `esm2_t33_650m` | 650M | 1280 | ESM-2 large |
| `esm2_t36_3b` | 3B | 2560 | ESM-2 extra large |
| `esm2_t48_15b` | 15B | 5120 | ESM-2 huge |
| `ankh_base` | - | 768 | Ankh base |
| `ankh_large` | - | 1536 | Ankh large |

---

## Embedding Types

| Type | Description | Output Shape |
|------|-------------|--------------|
| `mean` | Mean pooling over sequence | `(embedding_dim,)` |
| `sum` | Sum pooling over sequence | `(embedding_dim,)` |
| `cls` | CLS token embedding | `(embedding_dim,)` |
| `per_residue` | Per-residue embeddings | `(seq_len, embedding_dim)` |

---

## Usage Examples

### Basic Embedding Generation

```python
from protos.processing.embedding import EmbeddingProcessor

proc = EmbeddingProcessor(
    model_name="esm2_t30_150m",
    device="cuda",
    batch_size=8
)

# Check available models
models = proc.list_available_models()
for name, info in models.items():
    print(f"{name}: dim={info['embedding_dim']}")

# Generate single embedding
sequences = {"protein_1": "MKWVTFISLLLLFSSAYS..."}
results = proc.embed_sequences(sequences, embedding_type="mean")

# Access embedding
embedding = results["protein_1"]
print(f"Shape: {embedding.shape}")  # (640,)
```

### Batch Processing

```python
# Process multiple sequences
sequences = {
    "seq_1": "MKWVTFIS...",
    "seq_2": "MLKFTISA...",
    "seq_3": "MKVLTDIS..."
}

embeddings = proc.embed_sequences(
    sequences,
    embedding_type="mean",
    batch_size=16,
    dataset_name="my_embeddings"
)

# Load from FASTA
embeddings = proc.embed_from_fasta(
    "proteins.fasta",
    embedding_type="mean",
    dataset_name="fasta_embeddings"
)
```

### Per-Residue Embeddings

```python
# Get per-residue embeddings
results = proc.embed_sequences(
    {"protein_1": "MKWVTFIS..."},
    embedding_type="per_residue"
)

per_res = results["protein_1"]
print(f"Shape: {per_res.shape}")  # (seq_len, 640)

# Collapse to single vector
collapsed = proc.collapse_per_residue(per_res, method="mean")
print(f"Collapsed shape: {collapsed.shape}")  # (640,)
```

### Saving and Loading

```python
# Embeddings are auto-saved when dataset_name provided
embeddings = proc.embed_sequences(
    sequences,
    dataset_name="my_study"
)

# Load saved embeddings
loaded = proc.load_embeddings("my_study")

# Load single entity
emb = proc.load_entity("protein_1")
```

### Import from Artifacts

```python
# Import from NPZ file
proc.ingest_from_artifact(
    "embeddings.npz",
    dataset_name="imported"
)

# Import from ModelManager invocation
proc.ingest_from_invocation(invocation, dataset_name="model_output")
```

### GPU Memory Management

```python
import torch

# Check GPU
if torch.cuda.is_available():
    proc = EmbeddingProcessor(
        model_name="esm2_t33_650m",
        device="cuda",
        batch_size=4  # Smaller for large models
    )
else:
    proc = EmbeddingProcessor(
        model_name="esm2_t12_35m",
        device="cpu"
    )

# Clear model from memory
proc.clear_cache()
torch.cuda.empty_cache()
```

### Dependency Handling

```python
proc = EmbeddingProcessor()

# Check dependencies
deps = proc.check_dependencies()
print(f"PyTorch: {deps['torch']}")
print(f"Transformers: {deps['transformers']}")

if proc.dependencies_available:
    # Full functionality
    embeddings = proc.embed_sequences(sequences)
else:
    # Limited - can still load pre-computed
    embeddings = proc.load_embeddings("precomputed")
```

### Processor Configuration

```python
proc = EmbeddingProcessor(
    name="my_processor",
    model_name="esm2_t30_150m",
    device="cuda",
    batch_size=8,
    max_seq_length=1022
)

# View configuration
print(proc.metadata)
# {
#     'model_name': 'esm2_t30_150m',
#     'device': 'cuda',
#     'batch_size': 8,
#     'embedding_dim': 640,
#     ...
# }

# Get embedding dimension
dim = proc.get_embedding_dim()
print(f"Embedding dim: {dim}")  # 640
```
