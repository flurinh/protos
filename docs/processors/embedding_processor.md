# EmbeddingProcessor

Import: `from protos.processing.embedding import EmbeddingProcessor`

`EmbeddingProcessor` can generate protein-language-model embeddings or ingest
existing NPZ arrays. Model weights are loaded lazily from Hugging Face when
inference is first requested.

## Dependencies

The `embedding` extra installs Transformers but does not install PyTorch.
Install a PyTorch build suitable for the machine, then install ProtOS with the
extra:

```bash
python -m pip install torch
python -m pip install -e ".[embedding]"
```

The exact PyTorch installation command may differ for CUDA builds.

## Inspect support without loading a model

```python
from protos.processing.embedding import EmbeddingProcessor

models = EmbeddingProcessor.available_models()
assert models["esm2_t6_8m"]["embedding_dim"] == 320

processor = EmbeddingProcessor(model_name="esm2_t6_8m", device="cpu")
status = processor.check_dependencies()
assert set(status) == {"torch", "transformers", "ready"}
assert processor.get_embedding_dim() == 320
```

Construction does not download weights. Accessing `model`/`tokenizer` or
calling `embed_sequences()` does.

## Inference and storage

`embed_sequences()` accepts one sequence, a list, or an ID-to-sequence mapping.
Supported output types are `mean`, `sum`, `cls`, and `per_residue`. Passing
`save_dataset="name"` persists per-entity pickle files and a dataset definition.

`load_embeddings(dataset_name)` loads a persisted dataset.
`ingest_from_artifact(dataset_name, artifact_path, ...)` imports every array in
an NPZ file without running a model. `collapse_per_residue(array,
reduction="mean"|"sum")` aggregates an existing per-residue array.

The listed model registry currently contains ESM-2 and Ankh identifiers. A
listed model is a configuration entry, not a guarantee that its remote weights
are already cached or that the machine has enough memory to load it.
