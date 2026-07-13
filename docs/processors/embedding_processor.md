# EmbeddingProcessor

Import: `from protos.processing.embedding import EmbeddingProcessor`

`EmbeddingProcessor` can generate protein-language-model embeddings or ingest
existing NPZ arrays. Model weights are loaded lazily from Hugging Face when
inference is first requested.

## Dependencies

The generic `embedding` extra installs Transformers for ESM-2 and Ankh but does
not install PyTorch. ESMC has a separate `esmc` extra with the exact Biohub
fork and `accelerate`.
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

The listed model registry contains ESM-2, Ankh, and the pinned ESMC-6B
identifier. A listed model is a configuration entry, not a guarantee that its
remote weights are already cached or that the machine has enough memory to load
it.

## Pinned Biohub ESMC-6B residue embeddings

`esmc_6b` is a model-specific `EmbeddingProcessor` path; it does not use
`ModelManager`, jobs, or the generic ESM-2/Ankh tokenizer path. Install a
machine-appropriate PyTorch build, then the ESMC dependencies:

```bash
python -m pip install torch
python -m pip install -e ".[esmc]"
```

The public model configuration is intentionally pinned:

- model ID `biohub/ESMC-6B`;
- Hugging Face revision
  `45b0fa5d7fb06faefbd5e3b89bdcef35d564e79a`;
- Biohub source-lineage commit
  `ba4d7124864eed323a93bf3cfefcd958f573b75a`;
- 80 layers and `d_model=2560`;
- `output.last_hidden_state` from the final post-layernorm encoder state;
- bfloat16 inference, float32 artifact storage, and the exact PEP 508 VCS pin
  `transformers @ git+https://github.com/Biohub/transformers.git@ef32577f55da19a4989cd7b22e004dc43a4998cb`.

CUDA loading uses `device_map="auto"` with a 28 GiB GPU ceiling and 64 GiB
CPU allowance. This leaves headroom on a 32-GiB RTX 5090 and permits CPU
offload. CPU loading uses the same bfloat16 compute policy. The full repository
is roughly 25.4 GB, so construction remains download-free and live validation
is opt-in.

```python
from protos.processing.embedding import (
    EmbeddingProcessor,
    ESMCSourceLineage,
)

processor = EmbeddingProcessor(
    model_name="esmc_6b",
    device="cuda",
    batch_size=2,  # keep token batches conservative for the 6B model
)
results = processor.embed_esmc_sequences(
    {
        "protein-a": "ACDEFGHIK",
        "protein-b": "MNPQRSTVWY",
    },
    source_lineage={
        "protein-a": ESMCSourceLineage(
            source_kind="fasta",
            source_id="training.fasta",
            source_revision="release-v2",
            source_sha256="...",
        ),
        "protein-b": ESMCSourceLineage(
            source_kind="fasta",
            source_id="training.fasta",
            source_revision="release-v2",
            source_sha256="...",
        ),
    },
    resume=True,
)

array = results["protein-a"].embedding
assert array.shape == (9, 2560)
assert array.dtype.name == "float32"
```

The return value maps each input ID, in input order, to an
`ESMCEmbeddingResult` containing its array, artifact path, exact metadata, and
whether it was reused. `ESMCModelProvenance`, `ESMCLoadPolicy`,
`ESMCTokenPolicy`, `ESMCOutputProvenance`, and `ESMCSourceLineage` are immutable
typed configuration/provenance records. A narrow `ESMCBackend` protocol is
public for deterministic offline testing; production callers should omit the
`backend` argument to use the pinned direct Hugging Face loader.

Sequences are normalized by removing whitespace and uppercasing before their
SHA-256 is calculated. The adapter disables truncation, requests the special
token mask, verifies distinct BOS/EOS/pad/unknown IDs and exact token placement,
rejects unknown residue tokens, and retains exactly one row per normalized
residue. Missing or ambiguous token mappings, non-final outputs, wrong
dimensions/shapes, and non-finite values fail rather than being stored.

Artifacts are atomically written under the model-specific `esmc_6b` namespace.
Their identity includes the normalized sequence and hash, entity ID, model and
code revisions, layer/dimension, bfloat16/device-map load policy, token policy,
output type and float32 storage dtype, and caller source lineage. Resume loads
only an exact identity match and revalidates stored provenance, shape, dtype,
token mapping, and embedding content SHA-256. Conflicting sequences for one
entity are rejected. Ankh, other ESMC tiers, other adapter/revision, and other
policy artifacts cannot satisfy a 6B resume lookup.

Security and provenance assumptions: model and tokenizer files are retrieved
only from the pinned public Hugging Face revision with
`trust_remote_code=False`; the Biohub repository commit is both an installed
runtime dependency and typed provenance. Stock Transformers, another fork SHA,
ESMC-600M, Ankh, and the forbidden Biohub/esm package are incompatible.
Source-lineage values are caller assertions and
should themselves be derived from reviewed, hashed inputs. Cache integrity is
checked with SHA-256, but this local artifact store does not provide signatures,
access control, or protection from an attacker who can rewrite both arrays and
metadata.

The two-device live smoke is never part of the default suite. To opt in on a
host with the dependencies, network, memory, disk, and (for CUDA) at least
28 GiB free VRAM:

```bash
PROTOS_RUN_ESMC_6B_LIVE=1 pytest -q \
  tests/test_esmc_adapter.py -k live_esmc_6b_two_sequence_smoke
```

### Python 3.11 production preflight

Run `scripts/esmc_python311_preflight.sh` to create a disposable environment,
install the exact fork and `accelerate`, fetch only the pinned model config,
and prove AutoConfig plus AutoModelForMaskedLM registration and loader-argument
construction without allocating or downloading weights. On 2026-07-13 it
completed with Python 3.11.14, Transformers 4.57.6 (Biohub commit
`ef32577f55da19a4989cd7b22e004dc43a4998cb`), accelerate 1.14.0, model type
`esmc`, registered class `ESMCForMaskedLM`, and `weights: not loaded`.
