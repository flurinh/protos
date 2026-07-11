#!/usr/bin/env python3
"""Minimal EmbeddingProcessor test — matches draft code example."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
protos.set_data_path(str(REPO_ROOT / "data"))

from protos.processing.embedding import EmbeddingProcessor

emb_proc = EmbeddingProcessor(model_name="ankh_large")

# Embed a small set of sequences
sequences = {"test_seq": "MNGTEGPNFYVPFSNATGVVR"}
embeddings = emb_proc.embed_sequences(
    sequences,
    embedding_type="per_residue",
    save_dataset="test_ankh_large",
)
print(f"Embedding shape: {list(embeddings.values())[0].shape}")
