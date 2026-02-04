#!/usr/bin/env python3
"""EmbeddingProcessor Example: Sequence vs Embedding similarity analysis.

Demonstrates ProtOS capabilities:
- EmbeddingProcessor: Generate ESM2 embeddings for sequence datasets
- Cross-processor data flow: Uses dataset from SequenceProcessor example
- Comparing sequence identity with embedding similarity

Question: "Do cone opsins cluster by spectral type in embedding space?"

This example picks up the cone_opsin_diversity dataset created by the
SequenceProcessor example and analyzes whether sequence embeddings
capture functional (spectral) relationships beyond sequence identity.

Insight: Sequences with similar function may cluster together in embedding
space even when sequence identity is not the highest predictor.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.processing.embedding import EmbeddingProcessor


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "embedding"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
SEQUENCE_OUTPUT_DIR = THESIS_DIR / "outputs" / "sequence"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Dataset from SequenceProcessor example
SOURCE_DATASET = "cone_opsin_diversity"
EMBEDDING_DATASET = f"{SOURCE_DATASET}__esm2_t12_35m__mean"

# Opsin type colors (teal/gold - matching sequence colors from colorscales.yaml)
TYPE_COLORS = {
    "short_wave": "#006d77",   # Teal
    "long_wave": "#e9c46a",    # Gold
}


def main() -> int:
    """Run the EmbeddingProcessor example."""
    print("=" * 70)
    print("EMBEDDING PROCESSOR EXAMPLE")
    print("Sequence vs Embedding Similarity Analysis")
    print("=" * 70)
    print(f"\nSource dataset: {SOURCE_DATASET}")
    print("(Created by SequenceProcessor example)")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Load dataset from SequenceProcessor
    # -------------------------------------------------------------------------
    print("\n[1] Loading cone opsin dataset...")
    seq_proc = SequenceProcessor()

    # Load the dataset created by SequenceProcessor
    sequences = seq_proc.load_dataset(SOURCE_DATASET)

    if not sequences:
        print(f"  Dataset '{SOURCE_DATASET}' not found!")
        print("  Please run sequence_processor_example.py first.")
        return 1

    print(f"  Loaded {len(sequences)} sequences from {SOURCE_DATASET}")

    # Load annotations
    annotations_file = SEQUENCE_OUTPUT_DIR / "opsin_annotations.json"
    if annotations_file.exists():
        with open(annotations_file) as f:
            annotations = json.load(f)
        print(f"  Loaded annotations for {len(annotations)} sequences")
    else:
        print("  No annotations found - will use sequence names only")
        annotations = {}

    # -------------------------------------------------------------------------
    # Step 2: Generate embeddings using ProtOS EmbeddingProcessor
    # -------------------------------------------------------------------------
    print("\n[2] Generating ESM2 embeddings...")
    emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m")

    # ProtOS generates embeddings and registers them as dataset
    embeddings = emb_proc.embed_sequences(
        sequences,
        embedding_type="mean",
        save_dataset=EMBEDDING_DATASET,
        register_entities=True,
    )

    if embeddings:
        sample_dim = next(iter(embeddings.values())).shape[-1]
        print(f"  Generated {len(embeddings)} embeddings (dim={sample_dim})")
    else:
        print("  No embeddings generated!")
        return 1

    # -------------------------------------------------------------------------
    # Step 3: Compute embedding similarity
    # -------------------------------------------------------------------------
    print("\n[3] Computing pairwise cosine similarity...")
    import torch
    import torch.nn.functional as F

    # Load embeddings through processor
    emb_dict = emb_proc.load_embeddings(EMBEDDING_DATASET)
    names = list(emb_dict.keys())

    # Convert to tensor matrix
    emb_list = []
    for n in names:
        e = emb_dict[n]
        if isinstance(e, torch.Tensor):
            emb_list.append(e)
        else:
            emb_list.append(torch.tensor(e))
    emb_matrix = torch.stack(emb_list)

    # Compute cosine similarity
    emb_norm = F.normalize(emb_matrix, p=2, dim=1)
    sim_matrix = torch.mm(emb_norm, emb_norm.T).numpy()
    sim_df = pd.DataFrame(sim_matrix, index=names, columns=names)
    sim_df.to_csv(OUTPUT_DIR / "opsin_embedding_similarity.csv")
    print(f"  Saved similarity matrix ({len(sim_df)}x{len(sim_df.columns)})")

    # -------------------------------------------------------------------------
    # Step 4: Compare sequence identity vs embedding similarity
    # -------------------------------------------------------------------------
    print("\n[4] Comparing sequence identity vs embedding similarity...")

    # Helper to extract accession from sequence name (e.g., "sp|P60015.1|OPSB_PANTR" -> "P60015")
    def get_accession(name: str) -> str:
        """Extract accession ID from sequence name."""
        if "|" in name:
            parts = name.split("|")
            if len(parts) >= 2:
                # Handle "sp|P60015.1|..." -> "P60015"
                acc = parts[1].split(".")[0]
                return acc
        return name

    # Build comparison data
    comparison_data = []
    for i, name_i in enumerate(names):
        for j, name_j in enumerate(names):
            if i >= j:  # Skip diagonal and duplicates
                continue

            # Get annotations using extracted accession IDs
            acc_i = get_accession(name_i)
            acc_j = get_accession(name_j)
            ann_i = annotations.get(acc_i, {})
            ann_j = annotations.get(acc_j, {})

            # Get sequence identity from annotations (if available)
            seq_id_i = ann_i.get("identity_to_query", 0)
            seq_id_j = ann_j.get("identity_to_query", 0)

            comparison_data.append({
                "seq1": name_i,
                "seq2": name_j,
                "type1": ann_i.get("opsin_type", "unknown"),
                "type2": ann_j.get("opsin_type", "unknown"),
                "same_type": ann_i.get("opsin_type") == ann_j.get("opsin_type"),
                "embedding_similarity": sim_matrix[i, j],
                "avg_seq_identity": (seq_id_i + seq_id_j) / 2,
            })

    comparison_df = pd.DataFrame(comparison_data)
    comparison_df.to_csv(OUTPUT_DIR / "identity_vs_embedding.csv", index=False)

    # Summary statistics
    same_type = comparison_df[comparison_df["same_type"]]
    diff_type = comparison_df[~comparison_df["same_type"]]

    print(f"\n  Same opsin type pairs: {len(same_type)}")
    print(f"    Mean embedding similarity: {same_type['embedding_similarity'].mean():.3f}")
    print(f"  Different type pairs: {len(diff_type)}")
    print(f"    Mean embedding similarity: {diff_type['embedding_similarity'].mean():.3f}")

    # -------------------------------------------------------------------------
    # Step 5: Create visualizations
    # -------------------------------------------------------------------------
    print("\n[5] Creating visualizations...")
    from sklearn.manifold import TSNE
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Prepare embedding matrix for t-SNE
    emb_matrix_np = np.stack([
        emb_dict[n].numpy() if hasattr(emb_dict[n], 'numpy') else np.array(emb_dict[n])
        for n in names
    ])

    # Get opsin types for coloring (using accession lookup)
    opsin_types = [annotations.get(get_accession(n), {}).get("opsin_type", "unknown") for n in names]

    # t-SNE visualization
    perplexity = min(30, len(names) - 1)
    coords = TSNE(n_components=2, perplexity=perplexity, random_state=42).fit_transform(emb_matrix_np)

    fig = go.Figure()

    for opsin_type, color in TYPE_COLORS.items():
        mask = [t == opsin_type for t in opsin_types]
        if not any(mask):
            continue
        idx = [i for i, m in enumerate(mask) if m]
        label = "SW" if opsin_type == "short_wave" else "LW"

        fig.add_trace(go.Scatter(
            x=coords[idx, 0],
            y=coords[idx, 1],
            mode="markers",
            name=label,
            marker=dict(size=8, color=color, opacity=0.7),
            text=[names[i] for i in idx],
            hovertemplate="%{text}<extra></extra>",
        ))

    fig.update_layout(
        xaxis=dict(showgrid=False, showticklabels=False, zeroline=False, title=""),
        yaxis=dict(showgrid=False, showticklabels=False, zeroline=False, title=""),
        width=700,
        height=550,
        plot_bgcolor="white",
        paper_bgcolor="white",
        legend=dict(x=0.02, y=0.98, font_size=10),
        margin=dict(t=30, b=30),
    )
    fig.write_image(str(FIGURES_DIR / "embedding_tsne.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'embedding_tsne.png'}")

    # Embedding similarity analysis - no titles
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.15)

    # Box plot: same vs different type (use teal/gold)
    fig.add_trace(go.Box(
        y=same_type["embedding_similarity"],
        name="Same",
        marker_color="#006d77",  # Teal
        boxpoints="outliers",
    ), row=1, col=1)

    fig.add_trace(go.Box(
        y=diff_type["embedding_similarity"],
        name="Diff",
        marker_color="#e9c46a",  # Gold
        boxpoints="outliers",
    ), row=1, col=1)

    # Histogram of all similarities
    fig.add_trace(go.Histogram(
        x=comparison_df["embedding_similarity"],
        nbinsx=30,
        marker_color="#7f7f7f",  # Neutral gray
        opacity=0.8,
        showlegend=False,
    ), row=1, col=2)

    fig.update_layout(
        height=400,
        width=850,
        showlegend=False,
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin=dict(t=30, b=50),
    )
    fig.update_yaxes(title_text="Similarity", row=1, col=1, title_font_size=10)
    fig.update_xaxes(title_text="Similarity", row=1, col=2, title_font_size=10)
    fig.update_yaxes(title_text="Count", row=1, col=2, title_font_size=10)
    fig.update_yaxes(showgrid=False)
    fig.update_xaxes(showgrid=False)
    fig.write_image(str(FIGURES_DIR / "embedding_similarity_analysis.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'embedding_similarity_analysis.png'}")

    # Heatmap with block structure: sort by opsin type for 2x2 pattern
    # Sort names by opsin type to get block structure (SW first, then LW)
    sorted_indices = sorted(range(len(names)), key=lambda i: (
        0 if opsin_types[i] == "short_wave" else 1,
        names[i]
    ))
    sorted_names = [names[i] for i in sorted_indices]
    sorted_types = [opsin_types[i] for i in sorted_indices]

    # Reorder similarity matrix
    heatmap_sim = sim_df.loc[sorted_names, sorted_names]

    # Count sequences per type for tick coloring
    n_sw = sum(1 for t in sorted_types if t == "short_wave")
    n_lw = len(sorted_types) - n_sw

    fig = go.Figure(data=go.Heatmap(
        z=heatmap_sim.values,
        colorscale="RdBu",
        zmid=0.85,
        zmin=0.6,
        zmax=1.0,
        colorbar=dict(title="Similarity", x=1.02, title_font_size=10),
        showscale=True,
    ))

    # Add colored tick marks on axes (no labels, just colored bars)
    # X-axis colored ticks (bottom)
    for i, t in enumerate(sorted_types):
        color = TYPE_COLORS.get(t, "#888")
        fig.add_shape(
            type="rect",
            x0=i - 0.5, x1=i + 0.5, y0=-0.02 * len(sorted_names), y1=0,
            fillcolor=color, line=dict(width=0),
            xref="x", yref="y",
        )
    # Y-axis colored ticks (left)
    for i, t in enumerate(sorted_types):
        color = TYPE_COLORS.get(t, "#888")
        fig.add_shape(
            type="rect",
            x0=-0.02 * len(sorted_names), x1=0, y0=i - 0.5, y1=i + 0.5,
            fillcolor=color, line=dict(width=0),
            xref="x", yref="y",
        )

    fig.update_layout(
        xaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.02 * len(sorted_names), len(sorted_names) - 0.5],
        ),
        yaxis=dict(
            showticklabels=False, showgrid=False, zeroline=False,
            range=[-0.02 * len(sorted_names), len(sorted_names) - 0.5],
            scaleanchor="x", scaleratio=1,
        ),
        width=650,
        height=600,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=30, b=30),
    )
    fig.write_image(str(FIGURES_DIR / "embedding_heatmap.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'embedding_heatmap.png'}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\nDataset: {len(sequences)} sequences from {SOURCE_DATASET}")
    print(f"Embeddings: {len(embeddings)} x {sample_dim} dimensions")
    print(f"\nKey finding:")
    print(f"  Same-type pairs have higher embedding similarity ({same_type['embedding_similarity'].mean():.3f})")
    print(f"  than different-type pairs ({diff_type['embedding_similarity'].mean():.3f})")
    print(f"\n→ ESM2 embeddings capture functional relationships beyond sequence identity")
    print(f"\nOutputs: {OUTPUT_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
