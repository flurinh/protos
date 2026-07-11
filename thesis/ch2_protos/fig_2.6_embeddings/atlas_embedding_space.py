#!/usr/bin/env python3
"""Atlas-scale Figure 2.6: Type II Opsin Embedding Space.

Produces fig_2.6a_umap_embedding.png — UMAP 2D projection of 27,639 Ankh Large
embeddings across 13 Type II opsin subfamilies.  The 9 query/reference sequences
used to build the atlas are highlighted as star markers.

Data sources:
  Property CSV:  .../lambda/data/property/tables/type_ii_opsin_atlas_properties.csv
  Embeddings:    .../lambda/data/embedding/embeddings/ankh_large/type_ii_opsin_atlas_ankh_large/*.pkl

On first run, mean-pools all 27,639 per-residue embeddings and caches the result
as an NPZ file (~160 MB).  Subsequent runs load from cache in ~1 s.
"""
from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
THESIS_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(THESIS_DIR / "shared"))
from thesis_style import (
    DARK_GRAY,
    SUBFAMILY_COLORS,
    SUBFAMILY_LABELS,
    SUBFAMILY_ORDER,
    GRAY,
    plotly_layout_defaults,
)

REPO_ROOT = THESIS_DIR.parent
LAMBDA_DATA = Path(
    os.environ.get("PROTOS_LAMBDA_DATA", REPO_ROOT.parent / "lambda" / "data")
)
PROPERTY_CSV = (
    LAMBDA_DATA / "property" / "tables" / "type_ii_opsin_atlas_properties.csv"
)
EMBEDDING_DIR = (
    LAMBDA_DATA
    / "embedding"
    / "embeddings"
    / "ankh_large"
    / "type_ii_opsin_atlas_ankh_large"
)
FIGURES_DIR = Path(__file__).resolve().parent
CACHE_DIR = THESIS_DIR / "outputs" / "atlas"
CACHE_FILE = CACHE_DIR / "type_ii_atlas_ankh_large_mean_pooled.npz"


# ---------------------------------------------------------------------------
# Load & cache embeddings
# ---------------------------------------------------------------------------
def load_embeddings() -> tuple[np.ndarray, list[str]]:
    """Return (N, 1536) float32 matrix and list of entity_ids."""
    if CACHE_FILE.exists():
        print(f"  Loading cached embeddings from {CACHE_FILE}")
        data = np.load(CACHE_FILE, allow_pickle=True)
        return data["embeddings"], data["entity_ids"].tolist()

    from tqdm import tqdm

    pkl_files = sorted(EMBEDDING_DIR.glob("*.pkl"))
    print(f"  Mean-pooling {len(pkl_files):,} per-residue embeddings...")

    entity_ids: list[str] = []
    vectors: list[np.ndarray] = []

    for fp in tqdm(pkl_files, desc="  Pooling", unit="emb"):
        # Extract entity_id from filename: {entity_id}|{family}|...
        entity_id = fp.stem.split("|")[0]
        with open(fp, "rb") as f:
            obj = pickle.load(f)
        emb = obj["embedding"]  # (L, 1536)
        vec = np.mean(emb, axis=0).astype(np.float32)
        entity_ids.append(entity_id)
        vectors.append(vec)

    mat = np.stack(vectors)  # (N, 1536)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE_FILE, embeddings=mat, entity_ids=np.array(entity_ids))
    print(f"  Cached to {CACHE_FILE} ({mat.nbytes / 1e6:.0f} MB)")
    return mat, entity_ids


def match_subfamilies(
    entity_ids: list[str], prop_df: pd.DataFrame
) -> np.ndarray:
    """Return array of subfamily strings aligned with entity_ids."""
    lookup = dict(zip(prop_df["entity_id"].astype(str), prop_df["subfamily"]))
    return np.array([lookup.get(eid, "unknown") for eid in entity_ids])


# ---------------------------------------------------------------------------
# Panel 1: UMAP 2D projection
# ---------------------------------------------------------------------------
def make_umap(
    mat: np.ndarray,
    subfamilies: np.ndarray,
    entity_ids: list[str],
    is_query: np.ndarray,
) -> None:
    import umap
    import plotly.graph_objects as go

    layout = plotly_layout_defaults()

    print("  Computing UMAP (n=30, min_dist=0.3, cosine)...")
    reducer = umap.UMAP(
        n_neighbors=30, min_dist=0.3, metric="cosine", random_state=42
    )
    coords = reducer.fit_transform(mat)

    fig = go.Figure()

    # --- Subfamily clouds (all non-query points) ---
    for sf in SUBFAMILY_ORDER:
        mask = (subfamilies == sf) & ~is_query
        if not mask.any():
            continue
        label = SUBFAMILY_LABELS.get(sf, sf)
        n_total = int((subfamilies == sf).sum())
        fig.add_trace(go.Scattergl(
            x=coords[mask, 0],
            y=coords[mask, 1],
            mode="markers",
            name=f"{label} ({n_total:,})",
            marker=dict(
                size=3,
                color=SUBFAMILY_COLORS.get(sf, GRAY),
                opacity=0.5,
            ),
        ))

    # --- Query / reference sequences as stars ---
    query_mask = is_query
    if query_mask.any():
        query_sfs = subfamilies[query_mask]
        query_colors = [SUBFAMILY_COLORS.get(sf, GRAY) for sf in query_sfs]
        query_names = [eid for eid, q in zip(entity_ids, is_query) if q]

        fig.add_trace(go.Scattergl(
            x=coords[query_mask, 0],
            y=coords[query_mask, 1],
            mode="markers",
            name=f"Query references ({int(query_mask.sum())})",
            marker=dict(
                size=12,
                color=query_colors,
                symbol="star",
                opacity=1.0,
                line=dict(color=DARK_GRAY, width=1),
            ),
            text=query_names,
            hovertemplate="%{text}<extra>Query</extra>",
        ))
        print(f"  Highlighted {int(query_mask.sum())} query sequences as stars")

    fig.update_layout(
        xaxis=dict(showgrid=False, showticklabels=False, zeroline=False, title=""),
        yaxis=dict(showgrid=False, showticklabels=False, zeroline=False, title=""),
        width=900,
        height=700,
        legend=dict(
            font_size=9,
            itemsizing="constant",
            x=1.02,
            y=1,
        ),
        **{**layout, "margin": dict(t=20, b=20, l=20, r=180)},
    )

    out = FIGURES_DIR / "fig_2.6a_umap_embedding.png"
    fig.write_image(str(out), scale=2)
    print(f"  Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 70)
    print("ATLAS FIGURE 2.6 — Opsin Embedding Space")
    print("=" * 70)

    # Load property table
    prop_df = pd.read_csv(PROPERTY_CSV)
    print(f"Property table: {len(prop_df):,} rows")

    # Load / cache embeddings
    mat, entity_ids = load_embeddings()
    print(f"Embedding matrix: {mat.shape}")

    # Match subfamilies and query status
    subfamilies = match_subfamilies(entity_ids, prop_df)
    known = subfamilies != "unknown"
    print(f"Matched subfamilies: {known.sum():,}/{len(subfamilies):,}")

    # Build is_query lookup from property table
    query_ids = set(
        prop_df.loc[prop_df["is_query"] == True, "entity_id"].astype(str)
    )
    is_query_all = np.array([eid in query_ids for eid in entity_ids])
    print(f"Query sequences in embeddings: {is_query_all.sum()}")

    # Filter to known subfamilies only
    mat = mat[known]
    entity_ids = [eid for eid, k in zip(entity_ids, known) if k]
    subfamilies = subfamilies[known]
    is_query_arr = is_query_all[known]

    make_umap(mat, subfamilies, entity_ids, is_query_arr)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
