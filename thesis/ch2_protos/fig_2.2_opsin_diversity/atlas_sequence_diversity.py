#!/usr/bin/env python3
"""Atlas-scale Figure 2.2: Type II Opsin Sequence Diversity.

Produces fig_2.2a_identity_ridge.png — per-subfamily sequence identity distributions
from the LAMBDA opsin atlas (27,640 sequences, 13 subfamilies).

Data source (no ProtOS API needed, plain pandas):
  $PROTOS_LAMBDA_DATA/property/tables/type_ii_opsin_atlas_properties.csv
"""
from __future__ import annotations

import os
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
    SUBFAMILY_COLORS,
    SUBFAMILY_LABELS,
    SUBFAMILY_ORDER,
    DARK_GRAY,
    GRAY,
    apply_style,
)

REPO_ROOT = THESIS_DIR.parent
LAMBDA_DATA = Path(
    os.environ.get("PROTOS_LAMBDA_DATA", REPO_ROOT.parent / "lambda" / "data")
)
PROPERTY_CSV = (
    LAMBDA_DATA / "property" / "tables" / "type_ii_opsin_atlas_properties.csv"
)
FIGURES_DIR = Path(__file__).resolve().parent
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# The 9 subfamilies that have query references (matching GRN table & UMAP stars)
QUERY_SUBFAMILIES = [
    "rod_opsin",
    "cone_SWS1",
    "cone_MWS",
    "melanopsin",
    "neuropsin",
    "encephalopsin",
    "rgr",
    "peropsin",
    "parapinopsin",
]


def load_atlas() -> pd.DataFrame:
    df = pd.read_csv(PROPERTY_CSV)
    # Keep only the 9 subfamilies with atlas query references
    df = df[df["subfamily"].isin(QUERY_SUBFAMILIES)].copy()
    return df


# ---------------------------------------------------------------------------
# Ridge / joy plot — per-subfamily identity KDE
# ---------------------------------------------------------------------------
def make_ridge_plot(df: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt
    from scipy.stats import gaussian_kde

    apply_style()

    # Use query subfamilies in canonical order, skip any with <=1 data point
    order = [s for s in SUBFAMILY_ORDER if s in QUERY_SUBFAMILIES and df[df["subfamily"] == s].shape[0] > 1]
    n = len(order)

    fig, axes = plt.subplots(n, 1, figsize=(7, 0.7 * n), sharex=True)
    if n == 1:
        axes = [axes]

    x_grid = np.linspace(0, 100, 500)
    overlap = 0.6  # vertical overlap factor

    for i, sf in enumerate(order):
        ax = axes[i]
        vals = df.loc[df["subfamily"] == sf, "identity"].dropna().values
        if len(vals) < 2:
            continue

        kde = gaussian_kde(vals, bw_method=0.3)
        density = kde(x_grid)
        # Normalise peak to 1 for uniform visual height
        density = density / density.max()

        color = SUBFAMILY_COLORS.get(sf, GRAY)
        ax.fill_between(x_grid, density, alpha=0.40, color=color)
        ax.plot(x_grid, density, color=color, lw=1.2)

        label = SUBFAMILY_LABELS.get(sf, sf)
        count = len(vals)
        ax.text(
            -0.01, 0.5, f"{label} ({count:,})",
            transform=ax.transAxes, ha="right", va="center",
            fontsize=8, color=DARK_GRAY,
        )

        ax.set_yticks([])
        ax.set_ylim(0, 1 + overlap)
        for spine in ax.spines.values():
            spine.set_visible(False)

    axes[-1].set_xlabel("Sequence Identity to Query (%)", fontsize=10)
    axes[-1].spines["bottom"].set_visible(True)
    plt.subplots_adjust(hspace=-overlap * 0.4, left=0.25)

    out = FIGURES_DIR / "fig_2.2a_identity_ridge.png"
    fig.savefig(str(out), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 70)
    print("ATLAS FIGURE 2.2 — Opsin Sequence Diversity")
    print("=" * 70)

    df = load_atlas()
    print(f"Loaded {len(df):,} sequences across {df['subfamily'].nunique()} subfamilies")

    make_ridge_plot(df)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
