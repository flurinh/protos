"""Thesis figure styling — shared across all chapters.

Loads colorscales.yaml and provides:
  - COLORS dict: all color hex values by category
  - apply_style(): sets matplotlib rcParams for consistent figures
  - Plotly layout helpers

Usage:
    import sys; sys.path.insert(0, str(Path(__file__).parent))
    from thesis_style import COLORS, apply_style
    apply_style()
"""
from __future__ import annotations

from pathlib import Path

import yaml

# =============================================================================
# Load colorscales.yaml
# =============================================================================
_THESIS_DIR = Path(__file__).resolve().parent.parent
_COLORSCALES = _THESIS_DIR / "colorscales.yaml"

with open(_COLORSCALES) as _f:
    COLORS = yaml.safe_load(_f)

# =============================================================================
# Convenience accessors
# =============================================================================
# Structures
SLATE = COLORS["structures"]["type_i"]["hex"]           # #3d5a80
SLATE_LIGHT = COLORS["structures"]["type_i_light"]["hex"]  # #98c1d9
TERRACOTTA = COLORS["structures"]["type_ii"]["hex"]      # #c1666b
TERRACOTTA_LIGHT = COLORS["structures"]["type_ii_light"]["hex"]  # #e4c1c1

# Sequences / Embeddings
SAGE = COLORS["sequences"]["short_wave"]["hex"]          # #457b6b
OCHRE = COLORS["sequences"]["long_wave"]["hex"]          # #d4a03c

# Molecules
RUST = COLORS["molecules"]["retinal"]["hex"]             # #bc6c25
MAUVE = COLORS["molecules"]["carazolol"]["hex"]          # #7a5980

# Neutral
GRAY = COLORS["neutral"]["gray"]["hex"]                  # #6c757d
LIGHT_GRAY = COLORS["neutral"]["light_gray"]["hex"]      # #adb5bd
DARK_GRAY = COLORS["neutral"]["dark_gray"]["hex"]        # #343a40
WARM_WHITE = COLORS["neutral"]["warm_white"]["hex"]      # #f7f5f3

# Figure elements
GRID_COLOR = COLORS["figure_elements"]["grid"]["hex"]    # #e9ecef

# Per-structure colors
STRUCTURE_COLORS = {
    "1C3W": COLORS["structures"]["1C3W"]["hex"],
    "3UG9": COLORS["structures"]["3UG9"]["hex"],
    "1U19": COLORS["structures"]["1U19"]["hex"],
    "2Z73": COLORS["structures"]["2Z73"]["hex"],
}

# Sequence type colors
TYPE_COLORS = {
    "short_wave": SAGE,
    "long_wave": OCHRE,
}

# Subfamily colors (Type II Opsin Atlas — 13 subfamilies)
SUBFAMILY_COLORS = {k: v["hex"] for k, v in COLORS["subfamilies"].items()}

# Canonical display order (descending by typical count in the atlas)
SUBFAMILY_ORDER = [
    "rod_opsin",
    "cone_LWS",
    "cone_SWS1",
    "cone_MWS",
    "melanopsin",
    "neuropsin",
    "encephalopsin",
    "rgr",
    "peropsin",
    "parapinopsin",
    "pinopsin",
    "vertebrate_ancient",
    "parietopsin",
]

# Display-friendly labels
SUBFAMILY_LABELS = {k: v["name"] for k, v in COLORS["subfamilies"].items()}


def apply_style():
    """Apply thesis-consistent matplotlib rcParams."""
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        # Figure
        "figure.facecolor": "white",
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",

        # Axes
        "axes.facecolor": "white",
        "axes.edgecolor": DARK_GRAY,
        "axes.labelcolor": DARK_GRAY,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,

        # Ticks
        "xtick.color": DARK_GRAY,
        "ytick.color": DARK_GRAY,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,

        # Grid
        "axes.grid": False,
        "grid.color": GRID_COLOR,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.8,

        # Font
        "font.size": 10,
        "font.family": "sans-serif",

        # Legend
        "legend.fontsize": 9,
        "legend.frameon": True,
        "legend.edgecolor": LIGHT_GRAY,
        "legend.fancybox": False,
    })


def plotly_layout_defaults() -> dict:
    """Return common plotly layout kwargs for thesis figures."""
    return {
        "paper_bgcolor": "white",
        "plot_bgcolor": "white",
        "font": {"color": DARK_GRAY, "size": 11},
        "margin": dict(t=30, b=50, l=60, r=30),
    }
