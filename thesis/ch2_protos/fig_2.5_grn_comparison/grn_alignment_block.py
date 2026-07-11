#!/usr/bin/env python3
"""Figure 2.5a: GRN Microswitch & Spectral Tuning Table.

Produces fig_2.5a_alignment_block.png — a table of amino acid identities at key
GRN positions across the 9 Type II opsin reference sequences (atlas queries)
used to build the LAMBDA atlas.

Rows   = reference opsins (one per subfamily query)
Columns = functionally important GRN positions (microswitches + spectral tuning)
Cells  = amino acid letter, background colored by physicochemical property

Data source:
  .../lambda/data/raw_data/opsin_phylogeny/type_ii_grns.csv
  (GRN annotations generated during atlas construction by step1_collect_opsins.py)
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
    COLORS,
    DARK_GRAY,
    GRAY,
    LIGHT_GRAY,
    SUBFAMILY_COLORS,
    apply_style,
)

REPO_ROOT = THESIS_DIR.parent
LAMBDA_DATA = Path(
    os.environ.get("PROTOS_LAMBDA_DATA", REPO_ROOT.parent / "lambda" / "data")
)
GRN_CSV = LAMBDA_DATA / "raw_data" / "opsin_phylogeny" / "type_ii_grns.csv"
FIGURES_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Reference sequences (atlas queries, UniProt accession → display name)
# Order: visual opsins first, then non-visual, matching SUBFAMILY_ORDER
# ---------------------------------------------------------------------------
REFERENCES = [
    ("P02699", "rod_opsin",      "Rhodopsin",       "B. taurus"),
    ("P04001", "cone_MWS",       "MWS Opsin",       "H. sapiens"),
    ("P03999", "cone_SWS1",      "SWS1 Opsin",      "H. sapiens"),
    ("O42266", "parapinopsin",   "Parapinopsin",    "I. punctatus"),
    ("Q9H1Y3", "encephalopsin",  "Encephalopsin",   "H. sapiens"),
    ("Q9UHM6", "melanopsin",     "Melanopsin",      "H. sapiens"),
    ("Q6U736", "neuropsin",      "Neuropsin",       "H. sapiens"),
    ("O14718", "peropsin",       "Peropsin",        "H. sapiens"),
    ("P47804", "rgr",            "RGR",             "H. sapiens"),
]


# ---------------------------------------------------------------------------
# GRN positions: (grn, label, category)
# category: "tuning" = spectral tuning, "switch" = microswitch
# ---------------------------------------------------------------------------
GRN_POSITIONS = [
    # TM1
    ("1.50",  "Na\u207a/water", "switch"),
    # TM2
    ("2.50",  "Na\u207a/water", "switch"),
    # TM3 — spectral tuning + E/DRY + PIF
    ("3.28",  "Counterion",     "tuning"),
    ("3.32",  "Tuning",         "tuning"),
    ("3.40",  "PIF",            "switch"),
    ("3.49",  "E/DRY",          "switch"),
    ("3.50",  "E/DRY",          "switch"),
    ("3.51",  "E/DRY",          "switch"),
    # TM4
    ("4.50",  "Core lock",      "switch"),
    # TM5
    ("5.50",  "PIF",            "switch"),
    # TM6 — transmission + CWxP
    ("6.44",  "F\u2011switch",  "switch"),
    ("6.47",  "CWxP",           "switch"),
    ("6.48",  "CWxP",           "switch"),
    ("6.50",  "CWxP",           "switch"),
    # TM7 — Schiff base + NPxxY
    ("7.43",  "Schiff K",       "switch"),
    ("7.49",  "NPxxY",          "switch"),
    ("7.50",  "NPxxY",          "switch"),
    ("7.53",  "NPxxY",          "switch"),
]


# ---------------------------------------------------------------------------
# Amino acid property coloring (mapped to thesis palette)
# ---------------------------------------------------------------------------
AA_PROPERTY = {}
for aa in "AVILM":
    AA_PROPERTY[aa] = ("nonpolar",  "#e9ecef",  DARK_GRAY)  # Very light gray
for aa in "GP":
    AA_PROPERTY[aa] = ("nonpolar",  "#dde1e5",  DARK_GRAY)  # Slightly darker
for aa in "FWY":
    AA_PROPERTY[aa] = ("aromatic",  "#d4a03c",  "white")    # Ochre
for aa in "STNQ":
    AA_PROPERTY[aa] = ("polar",     "#a3c4b8",  DARK_GRAY)  # Sage light
for aa in "KRH":
    AA_PROPERTY[aa] = ("positive",  "#3d5a80",  "white")    # Slate
for aa in "DE":
    AA_PROPERTY[aa] = ("negative",  "#c1666b",  "white")    # Terracotta
AA_PROPERTY["C"] = ("special",  "#c4b3c7",  DARK_GRAY)      # Mauve light
AA_PROPERTY["-"] = ("gap",      "white",    LIGHT_GRAY)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_query_grns() -> list[tuple[str, str, str, str, dict]]:
    """Load GRN data for reference sequences.

    Returns list of (accession, subfamily, display_name, species, grn_dict).
    """
    grn = pd.read_csv(GRN_CSV, index_col=0)

    rows = []
    for acc, subfamily, name, species in REFERENCES:
        matches = [i for i in grn.index if str(i).startswith(acc + "|")]
        if not matches:
            print(f"  WARNING: {acc} ({name}) not found in GRN table")
            continue

        eid = matches[0]
        grn_row = grn.loc[eid]

        aa_dict = {}
        for pos, _, _ in GRN_POSITIONS:
            v = str(grn_row.get(pos, "-"))
            aa_dict[pos] = v[0] if v not in ("nan", "-", "") else "-"

        rows.append((acc, subfamily, name, species, aa_dict))

    return rows


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_table(rows: list) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    apply_style()

    n_rows = len(rows)
    n_cols = len(GRN_POSITIONS)

    cell_w = 0.52
    cell_h = 0.50
    gap_between_groups = 0.15  # Extra gap between TM helix groups

    # Compute column x-positions with gaps between helices
    col_x = []
    current_helix = None
    x = 0.0
    for pos, _, _ in GRN_POSITIONS:
        helix = pos.split(".")[0]
        if current_helix is not None and helix != current_helix:
            x += gap_between_groups
        col_x.append(x)
        x += cell_w
        current_helix = helix
    total_w = x

    fig_w = total_w + 3.5
    fig_h = n_rows * cell_h + 2.2

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Draw cells
    for i, (acc, subfamily, name, species, aa_dict) in enumerate(rows):
        y = (n_rows - 1 - i) * cell_h

        for j, (pos, label, cat) in enumerate(GRN_POSITIONS):
            aa = aa_dict.get(pos, "-")
            prop_name, bg, fg = AA_PROPERTY.get(aa, ("unknown", "#f0f0f0", DARK_GRAY))

            cx = col_x[j]
            rect = mpatches.FancyBboxPatch(
                (cx + 0.02, y + 0.02),
                cell_w - 0.04,
                cell_h - 0.04,
                boxstyle="round,pad=0.02",
                facecolor=bg,
                edgecolor="white",
                linewidth=0.8,
            )
            ax.add_patch(rect)

            ax.text(
                cx + cell_w / 2,
                y + cell_h / 2,
                aa,
                ha="center",
                va="center",
                fontsize=10,
                fontfamily="monospace",
                fontweight="bold",
                color=fg,
            )

    # Row labels — subfamily-colored dot + name + species
    for i, (acc, subfamily, name, species, _) in enumerate(rows):
        y = (n_rows - 1 - i) * cell_h + cell_h / 2
        sf_color = SUBFAMILY_COLORS.get(subfamily, GRAY)

        # Colored dot
        ax.plot(-0.45, y, "o", color=sf_color, markersize=6, clip_on=False)

        # Name
        ax.text(
            -0.6, y + 0.02,
            name,
            ha="right", va="center",
            fontsize=8.5, fontweight="bold", color=DARK_GRAY,
        )
        # Species
        ax.text(
            -0.6, y - 0.17,
            species,
            ha="right", va="center",
            fontsize=6.5, fontstyle="italic", color=GRAY,
        )

    # Column headers: GRN position + functional label
    for j, (pos, label, cat) in enumerate(GRN_POSITIONS):
        cx = col_x[j] + cell_w / 2
        y_top = n_rows * cell_h + 0.20

        # GRN position number
        lbl_color = "#bc6c25" if cat == "tuning" else DARK_GRAY
        ax.text(
            cx, y_top, pos,
            ha="center", va="bottom",
            fontsize=7, color=lbl_color,
            rotation=50, rotation_mode="anchor",
        )

        # Functional label above
        if label:
            ax.text(
                cx, y_top + 0.60, label,
                ha="center", va="bottom",
                fontsize=5.5, fontweight="bold", color=lbl_color,
                rotation=50, rotation_mode="anchor",
            )

    # Helix group brackets below
    bracket_y = -0.25
    current_helix = None
    group_start = None
    for j, (pos, _, _) in enumerate(GRN_POSITIONS):
        helix = pos.split(".")[0]
        if helix != current_helix:
            if current_helix is not None:
                # Draw bracket for previous group
                xs = col_x[group_start] + 0.05
                xe = col_x[j - 1] + cell_w - 0.05
                xm = (xs + xe) / 2
                ax.plot([xs, xe], [bracket_y, bracket_y],
                        color=LIGHT_GRAY, lw=1.5, solid_capstyle="round")
                ax.text(xm, bracket_y - 0.18, f"TM{current_helix}",
                        ha="center", va="top", fontsize=6.5,
                        color=GRAY, fontweight="bold")
            current_helix = helix
            group_start = j
    # Last group
    xs = col_x[group_start] + 0.05
    xe = col_x[-1] + cell_w - 0.05
    xm = (xs + xe) / 2
    ax.plot([xs, xe], [bracket_y, bracket_y],
            color=LIGHT_GRAY, lw=1.5, solid_capstyle="round")
    ax.text(xm, bracket_y - 0.18, f"TM{current_helix}",
            ha="center", va="top", fontsize=6.5,
            color=GRAY, fontweight="bold")

    # Legend — AA property colors
    legend_y = -0.85
    legend_items = [
        ("#3d5a80", "white", "Positive (K, R, H)"),
        ("#c1666b", "white", "Negative (D, E)"),
        ("#d4a03c", "white", "Aromatic (F, W, Y)"),
        ("#a3c4b8", DARK_GRAY, "Polar (S, T, N, Q)"),
        ("#e9ecef", DARK_GRAY, "Nonpolar"),
    ]
    lx = 0
    for bg, fg, label in legend_items:
        rect = mpatches.FancyBboxPatch(
            (lx, legend_y), 0.32, 0.25,
            boxstyle="round,pad=0.02", facecolor=bg, edgecolor="white",
        )
        ax.add_patch(rect)
        ax.text(lx + 0.45, legend_y + 0.12, label,
                ha="left", va="center", fontsize=6.5, color=DARK_GRAY)
        lx += 2.0

    ax.set_xlim(-0.6, total_w + 0.3)
    ax.set_ylim(legend_y - 0.4, n_rows * cell_h + 1.9)
    ax.set_aspect("equal")
    ax.axis("off")

    out = FIGURES_DIR / "fig_2.5a_alignment_block.png"
    fig.savefig(str(out), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 70)
    print("FIGURE 2.5a — GRN Microswitch & Spectral Tuning Table")
    print(f"Using {len(REFERENCES)} atlas reference sequences")
    print("=" * 70)

    rows = load_query_grns()
    print(f"Loaded {len(rows)} sequences with GRN annotations\n")

    # Print summary table
    positions = [p for p, _, _ in GRN_POSITIONS]
    header = f"{'':18s} " + " ".join(f"{p:>5s}" for p in positions)
    print(header)
    print("-" * len(header))
    for acc, sf, name, sp, aa in rows:
        vals = " ".join(f"{aa.get(p, '-'):>5s}" for p in positions)
        print(f"{name:18s} {vals}")

    print()
    make_table(rows)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
