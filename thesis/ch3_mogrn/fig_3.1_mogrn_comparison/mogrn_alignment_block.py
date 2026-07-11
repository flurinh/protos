#!/usr/bin/env python3
"""Figure 3.1a: MO-GRN Alignment Block for Microbial Rhodopsins.

Produces fig_3.1a_mogrn_alignment_block.png — a table of amino acid identities
at 20 key MO-GRN positions across 8 representative microbial rhodopsins spanning
the major functional classes.

Rows    = representative microbial rhodopsins (one per functional class)
Columns = functionally important MO-GRN positions (tuning, functional, pocket, channel)
Cells   = amino acid letter, background colored by physicochemical property

Data source:
  data/grn/reference/mo_ref.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
THESIS_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(THESIS_DIR / "shared"))
from thesis_style import DARK_GRAY, GRAY, LIGHT_GRAY, apply_style

PROJECT_ROOT = THESIS_DIR.parent
GRN_CSV = PROJECT_ROOT / "data" / "grn" / "reference" / "mo_ref.csv"
FIGURES_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# 8 representative microbial rhodopsins
# (entity_id, functional_class, display_name, species)
# ---------------------------------------------------------------------------
REFERENCES = [
    ("1C3W",          "H\u207a pump",        "Bacteriorhodopsin", "H. salinarum"),
    ("3A7K",          "Cl\u207b pump",       "Halorhodopsin",     "N. pharaonis"),
    ("3X3C",          "Na\u207a pump",       "KR2",               "K. eikastus"),
    ("3UG9",          "Cation channel",      "C1C2",              "C. reinhardtii"),
    ("6CSN",          "Anion channel",       "GtACR1",            "G. theta"),
    ("1JGJ",          "Sensory",             "NpSRII",            "N. pharaonis"),
    ("NeoR_model_0",  "Enzyme-fused",        "NeoR",              "R. globosum"),
    ("6IS6",          "Heliorhodopsin",      "TaHeR",             "Thermoplasm."),
]


# ---------------------------------------------------------------------------
# 20 key MO-GRN positions: (position, label, category)
# category: "tuning" = spectral tuning, "functional" = microswitch/conserved,
#           "pocket" = binding pocket, "channel" = channel-specific
# ---------------------------------------------------------------------------
MO_GRN_POSITIONS = [
    # TM2
    ("2.46", "E1 (ChR)",       "channel"),
    ("2.47", "E2 (ChR)",       "channel"),
    # TM3
    ("3.42", "Conserved R",    "functional"),
    ("3.45", "Counterion",     "functional"),
    ("3.46", "Pocket",         "pocket"),
    ("3.49", "TM3 motif",      "functional"),
    ("3.50", "DC gate",        "functional"),
    ("3.53", "L/Q switch",     "tuning"),
    ("3.56", "TM3 motif",      "functional"),
    # TM4
    ("4.47", "DC gate",        "functional"),
    ("4.51", "N/LI switch",    "tuning"),
    ("4.54", "Non-G rule",     "tuning"),
    # TM5
    ("5.44", "Retinal uptake", "functional"),
    ("5.47", "Carotenoid",     "functional"),
    # TM6
    ("6.50", "Pocket",         "pocket"),
    ("6.53", "Pocket",         "pocket"),
    ("6.54", "G/P switch",     "tuning"),
    # TM7
    ("7.46", "Counterion",     "functional"),
    ("7.49", "A/TS switch",    "tuning"),
    ("7.50", "Schiff base K",  "functional"),
]


# ---------------------------------------------------------------------------
# Functional class colors (muted tones from the thesis palette)
# ---------------------------------------------------------------------------
FUNCTIONAL_COLORS = {
    "H\u207a pump":        "#3d5a80",  # Slate
    "Cl\u207b pump":       "#457b6b",  # Sage
    "Na\u207a pump":       "#d4a03c",  # Ochre
    "Cation channel":      "#6a9bc3",  # Light slate
    "Anion channel":       "#8ab0a1",  # Light sage
    "Sensory":             "#c1666b",  # Terracotta
    "Enzyme-fused":        "#7a5980",  # Mauve
    "Heliorhodopsin":      "#bc6c25",  # Rust
}


# ---------------------------------------------------------------------------
# Amino acid property coloring (same scheme as Type II alignment block)
# ---------------------------------------------------------------------------
AA_PROPERTY = {}
for _aa in "AVILM":
    AA_PROPERTY[_aa] = ("nonpolar", "#e9ecef", DARK_GRAY)   # Very light gray
for _aa in "GP":
    AA_PROPERTY[_aa] = ("nonpolar", "#dde1e5", DARK_GRAY)   # Slightly darker
for _aa in "FWY":
    AA_PROPERTY[_aa] = ("aromatic", "#d4a03c", "white")      # Ochre
for _aa in "STNQ":
    AA_PROPERTY[_aa] = ("polar",    "#a3c4b8", DARK_GRAY)    # Sage light
for _aa in "KRH":
    AA_PROPERTY[_aa] = ("positive", "#3d5a80", "white")      # Slate
for _aa in "DE":
    AA_PROPERTY[_aa] = ("negative", "#c1666b", "white")      # Terracotta
AA_PROPERTY["C"] = ("special", "#c4b3c7", DARK_GRAY)         # Mauve light
AA_PROPERTY["-"] = ("gap",     "white",   LIGHT_GRAY)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def _col_name(pos: str, columns: list[str]) -> str:
    """Find matching CSV column, handling trailing-zero removal (3.50 → 3.5)."""
    if pos in columns:
        return pos
    normalized = str(float(pos))
    if normalized in columns:
        return normalized
    return pos


def load_mo_grns() -> list[tuple[str, str, str, str, dict]]:
    """Load MO-GRN data for the 8 reference microbial rhodopsins.

    Returns list of (entity_id, func_class, display_name, species, aa_dict).
    """
    df = pd.read_csv(GRN_CSV, index_col=0)
    df.columns = df.columns.astype(str)
    cols = list(df.columns)

    rows = []
    for entity_id, func_class, display_name, species in REFERENCES:
        if entity_id not in df.index:
            print(f"  WARNING: {entity_id} ({display_name}) not found in GRN table")
            continue

        grn_row = df.loc[entity_id]

        aa_dict = {}
        for pos, _, _ in MO_GRN_POSITIONS:
            col = _col_name(pos, cols)
            v = str(grn_row.get(col, "-"))
            if v in ("nan", "-", "", "None"):
                aa_dict[pos] = "-"
            else:
                aa_dict[pos] = v[0]

        rows.append((entity_id, func_class, display_name, species, aa_dict))

    return rows


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_table(rows: list) -> None:
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    apply_style()

    n_rows = len(rows)
    cell_w = 0.52
    cell_h = 0.50
    gap_between_groups = 0.15  # Extra gap between TM helix groups

    # Compute column x-positions with gaps between helices
    col_x = []
    current_helix = None
    x = 0.0
    for pos, _, _ in MO_GRN_POSITIONS:
        helix = pos.split(".")[0]
        if current_helix is not None and helix != current_helix:
            x += gap_between_groups
        col_x.append(x)
        x += cell_w
        current_helix = helix
    total_w = x

    fig_w = total_w + 4.0
    fig_h = n_rows * cell_h + 2.2

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Draw cells
    for i, (eid, func_class, name, species, aa_dict) in enumerate(rows):
        y = (n_rows - 1 - i) * cell_h

        for j, (pos, label, cat) in enumerate(MO_GRN_POSITIONS):
            aa = aa_dict.get(pos, "-")
            _, bg, fg = AA_PROPERTY.get(aa, ("unknown", "#f0f0f0", DARK_GRAY))

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

    # Row labels — name + species
    for i, (eid, func_class, name, species, _) in enumerate(rows):
        y = (n_rows - 1 - i) * cell_h + cell_h / 2

        # Name
        ax.text(
            -0.25, y + 0.02,
            name,
            ha="right", va="center",
            fontsize=8.5, fontweight="bold", color=DARK_GRAY,
        )
        # Species
        ax.text(
            -0.25, y - 0.17,
            species,
            ha="right", va="center",
            fontsize=6.5, fontstyle="italic", color=GRAY,
        )

    # Column headers: GRN position + functional label
    for j, (pos, label, cat) in enumerate(MO_GRN_POSITIONS):
        cx = col_x[j] + cell_w / 2
        y_top = n_rows * cell_h + 0.20

        # Color: tuning positions in ochre, others in dark gray
        lbl_color = "#bc6c25" if cat == "tuning" else DARK_GRAY

        # GRN position number
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
    for j, (pos, _, _) in enumerate(MO_GRN_POSITIONS):
        helix = pos.split(".")[0]
        if helix != current_helix:
            if current_helix is not None:
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

    out = FIGURES_DIR / "fig_3.1a_mogrn_alignment_block.png"
    fig.savefig(str(out), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 70)
    print("FIGURE 3.1a \u2014 MO-GRN Alignment Block (Microbial Rhodopsins)")
    print(f"Using {len(REFERENCES)} representative microbial rhodopsins")
    print("=" * 70)

    rows = load_mo_grns()
    print(f"Loaded {len(rows)} sequences with MO-GRN annotations\n")

    # Print summary table
    positions = [p for p, _, _ in MO_GRN_POSITIONS]
    header = f"{'':20s} " + " ".join(f"{p:>5s}" for p in positions)
    print(header)
    print("-" * len(header))
    for eid, fc, name, sp, aa in rows:
        vals = " ".join(f"{aa.get(p, '-'):>5s}" for p in positions)
        print(f"{name:20s} {vals}")

    print()
    make_table(rows)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
