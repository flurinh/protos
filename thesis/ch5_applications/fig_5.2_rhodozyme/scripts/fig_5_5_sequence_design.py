#!/usr/bin/env python3
"""Figure 5.5 — Sequence Design.

Panel A: Alignment of designed sequence vs wild-type rhodopsin (3PQR).
Panel B: Sequence logo across 8 sampled sequences for one backbone.

Reads LigandMPNN FASTA from rhodozyme_config.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rhodozyme_colors import C, hex, rgb
from rhodozyme_config import CFG

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# =========================================================================
# Parse FASTA
# =========================================================================

def parse_fasta(fasta_path):
    """Parse multi-sequence FASTA. Return list of (header, sequence)."""
    seqs = []
    header, seq = None, []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        seqs.append((header, "".join(seq)))
    return seqs


# =========================================================================
# Load sequences
# =========================================================================

fa_path = CFG["design_fa"]
if not fa_path.exists():
    print(f"ERROR: FASTA not found: {fa_path}")
    exit(1)

seqs = parse_fasta(fa_path)
wt_seq = seqs[0][1]  # First entry is the template/WT
designed_seqs = [(h, s) for h, s in seqs[1:]]  # Rest are designed

print(f"Loaded {len(designed_seqs)} designed sequences from {fa_path.name}")
print(f"WT length: {len(wt_seq)}, Design length: {len(designed_seqs[0][1])}")

# Use first designed sequence for alignment panel
design_seq = designed_seqs[0][1]

# Theozyme positions (0-indexed)
tz_positions = [CFG["theozyme_his"] - 1, CFG["theozyme_asp"] - 1, CFG["theozyme_ser"] - 1]

# =========================================================================
# Panel A — Sequence alignment strip
# =========================================================================

fig, ax = plt.subplots(figsize=(18, 2.5))

n = min(len(wt_seq), len(design_seq))
colors = []
for i in range(n):
    if i in tz_positions:
        colors.append(hex("theozyme"))
    elif wt_seq[i] == design_seq[i]:
        colors.append(hex("seq_identity"))
    else:
        colors.append(hex("seq_mutation"))

# Draw colored bars
for i in range(n):
    ax.bar(i, 1, width=1.0, color=colors[i], edgecolor="none")

# Count mutations
n_mut = sum(1 for i in range(n) if wt_seq[i] != design_seq[i] and i not in tz_positions)
n_ident = sum(1 for i in range(n) if wt_seq[i] == design_seq[i])

ax.set_xlim(0, n)
ax.set_ylim(0, 1)
ax.set_yticks([])
ax.set_xlabel("Residue position", fontsize=11, color=hex("label"))
ax.set_title(
    f"Designed vs WT — {n_ident}/{n} identical ({100*n_ident/n:.0f}%), "
    f"{n_mut} mutations, {len(tz_positions)} theozyme",
    fontsize=12, color=hex("label")
)
ax.tick_params(colors=hex("label"))

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=hex("seq_identity"), label="Identical"),
    Patch(facecolor=hex("seq_mutation"), label="Mutated"),
    Patch(facecolor=hex("theozyme"), label="Theozyme"),
]
ax.legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.9)

plt.tight_layout()
out_a = FIGURES_DIR / "fig_5_5a_sequence_alignment.png"
plt.savefig(out_a, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved {out_a}")

# =========================================================================
# Panel B — Sequence logo (conservation across 8 sequences)
# =========================================================================

# Compute per-position conservation
aa_order = "ACDEFGHIKLMNPQRSTVWY"
n_seqs = len(designed_seqs)
n_pos = min(len(s) for _, s in designed_seqs)

# Only show the designed region (positions that differ from WT)
designed_positions = [i for i in range(n_pos) if any(s[i] != wt_seq[i] for _, s in designed_seqs)]

if designed_positions:
    # Entropy-based information content
    conservation = []
    for i in designed_positions:
        counts = {}
        for _, s in designed_seqs:
            aa = s[i]
            counts[aa] = counts.get(aa, 0) + 1
        # Shannon entropy
        H = 0
        for c in counts.values():
            p = c / n_seqs
            if p > 0:
                H -= p * np.log2(p)
        max_H = np.log2(20)
        info = max_H - H  # information content
        conservation.append(info)

    fig, ax = plt.subplots(figsize=(18, 3))

    bar_colors = []
    for i, pos in enumerate(designed_positions):
        if pos in tz_positions:
            bar_colors.append(hex("theozyme"))
        elif conservation[i] > 3.0:
            bar_colors.append(hex("candidate_region"))  # sage = conserved
        else:
            bar_colors.append(hex("designed"))  # terracotta = variable

    ax.bar(range(len(designed_positions)), conservation, color=bar_colors, edgecolor="none")
    ax.set_xlabel("Designed region position", fontsize=11, color=hex("label"))
    ax.set_ylabel("Information content (bits)", fontsize=11, color=hex("label"))
    ax.set_title(f"Conservation across {n_seqs} sequences for design {CFG['design_num']}",
                 fontsize=12, color=hex("label"))
    ax.tick_params(colors=hex("label"))
    ax.set_xlim(-0.5, len(designed_positions) - 0.5)

    legend_elements = [
        Patch(facecolor=hex("candidate_region"), label="Conserved"),
        Patch(facecolor=hex("designed"), label="Variable"),
        Patch(facecolor=hex("theozyme"), label="Theozyme (fixed)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=9, framealpha=0.9)

    plt.tight_layout()
    out_b = FIGURES_DIR / "fig_5_5b_sequence_logo.png"
    plt.savefig(out_b, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"Saved {out_b}")
else:
    print("No variable positions found — skipping Panel B")

print(f"\nConfig: placement={CFG['placement_num']}, design={CFG['design_num']}")
