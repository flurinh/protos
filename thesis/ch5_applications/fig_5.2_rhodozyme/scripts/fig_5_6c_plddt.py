#!/usr/bin/env python3
"""Figure 5.6c — pLDDT per-residue plot.

Parses Boltz2 CIF for B_iso_or_equiv (pLDDT) and plots per-residue
confidence colored by region: locked (gray), designed (terracotta),
theozyme (green).

Reads Boltz2 CIF path from rhodozyme_config.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from rhodozyme_colors import C, hex, rgb
from rhodozyme_config import CFG

FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# =========================================================================
# Region definitions (must match fig_5_4 / fig_5_6)
# =========================================================================

# Locked TM core residues
TM_CORE = set()
for rng in ["1-55", "71-130", "156-214", "266-302"]:
    a, b = map(int, rng.split("-"))
    TM_CORE.update(range(a, b + 1))

# Designable ICL loop residues
ICL_LOOPS = set()
for rng in ["56-70", "131-155", "215-265", "303-326"]:
    a, b = map(int, rng.split("-"))
    ICL_LOOPS.update(range(a, b + 1))

# Theozyme positions
TZ_POS = set(CFG["theozyme_resi"])

# =========================================================================
# Parse pLDDT from Boltz2 CIF
# =========================================================================

cif_path = CFG["boltz_cif"]
if not cif_path.exists():
    print(f"ERROR: Boltz2 CIF not found: {cif_path}")
    exit(1)

plddt = {}  # resi -> pLDDT
with open(cif_path) as f:
    in_atom_site = False
    for line in f:
        if line.startswith("_atom_site."):
            in_atom_site = True
            continue
        if in_atom_site and line.startswith("ATOM"):
            fields = line.split()
            # Column indices (0-based): 6=label_seq_id, 17=B_iso_or_equiv
            resi = int(fields[6])
            if resi not in plddt:
                plddt[resi] = float(fields[17])
        elif in_atom_site and not line.startswith(("ATOM", "HETATM", "_", "#", "loop")):
            if line.strip() == "" or line.strip() == "#":
                in_atom_site = False

if not plddt:
    print("ERROR: No pLDDT values parsed from CIF")
    exit(1)

residues = sorted(plddt.keys())
values = np.array([plddt[r] for r in residues])

print(f"Parsed {len(residues)} residues from {cif_path.name}")
print(f"pLDDT range: {values.min():.1f} – {values.max():.1f}, mean: {values.mean():.1f}")

# =========================================================================
# Assign colors per residue
# =========================================================================

colors = []
for r in residues:
    if r in TZ_POS:
        colors.append(hex("theozyme"))
    elif r in ICL_LOOPS:
        colors.append(hex("designed"))
    else:
        colors.append(hex("backbone"))

# =========================================================================
# Plot
# =========================================================================

fig, ax = plt.subplots(figsize=(14, 3.5))

ax.bar(residues, values, width=1.0, color=colors, edgecolor="none")

# Confidence thresholds (AlphaFold convention)
ax.axhline(90, color="#444444", lw=0.6, ls="--", alpha=0.5)
ax.axhline(70, color="#444444", lw=0.6, ls="--", alpha=0.5)
ax.text(residues[-1] + 2, 91, "90", fontsize=7, color="#666666", va="bottom")
ax.text(residues[-1] + 2, 71, "70", fontsize=7, color="#666666", va="bottom")

ax.set_xlim(residues[0] - 1, residues[-1] + 1)
ax.set_ylim(0, 100)
ax.set_xlabel("Residue number", fontsize=11, color=hex("label"))
ax.set_ylabel("pLDDT", fontsize=11, color=hex("label"))
ax.set_title(
    f"Boltz2 confidence — design {CFG['design_num']}, seq {CFG['seq_num']} "
    f"(mean pLDDT {values.mean():.1f})",
    fontsize=12, color=hex("label"),
)
ax.tick_params(colors=hex("label"))

# Legend
from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor=hex("backbone"), label="Locked (TM core)"),
    Patch(facecolor=hex("designed"), label="Designed (ICL loops)"),
    Patch(facecolor=hex("theozyme"), label="Theozyme"),
]
ax.legend(handles=legend_elements, loc="lower left", fontsize=9, framealpha=0.9)

plt.tight_layout()
out = FIGURES_DIR / "fig_5_6c_plddt.png"
plt.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved {out}")

# Summary stats by region
locked_vals = [plddt[r] for r in residues if r in TM_CORE]
loop_vals = [plddt[r] for r in residues if r in ICL_LOOPS]
tz_vals = [plddt[r] for r in residues if r in TZ_POS]

print(f"\nRegion means:")
if locked_vals:
    print(f"  TM core (locked):  {np.mean(locked_vals):.1f} ({len(locked_vals)} residues)")
if loop_vals:
    print(f"  ICL loops (designed): {np.mean(loop_vals):.1f} ({len(loop_vals)} residues)")
if tz_vals:
    print(f"  Theozyme:          {np.mean(tz_vals):.1f} ({len(tz_vals)} residues)")
