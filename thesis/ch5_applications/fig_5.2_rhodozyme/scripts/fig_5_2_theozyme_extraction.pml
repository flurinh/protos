# Figure 5.2 — Theozyme Extraction
# ===========================================================================
# Single panel: 2AGE trypsin with catalytic triad + substrate highlighted,
# and the geometric abstraction overlaid (Cα spheres, Cα→Cβ vectors,
# pairwise Cα–Cα distances).
#
# Trypsin backbone in gray (experimental). Only theozyme + substrate colored.
# Style matches fig 4.2.
#
# Run: pymol -cq fig_5_2_theozyme_extraction.pml
# ===========================================================================

# Load color scheme + render settings
@rhodozyme_colors.pml

# =============================================================================
# Fetch and prepare
# =============================================================================
fetch 2AGE, name=trypsin, type=pdb
# 2AGE: chain X = trypsin protein, chain A = succinyl-AAPR substrate
remove trypsin and not (chain X or chain A)
remove resn HOH

# =============================================================================
# Backbone — gray, heavily transparent (context only)
# =============================================================================
hide everything, trypsin
show cartoon, trypsin
color backbone, trypsin
set cartoon_transparency, 0.7, trypsin

# =============================================================================
# Catalytic triad — theozyme green sticks
# =============================================================================
select triad, trypsin and resi 57+102+195 and chain X
show sticks, triad
color theozyme, triad

# =============================================================================
# Substrate — ochre, slightly transparent (context, not the focus)
# =============================================================================
# 2AGE: chain A = succinyl-AAPR substrate peptide
show sticks, trypsin and chain A
color substrate, trypsin and chain A
set stick_transparency, 0.3, trypsin and chain A

# =============================================================================
# Cα spheres — geometric anchor points
# =============================================================================
select triad_ca, trypsin and resi 57+102+195 and name CA and chain X
show spheres, triad_ca
color geo, triad_ca
set sphere_scale, 0.45, triad_ca

# =============================================================================
# Cα→Cβ direction vectors — CGO arrows
# =============================================================================
python

from pymol import cmd
from pymol.cgo import *
import numpy as np

def get_coord(selection):
    model = cmd.get_model(selection)
    if model.atom:
        a = model.atom[0]
        return np.array([a.coord[0], a.coord[1], a.coord[2]])
    return None

residues = [
    ("Ser195", "resi 195"),
    ("His57",  "resi 57"),
    ("Asp102", "resi 102"),
]

shaft_r = 0.10
head_r = 0.22
head_len = 0.5
vec_scale = 3.0
c = [0.204, 0.227, 0.251]  # geo color

cgo = []
for name, rsel in residues:
    ca = get_coord(f"trypsin and {rsel} and name CA and chain X")
    cb = get_coord(f"trypsin and {rsel} and name CB and chain X")
    if ca is not None and cb is not None:
        d = cb - ca
        d = d / np.linalg.norm(d)
        tip = ca + d * vec_scale
        cone_tip = tip + d * head_len
        # Shaft
        cgo.extend([CYLINDER, ca[0], ca[1], ca[2], tip[0], tip[1], tip[2],
                     shaft_r, c[0], c[1], c[2], c[0], c[1], c[2]])
        # Head
        cgo.extend([CONE, tip[0], tip[1], tip[2], cone_tip[0], cone_tip[1], cone_tip[2],
                     head_r, 0.0, c[0], c[1], c[2], c[0], c[1], c[2], 1.0, 0.0])

if cgo:
    cmd.load_cgo(cgo, "cb_vectors")

python end

# =============================================================================
# Pairwise Cα–Cα distances — dashed lines with Å labels
# =============================================================================
distance d_ser_his, trypsin and resi 195 and name CA and chain X, trypsin and resi 57 and name CA and chain X
distance d_his_asp, trypsin and resi 57 and name CA and chain X, trypsin and resi 102 and name CA and chain X
distance d_ser_asp, trypsin and resi 195 and name CA and chain X, trypsin and resi 102 and name CA and chain X

set dash_color, red, d_ser_his
set dash_color, red, d_his_asp
set dash_color, red, d_ser_asp
set label_color, black
set dash_width, 2.5
set dash_gap, 0.4
set label_size, 16


# =============================================================================
# View and export
# =============================================================================
zoom triad, buffer=10
orient triad

scene theozyme, store, message=Theozyme extraction

# Ray and save
ray 2400, 1800
png ../figures/fig_5_2_theozyme_extraction.png, dpi=300
save ../figures/fig_5_2_theozyme_extraction.pse

print("")
print("Figure 5.2 complete.")
print("  fig_5_2_theozyme_extraction.png")
print("  fig_5_2_theozyme_extraction.pse")
print("")
print("Distance values (copy to application_draft.md Figure 5.2 caption):")
print("  Ser195-His57:  see d_ser_his label")
print("  His57-Asp102:  see d_his_asp label")
print("  Ser195-Asp102: see d_ser_asp label")
print("")
print("Color key:")
print("  Gray (backbone) = experimental crystal structure")
print("  Green (theozyme) = catalytic triad sidechains")
print("  Ochre (substrate) = succinyl-AAPR")
print("  Dark gray (geo) = Ca spheres + Cb vectors")
print("  Mid gray (hbond) = distance dashes")
