# Figure 5.6 — Boltz2 Evaluation
# ===========================================================================
# Panel A: Predicted structure overlaid on 3PQR
# Panel B: Theozyme geometry comparison (predicted vs placement)
# Panel C: pLDDT plot (separate matplotlib script)
#
# Uses Boltz2 CIF from rhodozyme_config.pml
# Run: pymol -cq fig_5_6_boltz_evaluation.pml
# ===========================================================================

@rhodozyme_colors.pml
@rhodozyme_config.pml

python

from pymol import cmd
import builtins
import os
import numpy as np

C = builtins._rz_config
PLACEMENT_PDB = C["placement_pdb"]
BOLTZ_CIF = C["boltz_cif"]
THEOZYME_HIS = C["theozyme_his"]
THEOZYME_ASP = C["theozyme_asp"]
THEOZYME_SER = C["theozyme_ser"]
TZ = f"{THEOZYME_HIS}+{THEOZYME_ASP}+{THEOZYME_SER}"

# TM core vs designed loops
TM_CORE = "resi 1-55+71-130+156-214+266-302"
ICL_LOOPS = "resi 56-70+131-155+215-265+303-326"

# Extended regions (+4 residues at each boundary for continuous cartoon)
TM_CORE_EXT  = "resi 1-59+67-134+152-218+262-306"
ICL_LOOPS_EXT = "resi 52-74+127-159+211-269+299-326"

# =========================================================================
# Load structures
# =========================================================================

cmd.load(PLACEMENT_PDB, "reference")
cmd.remove("reference and not (chain A or chain B)")
cmd.remove("reference and hetatm and not resn RET+LYR")

if os.path.exists(BOLTZ_CIF):
    cmd.load(BOLTZ_CIF, "predicted")
    cmd.align("predicted and name CA", "reference and name CA")
    print(f"Loaded Boltz2 prediction: {BOLTZ_CIF}")
else:
    print(f"WARNING: Boltz2 CIF not found: {BOLTZ_CIF}")

# Compute overall RMSD
try:
    rmsd = cmd.align("predicted and name CA", "reference and name CA")[0]
    print(f"Overall backbone RMSD: {rmsd:.2f} A")
except:
    pass

# =========================================================================
# Create theozyme copies for Panel B (before any scene setup)
# =========================================================================
cmd.create("ref_tz", f"reference and resi {TZ} and chain A")
cmd.create("pred_tz", f"predicted and resi {TZ}")

# Align reference theozyme ONTO predicted (validation) theozyme
try:
    tz_rmsd = cmd.align("ref_tz and name CA", "pred_tz and name CA")[0]
    print(f"Theozyme Ca RMSD: {tz_rmsd:.2f} A")
except:
    pass

# Style theozyme copies
cmd.color("backbone", "ref_tz")
cmd.color("backbone_dark", "ref_tz and name CA")
cmd.color("theozyme", "pred_tz")

# Hide copies for now
cmd.disable("ref_tz")
cmd.disable("pred_tz")

# =========================================================================
# Panel A — Structural overlay
# =========================================================================

cmd.hide("everything")

# Reference: TM core extended (gray)
cmd.show("cartoon", f"reference and ({TM_CORE_EXT})")
cmd.color("backbone", f"reference and ({TM_CORE_EXT})")

# Predicted: designed loops extended (terracotta)
cmd.show("cartoon", f"predicted and ({ICL_LOOPS_EXT})")
cmd.color("designed", f"predicted and ({ICL_LOOPS_EXT})")

# Theozyme residues in predicted structure
cmd.show("sticks", f"predicted and resi {TZ}")
cmd.color("theozyme", f"predicted and resi {TZ}")
cmd.show("spheres", f"predicted and resi {TZ} and name CA")
cmd.color("theozyme", f"predicted and resi {TZ} and name CA")
cmd.set("sphere_scale", 0.35, f"predicted and resi {TZ} and name CA")

# Retinal (reference)
cmd.show("sticks", "reference and resn RET")
cmd.set("stick_radius", 0.25, "reference and resn RET")
cmd.color("retinal", "reference and resn RET")

# Schiff base Lys-296 sidechain (reference)
cmd.show("sticks", "reference and resi 296 and not name C+N+O")
cmd.color("retinal", "reference and resi 296 and not name C+N+O")

cmd.set_view((\
    -0.224436864,   -0.957964301,    0.178679466,\
     0.908177614,   -0.272100776,   -0.318075329,\
     0.353325307,    0.090886854,    0.931067526,\
    -0.000118424,   -0.000183463, -208.517562866,\
   -40.395961761,   -7.634953499,   38.179256439,\
   160.063751221,  256.901672363,   20.000000000))
cmd.scene("panel_A", "store", message="Boltz2 overlay")

# Export Panel A
cmd.ray(2400, 1800)
cmd.png("../figures/fig_5_6a_boltz_overlay.png", dpi=300)

# =========================================================================
# Panel B — Theozyme geometry comparison (separate objects)
# =========================================================================

cmd.hide("everything")
cmd.disable("reference")
cmd.disable("predicted")
cmd.enable("ref_tz")
cmd.enable("pred_tz")

# ref_tz: gray sticks + dark spheres
cmd.show("sticks", "ref_tz")
cmd.show("spheres", "ref_tz and name CA")
cmd.set("sphere_scale", 0.4, "ref_tz and name CA")

# pred_tz: green sticks + spheres
cmd.show("sticks", "pred_tz")
cmd.show("spheres", "pred_tz and name CA")
cmd.set("sphere_scale", 0.4, "pred_tz and name CA")

# Catalytic triad interaction distances (Ser-His, His-Asp) in both structures
sidechain_atoms = {
    "SER": ["OG"],
    "HIS": ["ND1", "NE2", "CD2", "CE1", "CG"],
    "ASP": ["OD1", "OD2", "CG"],
}
interactions = [
    ("Ser-His", THEOZYME_SER, "SER", THEOZYME_HIS, "HIS"),
    ("His-Asp", THEOZYME_HIS, "HIS", THEOZYME_ASP, "ASP"),
]

for obj, obj_label, dash_color in [("ref_tz", "ref", "backbone_dark"), ("pred_tz", "pred", "theozyme")]:
    for iname, resi1, resn1, resi2, resn2 in interactions:
        best_dist = 999.0
        best_a1 = best_a2 = None
        for a1 in sidechain_atoms[resn1]:
            for a2 in sidechain_atoms[resn2]:
                sel1 = f"{obj} and resi {resi1} and name {a1}"
                sel2 = f"{obj} and resi {resi2} and name {a2}"
                c1 = cmd.get_coords(sel1)
                c2 = cmd.get_coords(sel2)
                if c1 is not None and c2 is not None and len(c1) > 0 and len(c2) > 0:
                    d = np.linalg.norm(c1[0] - c2[0])
                    if d < best_dist:
                        best_dist = d
                        best_a1 = sel1
                        best_a2 = sel2
        if best_a1:
            dname = f"dist_{obj_label}_{iname.replace('-','_')}"
            cmd.distance(dname, best_a1, best_a2)
            cmd.set("dash_color", dash_color, dname)
            cmd.set("dash_gap", 0.3, dname)
            cmd.set("dash_length", 0.2, dname)
            cmd.set("dash_width", 2.0, dname)
            cmd.set("label_color", "black", dname)
            cmd.set("label_size", 16, dname)
            print(f"  {obj_label} {iname}: {best_dist:.2f} A")

cmd.set_view((\
     0.151739627,   -0.380330384,    0.912313700,\
     0.360970914,   -0.837928653,   -0.409356683,\
     0.920141935,    0.391435802,    0.010139735,\
    -0.000512086,   -0.000161117,  -29.542037964,\
   -13.538523674,  -12.698078156,   41.678886414,\
    14.325818062,   44.754657745,   20.000000000))
cmd.scene("panel_B", "store", message="Theozyme comparison")

# Export Panel B
cmd.ray(2400, 1800)
cmd.png("../figures/fig_5_6b_theozyme_comparison.png", dpi=300)

cmd.save("../figures/fig_5_6_boltz_evaluation.pse")

print("")
print(f"Figure 5.6 complete. Design={C['design_num']}, Seq={C['seq_num']}")
print("Panel C (pLDDT plot): run fig_5_6c_plddt.py separately")

python end
