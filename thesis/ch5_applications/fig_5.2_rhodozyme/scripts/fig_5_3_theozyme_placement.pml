# Figure 5.3 — Theozyme Placement
# ===========================================================================
# Panel A: Rhodopsin intracellular face with candidate region (GRN)
# Panel B: Winning triplet placed in the cavity
#
# Run interactively: pymol fig_5_3_theozyme_placement.pml
# Run headless:      pymol -cq fig_5_3_theozyme_placement.pml
# ===========================================================================

@rhodozyme_colors.pml
@rhodozyme_config.pml

# =========================================================================
# Load structures
# =========================================================================
python

from pymol import cmd
import builtins

C = builtins._rz_config
PLACEMENT_PDB = C["placement_pdb"]
THEOZYME_HIS = C["theozyme_his"]
THEOZYME_ASP = C["theozyme_asp"]
THEOZYME_SER = C["theozyme_ser"]
TZ = f"{THEOZYME_HIS}+{THEOZYME_ASP}+{THEOZYME_SER}"

cmd.load(PLACEMENT_PDB, "placement")
cmd.remove("placement and not (chain A or chain B)")
cmd.remove("placement and hetatm and not resn RET+LYR+SB1+SB2+SB3+SB4")


print(f"Theozyme: HIS-{THEOZYME_HIS}, ASP-{THEOZYME_ASP}, SER-{THEOZYME_SER}")

# --- Panel A setup ---
# Entire placement is experimental structure -> all gray
cmd.hide("everything", "placement")
cmd.show("cartoon", "placement")
cmd.color("backbone", "placement")


# Retinal
cmd.show("sticks", "placement and resn RET+LYR")
cmd.set("stick_radius", 0.25, "placement and resn RET+LYR")
cmd.color("retinal", "placement and resn RET+LYR")

# Substrate (SB1-SB4 tetrapeptide)
cmd.show("sticks", "placement and resn SB1+SB2+SB3+SB4")
cmd.set("stick_radius", 0.2, "placement and resn SB1+SB2+SB3+SB4")
cmd.color("substrate", "placement and resn SB1+SB2+SB3+SB4")

# Theozyme sidechains + Cα spheres
cmd.show("sticks", f"placement and resi {TZ} and chain A")
cmd.color("theozyme", f"placement and resi {TZ} and chain A")
cmd.show("spheres", f"placement and resi {TZ} and name CA and chain A")
cmd.color("theozyme", f"placement and resi {TZ} and name CA and chain A")
cmd.set("sphere_scale", 0.45, f"placement and resi {TZ} and name CA and chain A")

cmd.set_view((\
    -0.213809878,   -0.452497602,    0.865755320,\
     0.444764674,    0.743972659,    0.498686820,\
    -0.869752586,    0.491682291,    0.042184502,\
    -0.000148175,    0.000032455, -128.739257812,\
   -20.875499725,   -5.523052216,   42.475288391,\
   107.992958069,  149.482116699,   20.000000000))
cmd.scene("panel_A", "store", message="Candidate region")

# --- Panel B setup — same scene, different view ---
cmd.set_view((\
    -0.210006699,   -0.969128430,   -0.129196271,\
     0.973501384,   -0.219502002,    0.064115025,\
    -0.090493947,   -0.112307325,    0.989542186,\
     0.000070643,    0.000207264, -208.219360352,\
   -38.788764954,   -6.597043514,   42.305454254,\
   182.605590820,  233.826721191,   20.000000000))
cmd.scene("panel_B", "store", message="Placed theozyme")

# Start on panel A for interactive viewing
cmd.scene("panel_A", "recall")

print("")
print("Figure 5.3 ready.")
print("  Use scene panel_A / scene panel_B to switch views")
print("  Set your view, then run: @fig_5_3_export.pml")

python end

# =========================================================================
# Export — uncomment or run manually after setting views
# To export: scene panel_A, recall  -> ray 2400, 1800 -> png ...
# Or run:  @fig_5_3_export.pml
# =========================================================================
