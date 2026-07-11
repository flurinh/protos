# Figure 5.4 — RFdiffusion Mask and Output
# ===========================================================================
# Panel A: The mask (locked gray vs designable terracotta)
# Panel B: Three designs superimposed on locked scaffold
#
# Uses design PDB from rhodozyme_config.pml
# Run: pymol -cq fig_5_4_rfdiffusion_mask.pml
# ===========================================================================

@rhodozyme_colors.pml
@rhodozyme_config.pml

python

from pymol import cmd
import builtins
import os

C = builtins._rz_config
PLACEMENT_PDB = C["placement_pdb"]
DESIGN_PDB = C["design_pdb"]
RUN_DIR = C["run_dir"]
PLACEMENT_NUM = C["placement_num"]
DESIGN_NUM = C["design_num"]
THEOZYME_RESI = C["theozyme_resi"]

# Load placement as reference scaffold
cmd.load(PLACEMENT_PDB, "scaffold")
cmd.remove("scaffold and not (chain A or chain B)")
cmd.remove("scaffold and hetatm and not resn RET+LYR+SB1+SB2+SB3+SB4")

# Load the selected design
cmd.load(DESIGN_PDB, "design_0")
cmd.align("design_0 and name CA", "scaffold and name CA")

# Load two more designs for diversity panel (pick next available)
extra_designs = []
for i in [2, 5, 10, 15, 20]:
    if str(i) == DESIGN_NUM:
        continue
    pdb = f"{RUN_DIR}/rhodozyme_{PLACEMENT_NUM}_{i}-atomized-bb-False.pdb"
    if os.path.exists(pdb):
        name = f"design_{i}"
        cmd.load(pdb, name)
        cmd.align(f"{name} and name CA", "scaffold and name CA")
        extra_designs.append(name)
        if len(extra_designs) == 2:
            break

# =========================================================================
# Define locked vs designable regions
# =========================================================================
# Locked: TM helices (core) + theozyme positions
# Designable: ICL loops between locked regions
# This is approximate — the actual mask from RFdiffusion may differ slightly

# TM core residues (locked)
tm_core = "resi 1-55+71-130+156-214+266-302"
# ICL loops (designable)
icl_loops = "resi 56-70+131-155+215-265+303-326"
# Theozyme
tz = f"resi {C['theozyme_his']}+{C['theozyme_asp']}+{C['theozyme_ser']}"

# =========================================================================
# Panel A — The mask
# =========================================================================

cmd.hide("everything")
cmd.show("cartoon", "scaffold")

# Locked regions: gray
cmd.color("backbone", f"scaffold and ({tm_core})")
# Designable regions: terracotta
cmd.color("designed", f"scaffold and ({icl_loops})")
# Theozyme: green spheres
cmd.show("spheres", f"scaffold and ({tz}) and name CA")
cmd.color("theozyme", f"scaffold and ({tz}) and name CA")
cmd.set("sphere_scale", 0.45, f"scaffold and ({tz}) and name CA")
cmd.show("sticks", f"scaffold and ({tz})")
cmd.color("theozyme", f"scaffold and ({tz})")

# Retinal
cmd.show("sticks", "scaffold and resn RET+LYR")
cmd.set("stick_radius", 0.25, "scaffold and resn RET+LYR")
cmd.color("retinal", "scaffold and resn RET+LYR")
cmd.set("specular", 0.8, "scaffold and resn RET+LYR")

# Substrate (SB1-SB4 tetrapeptide)
cmd.show("sticks", "scaffold and resn SB1+SB2+SB3+SB4")
cmd.set("stick_radius", 0.2, "scaffold and resn SB1+SB2+SB3+SB4")
cmd.color("substrate", "scaffold and resn SB1+SB2+SB3+SB4")

cmd.set_view((\
    -0.210006699,   -0.969128430,   -0.129196271,\
     0.973501384,   -0.219502002,    0.064115025,\
    -0.090493947,   -0.112307325,    0.989542186,\
     0.000070643,    0.000207264, -208.219360352,\
   -38.788764954,   -6.597043514,   42.305454254,\
   182.605590820,  233.826721191,   20.000000000))
cmd.scene("panel_A", "store", message="RFdiffusion mask")

# =========================================================================
# Panel B — Three designs superimposed
# =========================================================================

# Reference scaffold: full gray
cmd.hide("everything")
cmd.show("cartoon", "scaffold")
cmd.color("backbone", "scaffold")

# Design 0: only show designed loop regions (locked core is from scaffold)
cmd.show("cartoon", f"design_0 and ({icl_loops})")
cmd.color("designed", f"design_0 and ({icl_loops})")

# Extra designs: only loop regions, lighter terracotta
for dname in extra_designs:
    cmd.show("cartoon", f"{dname} and ({icl_loops})")
    cmd.color("designed_light", f"{dname} and ({icl_loops})")

# Theozyme in all designs
for dname in ["design_0"] + extra_designs:
    cmd.show("spheres", f"{dname} and ({tz}) and name CA")
    cmd.color("theozyme", f"{dname} and ({tz}) and name CA")
    cmd.set("sphere_scale", 0.45, f"{dname} and ({tz}) and name CA")

# Retinal + substrate from scaffold
cmd.show("sticks", "scaffold and resn RET+LYR")
cmd.color("retinal", "scaffold and resn RET+LYR")
cmd.show("sticks", "scaffold and resn SB1+SB2+SB3+SB4")
cmd.color("substrate", "scaffold and resn SB1+SB2+SB3+SB4")

cmd.set_view((\
    -0.210006699,   -0.969128430,   -0.129196271,\
     0.973501384,   -0.219502002,    0.064115025,\
    -0.090493947,   -0.112307325,    0.989542186,\
     0.000070643,    0.000207264, -208.219360352,\
   -38.788764954,   -6.597043514,   42.305454254,\
   182.605590820,  233.826721191,   20.000000000))
cmd.scene("panel_B", "store", message="Three backbone designs")

# =========================================================================
# Export
# =========================================================================

cmd.scene("panel_A", "recall")
cmd.ray(2400, 1800)
cmd.png("../figures/fig_5_4a_mask.png", dpi=300)

cmd.scene("panel_B", "recall")
cmd.ray(2400, 1800)
cmd.png("../figures/fig_5_4b_designs.png", dpi=300)

cmd.save("../figures/fig_5_4_rfdiffusion_mask.pse")

print("")
print(f"Figure 5.4 complete. Design={DESIGN_NUM}, extras={[d for d in extra_designs]}")

python end
