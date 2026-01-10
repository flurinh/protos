# Render all GPCR pocket views
# Run with: pymol -cq scripts/gpcr_visualization.pml scripts/gpcr_render_all.pml

python
from pymol import cmd
import os

output_dir = "/tmp/gpcr_output"
os.makedirs(output_dir, exist_ok=True)

cmd.bg_color("white")
cmd.set("ray_shadows", 0)
cmd.set("antialias", 2)

def render_view(name):
    path = output_dir + "/" + name + ".png"
    cmd.png(path, width=1200, height=1200, dpi=150, ray=1)
    print("Rendered: " + path)

# Pocket 1: ADRB1
cmd.disable("all")
cmd.enable("2VT4_gpcr")
cmd.enable("2Y02_gpcr")
cmd.enable("2Y04_gpcr")
cmd.enable("2Y00_gpcr")
cmd.show("spheres", "(2VT4_gpcr or 2Y02_gpcr or 2Y04_gpcr or 2Y00_gpcr) and resn HOH")
cmd.view("pocket_view", "recall")
render_view("pocket1_adrb1")

# Pocket 2: ADRB2
cmd.disable("all")
cmd.enable("2RH1_gpcr")
cmd.enable("3NY9_gpcr")
cmd.enable("3SN6_gpcr")
cmd.enable("4LDO_gpcr")
cmd.show("spheres", "(2RH1_gpcr or 3NY9_gpcr or 3SN6_gpcr or 4LDO_gpcr) and resn HOH")
cmd.view("pocket_view", "recall")
render_view("pocket2_adrb2")

# Pocket 3: Active
cmd.disable("all")
cmd.enable("3SN6_gpcr")
cmd.enable("4LDO_gpcr")
cmd.enable("2Y02_gpcr")
cmd.show("spheres", "(3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and resn HOH")
cmd.view("pocket_view", "recall")
render_view("pocket3_active")

# Pocket 4: Inactive
cmd.disable("all")
cmd.enable("2RH1_gpcr")
cmd.enable("3NY9_gpcr")
cmd.enable("2VT4_gpcr")
cmd.show("spheres", "(2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and resn HOH")
cmd.view("pocket_view", "recall")
render_view("pocket4_inactive")

# Pocket 5: Agonist
cmd.disable("all")
cmd.enable("3SN6_gpcr")
cmd.enable("4LDO_gpcr")
cmd.enable("2Y02_gpcr")
cmd.enable("2Y04_gpcr")
cmd.enable("2Y00_gpcr")
cmd.show("spheres", "(3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr or 2Y04_gpcr or 2Y00_gpcr) and resn HOH")
cmd.view("pocket_view", "recall")
render_view("pocket5_agonist")

# Pocket 6: Inverse Agonist
cmd.disable("all")
cmd.enable("2RH1_gpcr")
cmd.enable("3NY9_gpcr")
cmd.show("spheres", "(2RH1_gpcr or 3NY9_gpcr) and resn HOH")
cmd.view("pocket_view", "recall")
render_view("pocket6_inverse_agonist")

# Overall view (no waters)
cmd.enable("all")
cmd.hide("spheres", "resn HOH")
cmd.view("overall_view", "recall")
render_view("overall_structures")

print("")
print("All 7 images saved to: " + output_dir)
python end
