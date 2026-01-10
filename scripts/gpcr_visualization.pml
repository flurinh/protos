# GPCR Binding Pocket Comparison - Main Script
# =============================================
#
# This script loads, extracts, and aligns 8 adrenergic receptor structures.
# Run pocket-specific scripts after this to visualize different subsets.
#
# Usage: pymol scripts/gpcr_visualization.pml
#
# Then run one of:
#   @scripts/gpcr_visualization_pocket1.pml  (ADRB1)
#   @scripts/gpcr_visualization_pocket2.pml  (ADRB2)
#   @scripts/gpcr_visualization_pocket3.pml  (Active)
#   @scripts/gpcr_visualization_pocket4.pml  (Inactive)
#   @scripts/gpcr_visualization_pocket5.pml  (Agonist)
#   @scripts/gpcr_visualization_pocket6.pml  (Inverse Agonist)

# Clear any existing objects
reinitialize

# Settings for better visualization
set fetch_path, /tmp/pdb_cache
set cif_use_auth, on
bg_color white
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set depth_cue, 0

# =============================================================================
# STEP 1: Fetch structures from PDB
# =============================================================================
print("Fetching structures from PDB...")

fetch 2RH1, type=cif, async=0
fetch 3NY9, type=cif, async=0
fetch 3SN6, type=cif, async=0
fetch 4LDO, type=cif, async=0
fetch 2VT4, type=cif, async=0
fetch 2Y02, type=cif, async=0
fetch 2Y04, type=cif, async=0
fetch 2Y00, type=cif, async=0

print("All structures loaded.")

# =============================================================================
# STEP 2: Extract GPCR chains with waters and ligands
# =============================================================================
print("Extracting GPCR chains...")

# Residue ranges determined from GRN annotations (7TM receptor only, no T4L/fusions)

# 2RH1: Chain A, inverse agonist (carazolol=CAU), inactive, ADRB2
select prot_2RH1, 2RH1 and chain A and polymer.protein and resi 29-341
select wat_2RH1, 2RH1 and resn HOH within 5 of prot_2RH1
select lig_2RH1, 2RH1 and resn CAU
create 2RH1_gpcr, prot_2RH1 or wat_2RH1 or lig_2RH1

# 3NY9: Chain A, inverse agonist (ICI-118551=JSZ), inactive, ADRB2
select prot_3NY9, 3NY9 and chain A and polymer.protein and resi 32-341
select wat_3NY9, 3NY9 and resn HOH within 5 of prot_3NY9
select lig_3NY9, 3NY9 and resn JSZ
create 3NY9_gpcr, prot_3NY9 or wat_3NY9 or lig_3NY9

# 3SN6: Chain R, full agonist (BI-167107=P0G), active, ADRB2
select prot_3SN6, 3SN6 and chain R and polymer.protein and resi 30-341
select wat_3SN6, 3SN6 and resn HOH within 5 of prot_3SN6
select lig_3SN6, 3SN6 and resn P0G
create 3SN6_gpcr, prot_3SN6 or wat_3SN6 or lig_3SN6

# 4LDO: Chain A, full agonist (adrenaline=ALE), active, ADRB2
# Note: residue numbering offset by ~1000, N-term starts at 1029
select prot_4LDO, 4LDO and chain A and polymer.protein and resi 1029-1341
select wat_4LDO, 4LDO and resn HOH within 5 of prot_4LDO
select lig_4LDO, 4LDO and resn ALE
create 4LDO_gpcr, prot_4LDO or wat_4LDO or lig_4LDO

# 2VT4: Chain A, antagonist (cyanopindolol=P32), inactive, ADRB1
select prot_2VT4, 2VT4 and chain A and polymer.protein and resi 40-358
select wat_2VT4, 2VT4 and resn HOH within 5 of prot_2VT4
select lig_2VT4, 2VT4 and resn P32 and chain A
create 2VT4_gpcr, prot_2VT4 or wat_2VT4 or lig_2VT4

# 2Y02: Chain A, full agonist (isoprenaline=WHJ), active-like, ADRB1
select prot_2Y02, 2Y02 and chain A and polymer.protein and resi 33-358
select wat_2Y02, 2Y02 and resn HOH within 5 of prot_2Y02
select lig_2Y02, 2Y02 and resn WHJ and chain A
create 2Y02_gpcr, prot_2Y02 or wat_2Y02 or lig_2Y02

# 2Y04: Chain A, partial agonist (salbutamol=68H), intermediate, ADRB1
select prot_2Y04, 2Y04 and chain A and polymer.protein and resi 33-358
select wat_2Y04, 2Y04 and resn HOH within 5 of prot_2Y04
select lig_2Y04, 2Y04 and resn 68H and chain A
create 2Y04_gpcr, prot_2Y04 or wat_2Y04 or lig_2Y04

# 2Y00: Chain A, partial agonist (dobutamine=Y00), intermediate, ADRB1
select prot_2Y00, 2Y00 and chain A and polymer.protein and resi 33-358
select wat_2Y00, 2Y00 and resn HOH within 5 of prot_2Y00
select lig_2Y00, 2Y00 and resn Y00 and chain A
create 2Y00_gpcr, prot_2Y00 or wat_2Y00 or lig_2Y00

# Delete original structures and temporary selections
delete 2RH1
delete 3NY9
delete 3SN6
delete 4LDO
delete 2VT4
delete 2Y02
delete 2Y04
delete 2Y00
delete prot_*
delete wat_*
delete lig_*

print("GPCR chains extracted.")

# =============================================================================
# STEP 3: Align all structures to reference (2RH1)
# =============================================================================
print("Aligning structures to 2RH1...")

# Use align for sequence-based alignment
align 3NY9_gpcr and name CA, 2RH1_gpcr and name CA
align 4LDO_gpcr and name CA, 2RH1_gpcr and name CA
align 2VT4_gpcr and name CA, 2RH1_gpcr and name CA
align 2Y02_gpcr and name CA, 2RH1_gpcr and name CA
align 2Y04_gpcr and name CA, 2RH1_gpcr and name CA
align 2Y00_gpcr and name CA, 2RH1_gpcr and name CA

# 3SN6 has different residue numbering - use cealign (structure-based)
cealign 2RH1_gpcr and polymer.protein, 3SN6_gpcr and polymer.protein

print("Alignment complete.")

# =============================================================================
# STEP 4: Set up visual representations
# =============================================================================
print("Setting up representations...")

hide everything

# Protein: elegant fancy helices
show cartoon, polymer.protein
set cartoon_fancy_helices, 1
set cartoon_dumbbell_width, 0.06
set cartoon_dumbbell_length, 0.8
set cartoon_dumbbell_radius, 0.05
set cartoon_loop_radius, 0.08
set cartoon_tube_radius, 0.08
set cartoon_transparency, 0.5

# Color by receptor state
# Inactive states: blue shades
color marine, 2RH1_gpcr and polymer.protein
color slate, 3NY9_gpcr and polymer.protein
color deepteal, 2VT4_gpcr and polymer.protein

# Active states: green shades
color forest, 3SN6_gpcr and polymer.protein
color lime, 4LDO_gpcr and polymer.protein
color chartreuse, 2Y02_gpcr and polymer.protein

# Intermediate states: yellow/orange
color yellow, 2Y04_gpcr and polymer.protein
color orange, 2Y00_gpcr and polymer.protein

# Ligands: sticks, colored by ligand category
# Two separate color scales: Agonist (greens) vs Inverse Agonist/Antagonist (magentas)
show sticks, resn CAU+JSZ+P0G+ALE+P32+WHJ+68H+Y00
set stick_radius, 0.25

# Inverse agonist/antagonist scale (magenta-pink-red)
color magenta, 2RH1_gpcr and resn CAU
color hotpink, 3NY9_gpcr and resn JSZ
color salmon, 2VT4_gpcr and resn P32

# Agonist scale (greens)
color forest, 3SN6_gpcr and resn P0G
color lime, 4LDO_gpcr and resn ALE
color chartreuse, 2Y02_gpcr and resn WHJ
color splitpea, 2Y04_gpcr and resn 68H
color olive, 2Y00_gpcr and resn Y00

# Waters: hide initially (pocket scripts will show them)
hide spheres, resn HOH
color cyan, resn HOH
set sphere_scale, 0.2, resn HOH

print("Representations set.")

# =============================================================================
# STEP 5: Create groups for easy toggling
# =============================================================================
print("Creating groups...")

# Groups by receptor subtype (each object can only be in ONE group)
group adrb1, 2VT4_gpcr, add
group adrb1, 2Y02_gpcr, add
group adrb1, 2Y04_gpcr, add
group adrb1, 2Y00_gpcr, add

group adrb2, 2RH1_gpcr, add
group adrb2, 3NY9_gpcr, add
group adrb2, 3SN6_gpcr, add
group adrb2, 4LDO_gpcr, add

# Selections for other categorizations (objects can be in multiple selections)
# By receptor state
select active_states, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr
select inactive_states, 2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr
select intermediate_states, 2Y04_gpcr or 2Y00_gpcr

# By ligand type
select agonist_bound, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr or 2Y04_gpcr or 2Y00_gpcr
select inverse_agonist_bound, 2RH1_gpcr or 3NY9_gpcr

# Other selections
select all_ligands, resn CAU+JSZ+P0G+ALE+P32+WHJ+68H+Y00
select all_waters, resn HOH

print("Groups created.")

# =============================================================================
# STEP 6: Store views
# =============================================================================
print("Storing views...")

# Store the zoomed-in pocket view (used by pocket scripts)
set_view (\
    -0.953111351,    0.207114682,   -0.220622301,\
    -0.187508062,   -0.976453245,   -0.106624432,\
    -0.237513170,   -0.060255140,    0.969507515,\
     0.000000000,    0.000000000,  -73.409294128,\
   -30.602531433,    8.374164581,    7.177645683,\
    57.731933594,   89.086647034,  -20.000000000 )
view pocket_view, store

# Set and store the zoomed-out overall view (default for main script)
set_view (\
    -0.953111351,    0.207114682,   -0.220622301,\
    -0.187508062,   -0.976453245,   -0.106624432,\
    -0.237513170,   -0.060255140,    0.969507515,\
     0.000000000,    0.000000000, -269.608856201,\
   -30.602531433,    8.374164581,    7.177645683,\
   211.258865356,  327.958648682,  -20.000000000 )
view overall_view, store

print("Views stored: pocket_view, overall_view")

# =============================================================================
# STEP 7: Define helper commands
# =============================================================================

python

from pymol import cmd

def save_session():
    cmd.save("/tmp/gpcr_comparison.pse")
    print("Session saved to /tmp/gpcr_comparison.pse")

def render_current(filename="gpcr_view"):
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    path = "/tmp/gpcr_output/" + filename + ".png"
    import os
    os.makedirs("/tmp/gpcr_output", exist_ok=True)
    cmd.png(path, width=1200, height=1200, dpi=150, ray=1)
    print("Rendered: " + path)

cmd.extend("save_session", save_session)
cmd.extend("render_current", render_current)

python end

# =============================================================================
# DONE
# =============================================================================

python
print("")
print("=" * 60)
print("GPCR VISUALIZATION LOADED")
print("=" * 60)
print("")
print("Structures:")
print("  ADRB2: 2RH1, 3NY9, 3SN6, 4LDO")
print("  ADRB1: 2VT4, 2Y02, 2Y04, 2Y00")
print("")
print("Load pocket views with:")
print("  @scripts/gpcr_visualization_pocket1.pml  (ADRB1)")
print("  @scripts/gpcr_visualization_pocket2.pml  (ADRB2)")
print("  @scripts/gpcr_visualization_pocket3.pml  (Active)")
print("  @scripts/gpcr_visualization_pocket4.pml  (Inactive)")
print("  @scripts/gpcr_visualization_pocket5.pml  (Agonist)")
print("  @scripts/gpcr_visualization_pocket6.pml  (Inverse Agonist)")
print("")
print("Commands: save_session, render_current [filename]")
print("=" * 60)
python end

deselect
