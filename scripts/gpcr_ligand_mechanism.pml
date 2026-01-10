# GPCR Ligand Mechanism Visualization
# =====================================
#
# Complete visualization for testing 4 hypotheses about ligand-dependent
# structural features in adrenergic receptors (ADRB1, ADRB2).
#
# Hypotheses:
#   H1: Agonists bind CLOSER to S5.43 than inverse agonists
#   H2: Inverse agonists bind CLOSER to W6.48 (toggle switch) than agonists
#   H3: Water at N6.55 is EXCLUSIVE to agonist-bound active structures
#   H4: D2.50-W6.48 distance is SHORTEST in inverse agonists
#
# Usage:
#   pymol scripts/gpcr_ligand_mechanism.pml
#
# Navigation:
#   Use scene buttons (H1A, H1B, H2A, H2B, H3A, H3B, H4A, H4B)
#   Or: scene H1A, recall
#
# Associated workflow:
#   python examples/workflow_tests/gpcr_ligand_mechanism.py
#
# =============================================================================
# STRUCTURE REFERENCES (DOIs for verification)
# =============================================================================
#
# FULL AGONISTS (active state):
#   3SN6 - ADRB2 + BI-167107 + Gs protein (FULLY ACTIVE)
#          DOI: 10.1038/nature10361 (Rasmussen et al. 2011)
#   4LDO - ADRB2 + adrenaline + Nb6B9 nanobody (active)
#          DOI: 10.1038/nature12572 (Ring et al. 2013)
#   2Y02 - ADRB1 + isoprenaline (active-like, NO G protein)
#          DOI: 10.1038/nature09746 (Warne et al. 2011)
#
# PARTIAL AGONISTS (intermediate state):
#   2Y04 - ADRB1 + salbutamol
#          DOI: 10.1038/nature09746 (Warne et al. 2011)
#   2Y00 - ADRB1 + dobutamine
#          DOI: 10.1038/nature09746 (Warne et al. 2011)
#
# INVERSE AGONISTS (inactive state):
#   2RH1 - ADRB2 + carazolol (inactive)
#          DOI: 10.1126/science.1150577 (Cherezov et al. 2007)
#   3NY9 - ADRB2 + ICI 118,551 (inactive)
#          DOI: 10.1021/ja105108q (Wacker et al. 2010)
#
# ANTAGONIST (inactive state):
#   2VT4 - ADRB1 + cyanopindolol (inactive)
#          DOI: 10.1038/nature07101 (Warne et al. 2008)
#
# =============================================================================

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

# Enable scene buttons
set scene_buttons, 1

# =============================================================================
# STEP 1: Fetch structures from PDB
# =============================================================================
print "Fetching structures from PDB..."

fetch 2RH1, type=cif, async=0
fetch 3NY9, type=cif, async=0
fetch 3SN6, type=cif, async=0
fetch 4LDO, type=cif, async=0
fetch 2VT4, type=cif, async=0
fetch 2Y02, type=cif, async=0
fetch 2Y04, type=cif, async=0
fetch 2Y00, type=cif, async=0

print "All structures loaded."

# =============================================================================
# STEP 2: Extract GPCR chains with waters and ligands
# =============================================================================
print "Extracting GPCR chains..."

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

print "GPCR chains extracted."

# =============================================================================
# STEP 3: Align all structures to reference (2RH1)
# =============================================================================
print "Aligning structures to 2RH1..."

# Use align for sequence-based alignment
align 3NY9_gpcr and name CA, 2RH1_gpcr and name CA
align 4LDO_gpcr and name CA, 2RH1_gpcr and name CA
align 2VT4_gpcr and name CA, 2RH1_gpcr and name CA
align 2Y02_gpcr and name CA, 2RH1_gpcr and name CA
align 2Y04_gpcr and name CA, 2RH1_gpcr and name CA
align 2Y00_gpcr and name CA, 2RH1_gpcr and name CA

# 3SN6 has different residue numbering - use cealign (structure-based)
cealign 2RH1_gpcr and polymer.protein, 3SN6_gpcr and polymer.protein

print "Alignment complete."

# =============================================================================
# STEP 4: Set up visual representations
# =============================================================================
print "Setting up representations..."

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

# Waters: hide initially (scenes will show them)
hide spheres, resn HOH
color cyan, resn HOH
set sphere_scale, 0.2, resn HOH

print "Representations set."

# =============================================================================
# STEP 5: Create groups for easy toggling
# =============================================================================
print "Creating groups..."

# Groups by receptor subtype
group adrb1, 2VT4_gpcr, add
group adrb1, 2Y02_gpcr, add
group adrb1, 2Y04_gpcr, add
group adrb1, 2Y00_gpcr, add

group adrb2, 2RH1_gpcr, add
group adrb2, 3NY9_gpcr, add
group adrb2, 3SN6_gpcr, add
group adrb2, 4LDO_gpcr, add

# Selections for other categorizations
select active_states, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr
select inactive_states, 2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr
select intermediate_states, 2Y04_gpcr or 2Y00_gpcr

# By ligand type
select agonist_bound, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr or 2Y04_gpcr or 2Y00_gpcr
select inverse_agonist_bound, 2RH1_gpcr or 3NY9_gpcr

# Other selections
select all_ligands, resn CAU+JSZ+P0G+ALE+P32+WHJ+68H+Y00
select all_waters, resn HOH

print "Groups created."

# =============================================================================
# STEP 6: Store views
# =============================================================================
print "Storing views..."

# Store the zoomed-in pocket view
set_view (\
    -0.953111351,    0.207114682,   -0.220622301,\
    -0.187508062,   -0.976453245,   -0.106624432,\
    -0.237513170,   -0.060255140,    0.969507515,\
     0.000000000,    0.000000000,  -73.409294128,\
   -30.602531433,    8.374164581,    7.177645683,\
    57.731933594,   89.086647034,  -20.000000000 )
view pocket_view, store

# Set and store the zoomed-out overall view
set_view (\
    -0.953111351,    0.207114682,   -0.220622301,\
    -0.187508062,   -0.976453245,   -0.106624432,\
    -0.237513170,   -0.060255140,    0.969507515,\
     0.000000000,    0.000000000, -269.608856201,\
   -30.602531433,    8.374164581,    7.177645683,\
   211.258865356,  327.958648682,  -20.000000000 )
view overall_view, store

print "Views stored: pocket_view, overall_view"

# =============================================================================
# STEP 7: Define helper commands
# =============================================================================

python

from pymol import cmd
import os

def save_session():
    cmd.save("/tmp/gpcr_ligand_mechanism.pse")
    print("Session saved to /tmp/gpcr_ligand_mechanism.pse")

def render_current(filename="gpcr_view"):
    cmd.bg_color("white")
    cmd.set("ray_shadows", 0)
    path = "/tmp/gpcr_output/" + filename + ".png"
    os.makedirs("/tmp/gpcr_output", exist_ok=True)
    cmd.png(path, width=1200, height=1200, dpi=150, ray=1)
    print("Rendered: " + path)

def render_all_scenes():
    scenes = ["H1A", "H1B", "H2A", "H2B", "H3A", "H3B", "H4A", "H4B"]
    os.makedirs("/tmp/gpcr_output", exist_ok=True)
    for scene in scenes:
        cmd.scene(scene, action="recall")
        path = "/tmp/gpcr_output/" + scene + ".png"
        cmd.png(path, width=1200, height=1200, dpi=150, ray=1)
        print("Rendered: " + path)

cmd.extend("save_session", save_session)
cmd.extend("render_current", render_current)
cmd.extend("render_all_scenes", render_all_scenes)

python end

# =============================================================================
# HYPOTHESIS VISUALIZATION SCENES
# =============================================================================
print "Setting up Hypothesis Visualization Scenes..."

# Color definitions for hypothesis views
set_color agonist_green, [0.2, 0.8, 0.2]
set_color inverse_magenta, [0.9, 0.2, 0.6]
set_color highlight_yellow, [1.0, 1.0, 0.0]
set_color water_cyan, [0.0, 0.8, 1.0]
set_color contact_red, [1.0, 0.3, 0.3]

# Common settings for all hypothesis views
set label_size, 24
set label_color, black
set label_font_id, 7
set dash_width, 2.5
set dash_gap, 0.2

# =============================================================================
# HYPOTHESIS 1: Agonists bind CLOSER to S5.43 (Serine)
# Evidence: Agonist 2.82-3.08A vs Inverse agonist 3.32-3.58A
# =============================================================================

print "Creating H1 scenes: S5.43 (Serine) contact..."

# Create selections for 5.43 residues
select S543_3SN6, 3SN6_gpcr and chain R and resi 203
select S543_4LDO, 4LDO_gpcr and chain A and resi 1203
select S543_2Y02, 2Y02_gpcr and chain A and resi 211

select agonist_S543, S543_3SN6 or S543_4LDO or S543_2Y02

select S543_2RH1, 2RH1_gpcr and chain A and resi 203
select S543_3NY9, 3NY9_gpcr and chain A and resi 203

select inverse_S543, S543_2RH1 or S543_3NY9

# -----------------------------------------------------------------------------
# H1A: Agonists close to S5.43
# -----------------------------------------------------------------------------

print "Storing scene H1A..."
disable all
enable 3SN6_gpcr
enable 4LDO_gpcr
enable 2Y02_gpcr

show cartoon, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr

show sticks, agonist_S543
color highlight_yellow, agonist_S543
set stick_radius, 0.25, agonist_S543

show sticks, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic
color agonist_green, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic

# Distance measurements
distance H1A_dist_3SN6, S543_3SN6 and name OG, 3SN6_gpcr and resn P0G and name OAE
distance H1A_dist_4LDO, S543_4LDO and name OG, 4LDO_gpcr and resn ALE and name O1
distance H1A_dist_2Y02, S543_2Y02 and name OG, 2Y02_gpcr and resn WHJ and name N1

color contact_red, H1A_dist_*
set dash_color, contact_red, H1A_dist_*

label S543_3SN6 and name CA, "S5.43"

set_view (\
     0.997919381,   -0.009538483,    0.063835151,\
    -0.002839159,   -0.982710123,   -0.185714707,\
     0.057281435,    0.187984124,   -0.980179429,\
     0.000000000,    0.000000000, -152.501068115,\
   -31.200321198,    9.068758011,    7.718489647,\
    72.501121521,  232.501098633,  -20.000000000 )

group H1A_objects, H1A_dist_*
scene H1A, store, Agonists close to S5.43

# -----------------------------------------------------------------------------
# H1B: Inverse agonists far from S5.43
# -----------------------------------------------------------------------------

print "Storing scene H1B..."
disable all
enable 2RH1_gpcr
enable 3NY9_gpcr

show cartoon, 2RH1_gpcr or 3NY9_gpcr

show sticks, inverse_S543
color highlight_yellow, inverse_S543
set stick_radius, 0.25, inverse_S543

show sticks, (2RH1_gpcr or 3NY9_gpcr) and organic
color inverse_magenta, (2RH1_gpcr or 3NY9_gpcr) and organic

distance H1B_dist_2RH1, S543_2RH1 and name OG, 2RH1_gpcr and resn CAU and name N7
distance H1B_dist_3NY9, S543_3NY9 and name OG, 3NY9_gpcr and resn JSZ and name O3

color contact_red, H1B_dist_*
set dash_color, contact_red, H1B_dist_*

label S543_2RH1 and name CA, "S5.43"

set_view (\
     0.997919381,   -0.009538483,    0.063835151,\
    -0.002839159,   -0.982710123,   -0.185714707,\
     0.057281435,    0.187984124,   -0.980179429,\
     0.000000000,    0.000000000, -152.501068115,\
   -31.200321198,    9.068758011,    7.718489647,\
    72.501121521,  232.501098633,  -20.000000000 )

group H1B_objects, H1B_dist_*
scene H1B, store, Inverse agonists far from S5.43

# =============================================================================
# HYPOTHESIS 2: Inverse agonists bind CLOSER to W6.48 (Toggle switch)
# Evidence: Inverse 3.42-3.44A vs Agonist 3.96A or no contact
# =============================================================================

print "Creating H2 scenes: W6.48 (Toggle switch) contact..."

select W648_2RH1, 2RH1_gpcr and chain A and resi 286
select W648_3NY9, 3NY9_gpcr and chain A and resi 286
select inverse_W648, W648_2RH1 or W648_3NY9

select W648_3SN6, 3SN6_gpcr and chain R and resi 286
select W648_4LDO, 4LDO_gpcr and chain A and resi 1286
select W648_2Y02, 2Y02_gpcr and chain A and resi 303
select agonist_W648, W648_3SN6 or W648_4LDO or W648_2Y02

# -----------------------------------------------------------------------------
# H2A: Inverse agonists close to W6.48
# -----------------------------------------------------------------------------

print "Storing scene H2A..."
disable all
enable 2RH1_gpcr
enable 3NY9_gpcr

show cartoon, 2RH1_gpcr or 3NY9_gpcr

show sticks, inverse_W648
color highlight_yellow, inverse_W648
set stick_radius, 0.25, inverse_W648

show sticks, (2RH1_gpcr or 3NY9_gpcr) and organic
color inverse_magenta, (2RH1_gpcr or 3NY9_gpcr) and organic

distance H2A_dist_2RH1, W648_2RH1 and name CH2, 2RH1_gpcr and resn CAU and name O17
distance H2A_dist_3NY9, W648_3NY9 and name CH2, 3NY9_gpcr and resn JSZ and name O5

color contact_red, H2A_dist_*
set dash_color, contact_red, H2A_dist_*

label W648_2RH1 and name CA, "W6.48"

set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

group H2A_objects, H2A_dist_*
scene H2A, store, Inverse agonists close to W6.48

# -----------------------------------------------------------------------------
# H2B: Agonists far from W6.48
# -----------------------------------------------------------------------------

print "Storing scene H2B..."
disable all
enable 3SN6_gpcr
enable 4LDO_gpcr
enable 2Y02_gpcr

show cartoon, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr

show sticks, agonist_W648
color highlight_yellow, agonist_W648
set stick_radius, 0.25, agonist_W648

show sticks, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic
color agonist_green, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic

# Only 2Y02 has a weak contact
distance H2B_dist_2Y02, W648_2Y02 and name CZ3, 2Y02_gpcr and resn WHJ and name O4

color contact_red, H2B_dist_*
set dash_color, contact_red, H2B_dist_*

label W648_2Y02 and name CA, "W6.48"

set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

group H2B_objects, H2B_dist_*
scene H2B, store, Agonists far from W6.48

# =============================================================================
# HYPOTHESIS 3: Water at N6.55 is EXCLUSIVE to ACTIVE state
# Evidence: Active (4LDO, 2Y02) have water; Inactive have none
# =============================================================================

print "Creating H3 scenes: N6.55 water (Active state marker)..."

select N655_4LDO, 4LDO_gpcr and chain A and resi 1293
select N655_2Y02, 2Y02_gpcr and chain A and resi 310
select active_N655, N655_4LDO or N655_2Y02

# Waters at 6.55 (only in active structures)
select water_655_4LDO, 4LDO_gpcr and resn HOH and resi 1501
select water_655_2Y02, 2Y02_gpcr and resn HOH and resi 2010
select active_water_655, water_655_4LDO or water_655_2Y02

select N655_2RH1, 2RH1_gpcr and chain A and resi 293
select N655_3NY9, 3NY9_gpcr and chain A and resi 293
select N655_2VT4, 2VT4_gpcr and chain A and resi 310
select inactive_N655, N655_2RH1 or N655_3NY9 or N655_2VT4

# -----------------------------------------------------------------------------
# H3A: Active structures WITH water at N6.55
# -----------------------------------------------------------------------------

print "Storing scene H3A..."
disable all
enable 4LDO_gpcr
enable 2Y02_gpcr

show cartoon, 4LDO_gpcr or 2Y02_gpcr

show sticks, active_N655
color highlight_yellow, active_N655
set stick_radius, 0.25, active_N655

show spheres, active_water_655
color water_cyan, active_water_655
set sphere_scale, 0.4, active_water_655

show sticks, (4LDO_gpcr or 2Y02_gpcr) and organic
color agonist_green, (4LDO_gpcr or 2Y02_gpcr) and organic

distance H3A_wat_4LDO, N655_4LDO and name OD1, water_655_4LDO and name O
distance H3A_wat_2Y02, N655_2Y02 and name OD1, water_655_2Y02 and name O

color water_cyan, H3A_wat_*
set dash_color, water_cyan, H3A_wat_*

label N655_4LDO and name CA, "N6.55"

set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

group H3A_objects, H3A_wat_*
scene H3A, store, Active state: water at N6.55

# -----------------------------------------------------------------------------
# H3B: Inactive structures WITHOUT water at N6.55
# -----------------------------------------------------------------------------

print "Storing scene H3B..."
disable all
enable 2RH1_gpcr
enable 3NY9_gpcr
enable 2VT4_gpcr

show cartoon, 2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr

show sticks, inactive_N655
color highlight_yellow, inactive_N655
set stick_radius, 0.25, inactive_N655

# Show all waters to demonstrate ABSENCE at 6.55
show spheres, (2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and resn HOH
color gray50, (2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and resn HOH
set sphere_scale, 0.25, resn HOH

show sticks, (2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and organic
color inverse_magenta, (2RH1_gpcr or 3NY9_gpcr) and organic
color gray70, 2VT4_gpcr and organic

label N655_2RH1 and name CA, "N6.55"

set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

scene H3B, store, Inactive state: NO water at N6.55

# =============================================================================
# HYPOTHESIS 4: D2.50-W6.48 Distance is SHORTEST in Inverse Agonists
# Evidence: Inverse agonist 6.51-6.64A vs Agonist 6.82-6.99A vs Antagonist 7.22A
# =============================================================================

print "Creating H4 scenes: D2.50-W6.48 distance..."

select D250_2RH1, 2RH1_gpcr and chain A and resi 79
select D250_3NY9, 3NY9_gpcr and chain A and resi 79
select inverse_D250, D250_2RH1 or D250_3NY9

select D250_3SN6, 3SN6_gpcr and chain R and resi 79
select D250_4LDO, 4LDO_gpcr and chain A and resi 1079
select D250_2VT4, 2VT4_gpcr and chain A and resi 87
select W648_2VT4, 2VT4_gpcr and chain A and resi 303

select compare_D250, D250_3SN6 or D250_4LDO or D250_2VT4

# -----------------------------------------------------------------------------
# H4A: Inverse agonist (2RH1) - SHORT D2.50-W6.48 distance (6.64A)
# Shows water network between D2.50 and W6.48 in cyan
# -----------------------------------------------------------------------------

print "Storing scene H4A..."
disable all
enable 2RH1_gpcr

show cartoon, 2RH1_gpcr

# Show D2.50 as sticks (sodium site)
show sticks, D250_2RH1
color orange, D250_2RH1
set stick_radius, 0.25, D250_2RH1

# Show W6.48 as sticks (toggle switch)
show sticks, W648_2RH1
color highlight_yellow, W648_2RH1
set stick_radius, 0.25, W648_2RH1

# Additional pathway residues between D2.50 and W6.48
select H4A_I640_2RH1, 2RH1_gpcr and chain A and resi 278   # I6.40
select H4A_F644_2RH1, 2RH1_gpcr and chain A and resi 282   # F6.44
show sticks, H4A_I640_2RH1 or H4A_F644_2RH1
color palegreen, H4A_I640_2RH1 or H4A_F644_2RH1
set stick_radius, 0.2, H4A_I640_2RH1 or H4A_F644_2RH1

# Waters in the D2.50-W6.48 pathway (2RH1)
select H4A_waters_2RH1, 2RH1_gpcr and resn HOH and (resi 528 or resi 532 or resi 534 or resi 548)
show spheres, H4A_waters_2RH1
color water_cyan, H4A_waters_2RH1
set sphere_scale, 0.4, H4A_waters_2RH1

# Show ligand
show sticks, 2RH1_gpcr and organic
color inverse_magenta, 2RH1_gpcr and organic

# Show polar contacts (H-bonds) within 3.5A cutoff
distance H4A_hbonds, (D250_2RH1 or W648_2RH1 or H4A_I640_2RH1 or H4A_F644_2RH1 or H4A_waters_2RH1), (D250_2RH1 or W648_2RH1 or H4A_I640_2RH1 or H4A_F644_2RH1 or H4A_waters_2RH1), 3.5, mode=2
color water_cyan, H4A_hbonds
set dash_color, water_cyan, H4A_hbonds
hide labels, H4A_hbonds

# D2.50-W6.48 distance measurement (RED - key metric)
distance H4A_dist_main, D250_2RH1 and name OD1, W648_2RH1 and name NE1
color red, H4A_dist_main
set dash_color, red, H4A_dist_main
set dash_width, 4, H4A_dist_main

label D250_2RH1 and name CA, "D2.50"
label W648_2RH1 and name CA, "W6.48"

set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

group H4A_objects, H4A_* D250_2RH1 W648_2RH1
scene H4A, store, Inverse agonist: SHORT D2.50-W6.48 (6.64A)

# -----------------------------------------------------------------------------
# H4B: Agonist (2Y02) - LONGER D2.50-W6.48 distance (6.95A)
# Shows water network between D2.50 and W6.48 in cyan
# -----------------------------------------------------------------------------

print "Storing scene H4B..."
disable all
enable 2Y02_gpcr

show cartoon, 2Y02_gpcr

# Show D2.50 as sticks (sodium site)
select D250_2Y02, 2Y02_gpcr and chain A and resi 87
show sticks, D250_2Y02
color orange, D250_2Y02
set stick_radius, 0.25, D250_2Y02

# Show W6.48 as sticks (toggle switch)
select W648_2Y02_h4, 2Y02_gpcr and chain A and resi 303
show sticks, W648_2Y02_h4
color highlight_yellow, W648_2Y02_h4
set stick_radius, 0.25, W648_2Y02_h4

# Additional pathway residues between D2.50 and W6.48
select H4B_I640_2Y02, 2Y02_gpcr and chain A and resi 295   # I6.40
select H4B_F644_2Y02, 2Y02_gpcr and chain A and resi 299   # F6.44
show sticks, H4B_I640_2Y02 or H4B_F644_2Y02
color palegreen, H4B_I640_2Y02 or H4B_F644_2Y02
set stick_radius, 0.2, H4B_I640_2Y02 or H4B_F644_2Y02

# Waters in the D2.50-W6.48 pathway (2Y02)
select H4B_waters_2Y02, 2Y02_gpcr and resn HOH and (resi 2002 or resi 2003 or resi 2016 or resi 2017)
show spheres, H4B_waters_2Y02
color water_cyan, H4B_waters_2Y02
set sphere_scale, 0.4, H4B_waters_2Y02

# Show ligand
show sticks, 2Y02_gpcr and organic
color agonist_green, 2Y02_gpcr and organic

# Show polar contacts (H-bonds) within 3.5A cutoff
distance H4B_hbonds, (D250_2Y02 or W648_2Y02_h4 or H4B_I640_2Y02 or H4B_F644_2Y02 or H4B_waters_2Y02), (D250_2Y02 or W648_2Y02_h4 or H4B_I640_2Y02 or H4B_F644_2Y02 or H4B_waters_2Y02), 3.5, mode=2
color water_cyan, H4B_hbonds
set dash_color, water_cyan, H4B_hbonds
hide labels, H4B_hbonds

# D2.50-W6.48 distance measurement (RED - key metric)
distance H4B_dist_main, D250_2Y02 and name OD1, W648_2Y02_h4 and name NE1
color red, H4B_dist_main
set dash_color, red, H4B_dist_main
set dash_width, 4, H4B_dist_main

label D250_2Y02 and name CA, "D2.50"
label W648_2Y02_h4 and name CA, "W6.48"

set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

group H4B_objects, H4B_* D250_2Y02 W648_2Y02_h4
scene H4B, store, Agonist: LONGER D2.50-W6.48 (6.95A)

# =============================================================================
# CLEANUP AND SUMMARY
# =============================================================================

# Hide selections from menu
disable S543_*
disable W648_*
disable N655_*
disable D250_*
disable water_*
disable agonist_*
disable inverse_*
disable active_*
disable inactive_*
disable compare_*

# Ensure scene buttons are visible
set scene_buttons, 1

# Reset to first scene
scene H1A, recall

deselect

python
print("")
print("=" * 60)
print("GPCR LIGAND MECHANISM VISUALIZATION LOADED")
print("=" * 60)
print("")
print("Structures (8 adrenergic receptors):")
print("  Full agonists:    3SN6, 4LDO (ADRB2), 2Y02 (ADRB1)")
print("  Partial agonists: 2Y04, 2Y00 (ADRB1)")
print("  Inverse agonists: 2RH1, 3NY9 (ADRB2)")
print("  Antagonist:       2VT4 (ADRB1)")
print("")
print("Hypothesis scenes:")
print("  H1A - Agonists close to S5.43 (2.82-3.25A)")
print("  H1B - Inverse agonists far from S5.43 (3.32-3.58A)")
print("")
print("  H2A - Inverse agonists close to W6.48 (3.42-3.44A)")
print("  H2B - Agonists far/no contact with W6.48")
print("")
print("  H3A - Agonist: water at N6.55 (active marker)")
print("  H3B - Inactive: NO water at N6.55")
print("")
print("  H4A - Inverse agonist: SHORT D2.50-W6.48 (6.64A)")
print("  H4B - Agonist: LONGER D2.50-W6.48 (6.95A)")
print("")
print("Usage: scene <name>, recall")
print("Example: scene H1A, recall")
print("")
print("Commands: save_session, render_current [filename], render_all_scenes")
print("=" * 60)
python end
