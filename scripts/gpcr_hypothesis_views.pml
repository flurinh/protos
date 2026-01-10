# GPCR Mechanism Hypothesis Visualizations
# =========================================
# 8 scenes demonstrating 4 hypotheses (2 contrasting views each)
#
# Run after: pymol scripts/gpcr_visualization.pml
#
# Usage:
#   scene H1A, recall   # Agonists close to S5.43
#   scene H1B, recall   # Inverse agonists far from S5.43
#   etc.
#
# =============================================================================
# STRUCTURE REFERENCES (DOIs for verification)
# =============================================================================
#
# ACTIVE STATE REPRESENTATIVES:
#   3SN6 - ADRB2 + BI-167107 agonist + Gs protein (FULLY ACTIVE)
#          DOI: 10.1038/nature10361 (Rasmussen et al. 2011)
#   4LDO - ADRB2 + adrenaline + Nb6B9 nanobody (active)
#          DOI: 10.1038/nature12572 (Ring et al. 2013)
#   2Y02 - ADRB1 + isoprenaline (active-like, NO G protein)
#          DOI: 10.1038/nature09746 (Warne et al. 2011)
#
# INACTIVE STATE REPRESENTATIVES:
#   2RH1 - ADRB2 + carazolol inverse agonist (inactive)
#          DOI: 10.1126/science.1150577 (Cherezov et al. 2007)
#   3NY9 - ADRB2 + ICI 118,551 inverse agonist (inactive)
#          DOI: 10.1021/ja105108q (Wacker et al. 2010)
#   2VT4 - ADRB1 + cyanopindolol antagonist (inactive)
#          DOI: 10.1038/nature07101 (Warne et al. 2008)
#
# INTERMEDIATE STATE:
#   2Y04 - ADRB1 + salbutamol partial agonist
#          DOI: 10.1038/nature09746 (Warne et al. 2011)
#   2Y00 - ADRB1 + dobutamine partial agonist
#          DOI: 10.1038/nature09746 (Warne et al. 2011)
#
# NOTE on 2Y02: While bound to full agonist, this structure lacks the G protein
# and shows "agonist-induced conformational changes" but is NOT fully active
# like 3SN6. Hence labeled "active_like" rather than "active".
# =============================================================================

print("Setting up Hypothesis Visualization Scenes...")

# Enable scene buttons panel
set scene_buttons, 1

# =============================================================================
# Color definitions for hypothesis views
# =============================================================================

# Ligand colors by type (consistent with main script)
set_color agonist_green, [0.2, 0.8, 0.2]
set_color inverse_magenta, [0.9, 0.2, 0.6]

# Highlight colors
set_color highlight_yellow, [1.0, 1.0, 0.0]
set_color water_cyan, [0.0, 0.8, 1.0]
set_color contact_red, [1.0, 0.3, 0.3]

# =============================================================================
# HYPOTHESIS 1: Agonists bind CLOSER to S5.43 (Serine)
# Evidence: Agonist 2.82-3.08A vs Inverse agonist 3.32-3.58A
# =============================================================================

print("Creating H1 scenes: S5.43 (Serine) contact...")

# --- H1A: Agonists close to S5.43 ---
# Structures: 3SN6, 4LDO, 2Y02

# Create selections for 5.43 residues
select S543_3SN6, 3SN6_gpcr and chain R and resi 203
select S543_4LDO, 4LDO_gpcr and chain A and resi 1203
select S543_2Y02, 2Y02_gpcr and chain A and resi 211

# Group agonist structures
select agonist_S543, S543_3SN6 or S543_4LDO or S543_2Y02

# --- H1B: Inverse agonists far from S5.43 ---
# Structures: 2RH1, 3NY9

select S543_2RH1, 2RH1_gpcr and chain A and resi 203
select S543_3NY9, 3NY9_gpcr and chain A and resi 203

select inverse_S543, S543_2RH1 or S543_3NY9

# =============================================================================
# HYPOTHESIS 2: Inverse agonists bind CLOSER to W6.48 (Toggle switch)
# Evidence: Inverse 3.42-3.44A vs Agonist 3.96A or no contact
# =============================================================================

print("Creating H2 scenes: W6.48 (Toggle switch) contact...")

# --- H2A: Inverse agonists close to W6.48 ---
select W648_2RH1, 2RH1_gpcr and chain A and resi 286
select W648_3NY9, 3NY9_gpcr and chain A and resi 286

select inverse_W648, W648_2RH1 or W648_3NY9

# --- H2B: Agonists far from W6.48 ---
select W648_3SN6, 3SN6_gpcr and chain R and resi 286
select W648_4LDO, 4LDO_gpcr and chain A and resi 1286
select W648_2Y02, 2Y02_gpcr and chain A and resi 303

select agonist_W648, W648_3SN6 or W648_4LDO or W648_2Y02

# =============================================================================
# HYPOTHESIS 3: Water at N6.55 is EXCLUSIVE to ACTIVE state
# Evidence: Active (4LDO, 2Y02) have water; Inactive have none
# =============================================================================

print("Creating H3 scenes: N6.55 water (Active state marker)...")

# --- H3A: Active structures WITH water at 6.55 ---
select N655_4LDO, 4LDO_gpcr and chain A and resi 1293
select N655_2Y02, 2Y02_gpcr and chain A and resi 310

select active_N655, N655_4LDO or N655_2Y02

# Waters at 6.55 (only in active structures)
select water_655_4LDO, 4LDO_gpcr and resn HOH and resi 1501
select water_655_2Y02, 2Y02_gpcr and resn HOH and resi 2010

select active_water_655, water_655_4LDO or water_655_2Y02

# --- H3B: Inactive structures WITHOUT water at 6.55 ---
select N655_2RH1, 2RH1_gpcr and chain A and resi 293
select N655_3NY9, 3NY9_gpcr and chain A and resi 293
select N655_2VT4, 2VT4_gpcr and chain A and resi 310

select inactive_N655, N655_2RH1 or N655_3NY9 or N655_2VT4

# =============================================================================
# HYPOTHESIS 4: D2.50-W6.48 Distance is SHORTEST in Inverse Agonists
# Evidence: Inverse agonist 6.51-6.64A vs Agonist 6.82-6.99A vs Antagonist 7.22A
# =============================================================================

print("Creating H4 scenes: D2.50-W6.48 distance...")

# --- H4A: Inverse agonists with SHORT D2.50-W6.48 distance ---
# D2.50 (sodium site)
select D250_2RH1, 2RH1_gpcr and chain A and resi 79
select D250_3NY9, 3NY9_gpcr and chain A and resi 79

# W6.48 (toggle switch) - same residue as H2
# W648_2RH1 and W648_3NY9 already defined above

select inverse_D250, D250_2RH1 or D250_3NY9

# --- H4B: Agonist/Antagonist with LONGER D2.50-W6.48 distance ---
select D250_3SN6, 3SN6_gpcr and chain R and resi 79
select D250_4LDO, 4LDO_gpcr and chain A and resi 1079
select D250_2VT4, 2VT4_gpcr and chain A and resi 87

# W6.48 for comparison structures
select W648_2VT4, 2VT4_gpcr and chain A and resi 303

select compare_D250, D250_3SN6 or D250_4LDO or D250_2VT4

# =============================================================================
# SCENE SETUP FUNCTIONS
# =============================================================================

# Common settings for all hypothesis views
set label_size, 24
set label_color, black
set label_font_id, 7
set dash_width, 2.5
set dash_gap, 0.2

# -----------------------------------------------------------------------------
# H1A: Agonists close to S5.43
# -----------------------------------------------------------------------------

print("Storing scene H1A...")
disable all
enable 3SN6_gpcr
enable 4LDO_gpcr
enable 2Y02_gpcr

# Show protein as cartoon
show cartoon, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr

# Show S5.43 as sticks (highlighted)
show sticks, agonist_S543
color highlight_yellow, agonist_S543
set stick_radius, 0.25, agonist_S543

# Show ligands
show sticks, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic
color agonist_green, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic

# Distance measurements: closest atom pairs only
# 3SN6: SER OG <-> P0G OAE = 2.82A
distance H1A_dist_3SN6, S543_3SN6 and name OG, 3SN6_gpcr and resn P0G and name OAE
# 4LDO: SER OG <-> ALE O1 = 3.25A
distance H1A_dist_4LDO, S543_4LDO and name OG, 4LDO_gpcr and resn ALE and name O1
# 2Y02: SER OG <-> WHJ N1 = 2.88A
distance H1A_dist_2Y02, S543_2Y02 and name OG, 2Y02_gpcr and resn WHJ and name N1

color contact_red, H1A_dist_*
set dash_color, contact_red, H1A_dist_*

# Labels
label S543_3SN6 and name CA, "S5.43"

# Set view for H1A (S5.43 - agonists)
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

print("Storing scene H1B...")
disable all
enable 2RH1_gpcr
enable 3NY9_gpcr

show cartoon, 2RH1_gpcr or 3NY9_gpcr

# Show S5.43 as sticks
show sticks, inverse_S543
color highlight_yellow, inverse_S543
set stick_radius, 0.25, inverse_S543

# Show ligands
show sticks, (2RH1_gpcr or 3NY9_gpcr) and organic
color inverse_magenta, (2RH1_gpcr or 3NY9_gpcr) and organic

# Distance measurements: closest atom pairs only
# 2RH1: SER OG <-> CAU N7 = 3.32A
distance H1B_dist_2RH1, S543_2RH1 and name OG, 2RH1_gpcr and resn CAU and name N7
# 3NY9: SER OG <-> JSZ O3 = 3.58A
distance H1B_dist_3NY9, S543_3NY9 and name OG, 3NY9_gpcr and resn JSZ and name O3

color contact_red, H1B_dist_*
set dash_color, contact_red, H1B_dist_*

label S543_2RH1 and name CA, "S5.43"

# Set view for H1B (same as H1A - comparing same position)
set_view (\
     0.997919381,   -0.009538483,    0.063835151,\
    -0.002839159,   -0.982710123,   -0.185714707,\
     0.057281435,    0.187984124,   -0.980179429,\
     0.000000000,    0.000000000, -152.501068115,\
   -31.200321198,    9.068758011,    7.718489647,\
    72.501121521,  232.501098633,  -20.000000000 )

group H1B_objects, H1B_dist_*
scene H1B, store, Inverse agonists far from S5.43

# -----------------------------------------------------------------------------
# H2A: Inverse agonists close to W6.48 (toggle switch)
# -----------------------------------------------------------------------------

print("Storing scene H2A...")
disable all
enable 2RH1_gpcr
enable 3NY9_gpcr

show cartoon, 2RH1_gpcr or 3NY9_gpcr

# Show W6.48 as sticks
show sticks, inverse_W648
color highlight_yellow, inverse_W648
set stick_radius, 0.25, inverse_W648

# Show ligands
show sticks, (2RH1_gpcr or 3NY9_gpcr) and organic
color inverse_magenta, (2RH1_gpcr or 3NY9_gpcr) and organic

# Distance measurements: closest atom pairs only
# 2RH1: TRP CH2 <-> CAU O17 = 3.42A
distance H2A_dist_2RH1, W648_2RH1 and name CH2, 2RH1_gpcr and resn CAU and name O17
# 3NY9: TRP CH2 <-> JSZ O5 = 3.44A
distance H2A_dist_3NY9, W648_3NY9 and name CH2, 3NY9_gpcr and resn JSZ and name O5

color contact_red, H2A_dist_*
set dash_color, contact_red, H2A_dist_*

label W648_2RH1 and name CA, "W6.48"

# Set view for H2A (W6.48 - toggle switch)
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

print("Storing scene H2B...")
disable all
enable 3SN6_gpcr
enable 4LDO_gpcr
enable 2Y02_gpcr

show cartoon, 3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr

# Show W6.48 as sticks
show sticks, agonist_W648
color highlight_yellow, agonist_W648
set stick_radius, 0.25, agonist_W648

# Show ligands
show sticks, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic
color agonist_green, (3SN6_gpcr or 4LDO_gpcr or 2Y02_gpcr) and organic

# Distance measurements: closest atom pairs only
# Note: 3SN6 and 4LDO do NOT contact W6.48 - only 2Y02 has weak contact
# 2Y02: TRP CZ3 <-> WHJ O4 = 3.96A (weak/distant)
distance H2B_dist_2Y02, W648_2Y02 and name CZ3, 2Y02_gpcr and resn WHJ and name O4

color contact_red, H2B_dist_*
set dash_color, contact_red, H2B_dist_*

label W648_2Y02 and name CA, "W6.48"

# Set view for H2B (same as H2A - comparing same position)
set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

group H2B_objects, H2B_dist_*
scene H2B, store, Agonists far from W6.48

# -----------------------------------------------------------------------------
# H3A: Active structures WITH water at N6.55
# -----------------------------------------------------------------------------

print("Storing scene H3A...")
disable all
enable 4LDO_gpcr
enable 2Y02_gpcr

show cartoon, 4LDO_gpcr or 2Y02_gpcr

# Show N6.55 as sticks
show sticks, active_N655
color highlight_yellow, active_N655
set stick_radius, 0.25, active_N655

# Show waters at 6.55 as spheres
show spheres, active_water_655
color water_cyan, active_water_655
set sphere_scale, 0.4, active_water_655

# Show ligands
show sticks, (4LDO_gpcr or 2Y02_gpcr) and organic
color agonist_green, (4LDO_gpcr or 2Y02_gpcr) and organic

# Distance from N6.55 to water: closest atom pairs
# 4LDO: ASN OD1 <-> HOH O = 2.91A
distance H3A_wat_4LDO, N655_4LDO and name OD1, water_655_4LDO and name O
# 2Y02: ASN OD1 <-> HOH O = 3.36A
distance H3A_wat_2Y02, N655_2Y02 and name OD1, water_655_2Y02 and name O

color water_cyan, H3A_wat_*
set dash_color, water_cyan, H3A_wat_*

label N655_4LDO and name CA, "N6.55"

# Set view for H3A (same as H2)
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

print("Storing scene H3B...")
disable all
enable 2RH1_gpcr
enable 3NY9_gpcr
enable 2VT4_gpcr

show cartoon, 2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr

# Show N6.55 as sticks
show sticks, inactive_N655
color highlight_yellow, inactive_N655
set stick_radius, 0.25, inactive_N655

# Show all waters to demonstrate ABSENCE at 6.55
show spheres, (2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and resn HOH
color gray50, (2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and resn HOH
set sphere_scale, 0.25, resn HOH

# Show ligands
show sticks, (2RH1_gpcr or 3NY9_gpcr or 2VT4_gpcr) and organic
color inverse_magenta, (2RH1_gpcr or 3NY9_gpcr) and organic
color gray70, 2VT4_gpcr and organic

label N655_2RH1 and name CA, "N6.55"

# Set view for H3B (same as H3A)
set_view (\
     0.918646157,    0.090081483,    0.384683996,\
     0.076466568,   -0.996494412,    0.036925454,\
     0.379639924,   -0.002264185,   -0.924785256,\
     0.000000000,    0.000000000, -147.662918091,\
   -31.105007172,   17.511806488,   11.639530182,\
    67.662979126,  227.662918091,  -20.000000000 )

scene H3B, store, Inactive state: NO water at N6.55

# -----------------------------------------------------------------------------
# H4A: Inverse agonist (2RH1) - SHORT D2.50-W6.48 distance (6.64A)
# Shows water network between D2.50 and W6.48 in cyan
# -----------------------------------------------------------------------------

print("Storing scene H4A...")
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
# Key waters: 548 (near D2.50), 534 (near W6.48), 528, 532
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

# Labels
label D250_2RH1 and name CA, "D2.50"
label W648_2RH1 and name CA, "W6.48"

# Set view for H4A (D2.50-W6.48 view)
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

print("Storing scene H4B...")
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
# Key waters: 2002 (near D2.50), 2017 (near W6.48), 2003, 2016
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

# Labels
label D250_2Y02 and name CA, "D2.50"
label W648_2Y02_h4 and name CA, "W6.48"

# Set view for H4B (same as H4A)
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

# Hide all selections from menu
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

print("")
print("=" * 60)
print("HYPOTHESIS VISUALIZATION SCENES READY")
print("=" * 60)
print("")
print("Available scenes:")
print("  H1A - Agonists close to S5.43 (2.82-3.25A)")
print("  H1B - Inverse agonists far from S5.43 (3.32-3.58A)")
print("")
print("  H2A - Inverse agonists close to W6.48 (3.42-3.44A)")
print("  H2B - Agonists far/no contact with W6.48")
print("")
print("  H3A - Agonist: water at N6.55 (active marker)")
print("  H3B - Inactive: NO water at N6.55")
print("")
print("  H4A - Inverse agonist (2RH1): SHORT D2.50-W6.48 (6.64A) + water network")
print("  H4B - Agonist (2Y02): LONGER D2.50-W6.48 (6.95A) + water network")
print("")
print("Usage: scene <name>, recall")
print("Example: scene H1A, recall")
print("=" * 60)
