# PyMOL Script: Type I vs Type II Opsin Structural Comparison
# ============================================================
#
# This script visualizes the structural alignment of Type I (microbial)
# and Type II (animal) opsins around their retinal binding sites.
#
# KEY INSIGHT: Despite completely different 7TM topologies, both opsin types
# bind the same chromophore (retinal) - this is convergent evolution at work.
#
# Structures:
#   Type I (Microbial):
#     - 1C3W: Bacteriorhodopsin (proton pump)
#     - 3UG9: C1C2 Channelrhodopsin (cation channel)
#   Type II (Animal/GPCR):
#     - 1U19: Bovine rhodopsin (dim light vision)
#     - 2Z73: Squid rhodopsin (vision)
#
# Usage: Open in PyMOL and run: @structure_type1_vs_type2.pml

# Load colorscales
@colorscales_pymol.pml

# Ensure white background
bg_color white

# =============================================================================
# Fetch Structures from PDB
# =============================================================================
print("Fetching structures from PDB...")

# Type I (Microbial)
fetch 1C3W, async=0
fetch 3UG9, async=0

# Type II (Animal)
fetch 1U19, async=0
fetch 2Z73, async=0

# =============================================================================
# Extract Chain A from each structure (single chain + retinal)
# =============================================================================
print("Extracting chain A from each structure...")

# Create selections for chain A only
create bR, 1C3W and chain A
create ChR, 3UG9 and chain A
create bovine, 1U19 and chain A
create squid, 2Z73 and chain A

# Delete original multi-chain structures
delete 1C3W
delete 3UG9
delete 1U19
delete 2Z73

# =============================================================================
# Align Type I structures to bacteriorhodopsin (1C3W)
# =============================================================================
print("Aligning Type I structures...")
align ChR, bR

# =============================================================================
# Align Type II structures to bovine rhodopsin
# =============================================================================
print("Aligning Type II structures...")
align squid, bovine

# =============================================================================
# Align Type II to Type I on retinal position
# =============================================================================
print("Aligning Type II to Type I on retinal...")

# Select retinal in each structure
select ret_bR, bR and resn RET
select ret_bovine, bovine and resn RET

# Align bovine to bR using retinal
align bovine and resn RET, bR and resn RET
align squid and resn RET, bR and resn RET
align ChR and resn RET, bR and resn RET

# =============================================================================
# Display Settings
# =============================================================================
print("Setting up display...")

# Hide everything first
hide everything

# Show cartoon for protein backbone
show cartoon, all

# Make cartoon transparent so retinal is visible
set cartoon_transparency, 0.7

# Color by structure type
color color_1C3W, bR
color color_3UG9, ChR
color color_1U19, bovine
color color_2Z73, squid

# =============================================================================
# Retinal Visualization
# =============================================================================
print("Highlighting retinal chromophore...")

# Select all retinals
select all_retinal, resn RET

# Show retinal as prominent sticks (thicker, fully opaque)
show sticks, all_retinal
color retinal_orange, all_retinal
set stick_radius, 0.25, all_retinal
set stick_transparency, 0, all_retinal

# =============================================================================
# Binding Pocket Visualization (7A around retinal)
# =============================================================================
print("Highlighting binding pockets...")

# Select binding pocket residues (within 7A of retinal)
select pocket_bR, bR and byres (bR within 7 of (bR and resn RET)) and not resn RET
select pocket_ChR, ChR and byres (ChR within 7 of (ChR and resn RET)) and not resn RET
select pocket_bovine, bovine and byres (bovine within 7 of (bovine and resn RET)) and not resn RET
select pocket_squid, squid and byres (squid within 7 of (squid and resn RET)) and not resn RET

# Create combined selections by type
select pocket_type1, pocket_bR or pocket_ChR
select pocket_type2, pocket_bovine or pocket_squid
select all_pockets, pocket_type1 or pocket_type2

# Show binding pocket residues as lines
show lines, all_pockets
set line_width, 2

# =============================================================================
# Labels and Annotations
# =============================================================================
# Create pseudoatom labels for structure identification
pseudoatom label_bR, pos=[0, 0, 30], label="1C3W: Bacteriorhodopsin (Type I)"
pseudoatom label_ChR, pos=[0, 0, 25], label="3UG9: Channelrhodopsin (Type I)"
pseudoatom label_bovine, pos=[0, 0, 20], label="1U19: Bovine rhodopsin (Type II)"
pseudoatom label_squid, pos=[0, 0, 15], label="2Z73: Squid rhodopsin (Type II)"

color color_1C3W, label_bR
color color_3UG9, label_ChR
color color_1U19, label_bovine
color color_2Z73, label_squid

# =============================================================================
# Camera and Final View
# =============================================================================
print("Setting camera view...")

# Center on retinal
center all_retinal
zoom all_retinal, 15

# Set a good viewing angle
set_view (\
     0.866025388,    0.000000000,    0.500000000,\
     0.250000000,    0.866025388,   -0.433012694,\
    -0.433012694,    0.500000000,    0.750000000,\
     0.000000000,    0.000000000, -100.000000000,\
     0.000000000,    0.000000000,    0.000000000,\
    80.000000000,  120.000000000,  -20.000000000 )

# =============================================================================
# Output Information
# =============================================================================
print("")
print("=" * 60)
print("TYPE I vs TYPE II OPSIN COMPARISON")
print("=" * 60)
print("")
print("Structures loaded:")
print("  Type I (Microbial) - Blue tones:")
print("    - bR (1C3W): Bacteriorhodopsin")
print("    - ChR (3UG9): Channelrhodopsin")
print("  Type II (Animal/GPCR) - Red tones:")
print("    - bovine (1U19): Bovine rhodopsin")
print("    - squid (2Z73): Squid rhodopsin")
print("")
print("KEY INSIGHT:")
print("  Despite completely different 7TM folds,")
print("  both opsin types bind retinal in a similar pocket.")
print("  This enables LAMBDA's fold-agnostic analysis.")
print("")
print("To export figure:")
print("  ray 2400, 1800")
print("  png structure_type1_vs_type2.png, dpi=300")
print("=" * 60)

# =============================================================================
# Alternate Views (uncomment to use)
# =============================================================================
# View 1: Type I only
# hide everything
# show cartoon, bR or ChR
# show sticks, (bR or ChR) and resn RET
# color color_1C3W, bR
# color color_3UG9, ChR

# View 2: Type II only
# hide everything
# show cartoon, bovine or squid
# show sticks, (bovine or squid) and resn RET
# color color_1U19, bovine
# color color_2Z73, squid

# View 3: Side by side comparison
# set grid_mode, 1
# set grid_slot, 1, bR
# set grid_slot, 2, bovine
