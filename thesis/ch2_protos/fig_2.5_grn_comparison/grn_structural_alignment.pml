# PyMOL Script: GRN-Based Structural Comparison of Animal Opsins
# ============================================================
#
# Compares bovine rhodopsin (1U19) and squid rhodopsin (2Z73) using
# Generic Residue Numbers (GRN) to identify structurally equivalent positions.
#
# GRN enables systematic comparison of binding pockets across GPCRs
# by defining structurally equivalent positions in the 7TM bundle.
#
# KEY BINDING POCKET POSITIONS (spectral tuning sites):
#   3.28: Counterion position - Bovine E113, Squid Y111
#   3.32: Spectral tuning site - Bovine A117, Squid G115
#   6.48: Rotamer switch - Bovine W265, Squid W274
#   7.42: Schiff base lysine - Bovine K296, Squid K305
#
# NOTE: Type I (microbial) opsins do NOT have GRN - this is the MOGRN gap
# that LAMBDA addresses using graph-based representations.
#
# Usage: Open PyMOL and run: @grn_structural_alignment.pml

# Load colorscales
@colorscales_pymol.pml

# =============================================================================
# Fetch Structures from PDB
# =============================================================================
print("Fetching opsin structures from PDB...")

fetch 1U19, async=0
fetch 2Z73, async=0

# =============================================================================
# Extract Chain A from each structure
# =============================================================================
print("Extracting chain A...")

create bovine, 1U19 and chain A
create squid, 2Z73 and chain A

delete 1U19
delete 2Z73

# =============================================================================
# Align Structures on Conserved x.50 Positions
# =============================================================================
print("Aligning structures on conserved GRN x.50 positions...")

# x.50 positions (conserved anchors in each helix):
#   1.50: Bovine N55, Squid N44
#   2.50: Bovine D83, Squid D72
#   3.50: Bovine R135, Squid R125
#   4.50: Bovine W161, Squid W152
#   5.50: Bovine P215, Squid P204
#   6.50: Bovine P267, Squid P268
#   7.50: Bovine P303, Squid P304

# Create selections for x.50 positions (CA atoms only)
select x50_bovine, bovine and name CA and (resi 55 or resi 83 or resi 135 or resi 161 or resi 215 or resi 267 or resi 303)
select x50_squid, squid and name CA and (resi 44 or resi 72 or resi 125 or resi 152 or resi 204 or resi 268 or resi 304)

# Align squid to bovine on x.50 positions
align squid and (resi 44 or resi 72 or resi 125 or resi 152 or resi 204 or resi 268 or resi 304) and name CA, bovine and (resi 55 or resi 83 or resi 135 or resi 161 or resi 215 or resi 267 or resi 303) and name CA

# =============================================================================
# Display Settings
# =============================================================================
print("Setting up display...")

bg_color white
hide everything
show cartoon, all

# Color by structure
color color_1U19, bovine
color color_2Z73, squid

# =============================================================================
# Retinal Visualization
# =============================================================================
print("Highlighting retinal chromophore...")

select all_retinal, resn RET
show sticks, all_retinal
color retinal_orange, all_retinal
set stick_radius, 0.2

# =============================================================================
# Binding Pocket GRN Positions
# =============================================================================
print("Highlighting GRN binding pocket positions...")

# 3.28: Counterion
select bp_3_28_bovine, bovine and resi 113
select bp_3_28_squid, squid and resi 111

# 3.32: Spectral tuning site
select bp_3_32_bovine, bovine and resi 117
select bp_3_32_squid, squid and resi 115

# 6.48: Rotamer switch (toggle switch)
select bp_6_48_bovine, bovine and resi 265
select bp_6_48_squid, squid and resi 274

# 7.42: Schiff base lysine
select bp_7_42_bovine, bovine and resi 296
select bp_7_42_squid, squid and resi 305

# Combined binding pocket selection
select all_bp, bp_3_28_* or bp_3_32_* or bp_6_48_* or bp_7_42_*

# Show binding pocket as sticks
show sticks, all_bp
set stick_radius, 0.15

# Highlight key functional residues with spheres
show spheres, bp_7_42_*
show spheres, bp_3_28_*
set sphere_scale, 0.25

# =============================================================================
# Labels
# =============================================================================
set label_size, 14
set label_color, black
set label_position, [2, 2, 2]

label bovine and resi 113 and name CA, "3.28"
label bovine and resi 117 and name CA, "3.32"
label bovine and resi 265 and name CA, "6.48"
label bovine and resi 296 and name CA, "7.42"

# =============================================================================
# Camera View
# =============================================================================
print("Setting camera view...")

orient all_bp
zoom all_bp, 10

# =============================================================================
# Summary Output
# =============================================================================
print("")
print("=" * 60)
print("GRN-BASED OPSIN STRUCTURAL ALIGNMENT")
print("=" * 60)
print("")
print("Structures aligned:")
print("  Bovine rhodopsin (1U19) - Reference")
print("  Squid rhodopsin (2Z73) - Mobile")
print("")
print("BINDING POCKET POSITIONS (GRN):")
print("  3.28 (Counterion):    Bovine=E113, Squid=Y111")
print("  3.32 (Tuning site):   Bovine=A117, Squid=G115")
print("  6.48 (Rotamer switch): Bovine=W265, Squid=W274")
print("  7.42 (Schiff base K): Bovine=K296, Squid=K305")
print("")
print("KEY INSIGHT: GRN enables systematic comparison of binding pockets")
print("Same GRN position = Same structural role in spectral tuning")
print("")
print("NOTE: Type I opsins (microbial) do NOT have GRN - MOGRN gap!")
print("=" * 60)

# =============================================================================
# Export Commands (uncomment to use)
# =============================================================================
# ray 2400, 1800
# png grn_structural_alignment.png, dpi=300
