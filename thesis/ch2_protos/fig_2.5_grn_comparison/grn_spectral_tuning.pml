# PyMOL Script: Spectral Tuning Sites in Animal Opsins (GRN-Based)
# ==================================================================
#
# Focuses on GRN positions involved in wavelength regulation
# with transparent cartoon to emphasize retinal and binding pocket.
#
# Spectral tuning in rhodopsins is controlled by:
#   - Electrostatic environment around the Schiff base
#   - Counterion distance and character
#   - Retinal-protein interactions in the binding pocket
#
# KEY GRN POSITIONS:
#   3.28: Counterion - E113 (bovine) vs Y111 (squid)
#   3.32: Tuning site - A117 (bovine) vs G115 (squid)
#   6.48: Rotamer switch - W265 (bovine) vs W274 (squid)
#   7.42: Schiff base K - K296 (bovine) vs K305 (squid)
#
# Usage: Open PyMOL and run: @grn_spectral_tuning.pml

# Load colorscales
@colorscales_pymol.pml

# =============================================================================
# Fetch Structures from PDB
# =============================================================================
print("Fetching opsin structures from PDB...")

fetch 1U19, async=0
fetch 2Z73, async=0

# =============================================================================
# Extract Chain A
# =============================================================================
print("Extracting chain A...")

create bovine, 1U19 and chain A
create squid, 2Z73 and chain A

delete 1U19
delete 2Z73

# =============================================================================
# Align on x.50 Conserved Positions
# =============================================================================
print("Aligning on GRN x.50 positions...")

# Align squid to bovine on x.50 positions
align squid and (resi 44 or resi 72 or resi 125 or resi 152 or resi 204 or resi 268 or resi 304) and name CA, bovine and (resi 55 or resi 83 or resi 135 or resi 161 or resi 215 or resi 267 or resi 303) and name CA

# =============================================================================
# Display Settings - Transparent for Retinal Focus
# =============================================================================
print("Setting up transparent display...")

bg_color white
hide everything
show cartoon, all

# Make cartoon transparent so retinal and binding pocket are prominent
set cartoon_transparency, 0.75

# Color by structure
color color_1U19, bovine
color color_2Z73, squid

# =============================================================================
# Retinal - The Chromophore (Prominent Display)
# =============================================================================
print("Highlighting retinal chromophore...")

select all_retinal, resn RET
show sticks, all_retinal
color retinal_orange, all_retinal
set stick_radius, 0.25
set stick_transparency, 0, all_retinal

# =============================================================================
# Spectral Tuning Residues (GRN Positions)
# =============================================================================
print("Highlighting spectral tuning GRN positions...")

# 3.28: Counterion - critical for Schiff base pKa
show sticks, bovine and resi 113
color color_1U19, bovine and resi 113
show sticks, squid and resi 111
color color_2Z73, squid and resi 111

# 3.32: Spectral tuning site - modulates counterion
show sticks, bovine and resi 117
color color_1U19, bovine and resi 117
show sticks, squid and resi 115
color color_2Z73, squid and resi 115

# 6.48: Rotamer switch (toggle switch) - activation
show sticks, bovine and resi 265
color color_1U19, bovine and resi 265
show sticks, squid and resi 274
color color_2Z73, squid and resi 274

# 7.42: Schiff base lysine - covalent retinal attachment
show sticks, bovine and resi 296
color color_1U19, bovine and resi 296
show sticks, squid and resi 305
color color_2Z73, squid and resi 305

set stick_radius, 0.18

# =============================================================================
# Labels for GRN Positions
# =============================================================================
set label_size, 12
set label_color, black
set label_position, [1.5, 1.5, 1.5]

label bovine and resi 113 and name CA, "3.28"
label bovine and resi 296 and name CA, "7.42"

# =============================================================================
# Camera View - Focus on Retinal
# =============================================================================
print("Setting camera view on retinal...")

center all_retinal
zoom all_retinal, 15

# =============================================================================
# Summary Output
# =============================================================================
print("")
print("=" * 60)
print("SPECTRAL TUNING IN ANIMAL OPSINS")
print("=" * 60)
print("")
print("SPECTRAL TUNING POSITIONS (GRN):")
print("  3.28 (Counterion):     E113 (bovine) vs Y111 (squid)")
print("  3.32 (Tuning site):    A117 (bovine) vs G115 (squid)")
print("  6.48 (Rotamer switch): W265 (bovine) vs W274 (squid)")
print("  7.42 (Schiff base K):  K296 (bovine) vs K305 (squid)")
print("")
print("The counterion at 3.28 stabilizes the protonated Schiff base.")
print("Changes at these positions shift the absorption spectrum.")
print("=" * 60)

# =============================================================================
# Export Commands (uncomment to use)
# =============================================================================
# ray 1600, 1200
# png grn_spectral_tuning.png, dpi=300
