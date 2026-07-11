# PyMOL script: Spectral tuning sites in animal opsins
# Focus on GRN positions involved in wavelength regulation

# Load shared colorscales (colors + render settings)
@thesis/shared/colorscales_pymol.pml

# Load structures
load 1u19_grn_aligned.cif, bovine
load 2z73_grn_aligned.cif, squid

# Setup
hide everything
show cartoon, all
set cartoon_transparency, 0.7
color color_1U19, bovine
color color_2Z73, squid

# Retinal - the chromophore
show sticks, resn RET
color retinal_rust, resn RET
set stick_radius, 0.2


# 3.28: Counterion
show sticks, bovine and resi 113
color color_1U19, bovine and resi 113
show sticks, squid and resi 111
color color_2Z73, squid and resi 111

# 3.32: Tuning site
show sticks, bovine and resi 117
color color_1U19, bovine and resi 117
show sticks, squid and resi 115
color color_2Z73, squid and resi 115

# 6.48: Rotamer switch
show sticks, bovine and resi 265
color color_1U19, bovine and resi 265
show sticks, squid and resi 274
color color_2Z73, squid and resi 274

# 7.42: Schiff base K
show sticks, bovine and resi 296
color color_1U19, bovine and resi 296
show sticks, squid and resi 305
color color_2Z73, squid and resi 305

set stick_radius, 0.15

# Center on retinal
center resn RET
zoom resn RET, 12

# Ray trace
# ray 1600, 1200
# png grn_spectral_tuning.png, dpi=300
