# PyMOL Color Definitions from colorscales.yaml
# Source this file at the start of any PyMOL script for consistent colors
#
# Usage: @colorscales_pymol.pml
#
# Generated from thesis/colorscales.yaml

# =============================================================================
# Structure Type Colors
# =============================================================================
set_color type1_blue, [0.122, 0.467, 0.706]    # Type I (Microbial) - #1f77b4
set_color type2_red, [0.839, 0.153, 0.157]     # Type II (Animal/GPCR) - #d62728

# =============================================================================
# Individual Structure Colors
# =============================================================================
# Type I (Microbial Rhodopsins)
set_color color_1C3W, [0.122, 0.467, 0.706]    # Bacteriorhodopsin - Blue
set_color color_3UG9, [0.682, 0.780, 0.910]    # C1C2 Channelrhodopsin - Light Blue

# Type II (Animal Opsins / GPCRs)
set_color color_1U19, [0.839, 0.153, 0.157]    # Bovine Rhodopsin - Red
set_color color_2Z73, [1.0, 0.596, 0.588]      # Squid Rhodopsin - Light Red/Salmon

# Non-opsin GPCR
set_color color_2RH1, [0.173, 0.627, 0.173]    # Beta2 Adrenergic Receptor - Green

# =============================================================================
# Ligand Colors
# =============================================================================
set_color retinal_orange, [1.0, 0.498, 0.055]  # Retinal (RET) - #ff7f0e
set_color carazolol_purple, [0.580, 0.404, 0.741]  # Carazolol (CAU) - #9467bd

# =============================================================================
# Functional Element Colors
# =============================================================================
set_color grn_conserved, [0.737, 0.741, 0.133]   # Conserved x.50 positions - Yellow-green
set_color binding_pocket, [0.498, 0.498, 0.498]  # Binding pocket residues - Gray
set_color schiff_base, [0.549, 0.337, 0.294]     # Schiff base lysine - Brown
set_color counterion_color, [0.890, 0.467, 0.761]  # Counterion (Glu) - Pink

# =============================================================================
# Default Display Settings (ALWAYS white background)
# =============================================================================
bg_color white
set bg_rgb, [1.0, 1.0, 1.0]
set ray_shadows, 0
set antialias, 2
set orthoscopic, 1
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set stick_radius, 0.15
set sphere_scale, 0.25
set label_size, 14
set label_color, black

print("Colorscales loaded from colorscales_pymol.pml")
