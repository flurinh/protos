# =============================================================================
# THESIS PyMOL Settings - Scientific Elegance Palette
# =============================================================================
#
# Source this file at the start of any PyMOL script for consistent colors
# and rendering settings across all thesis figures.
#
# Usage: @thesis_pymol_settings.pml
#
# Reference: colorscales.yaml (Scientific Elegance palette)
#
# =============================================================================

# =============================================================================
# STRUCTURE COLORS - Slate (Type I) / Terracotta (Type II)
# =============================================================================

# Type I (Microbial Opsins) - Slate blue family
set_color type1_slate, [0.239, 0.353, 0.502]          # #3d5a80 - Primary
set_color type1_light, [0.596, 0.757, 0.851]          # #98c1d9 - Light variant

# Type II (Animal Opsins) - Terracotta family
set_color type2_terracotta, [0.757, 0.400, 0.420]     # #c1666b - Primary
set_color type2_light, [0.894, 0.757, 0.757]          # #e4c1c1 - Light variant

# =============================================================================
# INDIVIDUAL STRUCTURE COLORS
# =============================================================================

# Type I (Microbial Rhodopsins)
set_color color_1C3W, [0.239, 0.353, 0.502]           # Bacteriorhodopsin - Slate
set_color color_3UG9, [0.596, 0.757, 0.851]           # Channelrhodopsin - Light slate

# Type II (Animal Opsins)
set_color color_1U19, [0.757, 0.400, 0.420]           # Bovine Rhodopsin - Terracotta
set_color color_2Z73, [0.894, 0.757, 0.757]           # Squid Rhodopsin - Light terracotta

# =============================================================================
# MOLECULE COLORS - Rust (Retinal) / Mauve (Other ligands)
# =============================================================================

# Retinal states
set_color retinal_rust, [0.737, 0.424, 0.145]         # #bc6c25 - 11-cis (dark state)
set_color retinal_light, [0.867, 0.722, 0.573]        # #ddb892 - All-trans (activated)
set_color retinal_deprot, [0.478, 0.349, 0.502]       # #7a5980 - Deprotonated

# Other ligands
set_color carazolol_mauve, [0.478, 0.349, 0.502]      # #7a5980
set_color carazolol_light, [0.769, 0.702, 0.780]      # #c4b3c7

# =============================================================================
# FUNCTIONAL ELEMENT COLORS
# =============================================================================

# Conserved/functional positions - Olive green
set_color conserved_olive, [0.416, 0.600, 0.306]      # #6a994e
set_color conserved_light, [0.655, 0.769, 0.580]      # #a7c494

# Schiff base nitrogen (same as conserved - key functional)
set_color schiff_base, [0.416, 0.600, 0.306]          # #6a994e

# =============================================================================
# SEQUENCE COLORS - Sage (SW) / Ochre (LW)
# =============================================================================

set_color short_wave_sage, [0.271, 0.482, 0.420]      # #457b6b
set_color short_wave_light, [0.639, 0.769, 0.722]     # #a3c4b8
set_color long_wave_ochre, [0.831, 0.627, 0.235]      # #d4a03c
set_color long_wave_light, [0.929, 0.851, 0.659]      # #edd9a8

# =============================================================================
# NEUTRAL COLORS
# =============================================================================

set_color neutral_gray, [0.424, 0.459, 0.490]         # #6c757d
set_color neutral_light, [0.678, 0.710, 0.741]        # #adb5bd
set_color neutral_dark, [0.204, 0.227, 0.251]         # #343a40

# =============================================================================
# GLOBAL DISPLAY SETTINGS
# =============================================================================

# Background - ALWAYS white for publication
bg_color white
set bg_rgb, [1.0, 1.0, 1.0]

# Rendering quality
set antialias, 2
set ray_shadows, 0
set ray_opaque_background, 1
set specular, 0.15
set spec_reflect, 0.2
set shininess, 20

# Orthographic projection (no perspective distortion)
set orthoscopic, 1
set depth_cue, 0
set fog, 0

# =============================================================================
# CARTOON SETTINGS
# =============================================================================

set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set cartoon_side_chain_helper, 0
set cartoon_transparency, 0.0

# Helix/sheet appearance
set cartoon_oval_length, 1.2
set cartoon_oval_width, 0.25
set cartoon_loop_radius, 0.2
set cartoon_rect_length, 1.4
set cartoon_rect_width, 0.3

# =============================================================================
# STICK/SPHERE SETTINGS
# =============================================================================

set stick_radius, 0.15
set stick_transparency, 0.0
set stick_ball, 0
set stick_ball_ratio, 1.5

set sphere_scale, 0.25
set sphere_transparency, 0.0

# =============================================================================
# LINE SETTINGS
# =============================================================================

set line_width, 2
set line_smooth, 1

# =============================================================================
# LABEL SETTINGS
# =============================================================================

set label_size, 14
set label_color, neutral_dark
set label_font_id, 7
set label_outline_color, white
set label_position, [0, 0, 2]

# =============================================================================
# RAY TRACING PRESETS
# =============================================================================

# For quick preview
alias ray_preview, ray 800, 600

# For publication (300 DPI at 8 inch width = 2400 px)
alias ray_pub, ray 2400, 1800

# For high-res poster
alias ray_poster, ray 4800, 3600

# =============================================================================
# COMMON SELECTIONS (can be customized per script)
# =============================================================================

# Retinal selection
alias sel_retinal, select retinal, resn RET or resn RTL or resn LYR

# Binding pocket (7A from retinal)
alias sel_pocket, select pocket, byres (all within 7 of retinal) and not retinal

# Schiff base lysine
alias sel_schiff, select schiff_lys, (resn LYS and name NZ) within 5 of retinal

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Show Type I structure style
alias style_type1, color type1_slate; set cartoon_transparency, 0

# Show Type II structure style
alias style_type2, color type2_terracotta; set cartoon_transparency, 0

# Transparent cartoon (for overlays)
alias cartoon_trans, set cartoon_transparency, 0.7

# Highlight retinal
alias show_retinal, show sticks, retinal; color retinal_rust, retinal; set stick_radius, 0.25, retinal

# Highlight binding pocket
alias show_pocket, show lines, pocket; color neutral_gray, pocket

# Standard figure export
alias export_fig, ray_pub; png

# =============================================================================
# CONFIRM LOAD
# =============================================================================

print ""
print "============================================================"
print "THESIS PyMOL Settings Loaded (Scientific Elegance Palette)"
print "============================================================"
print ""
print "Structure Colors:"
print "  Type I:  type1_slate (#3d5a80), type1_light (#98c1d9)"
print "  Type II: type2_terracotta (#c1666b), type2_light (#e4c1c1)"
print ""
print "Molecule Colors:"
print "  Retinal: retinal_rust (#bc6c25), retinal_light (#ddb892)"
print ""
print "Functional:"
print "  Conserved: conserved_olive (#6a994e)"
print ""
print "Quick commands:"
print "  ray_preview  - Quick 800x600 render"
print "  ray_pub      - Publication 2400x1800 render"
print "  show_retinal - Highlight retinal as sticks"
print "  show_pocket  - Show binding pocket residues"
print "============================================================"
print ""
