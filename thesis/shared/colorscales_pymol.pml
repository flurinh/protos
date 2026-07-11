# Thesis PyMOL Color Definitions
# ===========================================================================
# All colors from colorscales.yaml — use in any chapter's PyMOL scripts.
#
# Source at the top of any .pml script:  @colorscales_pymol.pml
#   (or use full path: @thesis/shared/colorscales_pymol.pml)
#
# Naming convention:
#   color_<PDB>      — per-structure color
#   <category>_<name> — semantic name
# ===========================================================================

# --- Structures: Type I (Microbial) — Slate palette -------------------------
set_color color_1C3W,   [0.239, 0.353, 0.502]  # #3d5a80  Slate
set_color color_3UG9,   [0.596, 0.757, 0.851]  # #98c1d9  Slate Light
set_color type_i,       [0.239, 0.353, 0.502]  # #3d5a80  alias
set_color type_i_light, [0.596, 0.757, 0.851]  # #98c1d9  alias

# --- Structures: Type II (Animal) — Terracotta palette ----------------------
set_color color_1U19,     [0.757, 0.400, 0.420]  # #c1666b  Terracotta
set_color color_2Z73,     [0.894, 0.757, 0.757]  # #e4c1c1  Terracotta Light
set_color color_3PQR,     [0.757, 0.400, 0.420]  # #c1666b  (same as 1U19)
set_color type_ii,        [0.757, 0.400, 0.420]  # #c1666b  alias
set_color type_ii_light,  [0.894, 0.757, 0.757]  # #e4c1c1  alias

# --- Molecules: Ligands -----------------------------------------------------
set_color retinal_rust,   [0.737, 0.424, 0.145]  # #bc6c25  Retinal
set_color retinal_light,  [0.867, 0.722, 0.573]  # #ddb892  Retinal Light
set_color carazolol_mauve,[0.478, 0.349, 0.502]  # #7a5980  Carazolol
set_color schiff_base,    [0.416, 0.600, 0.306]  # #6a994e  Schiff Base N

# --- Functional elements ----------------------------------------------------
set_color conserved_green,[0.416, 0.600, 0.306]  # #6a994e  Conserved positions

# --- Neutral / backbone (gray palette) --------------------------------------
set_color backbone,       [0.678, 0.710, 0.741]  # #adb5bd  Light gray
set_color backbone_dark,  [0.424, 0.459, 0.490]  # #6c757d  Mid gray
set_color backbone_faint, [0.914, 0.926, 0.937]  # #e9ecef  Very faint
set_color dark_gray,      [0.204, 0.227, 0.251]  # #343a40  Annotation/label

# --- Sequences / Embeddings (Sage/Ochre) ------------------------------------
set_color sage,           [0.271, 0.482, 0.420]  # #457b6b  Short Wave
set_color ochre,          [0.831, 0.627, 0.235]  # #d4a03c  Long Wave

# =============================================================================
# Render settings (consistent across all chapters)
# =============================================================================
bg_color white
set orthoscopic, 1

# Lighting
set ray_shadow, 0
set ray_opaque_background, 1
set antialias, 2
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set ray_texture, 0
set direct, 0.7
set ambient, 0.3
set specular, 0.2

# Sticks
set stick_radius, 0.15
set stick_ball, on
set stick_ball_ratio, 1.5

# Cartoon
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_tube_radius, 0.2

# Spheres
set sphere_scale, 0.25

# Labels
set label_size, 14
set label_font_id, 7
set label_color, dark_gray

print("Thesis colorscales loaded (colorscales_pymol.pml)")
