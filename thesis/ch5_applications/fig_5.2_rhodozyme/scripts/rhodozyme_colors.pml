# Rhodozyme PyMOL Colors + Style
# ===========================================================================
# Derived from rhodozyme_colors.yaml — keep in sync.
# Consistent with fig 4.2 spectral TM helix colors and render settings.
#
# Convention: gray = experimental / crystal structure, color = highlight / model
#
# Source at the top of any Ch5 .pml script:  @rhodozyme_colors.pml
# ===========================================================================

# --- Backbone / experimental (gray hues) ------------------------------------
set_color backbone,       [0.678, 0.710, 0.741]  # #adb5bd  light gray
set_color backbone_dark,  [0.424, 0.459, 0.490]  # #6c757d  mid gray
set_color backbone_faint, [0.914, 0.926, 0.937]  # #e9ecef  very faint

# --- Designed backbone (model output) ---------------------------------------
set_color designed,       [0.757, 0.400, 0.420]  # #c1666b  terracotta
set_color designed_light, [0.894, 0.757, 0.757]  # #e4c1c1  terracotta light

# --- TM helix spectral colors (from fig 4.2) --------------------------------
set_color tm1_color, [0.361, 0.353, 0.478]  # #5c5a7a  400nm
set_color tm2_color, [0.239, 0.353, 0.502]  # #3d5a80  450nm
set_color tm3_color, [0.271, 0.482, 0.420]  # #457b6b  490nm
set_color tm4_color, [0.416, 0.600, 0.306]  # #6a994e  530nm
set_color tm5_color, [0.831, 0.627, 0.235]  # #d4a03c  570nm
set_color tm6_color, [0.737, 0.424, 0.145]  # #bc6c25  580nm
set_color tm7_color, [0.757, 0.400, 0.420]  # #c1666b  620nm
set_color loop_color, [0.678, 0.710, 0.741] # #adb5bd  same as backbone

# --- Functional roles -------------------------------------------------------
set_color retinal,    [0.737, 0.424, 0.145]  # #bc6c25  rust (type II opsin)
set_color theozyme,   [0.416, 0.600, 0.306]  # #6a994e  green (530nm / conserved)
set_color substrate,  [0.831, 0.627, 0.235]  # #d4a03c  ochre (570nm)
set_color candidate,  [0.271, 0.482, 0.420]  # #457b6b  sage (490nm)

# --- Geometry / annotation --------------------------------------------------
set_color geo,    [0.204, 0.227, 0.251]  # #343a40  dark gray
set_color hbond,  [0.424, 0.459, 0.490]  # #6c757d  mid gray

# =============================================================================
# Render settings (matching fig 4.2)
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
set label_color, geo

print("Rhodozyme colors + style loaded (rhodozyme_colors.pml)")
