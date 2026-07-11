# Figure 2.3 — Structure Processor: Type I vs Type II Alignment on Retinal
# ===========================================================================
# Three panels comparing bacteriorhodopsin and bovine rhodopsin:
#
#   A: 1C3W  Bacteriorhodopsin  (Type I)  — full structure
#   B: 1U19  Bovine rhodopsin   (Type II) — full structure
#   C: Both overlaid, aligned on retinal  — fold comparison
#
# All panels: slate (Type I) vs terracotta (Type II), rust retinal, green K.
#
# KEY INSIGHT: Different 7TM topologies converge on the same retinal pocket.
# The Structure Processor enables this comparison through retinal-based alignment.
#
# Usage:
#   pymol render_aligned.pml
#   Then: scene panel_A, recall  (B, C)
# ===========================================================================

# Load shared colorscales (colors + render settings)
@thesis/shared/colorscales_pymol.pml


# =============================================================================
# Fetch structures
# =============================================================================
print("Fetching structures...")

fetch 1C3W, name=bR_raw, type=pdb, async=0
fetch 1U19, name=bovine_raw, type=pdb, async=0

# Extract chain A only
create bR, bR_raw and chain A
create bovine, bovine_raw and chain A

delete bR_raw
delete bovine_raw

remove resn HOH


# =============================================================================
# Align on retinal (common reference frame)
# =============================================================================
print("Aligning on retinal...")

align bovine and resn RET, bR and resn RET


# =============================================================================
# Selections
# =============================================================================
print("Building selections...")

# Retinal
select ret_bR, bR and resn RET
select ret_bovine, bovine and resn RET

# Schiff base lysines
select schiff_bR, bR and chain A and resi 216
select schiff_bovine, bovine and chain A and resi 296

# Binding pocket residues (within 6A of retinal, protein only)
select pocket_bR, bR and byres (bR within 6 of (bR and resn RET)) and not resn RET and polymer
select pocket_bovine, bovine and byres (bovine within 6 of (bovine and resn RET)) and not resn RET and polymer


# =============================================================================
# Panel A — 1C3W full structure
# =============================================================================
hide everything

# Cartoon: slate (Type I), transparent
show cartoon, bR and polymer
color color_1C3W, bR
set cartoon_transparency, 0.7, bR

# Retinal: rust, prominent
show sticks, ret_bR
color retinal_rust, ret_bR
set stick_radius, 0.25, ret_bR

# Schiff base lysine: green sticks
show sticks, schiff_bR
color schiff_base, schiff_bR

set_view (\
     0.521108687,   -0.214197397,    0.826174915,\
    -0.851958215,   -0.072513409,    0.518566728,\
    -0.051166043,   -0.974096894,   -0.220273525,\
    -0.000103682,   -0.000014109, -168.123947144,\
    18.453006744,   40.861843109,    6.949202061,\
   135.654693604,  200.583694458,   20.000000000 )

scene panel_A, store, message=1C3W full structure


# =============================================================================
# Panel B — 1U19 full structure
# =============================================================================
hide everything

# Cartoon: terracotta (Type II), transparent
show cartoon, bovine and polymer
color color_1U19, bovine
set cartoon_transparency, 0.7, bovine

show sticks, ret_bovine
color retinal_rust, ret_bovine
set stick_radius, 0.25, ret_bovine

show sticks, schiff_bovine
color schiff_base, schiff_bovine

set_view (\
     0.778251946,   -0.559981763,    0.284156501,\
    -0.496124566,   -0.825713933,   -0.268439978,\
     0.384953201,    0.067935660,   -0.920435607,\
    -0.000135182,    0.000050765, -214.546783447,\
    21.781805038,   48.353328705,    7.465656757,\
   183.837738037,  245.230148315,   20.000000000 )

scene panel_B, store, message=1U19 full structure


# =============================================================================
# Panel C — Both overlaid, aligned on retinal
# =============================================================================
hide everything

# bR: slate cartoon
show cartoon, bR and polymer
color color_1C3W, bR
set cartoon_transparency, 0.5, bR

# Bovine: terracotta cartoon
show cartoon, bovine and polymer
color color_1U19, bovine
set cartoon_transparency, 0.5, bovine

# Both retinals: rust
show sticks, ret_bR or ret_bovine
color retinal_rust, ret_bR or ret_bovine
set stick_radius, 0.25, ret_bR or ret_bovine

# Both Schiff base lysines: green
show sticks, schiff_bR or schiff_bovine
color schiff_base, schiff_bR or schiff_bovine

set_view (\
    -0.027786685,   -0.487263709,    0.872813642,\
    -0.886746287,   -0.391030431,   -0.246529534,\
     0.461420357,   -0.780814409,   -0.421211869,\
     0.000011012,   -0.000062842, -177.921432495,\
    18.612159729,   45.512092590,    7.222129822,\
   145.457748413,  210.386764526,   20.000000000 )

scene panel_C, store, message=Overlay aligned on retinal


# =============================================================================
# Interactive navigation
# =============================================================================
scene panel_A, recall
deselect

print("")
print("=" * 60)
print("Figure 2.3 — Structure Processor: Retinal Alignment")
print("=" * 60)
print("")
print("Panels:")
print("  A: 1C3W  Bacteriorhodopsin  (Type I)  — full structure")
print("  B: 1U19  Bovine rhodopsin   (Type II) — full structure")
print("  C: Both overlaid, aligned on retinal")
print("")
print("Color key:")
print("  Slate      (color_1C3W) = Type I  (bR)")
print("  Terracotta (color_1U19) = Type II (bovine)")
print("  Rust   (retinal_rust)   = retinal chromophore")
print("  Green  (schiff_base)    = Schiff base lysine")
print("")
print("Navigate: scene panel_A, recall  (B, C)")
print("Get view: get_view")
print("=" * 60)


# =============================================================================
# Export
# =============================================================================
scene panel_A, recall
ray 2400, 1800
png thesis/ch2_protos/fig_2.3_br_binding_pocket/fig_2.3a_bR_structure.png, dpi=300

scene panel_B, recall
ray 2400, 1800
png thesis/ch2_protos/fig_2.3_br_binding_pocket/fig_2.3b_bovine_structure.png, dpi=300

scene panel_C, recall
ray 2400, 1800
png thesis/ch2_protos/fig_2.3_br_binding_pocket/fig_2.3c_overlay.png, dpi=300

save thesis/ch2_protos/fig_2.3_br_binding_pocket/fig_2.3_alignment.pse
