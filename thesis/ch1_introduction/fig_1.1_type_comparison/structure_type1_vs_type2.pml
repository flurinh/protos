# Figure 1.1 — Type I vs Type II Opsin Structural Comparison
# ===========================================================================
# Two structures, four panels:
#
#   A: 1C3W  Bacteriorhodopsin  (Type I)  — full structure
#   B: 1U19  Bovine rhodopsin   (Type II) — full structure
#   C: 1C3W  Bacteriorhodopsin  (Type I)  — binding pocket close-up
#   D: 1U19  Bovine rhodopsin   (Type II) — binding pocket close-up
#
# Full structure panels (A, B):
#   Gray transparent cartoon, retinal (rust), Schiff base lysine (green).
#   No other sidechains shown.
#
# Binding pocket panels (C, D):
#   Pocket residues within 6A of retinal as gray sticks.
#   Retinal (rust), Schiff base lysine (green). No cartoon.
#
# Style matches ch5 figures (gray backbone, rust retinal, green Schiff base).
#
# Usage:
#   pymol structure_type1_vs_type2.pml
#   Then: scene panel_A, recall  (B, C, D)
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

# Cartoon: gray, transparent
show cartoon, bR and polymer
color backbone, bR
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
# Panel A_nolig — 1C3W full structure, no ligand/lysine
# =============================================================================
hide everything

show cartoon, bR and polymer
color backbone, bR
set cartoon_transparency, 0.7, bR

set_view (\
     0.521108687,   -0.214197397,    0.826174915,\
    -0.851958215,   -0.072513409,    0.518566728,\
    -0.051166043,   -0.974096894,   -0.220273525,\
    -0.000103682,   -0.000014109, -168.123947144,\
    18.453006744,   40.861843109,    6.949202061,\
   135.654693604,  200.583694458,   20.000000000 )

scene panel_A_nolig, store, message=1C3W no ligand


# =============================================================================
# Panel B — 1U19 full structure
# =============================================================================
hide everything

show cartoon, bovine and polymer
color backbone, bovine
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
# Panel B_nolig — 1U19 full structure, no ligand/lysine
# =============================================================================
hide everything

show cartoon, bovine and polymer
color backbone, bovine
set cartoon_transparency, 0.7, bovine

set_view (\
     0.778251946,   -0.559981763,    0.284156501,\
    -0.496124566,   -0.825713933,   -0.268439978,\
     0.384953201,    0.067935660,   -0.920435607,\
    -0.000135182,    0.000050765, -214.546783447,\
    21.781805038,   48.353328705,    7.465656757,\
   183.837738037,  245.230148315,   20.000000000 )

scene panel_B_nolig, store, message=1U19 no ligand


# =============================================================================
# Panel C — 1C3W binding pocket
# =============================================================================
hide everything

# Pocket: gray sticks (no cartoon)
show sticks, pocket_bR
color backbone, pocket_bR

# Retinal: rust
show sticks, ret_bR
color retinal_rust, ret_bR
set stick_radius, 0.25, ret_bR

# Schiff base lysine: green
color schiff_base, schiff_bR

# Same view as Panel A
set_view (\
     0.521108687,   -0.214197397,    0.826174915,\
    -0.851958215,   -0.072513409,    0.518566728,\
    -0.051166043,   -0.974096894,   -0.220273525,\
    -0.000103682,   -0.000014109, -168.123947144,\
    18.453006744,   40.861843109,    6.949202061,\
   135.654693604,  200.583694458,   20.000000000 )

scene panel_C, store, message=1C3W binding pocket


# =============================================================================
# Panel D — 1U19 binding pocket
# =============================================================================
hide everything

# Pocket: gray sticks (no cartoon)
show sticks, pocket_bovine
color backbone, pocket_bovine

# Retinal: rust
show sticks, ret_bovine
color retinal_rust, ret_bovine
set stick_radius, 0.25, ret_bovine

# Schiff base lysine: green
color schiff_base, schiff_bovine

# Same view as Panel B
set_view (\
     0.778251946,   -0.559981763,    0.284156501,\
    -0.496124566,   -0.825713933,   -0.268439978,\
     0.384953201,    0.067935660,   -0.920435607,\
    -0.000135182,    0.000050765, -214.546783447,\
    21.781805038,   48.353328705,    7.465656757,\
   183.837738037,  245.230148315,   20.000000000 )

scene panel_D, store, message=1U19 binding pocket


# =============================================================================
# Interactive navigation
# =============================================================================
scene panel_A, recall
deselect

print("")
print("=" * 60)
print("Figure 1.1 — Type I vs Type II Opsin Comparison")
print("=" * 60)
print("")
print("Panels:")
print("  A: 1C3W  Bacteriorhodopsin  (Type I)  — full structure")
print("  B: 1U19  Bovine rhodopsin   (Type II) — full structure")
print("  C: 1C3W  Bacteriorhodopsin  (Type I)  — binding pocket")
print("  D: 1U19  Bovine rhodopsin   (Type II) — binding pocket")
print("")
print("Color key:")
print("  Gray   (backbone)     = cartoon / pocket sticks")
print("  Rust   (retinal_rust) = retinal chromophore")
print("  Green  (schiff_base)  = Schiff base lysine")
print("")
print("Navigate: scene panel_A, recall  (B, C, D)")
print("Get view: get_view")
print("=" * 60)


# =============================================================================
# Export
# =============================================================================
scene panel_A, recall
ray 2400, 1800
png thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1a_bR_structure.png, dpi=300

scene panel_A_nolig, recall
ray 2400, 1800
png thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1a_bR_structure_nolig.png, dpi=300

scene panel_B, recall
ray 2400, 1800
png thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1b_bovine_structure.png, dpi=300

scene panel_B_nolig, recall
ray 2400, 1800
png thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1b_bovine_structure_nolig.png, dpi=300

scene panel_C, recall
ray 2400, 1800
png thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1c_bR_pocket.png, dpi=300

scene panel_D, recall
ray 2400, 1800
png thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1d_bovine_pocket.png, dpi=300

save thesis/ch1_introduction/fig_1.1_type_comparison/fig_1.1_type_comparison.pse
