# Figure 2.5b — x.50 Anchor Positions on Bovine Rhodopsin (1U19)
# ===========================================================================
# The 7 x.50 reference positions, one per TM helix, colored by helix
# (spectral color scheme matching Fig 4.2).
#
# Usage:
#   pymol grn_microswitches.pml
# ===========================================================================

# Load shared colorscales (colors + render settings)
@thesis/shared/colorscales_pymol.pml

# --- TM helix colors (spectral scheme, matching Fig 4.2) ---
set_color tm1_color, [0.361, 0.353, 0.478]
set_color tm2_color, [0.239, 0.353, 0.502]
set_color tm3_color, [0.271, 0.482, 0.420]
set_color tm4_color, [0.416, 0.600, 0.306]
set_color tm5_color, [0.831, 0.627, 0.235]
set_color tm6_color, [0.737, 0.424, 0.145]
set_color tm7_color, [0.757, 0.400, 0.420]


# =============================================================================
# Load structure
# =============================================================================
fetch 1U19, rhodopsin, type=pdb, async=0
remove solvent
remove rhodopsin and not chain A


# =============================================================================
# Backbone: gray, transparent
# =============================================================================
hide everything
show cartoon, rhodopsin and chain A and polymer
color backbone, rhodopsin
set cartoon_transparency, 0.75, rhodopsin



# =============================================================================
# GRN .50 positions — colored by TM helix
# =============================================================================

# 1.50 — N55 (TM1)
select grn_1_50, rhodopsin and chain A and resi 55
show sticks, grn_1_50
color tm1_color, grn_1_50

# 2.50 — D83 (TM2)
select grn_2_50, rhodopsin and chain A and resi 83
show sticks, grn_2_50
color tm2_color, grn_2_50

# 3.50 — R135 (TM3)
select grn_3_50, rhodopsin and chain A and resi 135
show sticks, grn_3_50
color tm3_color, grn_3_50

# 4.50 — W161 (TM4)
select grn_4_50, rhodopsin and chain A and resi 161
show sticks, grn_4_50
color tm4_color, grn_4_50

# 5.50 — P215 (TM5)
select grn_5_50, rhodopsin and chain A and resi 215
show sticks, grn_5_50
color tm5_color, grn_5_50

# 6.50 — P267 (TM6)
select grn_6_50, rhodopsin and chain A and resi 267
show sticks, grn_6_50
color tm6_color, grn_6_50

# 7.50 — P303 (TM7)
select grn_7_50, rhodopsin and chain A and resi 303
show sticks, grn_7_50
color tm7_color, grn_7_50

# Combined selection
select all_grn, grn_1_50 or grn_2_50 or grn_3_50 or grn_4_50 or grn_5_50 or grn_6_50 or grn_7_50
set stick_radius, 0.15, all_grn


# =============================================================================
# Labels
# =============================================================================
set label_position, [2.5, 1.5, 1.5]
set label_size, 14
set label_outline_color, white

label resi 55 and name CA, "1.50"
label resi 83 and name CA, "2.50"
label resi 135 and name CA, "3.50"
label resi 161 and name CA, "4.50"
label resi 215 and name CA, "5.50"
label resi 267 and name CA, "6.50"
label resi 303 and name CA, "7.50"


# =============================================================================
# Camera and render
# =============================================================================
deselect

set_view (\
     0.776293337,    0.612253726,    0.150038302,\
     0.445581734,   -0.364597648,   -0.817633450,\
    -0.445895404,    0.701580584,   -0.555844843,\
    -0.000121266,   -0.000010204, -234.715469360,\
    37.545406342,   12.098978043,   11.371767044,\
   199.695877075,  269.732818604,   20.000000000 )

print("")
print("=" * 60)
print("Figure 2.5b — x.50 Anchor Positions on Bovine Rhodopsin")
print("=" * 60)
print("7 x.50 positions, colored by TM helix:")
print("")
print("  TM1: N55  (1.50)")
print("  TM2: D83  (2.50)")
print("  TM3: R135 (3.50)")
print("  TM4: W161 (4.50)")
print("  TM5: P215 (5.50)")
print("  TM6: P267 (6.50)")
print("  TM7: P303 (7.50)")
print("")
print("Get view: get_view")
print("Render:   ray 2400, 1800")
print("=" * 60)

# Render with labels
ray 2400, 1800
png thesis/ch2_protos/fig_2.5_grn_comparison/fig_2.5b_grn_x50.png, dpi=300

# Render without labels
hide labels
ray 2400, 1800
png thesis/ch2_protos/fig_2.5_grn_comparison/fig_2.5b_grn_x50_nolabel.png, dpi=300

save thesis/ch2_protos/fig_2.5_grn_comparison/fig_2.5b_grn_x50.pse
