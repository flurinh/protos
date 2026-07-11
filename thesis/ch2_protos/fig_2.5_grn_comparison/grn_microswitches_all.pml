# Figure 2.5c — GRN Microswitches on Bovine Rhodopsin (1U19)
# ===========================================================================
# All 18 GRN microswitch positions from the Fig 2.5a alignment table,
# mapped onto bovine rhodopsin. Sidechains colored by TM helix (spectral
# color scheme matching Fig 4.2).
#
# GRN→residue mapping from type_ii_grns.csv (P02699):
#
#   1.50→N55   2.50→D83   3.28→E113  3.32→A117  3.40→L125
#   3.49→E134  3.50→R135  3.51→Y136  4.50→W161  5.50→P215
#   6.44→F261  6.47→C264  6.48→W265  6.50→P267  7.43→K296
#   7.49→N302  7.50→P303  7.53→Y306
#
# Usage:
#   pymol grn_microswitches_all.pml
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
# GRN microswitches — PIF + CWxP, colored by TM helix
# =============================================================================

# PIF connector: 3.40 (L125, TM3), 5.50 (P215, TM5), 6.44 (F261, TM6)
select grn_3_40, rhodopsin and chain A and resi 125
select grn_5_50, rhodopsin and chain A and resi 215
select grn_6_44, rhodopsin and chain A and resi 261
select pif, grn_3_40 or grn_5_50 or grn_6_44
show sticks, pif
color tm3_color, grn_3_40
color tm5_color, grn_5_50
color tm6_color, grn_6_44

# CWxP toggle: 6.48 (W265), 6.50 (P267) — TM6
select grn_6_48, rhodopsin and chain A and resi 265
select grn_6_50, rhodopsin and chain A and resi 267
select cwxp, grn_6_48 or grn_6_50
show sticks, cwxp
color tm6_color, cwxp

# Combined selection
select all_grn, pif or cwxp
set stick_radius, 0.15, all_grn


# =============================================================================
# Labels (GRN numbers)
# =============================================================================
set label_position, [2.5, 1.5, 1.5]
set label_size, 14
set label_outline_color, white

label resi 125 and name CA, "3.40"
label resi 215 and name CA, "5.50"
label resi 261 and name CA, "6.44"
label resi 265 and name CA, "6.48"
label resi 267 and name CA, "6.50"


# =============================================================================
# Camera and render
# =============================================================================
deselect

set_view (\
     0.713131785,    0.668824077,    0.210029781,\
     0.503620923,   -0.280373693,   -0.817162871,\
    -0.487653702,    0.688527882,   -0.536779761,\
    -0.000121266,   -0.000010204, -224.691619873,\
    37.545406342,   12.098978043,   11.371767044,\
   189.672027588,  259.708953857,   20.000000000 )

print("")
print("=" * 60)
print("Figure 2.5c — PIF + CWxP Microswitches on Bovine Rhodopsin")
print("=" * 60)
print("PIF + CWxP microswitches, colored by TM helix:")
print("")
print("  PIF:  L125 (3.40, TM3), P215 (5.50, TM5), F261 (6.44, TM6)")
print("  CWxP: W265 (6.48), P267 (6.50) — TM6")
print("")
print("Get view: get_view")
print("Render:   ray 2400, 1800")
print("=" * 60)

# Render with labels
ray 2400, 1800
png thesis/ch2_protos/fig_2.5_grn_comparison/fig_2.5c_grn_microswitches.png, dpi=300

# Render without labels
hide labels
ray 2400, 1800
png thesis/ch2_protos/fig_2.5_grn_comparison/fig_2.5c_grn_microswitches_nolabel.png, dpi=300

save thesis/ch2_protos/fig_2.5_grn_comparison/fig_2.5c_grn_microswitches.pse
