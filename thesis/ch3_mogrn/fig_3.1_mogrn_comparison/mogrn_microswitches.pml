# Figure 3.1b — MO-GRN Functional Sets on Bacteriorhodopsin (1C3W)
# ===========================================================================
# All 20 key MO-GRN positions from Fig 3.1a, organized into 6 functional
# sets and rendered as separate panels. No labels. Retinal shown as context.
#
# Residue mapping (MO-GRN → 1C3W residue):
#   TM2: 2.46→I45   2.47→T46
#   TM3: 3.42→R82   3.45→D85  3.46→W86  3.49→T89  3.50→T90
#         3.53→L93   3.56→D96
#   TM4: 4.47→D115  4.51→I119  4.54→G122
#   TM5: 5.44→F135  5.47→W138
#   TM6: 6.50→W182  6.53→Y185  6.54→P186
#   TM7: 7.46→D212  7.49→A215  7.50→K216
#
# Sets:
#   1. Schiff base & counterions  (3.42, 3.45, 7.46, 7.50)
#   2. Proton pathway / TM3 motif (3.45, 3.49, 3.53, 3.56)
#   3. DC gate                    (3.50, 4.47)
#   4. Binding pocket             (3.46, 5.44, 5.47, 6.50, 6.53)
#   5. Spectral tuning switches   (3.53, 4.51, 4.54, 6.54, 7.49)
#   6. Channel gate               (2.46, 2.47)
#
# Usage:
#   pymol mogrn_microswitches.pml
# ===========================================================================

# Load shared colorscales (colors + render settings)
@thesis/shared/colorscales_pymol.pml

# --- TM helix colors (spectral scheme, matching Fig 4.2) ---
set_color tm2_color, [0.239, 0.353, 0.502]
set_color tm3_color, [0.271, 0.482, 0.420]
set_color tm4_color, [0.416, 0.600, 0.306]
set_color tm5_color, [0.831, 0.627, 0.235]
set_color tm6_color, [0.737, 0.424, 0.145]
set_color tm7_color, [0.757, 0.400, 0.420]



# =============================================================================
# Load structure
# =============================================================================
fetch 1C3W, br, type=pdb, async=0
remove solvent
remove br and not chain A


# =============================================================================
# Backbone: gray, transparent
# =============================================================================
hide everything
show cartoon, br and chain A and polymer
color backbone, br
set cartoon_transparency, 0.75, br


# =============================================================================
# Retinal chromophore (always visible)
# =============================================================================
select retinal, br and chain A and resn RET
show sticks, retinal
color retinal_rust, retinal
set stick_radius, 0.20, retinal


# =============================================================================
# Define ALL residue selections (by MO-GRN position)
# =============================================================================

# TM2
select grn_2_46, br and chain A and resi 45
select grn_2_47, br and chain A and resi 46

# TM3
select grn_3_42, br and chain A and resi 82
select grn_3_45, br and chain A and resi 85
select grn_3_46, br and chain A and resi 86
select grn_3_49, br and chain A and resi 89
select grn_3_50, br and chain A and resi 90
select grn_3_53, br and chain A and resi 93
select grn_3_56, br and chain A and resi 96

# TM4
select grn_4_47, br and chain A and resi 115
select grn_4_51, br and chain A and resi 119
select grn_4_54, br and chain A and resi 122

# TM5
select grn_5_44, br and chain A and resi 135
select grn_5_47, br and chain A and resi 138

# TM6
select grn_6_50, br and chain A and resi 182
select grn_6_53, br and chain A and resi 185
select grn_6_54, br and chain A and resi 186

# TM7
select grn_7_46, br and chain A and resi 212
select grn_7_49, br and chain A and resi 215
select grn_7_50, br and chain A and resi 216


# =============================================================================
# Define functional sets
# =============================================================================
select set1_schiff,  grn_3_42 or grn_3_45 or grn_7_46 or grn_7_50
select set2_proton,  grn_3_45 or grn_3_49 or grn_3_53 or grn_3_56
select set3_dcgate,  grn_3_50 or grn_4_47
select set4_pocket,  grn_3_46 or grn_5_44 or grn_5_47 or grn_6_50 or grn_6_53
select set5_tuning,  grn_3_53 or grn_4_51 or grn_4_54 or grn_6_54 or grn_7_49
select set6_channel, grn_2_46 or grn_2_47

# All 20 positions combined
select all_grn, set1_schiff or set2_proton or set3_dcgate or set4_pocket or set5_tuning or set6_channel

# Start with all sticks hidden
hide sticks, all_grn
deselect


# =============================================================================
# Camera — consistent view for all panels
# =============================================================================
set_view (\
    -0.460567087,   -0.024940038,   -0.887273133,\
     0.881179094,    0.107397772,   -0.460422933,\
     0.106775478,   -0.993902326,   -0.027486838,\
     0.000029527,   -0.000008181, -171.372528076,\
    16.382932663,   35.955245972,    5.651223660,\
   148.000717163,  194.747482300,   20.000000000 )


# =============================================================================
# SET 1 — Schiff base & counterions (R82, D85, D212, K216)
# =============================================================================
color backbone, all_grn
show sticks, set1_schiff
color tm3_color, grn_3_42 or grn_3_45
color tm7_color, grn_7_46 or grn_7_50
set stick_radius, 0.15, set1_schiff

# Dashed lines: counterion–Schiff base interactions
distance d_D85_K216, (resi 85 and name OD1), (resi 216 and name NZ), cutoff=5
distance d_D212_K216, (resi 212 and name OD1), (resi 216 and name NZ), cutoff=5
set dash_color, dark_gray, d_D85_K216
set dash_color, dark_gray, d_D212_K216
set dash_gap, 0.25, d_D85_K216
set dash_gap, 0.25, d_D212_K216
set dash_width, 1.5, d_D85_K216
set dash_width, 1.5, d_D212_K216
hide labels, d_D85_K216
hide labels, d_D212_K216

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_set1_schiff_base.png, dpi=300

hide sticks, set1_schiff
hide dashes, d_D85_K216
hide dashes, d_D212_K216


# =============================================================================
# SET 2 — Proton pathway / TM3 motif (D85, T89, L93, D96)
# =============================================================================
color backbone, all_grn
show sticks, set2_proton
color tm3_color, set2_proton
set stick_radius, 0.15, set2_proton

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_set2_proton_pathway.png, dpi=300

hide sticks, set2_proton


# =============================================================================
# SET 3 — DC gate (T90, D115)
# =============================================================================
color backbone, all_grn
show sticks, set3_dcgate
color tm3_color, grn_3_50
color tm4_color, grn_4_47
set stick_radius, 0.15, set3_dcgate

# Dashed line: cross-helix gate interaction
distance d_dcgate, (resi 90 and name OG1), (resi 115 and name OD1), cutoff=6
set dash_color, dark_gray, d_dcgate
set dash_gap, 0.25, d_dcgate
set dash_width, 1.5, d_dcgate
hide labels, d_dcgate

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_set3_dc_gate.png, dpi=300

hide sticks, set3_dcgate
hide dashes, d_dcgate


# =============================================================================
# SET 4 — Binding pocket (W86, F135, W138, W182, Y185)
# =============================================================================
color backbone, all_grn
show sticks, set4_pocket
color tm3_color, grn_3_46
color tm5_color, grn_5_44 or grn_5_47
color tm6_color, grn_6_50 or grn_6_53
set stick_radius, 0.15, set4_pocket

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_set4_binding_pocket.png, dpi=300

hide sticks, set4_pocket


# =============================================================================
# SET 5 — Spectral tuning switches (L93, I119, G122, P186, A215)
# =============================================================================
color backbone, all_grn
show sticks, set5_tuning
color tm3_color, grn_3_53
color tm4_color, grn_4_51 or grn_4_54
color tm6_color, grn_6_54
color tm7_color, grn_7_49
set stick_radius, 0.15, set5_tuning

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_set5_spectral_tuning.png, dpi=300

hide sticks, set5_tuning


# =============================================================================
# SET 6 — Channel gate (I45, T46)
# =============================================================================
color backbone, all_grn
show sticks, set6_channel
color tm2_color, set6_channel
set stick_radius, 0.15, set6_channel

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_set6_channel_gate.png, dpi=300

hide sticks, set6_channel


# =============================================================================
# ALL SETS combined — overview render
# =============================================================================
show sticks, all_grn
color tm2_color, grn_2_46 or grn_2_47
color tm3_color, grn_3_42 or grn_3_45 or grn_3_46 or grn_3_49 or grn_3_50 or grn_3_53 or grn_3_56
color tm4_color, grn_4_47 or grn_4_51 or grn_4_54
color tm5_color, grn_5_44 or grn_5_47
color tm6_color, grn_6_50 or grn_6_53 or grn_6_54
color tm7_color, grn_7_46 or grn_7_49 or grn_7_50
set stick_radius, 0.15, all_grn

ray 2400, 1800
png thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_all_positions.png, dpi=300


# =============================================================================
# Save session
# =============================================================================
deselect
save thesis/ch3_mogrn/fig_3.1_mogrn_comparison/fig_3.1b_mogrn_microswitches.pse

print("")
print("=" * 60)
print("Figure 3.1b — MO-GRN Functional Sets (6 panels + overview)")
print("=" * 60)
print("  Set 1: Schiff base & counterions  (R82, D85, D212, K216)")
print("  Set 2: Proton pathway / TM3 motif (D85, T89, L93, D96)")
print("  Set 3: DC gate                    (T90, D115)")
print("  Set 4: Binding pocket             (W86, F135, W138, W182, Y185)")
print("  Set 5: Spectral tuning switches   (L93, I119, G122, P186, A215)")
print("  Set 6: Channel gate               (I45, T46)")
print("=" * 60)
