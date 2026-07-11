# Figure 5.1 — Input Structures and the Design Premise
# ===========================================================================
# Panel A: 1U19 (dark) vs 3PQR (active) — TM5/6 displacement highlighted
# Panel B: 2AGE trypsin — catalytic triad with H-bond distances
#
# All experimental structures in gray. Only highlights use color.
# Style matches fig 4.2 (fancy helices, stick_ball, ray_trace_mode 1).
#
# Run: pymol -cq fig_5_1_input_structures.pml
# ===========================================================================

# Load color scheme + render settings
@rhodozyme_colors.pml

# =============================================================================
# Fetch structures
# =============================================================================
fetch 1U19, name=dark, type=pdb
fetch 3PQR, name=active, type=pdb
fetch 2AGE, name=trypsin, type=pdb

# Clean: keep chain A for rhodopsins, chains X+A for trypsin (X=protein, A=substrate)
remove dark and not chain A
remove active and not chain A
remove trypsin and not (chain X or chain A)
remove resn HOH

# =============================================================================
# Panel A — Dark vs Active rhodopsin
# =============================================================================
# Both structures are experimental → gray backbone.
# TM5 and TM6 are highlighted with their spectral colors to show displacement.

# --- Align dark onto active (Cα superposition) ---
align dark and name CA, active and name CA

# --- Dark state: mid-gray, slightly transparent ---
hide everything, dark
show cartoon, dark
color backbone_dark, dark
set cartoon_transparency, 0.5, dark

# --- Active state: light gray (slightly brighter than dark) ---
hide everything, active
show cartoon, active
color backbone, active

# --- Highlight TM5 in both states (ochre, spectral 570nm) ---
# TM5 approx residues 200-225
color tm5_color, dark and resi 200-225
color tm5_color, active and resi 200-225

# --- Highlight TM6 in both states (rust, spectral 580nm) ---
# TM6 approx residues 247-277
color tm6_color, dark and resi 247-277
color tm6_color, active and resi 247-277

# --- Retinal (active state, sticks) ---
show sticks, active and resn RET+LYR
set stick_radius, 0.25, active and resn RET+LYR
color retinal, active and resn RET+LYR
set specular, 0.8, active and resn RET+LYR
set spec_power, 200, active and resn RET+LYR

# --- Retinal (dark state, faint sticks) ---
show sticks, dark and resn RET+LYR
set stick_radius, 0.25, dark and resn RET+LYR
color retinal, dark and resn RET+LYR
set stick_transparency, 0.4, dark and resn RET+LYR



# --- Store Panel A ---
set_view (\
    -0.400470138,   -0.914588273,   -0.056143947,\
     0.907670140,   -0.404340893,    0.112441264,\
    -0.125537932,   -0.005929092,    0.992068827,\
     0.000143170,   -0.000104382, -229.661300659,\
   -37.423156738,  -11.309287071,   39.157096863,\
   179.299621582,  280.018035889,  -20.000000000 )
scene panel_A, store, message=Dark vs Active rhodopsin

# =============================================================================
# Panel B — Trypsin catalytic triad
# =============================================================================

# Hide rhodopsins
hide everything, dark
hide everything, active
# --- Trypsin backbone: gray, transparent (context) ---
show cartoon, trypsin
color backbone, trypsin
set cartoon_transparency, 0.7, trypsin

# --- Catalytic triad: theozyme green, sticks + sidechains ---
# 2AGE: chain X = trypsin protein, chain A = succinyl-AAPR substrate
select triad, trypsin and resi 57+102+195 and chain X
show sticks, triad
color theozyme, triad

# --- Substrate (chain A = succinyl-AAPR peptide) ---
show sticks, trypsin and chain A
color substrate, trypsin and chain A

# --- H-bond distances (sidechain atoms) ---
# Ser195 OG → His57 NE2
distance ser_his, trypsin and resi 195 and name OG and chain X, trypsin and resi 57 and name NE2 and chain X
set dash_color, red, ser_his
set dash_width, 2.5, ser_his
set dash_gap, 0.35, ser_his

# His57 ND1 → Asp102 OD2
distance his_asp, trypsin and resi 57 and name ND1 and chain X, trypsin and resi 102 and name OD2 and chain X
set dash_color, red, his_asp
set dash_width, 2.5, his_asp
set dash_gap, 0.35, his_asp

# Hide distance number labels
hide labels, ser_his
hide labels, his_asp

# --- Store Panel B ---
set_view (\
    -0.555999279,    0.488205671,   -0.672696590,\
    -0.716420591,   -0.691831648,    0.090046816,\
    -0.421430856,    0.531998515,    0.734418094,\
     0.000274185,    0.000007391, -122.816200256,\
    19.811824799,   69.337989807,   10.551140785,\
    26.688549042,  218.920486450,  -20.000000000 )
scene panel_B, store, message=Trypsin catalytic triad

# =============================================================================
# Export
# =============================================================================
# Panel A
scene panel_A, recall
ray 2400, 1800
png ../figures/fig_5_1a_dark_vs_active.png, dpi=300

# Panel B
scene panel_B, recall
ray 2400, 1800
png ../figures/fig_5_1b_trypsin_triad.png, dpi=300

# Save session
save ../figures/fig_5_1_input_structures.pse

print("")
print("Figure 5.1 complete.")
print("  Panel A: fig_5_1a_dark_vs_active.png")
print("  Panel B: fig_5_1b_trypsin_triad.png")
print("  Session: fig_5_1_input_structures.pse")
print("")
print("Color key:")
print("  Gray (backbone/backbone_dark) = experimental crystal structure")
print("  Ochre (TM5) = TM5 helix, spectral 570nm")
print("  Rust  (TM6) = TM6 helix, spectral 580nm")
print("  Rust  (retinal) = retinal cofactor")
print("  Green (theozyme) = catalytic triad")
print("  Ochre (substrate) = succinyl-AAPR")
