# Figure 1.2 — Charge Redistribution by Binding Pocket Sidechains
# ===========================================================================
# QM/MM (B3LYP/6-31G* + AMBER point charges) visualization of how
# binding pocket sidechains locally modulate the +1 formal charge
# distributed across the retinal conjugated system.
#
# B-factor = Δq × 100 (H-collapsed centi-electrons):
#   Positive (blue)  = protein pulled e⁻ away → more positive
#   Negative (red)   = protein pushed e⁻ here → more negative
#
# Surface = Gaussian electron cloud on conjugated backbone,
# colored by a spatially-interpolated Δq field (custom DX map).
# Retinal sticks shown in neutral gray.
#
# Prereq: run calculate_pocket_potential.py to generate the modified PDB.
#
# Usage:
#   pymol retinal_pocket_potential.pml
# ===========================================================================

@thesis/shared/colorscales_pymol.pml

# =============================================================================
# Load and prepare
# =============================================================================
fetch 1U19, rhodopsin, type=pdb, async=0
remove solvent
remove rhodopsin and not chain A

# Modified retinal with Δq as B-factor
load thesis/ch1_introduction/fig_1.2_spectral_tuning/1u19_retinal_potential.pdb, ret_potential

# =============================================================================
# Selections
# =============================================================================
select retinal, ret_potential and resn RET
select lys_nz, ret_potential and resn LYS and name NZ
select chromophore, retinal or lys_nz

# Conjugated backbone only (exclude sp3 methyl branches)
select conjugated, (ret_potential and resn RET and not name C1+C2+C3+C4+C16+C17+C18+C19+C20) or lys_nz

# All GRN microswitch positions (from grn_microswitches_all.pml)
# 1.50→N55  2.50→D83  3.28→E113 3.32→A117 3.40→L125
# 3.49→E134 3.50→R135 3.51→Y136 4.50→W161 5.50→P215
# 6.44→F261 6.47→C264 6.48→W265 6.50→P267 7.49→N302 7.50→P303 7.53→Y306
# (K296/7.43 excluded — its NZ is part of the chromophore)
select pocket_key, rhodopsin and chain A and resi 55+83+113+117+125+134+135+136+161+215+261+264+265+267+302+303+306 and not name N+CA+C+O

hide everything

# =============================================================================
# Retinal: neutral gray sticks
# =============================================================================
show sticks, chromophore
set stick_radius, 0.15, chromophore
color backbone_dark, chromophore
set stick_ball, off, lys_nz

# =============================================================================
# Build custom Gaussian maps in Python:
#   1. density_map  — Gaussian shape envelope (for isosurface geometry)
#   2. bfac_map     — Gaussian-smoothed Δq field (for surface coloring)
# =============================================================================
python
import numpy as np
from pymol import cmd, stored

# ---- Collect atom data ----
stored.data = []
cmd.iterate_state(1, "conjugated", "stored.data.append((x, y, z, b))")
atoms = np.array(stored.data)
pos = atoms[:, :3]
bfac = atoms[:, 3]

print(f"Conjugated atoms: {len(pos)}")
print(f"Δq×100 range: {bfac.min():.2f} to {bfac.max():.2f}")

# ---- Grid parameters ----
SIGMA = 1.2       # Gaussian width (Å) — wide enough for continuous surface
SPACING = 0.4     # grid resolution (Å)
BUFFER = 5.0      # padding around atoms (Å)

mins = pos.min(axis=0) - BUFFER
maxs = pos.max(axis=0) + BUFFER
ns = (np.ceil((maxs - mins) / SPACING) + 1).astype(int)

xs = mins[0] + np.arange(ns[0]) * SPACING
ys = mins[1] + np.arange(ns[1]) * SPACING
zs = mins[2] + np.arange(ns[2]) * SPACING

print(f"Grid: {ns[0]}×{ns[1]}×{ns[2]} = {ns.prod()} points")

# ---- Compute fields ----
density = np.zeros(ns, dtype=np.float64)
weighted = np.zeros(ns, dtype=np.float64)

X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

for i in range(len(pos)):
    d2 = (X - pos[i, 0])**2 + (Y - pos[i, 1])**2 + (Z - pos[i, 2])**2
    g = np.exp(-d2 / (2.0 * SIGMA**2))
    density += g
    weighted += bfac[i] * g

# B-factor field = weighted average where density is non-negligible
bfac_field = np.zeros_like(density)
mask = density > 1e-6
bfac_field[mask] = weighted[mask] / density[mask]

print(f"Density range: {density.min():.4f} to {density.max():.4f}")
print(f"B-factor field range: {bfac_field[mask].min():.2f} to {bfac_field[mask].max():.2f}")

# ---- Write OpenDX files ----
def write_dx(path, data, origin, spacing):
    nx, ny, nz = data.shape
    n = nx * ny * nz
    with open(path, 'w') as f:
        f.write(f"object 1 class gridpositions counts {nx} {ny} {nz}\n")
        f.write(f"origin {origin[0]:.6f} {origin[1]:.6f} {origin[2]:.6f}\n")
        f.write(f"delta {spacing:.6f} 0.000000 0.000000\n")
        f.write(f"delta 0.000000 {spacing:.6f} 0.000000\n")
        f.write(f"delta 0.000000 0.000000 {spacing:.6f}\n")
        f.write(f"object 2 class gridconnections counts {nx} {ny} {nz}\n")
        f.write(f"object 3 class array type double rank 0 items {n} data follows\n")
        flat = data.flatten(order='C')
        for idx in range(0, len(flat), 3):
            chunk = flat[idx:idx + 3]
            f.write(" ".join(f"{v:.6e}" for v in chunk) + "\n")
        f.write(f'\nobject "field" class field\n')

dx_dir = "thesis/ch1_introduction/fig_1.2_spectral_tuning"
density_path = f"{dx_dir}/density_map.dx"
bfac_path = f"{dx_dir}/bfac_map.dx"

write_dx(density_path, density, mins, SPACING)
write_dx(bfac_path, bfac_field, mins, SPACING)
print(f"Wrote {density_path}")
print(f"Wrote {bfac_path}")

# ---- Load maps ----
cmd.load(density_path, "density_map")
cmd.load(bfac_path, "bfac_map")

# ---- Isosurface from density map ----
# Level chosen so Gaussians from adjacent bonded atoms (1.4Å apart, σ=1.2)
# overlap into a continuous surface: at midpoint, density ≈ 2×exp(-0.34) ≈ 1.42
cmd.isosurface("ret_surf", "density_map", level=1.0)
cmd.set("transparency", 0.35, "ret_surf")

# ---- Color surface by Δq field ----
CLIM = 1.0  # ±1 centi-electron saturates color
cmd.ramp_new("bfac_ramp", "bfac_map", [-CLIM, 0.0, CLIM],
             ["color_1U19", "white", "color_1C3W"])
cmd.set("surface_color", "bfac_ramp", "ret_surf")
cmd.color("bfac_ramp", "ret_surf")

print(f"Color ramp: -{CLIM} (red) -> 0 (white) -> +{CLIM} (blue)")
python end

# =============================================================================
# Pocket residues: thin gray sticks with labels
# =============================================================================
show sticks, pocket_key
set stick_radius, 0.10, pocket_key
color backbone, pocket_key


# =============================================================================
# Camera and render
# =============================================================================
orient chromophore
zoom chromophore, 4

deselect

ray 2400, 1200
png thesis/ch1_introduction/fig_1.2_spectral_tuning/fig_1.2_pocket_potential.png, dpi=300

save thesis/ch1_introduction/fig_1.2_spectral_tuning/fig_1.2_pocket_potential.pse
