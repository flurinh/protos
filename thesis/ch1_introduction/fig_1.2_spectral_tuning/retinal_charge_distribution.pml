# Figure 1.2a — Retinal in Vacuum (No Sidechain Interactions)
# ===========================================================================
# The protonated Schiff base retinal chromophore from bovine rhodopsin (1U19)
# shown with a Gaussian electron cloud, colored on the same Δq scale as the
# pocket potential figure. Without binding pocket sidechains, Δq = 0
# everywhere → the surface is uniformly white.
#
# Paired with retinal_pocket_potential.pml to show that sidechains create
# the red/blue charge redistribution pattern.
#
# Usage:
#   pymol retinal_charge_distribution.pml
# ===========================================================================

@thesis/shared/colorscales_pymol.pml

# =============================================================================
# Load and prepare
# =============================================================================
# Full rhodopsin for pocket context
fetch 1U19, rhodopsin, type=pdb, async=0
remove solvent
remove rhodopsin and not chain A

# Use the same retinal coordinates as the pocket potential figure
load thesis/ch1_introduction/fig_1.2_spectral_tuning/1u19_retinal_potential.pdb, ret_vacuum

# =============================================================================
# Selections
# =============================================================================
select retinal, ret_vacuum and resn RET
select lys_nz, ret_vacuum and resn LYS and name NZ
select chromophore, retinal or lys_nz
select conjugated, (ret_vacuum and resn RET and not name C1+C2+C3+C4+C16+C17+C18+C19+C20) or lys_nz

# All GRN microswitch positions (same as pocket potential figure)
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
# Gaussian electron cloud — uniform white (Δq = 0, no sidechain interactions)
# Same DX map pipeline as pocket potential figure for consistency.
# =============================================================================
python
import numpy as np
from pymol import cmd, stored

# ---- Collect atom positions, set Δq = 0 (vacuum reference) ----
stored.data = []
cmd.iterate_state(1, "conjugated", "stored.data.append((x, y, z))")
pos = np.array(stored.data)

# All B-factors zero: no sidechain interaction
bfac = np.zeros(len(pos))

print(f"Conjugated atoms: {len(pos)}")
print(f"Δq = 0 everywhere (vacuum reference)")

# ---- Grid parameters (identical to pocket potential figure) ----
SIGMA = 1.2
SPACING = 0.4
BUFFER = 5.0

mins = pos.min(axis=0) - BUFFER
maxs = pos.max(axis=0) + BUFFER
ns = (np.ceil((maxs - mins) / SPACING) + 1).astype(int)

xs = mins[0] + np.arange(ns[0]) * SPACING
ys = mins[1] + np.arange(ns[1]) * SPACING
zs = mins[2] + np.arange(ns[2]) * SPACING

# ---- Compute Gaussian density field ----
density = np.zeros(ns, dtype=np.float64)
X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')

for i in range(len(pos)):
    d2 = (X - pos[i, 0])**2 + (Y - pos[i, 1])**2 + (Z - pos[i, 2])**2
    g = np.exp(-d2 / (2.0 * SIGMA**2))
    density += g

# B-factor field is all zeros (no sidechain interaction)
bfac_field = np.zeros_like(density)

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
density_path = f"{dx_dir}/density_map_vacuum.dx"
bfac_path = f"{dx_dir}/bfac_map_vacuum.dx"

write_dx(density_path, density, mins, SPACING)
write_dx(bfac_path, bfac_field, mins, SPACING)

# ---- Load maps ----
cmd.load(density_path, "density_map")
cmd.load(bfac_path, "bfac_map")

# ---- Isosurface from density map (same level as pocket figure) ----
cmd.isosurface("ret_surf", "density_map", level=1.0)
cmd.set("transparency", 0.35, "ret_surf")

# ---- Color surface: same ramp as pocket figure (±1.0 centi-electron) ----
CLIM = 1.0
cmd.ramp_new("bfac_ramp", "bfac_map", [-CLIM, 0.0, CLIM],
             ["color_1U19", "white", "color_1C3W"])
cmd.set("surface_color", "bfac_ramp", "ret_surf")
cmd.color("bfac_ramp", "ret_surf")

print(f"Color ramp: -{CLIM} (red) -> 0 (white) -> +{CLIM} (blue)")
print(f"Surface is uniformly white (Δq = 0)")
python end

# =============================================================================
# Pocket residues: thin gray sticks (same set as pocket potential figure)
# =============================================================================
show sticks, pocket_key
set stick_radius, 0.10, pocket_key
color backbone, pocket_key

# =============================================================================
# Camera and render (same framing as pocket potential figure)
# =============================================================================
orient chromophore
zoom chromophore, 4

deselect

ray 2400, 1200
png thesis/ch1_introduction/fig_1.2_spectral_tuning/fig_1.2_retinal_vacuum.png, dpi=300

save thesis/ch1_introduction/fig_1.2_spectral_tuning/fig_1.2_retinal_vacuum.pse
