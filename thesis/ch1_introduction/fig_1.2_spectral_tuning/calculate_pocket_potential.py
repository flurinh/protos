#!/usr/bin/env python3
"""QM/MM charge redistribution on retinal by rhodopsin binding pocket.

Computes how individual binding-pocket sidechains locally modulate the
distribution of the +1 formal charge along the retinal conjugated system.

Method:
    1. DFT (B3LYP/6-31G*) on retinal protonated Schiff base in vacuum
    2. DFT with electrostatic embedding: protein heavy atoms as AMBER point charges
    3. Δq = Mulliken charge (embedded) − Mulliken charge (vacuum)
    4. Hydrogen Δq collapsed onto bonded heavy atoms

Output: PDB with retinal B-factor = Δq × 100 (centi-electrons, H-collapsed).
Positive Δq = protein made this atom more positive (pulled e⁻ away).
Negative Δq = protein pushed e⁻ toward this atom.

Expected: NZ/C15 most positive (E113 counterion localizes + charge at SB),
with local modulations from E122, H211, W265 etc. along the chain.

Usage:
    python calculate_pocket_potential.py
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
FIG11_DIR = SCRIPT_DIR.parent / "fig_1.1_type_comparison"
PDB_PATH = FIG11_DIR / "1u19.pdb"
OUTPUT_PDB = SCRIPT_DIR / "1u19_retinal_potential.pdb"

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
CHAIN = "A"
MM_CUTOFF = 15.0     # Å from retinal centroid for MM charges
BASIS = "6-31g*"
XC = "b3lyp"

# ---------------------------------------------------------------------------
# AMBER ff14SB united-atom charges (H collapsed onto bonded heavy atom)
# ---------------------------------------------------------------------------
BACKBONE_CHARGES = {"N": -0.14, "CA": 0.12, "C": 0.60, "O": -0.57}

SIDECHAIN_CHARGES: dict[str, dict[str, float]] = {
    "ALA": {"CB": 0.0},
    "VAL": {"CB": 0.0, "CG1": 0.0, "CG2": 0.0},
    "LEU": {"CB": 0.0, "CG": 0.0, "CD1": 0.0, "CD2": 0.0},
    "ILE": {"CB": 0.0, "CG1": 0.0, "CG2": 0.0, "CD1": 0.0},
    "PRO": {"CB": 0.0, "CG": 0.0, "CD": 0.10},
    "MET": {"CB": 0.0, "CG": 0.0, "SD": -0.274, "CE": 0.034},
    "PHE": {"CB": 0.0, "CG": 0.0, "CD1": 0.0, "CD2": 0.0,
            "CE1": 0.0, "CE2": 0.0, "CZ": 0.0},
    "TYR": {"CB": 0.0, "CG": 0.0, "CD1": 0.0, "CD2": 0.0,
            "CE1": 0.0, "CE2": 0.0, "CZ": 0.265, "OH": -0.558},
    "TRP": {"CB": 0.0, "CG": -0.149, "CD1": 0.016, "CD2": 0.075,
            "NE1": -0.096, "CE2": 0.138, "CE3": -0.114,
            "CZ2": -0.114, "CZ3": -0.114, "CH2": -0.114},
    "SER": {"CB": 0.0, "OG": -0.530},
    "THR": {"CB": 0.0, "OG1": -0.530, "CG2": 0.0},
    "ASN": {"CB": 0.0, "CG": 0.550, "OD1": -0.555, "ND2": -0.295},
    "GLN": {"CB": 0.0, "CG": 0.0, "CD": 0.550, "OE1": -0.555, "NE2": -0.295},
    "CYS": {"CB": 0.0, "SG": -0.312},
    "HIS": {"CB": 0.0, "CG": -0.032, "ND1": -0.146, "CD2": 0.114,
            "CE1": 0.166, "NE2": -0.572},
    "ASP": {"CB": 0.0, "CG": 0.795, "OD1": -0.556, "OD2": -0.556},
    "GLU": {"CB": 0.0, "CG": 0.0, "CD": 0.805, "OE1": -0.556, "OE2": -0.556},
    "LYS": {"CB": 0.0, "CG": 0.0, "CD": 0.0, "CE": 0.0, "NZ": 0.0},
    "ARG": {"CB": 0.0, "CG": 0.0, "CD": 0.0,
            "NE": -0.101, "CZ": 0.807, "NH1": -0.153, "NH2": -0.153},
    "GLY": {},
}

# ---------------------------------------------------------------------------
# Retinal connectivity (11-cis RPSB)
# ---------------------------------------------------------------------------
# Bonds between heavy atoms in retinal + NZ
RETINAL_BONDS = {
    "C1":  ["C2", "C6", "C16", "C17"],
    "C2":  ["C1", "C3"],
    "C3":  ["C2", "C4"],
    "C4":  ["C3", "C5"],
    "C5":  ["C4", "C6", "C18"],
    "C6":  ["C5", "C1", "C7"],
    "C7":  ["C6", "C8"],
    "C8":  ["C7", "C9"],
    "C9":  ["C8", "C10", "C19"],
    "C10": ["C9", "C11"],
    "C11": ["C10", "C12"],
    "C12": ["C11", "C13"],
    "C13": ["C12", "C14", "C20"],
    "C14": ["C13", "C15"],
    "C15": ["C14", "NZ"],
    "C16": ["C1"],
    "C17": ["C1"],
    "C18": ["C5"],
    "C19": ["C9"],
    "C20": ["C13"],
    "NZ":  ["C15"],  # CE replaced by cap H, proton H added
}

# Hybridization: sp2 atoms in the conjugated system
SP2_ATOMS = {"C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12",
             "C13", "C14", "C15", "NZ"}
# sp3 with 1 neighbor = methyl
METHYL_ATOMS = {"C16", "C17", "C18", "C19", "C20"}
# sp3 with 2 neighbors = ring CH2
SP3_CH2_ATOMS = {"C2", "C3", "C4"}
# sp3 with 4 heavy atom neighbors = no H
QUATERNARY = {"C1"}


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------
def parse_pdb(path: Path) -> list[dict]:
    atoms = []
    with open(path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atoms.append({
                "record": line[:6].strip(),
                "serial": int(line[6:11]),
                "name":   line[12:16].strip(),
                "resn":   line[17:20].strip(),
                "chain":  line[21],
                "resi":   int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": line[76:78].strip() if len(line) > 76 else "",
            })
    return atoms


def is_heavy(atom: dict) -> bool:
    elem = atom.get("element", "")
    if elem:
        return elem not in ("H", "D")
    return not atom["name"].startswith("H")


# ---------------------------------------------------------------------------
# Hydrogen placement geometry
# ---------------------------------------------------------------------------
def _norm(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-10 else v


def add_sp2_H(center, n1, n2, bl=1.08):
    """One H for sp2 atom with 2 heavy-atom neighbors."""
    v1 = _norm(n1 - center)
    v2 = _norm(n2 - center)
    h_dir = -_norm(v1 + v2)
    return [center + bl * h_dir]


def add_sp3_2H(center, n1, n2, bl=1.09):
    """Two H's for sp3 CH2 with 2 heavy-atom neighbors (tetrahedral)."""
    v1 = _norm(n1 - center)
    v2 = _norm(n2 - center)
    c12 = np.dot(v1, v2)
    normal = np.cross(v1, v2)
    nn = np.linalg.norm(normal)
    if nn < 1e-10:
        # Degenerate: pick arbitrary perpendicular
        perp = np.array([1, 0, 0]) if abs(v1[0]) < 0.9 else np.array([0, 1, 0])
        normal = _norm(np.cross(v1, perp))
    else:
        normal = normal / nn
    # Exact tetrahedral solution: d·v1 = d·v2 = -1/3
    a = -1.0 / (3.0 * (1.0 + c12))
    c_sq = 1.0 - 2.0 / (9.0 * (1.0 + c12))
    c_val = np.sqrt(max(c_sq, 0.0))
    d1 = a * (v1 + v2) + c_val * normal
    d2 = a * (v1 + v2) - c_val * normal
    return [center + bl * d1, center + bl * d2]


def add_methyl_3H(center, parent, bl=1.09):
    """Three H's for methyl group (sp3, 1 heavy-atom neighbor)."""
    v = _norm(parent - center)
    perp = np.array([1, 0, 0]) if abs(v[0]) < 0.9 else np.array([0, 1, 0])
    p1 = _norm(np.cross(v, perp))
    p2 = _norm(np.cross(v, p1))
    # H-C-parent angle = 109.47°: cos = -1/3, sin = sqrt(8/9)
    cos_t = -1.0 / 3.0
    sin_t = np.sqrt(8.0 / 9.0)
    positions = []
    for i in range(3):
        phi = 2.0 * np.pi * i / 3.0
        h_dir = cos_t * (-v) + sin_t * (np.cos(phi) * p1 + np.sin(phi) * p2)
        positions.append(center + bl * _norm(h_dir))
    return positions


# ---------------------------------------------------------------------------
# Build QM region: retinal + NZ + hydrogens
# ---------------------------------------------------------------------------
def build_qm_region(all_atoms):
    """Build QM atom list: [(element, [x,y,z]), ...], index map, and H→heavy mapping."""
    # Collect retinal heavy atoms + NZ
    heavy = {}  # name -> (element, coord)
    pdb_atoms = {}  # name -> PDB atom dict
    for a in all_atoms:
        if a["chain"] != CHAIN:
            continue
        if a["resn"] == "RET" and is_heavy(a) and a["name"] in RETINAL_BONDS:
            heavy[a["name"]] = ("C", np.array([a["x"], a["y"], a["z"]]))
            pdb_atoms[a["name"]] = a
        if a["resi"] == 296 and a["resn"] == "LYS" and a["name"] == "NZ":
            heavy["NZ"] = ("N", np.array([a["x"], a["y"], a["z"]]))
            pdb_atoms["NZ"] = a

    # K296 CE position (for cap H direction)
    ce_pos = None
    for a in all_atoms:
        if a["chain"] == CHAIN and a["resi"] == 296 and a["name"] == "CE":
            ce_pos = np.array([a["x"], a["y"], a["z"]])

    if "NZ" not in heavy or "C15" not in heavy:
        raise RuntimeError("Missing NZ or C15 in PDB")

    # Build atom list: heavy atoms first, then hydrogens
    qm_atoms = []
    heavy_order = []  # track which heavy atom each index corresponds to

    # Add heavy atoms in a fixed order
    atom_names_ordered = [
        "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10",
        "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20",
        "NZ",
    ]
    for name in atom_names_ordered:
        if name not in heavy:
            raise RuntimeError(f"Missing atom {name}")
        elem, coord = heavy[name]
        qm_atoms.append((elem, coord.tolist()))
        heavy_order.append(name)

    n_heavy = len(qm_atoms)

    # Add hydrogens, tracking which heavy atom each H belongs to
    h_parent = []  # index into heavy_order for each H atom

    def coord(name):
        return heavy[name][1]

    for name in atom_names_ordered:
        if name in QUATERNARY:
            continue  # C1: 4 heavy neighbors, no H

        heavy_idx = heavy_order.index(name)
        neighbors = RETINAL_BONDS[name]
        n_coords = [coord(n) for n in neighbors if n in heavy]

        if name in SP2_ATOMS:
            if len(n_coords) == 2:
                # sp2 with 2 heavy neighbors → 1 H
                for h in add_sp2_H(coord(name), n_coords[0], n_coords[1]):
                    qm_atoms.append(("H", h.tolist()))
                    h_parent.append(heavy_idx)
            elif name == "NZ":
                # Schiff base: add proton H and cap H (replacing CE)
                nz_c = coord("NZ")
                c15_c = coord("C15")
                if ce_pos is not None:
                    v_c15 = _norm(c15_c - nz_c)
                    v_ce = _norm(ce_pos - nz_c)
                    # Proton: in plane, bisecting the open angle
                    h_dir = -_norm(v_c15 + v_ce)
                    qm_atoms.append(("H", (nz_c + 1.01 * h_dir).tolist()))
                    h_parent.append(heavy_idx)
                    # Cap H: along NZ→CE direction
                    qm_atoms.append(("H", (nz_c + 1.01 * v_ce).tolist()))
                    h_parent.append(heavy_idx)
                else:
                    # Fallback: place H along C15→NZ extended
                    h_dir = _norm(nz_c - c15_c)
                    qm_atoms.append(("H", (nz_c + 1.01 * h_dir).tolist()))
                    h_parent.append(heavy_idx)
            # sp2 with 3 heavy neighbors → 0 H (C5, C6, C9, C13)

        elif name in SP3_CH2_ATOMS:
            if len(n_coords) == 2:
                for h in add_sp3_2H(coord(name), n_coords[0], n_coords[1]):
                    qm_atoms.append(("H", h.tolist()))
                    h_parent.append(heavy_idx)

        elif name in METHYL_ATOMS:
            if len(n_coords) == 1:
                for h in add_methyl_3H(coord(name), n_coords[0]):
                    qm_atoms.append(("H", h.tolist()))
                    h_parent.append(heavy_idx)

    print(f"QM region: {n_heavy} heavy atoms + {len(qm_atoms) - n_heavy} H = {len(qm_atoms)} total")
    return qm_atoms, heavy_order, n_heavy, h_parent


# ---------------------------------------------------------------------------
# Build MM region: protein point charges
# ---------------------------------------------------------------------------
def build_mm_region(all_atoms, retinal_centroid):
    """Protein heavy atoms within MM_CUTOFF as AMBER point charges."""
    mm_coords = []
    mm_charges = []

    # Atoms to exclude from MM (they are in QM region)
    qm_resi_excl = set()  # RET residue
    for a in all_atoms:
        if a["chain"] == CHAIN and a["resn"] == "RET":
            qm_resi_excl.add(("RET", a["resi"]))

    for a in all_atoms:
        if a["chain"] != CHAIN or not is_heavy(a):
            continue
        if a["record"] != "ATOM":
            continue

        coord = np.array([a["x"], a["y"], a["z"]])
        dist = np.linalg.norm(coord - retinal_centroid)
        if dist > MM_CUTOFF:
            continue

        # Exclude K296 NZ and CE (in QM region or capped)
        if a["resi"] == 296 and a["name"] in ("NZ", "CE"):
            continue

        # Assign charge
        name = a["name"]
        resn = a["resn"]
        q = 0.0

        if name in BACKBONE_CHARGES:
            q = BACKBONE_CHARGES[name]
        else:
            sc = SIDECHAIN_CHARGES.get(resn, {})
            q = sc.get(name, 0.0)

        if abs(q) < 1e-6:
            continue  # skip uncharged atoms

        mm_coords.append(coord)
        mm_charges.append(q)

    mm_coords = np.array(mm_coords) if mm_coords else np.zeros((0, 3))
    mm_charges = np.array(mm_charges)

    print(f"MM region: {len(mm_charges)} charged atoms within {MM_CUTOFF}Å")
    return mm_coords, mm_charges


# ---------------------------------------------------------------------------
# PySCF QM/MM calculation
# ---------------------------------------------------------------------------
def run_qm(qm_atoms, mm_coords, mm_charges):
    """Run vacuum and embedded DFT, return Mulliken charges for both."""
    from pyscf import gto, dft, qmmm

    # Build molecule
    mol = gto.Mole()
    mol.atom = qm_atoms
    mol.basis = BASIS
    mol.charge = 1   # protonated Schiff base
    mol.spin = 0     # closed shell singlet
    mol.unit = "Angstrom"
    mol.verbose = 3
    mol.build()

    print(f"\nBasis functions: {mol.nao}")
    print(f"Electrons: {mol.nelectron}")

    # --- Vacuum DFT ---
    print("\n" + "=" * 60)
    print("Running vacuum DFT...")
    print("=" * 60)
    t0 = time.time()
    mf_vac = dft.RKS(mol)
    mf_vac.xc = XC
    mf_vac.kernel()
    t_vac = time.time() - t0
    print(f"Vacuum DFT converged in {t_vac:.1f}s, E = {mf_vac.e_tot:.6f} Ha")

    # Mulliken charges: mulliken_pop returns (AO_pop, atomic_charges)
    # where atomic_charges = nuclear_charge - electron_population
    _, q_vac = mf_vac.mulliken_pop(verbose=0)
    q_vac = np.array(q_vac)

    # --- Embedded DFT (protein point charges) ---
    print("\n" + "=" * 60)
    print("Running embedded DFT (protein point charges)...")
    print("=" * 60)
    t0 = time.time()

    # PySCF qmmm expects coordinates in Bohr
    ANG2BOHR = 1.8897259886
    mm_coords_bohr = mm_coords * ANG2BOHR

    mf_emb = qmmm.mm_charge(dft.RKS(mol), mm_coords_bohr, mm_charges)
    mf_emb.xc = XC
    mf_emb.kernel()
    t_emb = time.time() - t0
    print(f"Embedded DFT converged in {t_emb:.1f}s, E = {mf_emb.e_tot:.6f} Ha")

    _, q_emb = mf_emb.mulliken_pop(verbose=0)
    q_emb = np.array(q_emb)

    return q_vac, q_emb


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    if not PDB_PATH.exists():
        print(f"ERROR: PDB file not found: {PDB_PATH}")
        return 1

    all_atoms = parse_pdb(PDB_PATH)
    print(f"Parsed {len(all_atoms)} atoms from {PDB_PATH.name}")

    # Build QM region
    qm_atoms, heavy_order, n_heavy, h_parent = build_qm_region(all_atoms)

    # Retinal centroid for MM cutoff
    ret_coords = np.array([c for _, c in qm_atoms[:n_heavy]])
    centroid = ret_coords.mean(axis=0)

    # Build MM region
    mm_coords, mm_charges = build_mm_region(all_atoms, centroid)

    # Run QM
    q_vac, q_emb = run_qm(qm_atoms, mm_coords, mm_charges)

    # Collapse H charges onto bonded heavy atoms
    q_vac_col = q_vac[:n_heavy].copy()
    q_emb_col = q_emb[:n_heavy].copy()
    for h_idx, parent_idx in enumerate(h_parent):
        q_vac_col[parent_idx] += q_vac[n_heavy + h_idx]
        q_emb_col[parent_idx] += q_emb[n_heavy + h_idx]

    delta_col = q_emb_col - q_vac_col

    # Summary table
    print("\n" + "=" * 80)
    print(f"{'Atom':<6s} {'q_vac(col)':>10s} {'q_emb(col)':>10s} {'Δq':>10s} {'Δq×100':>8s}  Interpretation")
    print("-" * 80)
    for i, name in enumerate(heavy_order):
        dq = delta_col[i]
        interp = ""
        if dq > 0.003:
            interp = "← protein pulled e⁻ away (more +)"
        elif dq < -0.003:
            interp = "← protein pushed e⁻ here (more −)"
        print(f"{name:<6s} {q_vac_col[i]:>10.4f} {q_emb_col[i]:>10.4f} {dq:>+10.4f} {dq*100:>+8.2f}  {interp}")
    print("=" * 80)
    print(f"\nΔq range: {delta_col.min():+.4f} to {delta_col.max():+.4f}")
    print(f"Sum Δq (should be ~0): {delta_col.sum():+.6f}")

    # Write output PDB: B-factor = Δq × 100 (centi-electrons, H-collapsed)
    ret_pdb = {}
    for a in all_atoms:
        if a["chain"] == CHAIN and a["resn"] == "RET" and is_heavy(a):
            ret_pdb[a["name"]] = a
        if a["chain"] == CHAIN and a["resi"] == 296 and a["name"] == "NZ":
            ret_pdb["NZ"] = a

    with open(OUTPUT_PDB, "w") as f:
        f.write("REMARK   QM/MM charge redistribution: binding pocket effect on retinal\n")
        f.write(f"REMARK   Method: {XC.upper()}/{BASIS} + AMBER point charges\n")
        f.write(f"REMARK   B-factor = Δq × 100 (H-collapsed, centi-electrons)\n")
        f.write(f"REMARK   Δq = Mulliken charge (embedded) − Mulliken charge (vacuum)\n")
        f.write(f"REMARK   Positive = protein made atom more positive (pulled e⁻ away)\n")
        f.write(f"REMARK   Negative = protein pushed e⁻ toward this atom\n")

        serial = 1
        for i, name in enumerate(heavy_order):
            a = ret_pdb[name]
            record = "HETATM" if a["resn"] == "RET" else "ATOM  "
            resn = a["resn"]
            resi = a["resi"]
            elem = a.get("element", "C") if a["resn"] == "RET" else "N"

            bfac = delta_col[i] * 100.0

            f.write(
                f"{record}{serial:5d}  {name:<3s} {resn:>3s} {CHAIN}"
                f"{resi:4d}    "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                f"{1.00:6.2f}{bfac:6.2f}"
                f"          {elem:>2s}\n"
            )
            serial += 1

        f.write("END\n")

    print(f"\nWrote {OUTPUT_PDB}")
    print(f"B-factor = Δq × 100 (centi-electrons, H-collapsed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
