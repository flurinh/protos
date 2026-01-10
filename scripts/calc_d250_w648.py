#!/usr/bin/env python3
"""Calculate D2.50 to W6.48 closest atom distance for all structures."""

import sys
from pathlib import Path
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.structure import StructureProcessor

protos.set_data_path(str(REPO_ROOT / "data"))
sp = StructureProcessor()

STRUCTURES = {
    "2RH1": ("A", "inverse_agonist", "inactive"),
    "3NY9": ("A", "inverse_agonist", "inactive"),
    "3SN6": ("R", "full_agonist", "active"),
    "4LDO": ("A", "full_agonist", "active"),
    "2VT4": ("A", "antagonist", "inactive"),
    "2Y02": ("A", "full_agonist", "active_like"),
    "2Y04": ("A", "partial_agonist", "intermediate"),
    "2Y00": ("A", "partial_agonist", "intermediate"),
}

print("D2.50 to W6.48 - Closest Atom Distance")
print("=" * 85)
print(f"{'PDB':<6} {'Ligand Type':<18} {'State':<12} {'Distance':<10} {'D2.50 atom':<12} {'W6.48 atom'}")
print("-" * 85)

results = []

for pdb_id, (chain, lig_type, state) in STRUCTURES.items():
    df = sp.load_entity(pdb_id)
    if df is None:
        continue
    df = df.reset_index()

    chain_df = df[(df["auth_chain_id"] == chain) & (df["grn"].notna())]

    # Get D2.50 heavy atoms
    d250_rows = chain_df[(chain_df["grn"] == "2.50") & (~chain_df["element"].isin(["H"]))]
    # Get W6.48 heavy atoms
    w648_rows = chain_df[(chain_df["grn"] == "6.48") & (~chain_df["element"].isin(["H"]))]

    if d250_rows.empty or w648_rows.empty:
        print(f"{pdb_id:<6} {lig_type:<18} {state:<12} Missing residues")
        continue

    d250_coords = d250_rows[["x", "y", "z"]].values
    d250_atoms = d250_rows["atom_name"].values

    w648_coords = w648_rows[["x", "y", "z"]].values
    w648_atoms = w648_rows["atom_name"].values

    # Find closest atom pair
    min_dist = float("inf")
    closest_d250_atom = None
    closest_w648_atom = None

    for i, d_coord in enumerate(d250_coords):
        for j, w_coord in enumerate(w648_coords):
            dist = np.linalg.norm(d_coord - w_coord)
            if dist < min_dist:
                min_dist = dist
                closest_d250_atom = d250_atoms[i]
                closest_w648_atom = w648_atoms[j]

    print(f"{pdb_id:<6} {lig_type:<18} {state:<12} {min_dist:<10.2f} {closest_d250_atom:<12} {closest_w648_atom}")
    results.append((pdb_id, lig_type, state, min_dist))

print()
print("=" * 85)
print("SUMMARY BY LIGAND TYPE:")
print("=" * 85)

from collections import defaultdict
by_type = defaultdict(list)
for pdb, ltype, state, dist in results:
    by_type[ltype].append((pdb, dist))

for ltype in ["full_agonist", "partial_agonist", "antagonist", "inverse_agonist"]:
    if ltype in by_type:
        dists = [d for _, d in by_type[ltype]]
        pdbs = [p for p, _ in by_type[ltype]]
        pdb_str = ", ".join(pdbs)
        print(f"{ltype:<18}: mean={np.mean(dists):.2f}A  range={min(dists):.2f}-{max(dists):.2f}A  [{pdb_str}]")

print()
print("SUMMARY BY STATE:")
print("=" * 85)
by_state = defaultdict(list)
for pdb, ltype, state, dist in results:
    by_state[state].append((pdb, dist))

for state in ["active", "active_like", "intermediate", "inactive"]:
    if state in by_state:
        dists = [d for _, d in by_state[state]]
        pdbs = [p for p, _ in by_state[state]]
        pdb_str = ", ".join(pdbs)
        print(f"{state:<12}: mean={np.mean(dists):.2f}A  range={min(dists):.2f}-{max(dists):.2f}A  [{pdb_str}]")
