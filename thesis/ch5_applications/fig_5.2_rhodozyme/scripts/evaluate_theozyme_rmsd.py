#!/usr/bin/env python3
"""Evaluate theozyme all-atom RMSD for all Boltz2 predictions.

For each Boltz2 CIF, aligns ONLY the 3 theozyme residues (all atoms)
to the reference placement PDB, then computes the RMSD.

Outputs a ranked table: design, seq, RMSD, n_atoms, pLDDT(theozyme), pLDDT(global).
"""

import sys
import os
import glob
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration (from rhodozyme_config)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[5]
DATA = Path(os.environ.get("PROTOS_MODEL_DATA", REPO_ROOT / "data" / "models"))
PLACEMENT_NUM = "00"

THEOZYME_SER = 230
THEOZYME_HIS = 245
THEOZYME_ASP = 250
TZ_RESI = [THEOZYME_SER, THEOZYME_HIS, THEOZYME_ASP]
TZ_SEL = f"{THEOZYME_SER}+{THEOZYME_HIS}+{THEOZYME_ASP}"

PLACEMENT_PDB = DATA / "rfdiffusion2" / "input" / f"placement_{PLACEMENT_NUM}_triad_ori.pdb"
BOLTZ_DIR = DATA / "boltz2" / f"20260203_placement_{PLACEMENT_NUM}_production"

# ---------------------------------------------------------------------------
# PyMOL setup (headless)
# ---------------------------------------------------------------------------
import pymol
pymol.finish_launching(["pymol", "-cq"])
from pymol import cmd

cmd.feedback("disable", "all", "everything")

# Load reference
cmd.load(str(PLACEMENT_PDB), "ref")
cmd.remove("ref and not chain A")
cmd.remove("ref and hetatm")

ref_tz = f"ref and resi {TZ_SEL} and chain A"
ref_model = cmd.get_model(ref_tz)
ref_n_atoms = len(ref_model.atom)
print(f"Reference: {PLACEMENT_PDB.name}")
print(f"Theozyme residues: Ser{THEOZYME_SER}, His{THEOZYME_HIS}, Asp{THEOZYME_ASP}")
print(f"Reference theozyme atoms: {ref_n_atoms}")
print()

# ---------------------------------------------------------------------------
# Parse pLDDT from CIF (B_iso_or_equiv column)
# ---------------------------------------------------------------------------
def parse_plddt(cif_path, resi_list=None):
    """Return per-residue pLDDT dict and global mean."""
    plddt = {}
    with open(cif_path) as f:
        in_atom = False
        for line in f:
            if line.startswith("_atom_site."):
                in_atom = True
                continue
            if in_atom and line.startswith("ATOM"):
                fields = line.split()
                resi = int(fields[6])
                if resi not in plddt:
                    plddt[resi] = float(fields[17])
            elif in_atom and not line.startswith(("ATOM", "HETATM", "_", "#", "loop")):
                if line.strip() in ("", "#"):
                    in_atom = False
    all_vals = list(plddt.values())
    global_mean = np.mean(all_vals) if all_vals else 0.0
    if resi_list:
        tz_vals = [plddt.get(r, 0.0) for r in resi_list]
        tz_mean = np.mean(tz_vals) if tz_vals else 0.0
    else:
        tz_mean = global_mean
    return tz_mean, global_mean

# ---------------------------------------------------------------------------
# Evaluate all Boltz2 predictions
# ---------------------------------------------------------------------------
cif_pattern = str(BOLTZ_DIR / "rhodozyme_00_*_seq*" / "outputs" / "boltz_results_config" / "predictions" / "config" / "config_model_0.cif")
cif_files = sorted(glob.glob(cif_pattern))

print(f"Found {len(cif_files)} Boltz2 predictions")
print("=" * 80)

results = []
for i, cif_path in enumerate(cif_files):
    # Parse design and seq number from path
    dirname = Path(cif_path).parents[4].name  # rhodozyme_00_N-atomized-bb-False_seqM
    parts = dirname.replace("rhodozyme_00_", "").replace("-atomized-bb-False_seq", " ").split()
    design_num = parts[0]
    seq_num = parts[1]

    try:
        # Load prediction
        obj_name = f"pred_{i}"
        cmd.load(cif_path, obj_name)

        # Align ONLY theozyme residues (all atoms) to reference
        pred_tz = f"{obj_name} and resi {TZ_SEL}"

        # Check atom count
        pred_model = cmd.get_model(pred_tz)
        n_pred = len(pred_model.atom)

        if n_pred == 0:
            print(f"  SKIP design={design_num} seq={seq_num}: no theozyme atoms found")
            cmd.delete(obj_name)
            continue

        # align returns (RMSD, n_aligned, n_cycles, RMSD_before, n_aligned_before, score, n_residues)
        # Use cycles=0 to get pure RMSD without outlier rejection
        aln = cmd.align(pred_tz, ref_tz, cycles=5)
        rmsd = aln[0]
        n_aligned = aln[1]

        # pLDDT
        tz_plddt, global_plddt = parse_plddt(cif_path, TZ_RESI)

        results.append({
            "design": design_num,
            "seq": seq_num,
            "rmsd": rmsd,
            "n_atoms": n_aligned,
            "tz_plddt": tz_plddt,
            "global_plddt": global_plddt,
            "cif": cif_path,
        })

        cmd.delete(obj_name)

    except Exception as e:
        print(f"  ERROR design={design_num} seq={seq_num}: {e}")
        try:
            cmd.delete(obj_name)
        except:
            pass

    if (i + 1) % 50 == 0:
        print(f"  Processed {i + 1}/{len(cif_files)}...")

# ---------------------------------------------------------------------------
# Sort and display results
# ---------------------------------------------------------------------------
results.sort(key=lambda r: r["rmsd"])

print()
print(f"{'Rank':>4}  {'Design':>6}  {'Seq':>3}  {'Tz RMSD':>8}  {'Atoms':>5}  {'Tz pLDDT':>8}  {'Global pLDDT':>12}")
print("-" * 65)
for rank, r in enumerate(results, 1):
    print(f"{rank:>4}  {r['design']:>6}  {r['seq']:>3}  {r['rmsd']:>8.3f}  {r['n_atoms']:>5}  {r['tz_plddt']:>8.1f}  {r['global_plddt']:>12.1f}")
    if rank == 30:
        print(f"  ... ({len(results) - 30} more)")
        break

# Top 10 summary
print()
print("=" * 65)
print("TOP 10 candidates (lowest theozyme all-atom RMSD):")
print("=" * 65)
for rank, r in enumerate(results[:10], 1):
    print(f"  {rank}. design={r['design']}, seq={r['seq']}  "
          f"RMSD={r['rmsd']:.3f} A  "
          f"tz_pLDDT={r['tz_plddt']:.1f}  "
          f"global_pLDDT={r['global_plddt']:.1f}")

# Best candidate
if results:
    best = results[0]
    print()
    print(f"BEST: design={best['design']}, seq={best['seq']}")
    print(f"  Theozyme all-atom RMSD: {best['rmsd']:.3f} A ({best['n_atoms']} atoms)")
    print(f"  Theozyme pLDDT: {best['tz_plddt']:.1f}")
    print(f"  Global pLDDT: {best['global_plddt']:.1f}")
    print(f"  CIF: {best['cif']}")

# Save full results to CSV
out_csv = Path(__file__).resolve().parent.parent / "figures" / "theozyme_rmsd_ranking.csv"
with open(out_csv, "w") as f:
    f.write("rank,design,seq,theozyme_rmsd,n_atoms,tz_plddt,global_plddt,cif\n")
    for rank, r in enumerate(results, 1):
        f.write(f"{rank},{r['design']},{r['seq']},{r['rmsd']:.4f},{r['n_atoms']},{r['tz_plddt']:.1f},{r['global_plddt']:.1f},{r['cif']}\n")
print(f"\nFull ranking saved to: {out_csv}")
