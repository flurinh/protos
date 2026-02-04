#!/usr/bin/env python3
"""BoltzGen Refined Alignment - Fine-tune top rotations with smaller steps.

Takes the top placements from coarse search and refines with:
- Smaller rotation steps (5° instead of 30°)
- Local search around each top rotation (±45°)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime
from itertools import product

import numpy as np
import pandas as pd
import yaml
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

import protos
from protos.processing.structure import StructureProcessor
from protos.io.formats.cif_utils import write_cif_file
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Configuration
# =============================================================================

TRYPSIN_PDB = "2age"
TRYPSIN_CHAIN = "X"
SUBSTRATE_CHAIN = "A"

TRYPSIN_TRIAD = [
    {"res_num": 195, "res_name": "SER", "role": "nucleophile"},
    {"res_num": 57, "res_name": "HIS", "role": "base"},
    {"res_num": 102, "res_name": "ASP", "role": "electrostatic"},
]

RHODOPSIN_PDB = "3pqr"
RHODOPSIN_CHAIN = "A"
PEPTIDE_CHAIN = "B"

DESIGNABLE_REGIONS = "58..76,134..157,222..256,306..326"
RETINAL_CODES = ["RET", "LYR"]
RETINAL_SMILES = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C"
SUBSTRATE_SMILES = "NC(=N)c1ccccc1"

# Top rotations from coarse search to refine
TOP_ROTATIONS = [
    {"rotation": [150, 30, 150], "residues": [135, 72, 73], "rmsd": 1.63},
    {"rotation": [90, 30, 240], "residues": [138, 135, 226], "rmsd": 1.76},
    {"rotation": [180, 60, 210], "residues": [139, 138, 134], "rmsd": 2.00},
]

# Refinement parameters
REFINE_RANGE = 45  # Search ±45° around each rotation
REFINE_STEP = 5    # 5° steps for fine search


def parse_designable_regions(region_str: str) -> List[int]:
    residues = []
    for part in region_str.split(","):
        if ".." in part:
            start, end = map(int, part.split(".."))
            residues.extend(range(start, end + 1))
        else:
            residues.append(int(part))
    return residues


def get_ca_positions(df: pd.DataFrame, chain: str, residues: List[int] = None) -> Tuple[np.ndarray, List[int]]:
    mask = (df["auth_chain_id"] == chain) & (df["atom_name"] == "CA")
    if residues is not None:
        mask = mask & (df["auth_seq_id"].isin(residues))
    ca_df = df[mask].sort_values("auth_seq_id")
    coords = ca_df[["x", "y", "z"]].values
    resids = ca_df["auth_seq_id"].values.tolist()
    return coords, resids


def get_residue_atoms(df: pd.DataFrame, chain: str, res_num: int) -> pd.DataFrame:
    mask = (df["auth_chain_id"] == chain) & (df["auth_seq_id"] == res_num)
    return df[mask].copy()


def rotate_points(points: np.ndarray, center: np.ndarray, rotation: Rotation) -> np.ndarray:
    centered = points - center
    rotated = rotation.apply(centered)
    return rotated + center


def compute_placement_score(triad_ca: np.ndarray, target_ca: np.ndarray) -> Tuple[float, List[Tuple[int, float]]]:
    assignments = []
    used = set()
    for triad_pos in triad_ca:
        dists = np.linalg.norm(target_ca - triad_pos, axis=1)
        for idx in np.argsort(dists):
            if idx not in used:
                assignments.append((idx, dists[idx]))
                used.add(idx)
                break
    total_dist = sum(d for _, d in assignments)
    return total_dist, assignments


def refine_rotation(
    initial_rotation: List[int],
    triad_ca_original: np.ndarray,
    substrate_centroid: np.ndarray,
    peptide_centroid: np.ndarray,
    target_ca: np.ndarray,
    target_resids: List[int],
    search_range: int = 45,
    step: int = 5,
) -> Dict:
    """Refine a rotation with finer grid search."""

    # Translation
    translation = peptide_centroid - substrate_centroid
    triad_ca_translated = triad_ca_original + translation

    rx0, ry0, rz0 = initial_rotation

    best_result = None
    best_score = float('inf')

    # Grid search around initial rotation
    rx_range = range(rx0 - search_range, rx0 + search_range + 1, step)
    ry_range = range(ry0 - search_range, ry0 + search_range + 1, step)
    rz_range = range(rz0 - search_range, rz0 + search_range + 1, step)

    total_evals = 0

    for rx in rx_range:
        for ry in ry_range:
            for rz in rz_range:
                # Normalize angles to 0-360
                rx_norm = rx % 360
                ry_norm = ry % 360
                rz_norm = rz % 360

                rot = Rotation.from_euler('xyz', [rx_norm, ry_norm, rz_norm], degrees=True)
                triad_ca_rotated = rotate_points(triad_ca_translated, peptide_centroid, rot)

                score, assignments = compute_placement_score(triad_ca_rotated, target_ca)
                total_evals += 1

                if score < best_score:
                    best_score = score
                    assigned_resids = [target_resids[idx] for idx, _ in assignments]
                    distances = [d for _, d in assignments]

                    best_result = {
                        "rotation_xyz_deg": [rx_norm, ry_norm, rz_norm],
                        "score": score,
                        "assigned_residues": assigned_resids,
                        "distances": distances,
                        "triad_ca_positions": triad_ca_rotated.tolist(),
                        "rmsd": np.sqrt(np.mean(np.array(distances) ** 2)),
                    }

    best_result["evaluations"] = total_evals
    return best_result


def scipy_refine_rotation(
    initial_rotation: List[int],
    triad_ca_original: np.ndarray,
    substrate_centroid: np.ndarray,
    peptide_centroid: np.ndarray,
    target_ca: np.ndarray,
    target_resids: List[int],
) -> Dict:
    """Use scipy optimization for continuous refinement."""

    translation = peptide_centroid - substrate_centroid
    triad_ca_translated = triad_ca_original + translation

    def objective(angles):
        rx, ry, rz = angles
        rot = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True)
        triad_ca_rotated = rotate_points(triad_ca_translated, peptide_centroid, rot)
        score, _ = compute_placement_score(triad_ca_rotated, target_ca)
        return score

    # Optimize
    result = minimize(
        objective,
        initial_rotation,
        method='L-BFGS-B',
        options={'maxiter': 1000}
    )

    # Get final result
    rx, ry, rz = result.x
    rot = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True)
    triad_ca_rotated = rotate_points(triad_ca_translated, peptide_centroid, rot)
    score, assignments = compute_placement_score(triad_ca_rotated, target_ca)

    assigned_resids = [target_resids[idx] for idx, _ in assignments]
    distances = [d for _, d in assignments]

    return {
        "rotation_xyz_deg": [rx % 360, ry % 360, rz % 360],
        "score": score,
        "assigned_residues": assigned_resids,
        "distances": distances,
        "triad_ca_positions": triad_ca_rotated.tolist(),
        "rmsd": np.sqrt(np.mean(np.array(distances) ** 2)),
        "scipy_success": result.success,
    }


def transform_residue_atoms(
    atoms_df: pd.DataFrame,
    original_centroid: np.ndarray,
    target_centroid: np.ndarray,
    rotation: Rotation,
) -> pd.DataFrame:
    atoms_df = atoms_df.copy()
    coords = atoms_df[["x", "y", "z"]].values
    translation = target_centroid - original_centroid
    coords_translated = coords + translation
    coords_rotated = rotate_points(coords_translated, target_centroid, rotation)
    atoms_df["x"] = coords_rotated[:, 0]
    atoms_df["y"] = coords_rotated[:, 1]
    atoms_df["z"] = coords_rotated[:, 2]
    return atoms_df


def create_grafted_structure(
    rhodopsin_df: pd.DataFrame,
    transformed_triad: List[pd.DataFrame],
    triad_info: List[Dict],
    target_residues: List[int],
    output_path: Path,
) -> pd.DataFrame:
    grafted_df = rhodopsin_df.copy()
    grafted_df = grafted_df[grafted_df["auth_chain_id"] != PEPTIDE_CHAIN].copy()

    existing_mapping = grafted_df[["auth_seq_id", "label_seq_id"]].drop_duplicates()
    auth_to_label = dict(zip(existing_mapping["auth_seq_id"], existing_mapping["label_seq_id"]))

    for i, (triad_atoms, target_res) in enumerate(zip(transformed_triad, target_residues)):
        mask = ~((grafted_df["auth_seq_id"] == target_res) &
                 (grafted_df["auth_chain_id"] == RHODOPSIN_CHAIN))
        grafted_df = grafted_df[mask].copy()

        triad_atoms = triad_atoms.copy()
        triad_atoms["auth_seq_id"] = target_res
        triad_atoms["auth_chain_id"] = RHODOPSIN_CHAIN
        triad_atoms["label_chain_id"] = RHODOPSIN_CHAIN
        triad_atoms["entity_id"] = "1"
        triad_atoms["gen_seq_id"] = target_res

        if target_res in auth_to_label:
            triad_atoms["label_seq_id"] = auth_to_label[target_res]
        else:
            triad_atoms["label_seq_id"] = target_res

        triad_atoms["res_name"] = triad_info[i]["res_name"]
        grafted_df = pd.concat([grafted_df, triad_atoms], ignore_index=True)

    grafted_df = grafted_df.sort_values(["auth_chain_id", "auth_seq_id", "atom_id"])

    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    ]
    keep_codes = amino_acids + RETINAL_CODES
    grafted_df = grafted_df[grafted_df["res_name"].isin(keep_codes)].copy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_cif_file(str(output_path), grafted_df, force_overwrite=True)

    return grafted_df


def main() -> int:
    print("=" * 70)
    print("REFINED ALIGNMENT OPTIMIZATION")
    print("Fine-tuning top rotations with smaller steps")
    print("=" * 70)

    # Initialize
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))

    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)

    # ==========================================================================
    # Step 1: Load structures
    # ==========================================================================
    print("\n[1] Loading structures...")

    loader.download_and_register(TRYPSIN_PDB, name="trypsin_2age")
    loader.download_and_register(RHODOPSIN_PDB, name=RHODOPSIN_PDB)

    trypsin_df = sp.load_entity("trypsin_2age").reset_index()
    rhodopsin_df = sp.load_entity(RHODOPSIN_PDB).reset_index()

    # ==========================================================================
    # Step 2: Extract positions
    # ==========================================================================
    print("\n[2] Extracting key positions...")

    triad_ca = []
    triad_atoms_list = []
    for res_info in TRYPSIN_TRIAD:
        atoms = get_residue_atoms(trypsin_df, TRYPSIN_CHAIN, res_info["res_num"])
        triad_atoms_list.append(atoms)
        ca = atoms[atoms["atom_name"] == "CA"][["x", "y", "z"]].values[0]
        triad_ca.append(ca)
    triad_ca = np.array(triad_ca)

    substrate_ca, _ = get_ca_positions(trypsin_df, SUBSTRATE_CHAIN)
    substrate_centroid = substrate_ca.mean(axis=0)

    peptide_ca, _ = get_ca_positions(rhodopsin_df, PEPTIDE_CHAIN)
    peptide_centroid = peptide_ca.mean(axis=0)

    designable_resids = parse_designable_regions(DESIGNABLE_REGIONS)
    target_ca, target_resids = get_ca_positions(rhodopsin_df, RHODOPSIN_CHAIN, designable_resids)

    print(f"  Substrate centroid: [{substrate_centroid[0]:.1f}, {substrate_centroid[1]:.1f}, {substrate_centroid[2]:.1f}]")
    print(f"  Peptide centroid: [{peptide_centroid[0]:.1f}, {peptide_centroid[1]:.1f}, {peptide_centroid[2]:.1f}]")

    # ==========================================================================
    # Step 3: Refine each top rotation
    # ==========================================================================
    print("\n[3] Refining top rotations...")
    print(f"  Search range: ±{REFINE_RANGE}°, Step: {REFINE_STEP}°")

    refined_placements = []

    for i, top in enumerate(TOP_ROTATIONS):
        print(f"\n  Refining rotation {i+1}: {top['rotation']} (initial RMSD: {top['rmsd']:.2f}Å)")

        # Grid refinement
        grid_result = refine_rotation(
            top["rotation"],
            triad_ca,
            substrate_centroid,
            peptide_centroid,
            target_ca,
            target_resids,
            search_range=REFINE_RANGE,
            step=REFINE_STEP,
        )

        print(f"    Grid search: {grid_result['evaluations']} evaluations")
        print(f"    Grid result: rotation={grid_result['rotation_xyz_deg']}, RMSD={grid_result['rmsd']:.3f}Å")

        # Scipy continuous refinement
        scipy_result = scipy_refine_rotation(
            grid_result["rotation_xyz_deg"],
            triad_ca,
            substrate_centroid,
            peptide_centroid,
            target_ca,
            target_resids,
        )

        print(f"    Scipy result: rotation=[{scipy_result['rotation_xyz_deg'][0]:.1f}, {scipy_result['rotation_xyz_deg'][1]:.1f}, {scipy_result['rotation_xyz_deg'][2]:.1f}], RMSD={scipy_result['rmsd']:.3f}Å")

        # Use the better result
        if scipy_result["rmsd"] < grid_result["rmsd"]:
            best_result = scipy_result
            best_result["method"] = "scipy"
        else:
            best_result = grid_result
            best_result["method"] = "grid"

        best_result["original_rotation"] = top["rotation"]
        best_result["original_rmsd"] = top["rmsd"]
        best_result["improvement"] = top["rmsd"] - best_result["rmsd"]

        refined_placements.append(best_result)

        print(f"    ✓ Best: RMSD={best_result['rmsd']:.3f}Å (improved by {best_result['improvement']:.3f}Å)")

    # Sort by RMSD
    refined_placements.sort(key=lambda x: x["rmsd"])

    # ==========================================================================
    # Step 4: Summary of refinements
    # ==========================================================================
    print("\n" + "=" * 70)
    print("REFINEMENT RESULTS")
    print("=" * 70)
    print(f"\n{'Rank':<5} {'Original RMSD':<15} {'Refined RMSD':<15} {'Improvement':<12} {'Residues'}")
    print("-" * 70)

    for i, p in enumerate(refined_placements):
        res_str = ", ".join(str(r) for r in p['assigned_residues'])
        print(f"{i+1:<5} {p['original_rmsd']:<15.3f} {p['rmsd']:<15.3f} {p['improvement']:<12.3f} {res_str}")

    # ==========================================================================
    # Step 5: Create grafted structures
    # ==========================================================================
    print("\n[4] Creating refined grafted structures...")

    grafted_dir = OUTPUT_DIR / "grafted_structures_refined"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    yaml_dir = OUTPUT_DIR / "boltzgen_configs_refined"
    yaml_dir.mkdir(parents=True, exist_ok=True)

    # Prepare rhodopsin
    rhodopsin_chain_a = rhodopsin_df[rhodopsin_df["auth_chain_id"] == RHODOPSIN_CHAIN].copy()
    unique_residues = rhodopsin_chain_a[["auth_seq_id", "res_name"]].drop_duplicates()
    unique_residues = unique_residues.sort_values("auth_seq_id")
    auth_to_label = {auth_seq: i + 1 for i, auth_seq in enumerate(unique_residues["auth_seq_id"].values)}
    rhodopsin_chain_a["label_seq_id"] = rhodopsin_chain_a["auth_seq_id"].map(auth_to_label)
    rhodopsin_chain_a["label_chain_id"] = rhodopsin_chain_a["auth_chain_id"]
    rhodopsin_chain_a["entity_id"] = "1"

    all_configs = []

    for rank, placement in enumerate(refined_placements, 1):
        rot_angles = placement["rotation_xyz_deg"]
        rotation = Rotation.from_euler('xyz', rot_angles, degrees=True)

        # Transform triad
        transformed_triad = []
        for atoms_df in triad_atoms_list:
            transformed = transform_residue_atoms(
                atoms_df,
                substrate_centroid,
                peptide_centroid,
                rotation,
            )
            transformed_triad.append(transformed)

        # Create structure
        config_name = f"trypsin_refined{rank:02d}_rmsd{placement['rmsd']:.2f}"
        output_path = grafted_dir / f"{config_name}.cif"

        grafted_df = create_grafted_structure(
            rhodopsin_chain_a,
            transformed_triad,
            TRYPSIN_TRIAD,
            placement["assigned_residues"],
            output_path,
        )

        print(f"  Rank {rank}: {config_name}")
        print(f"    Residues: {placement['assigned_residues']}")
        print(f"    RMSD: {placement['rmsd']:.3f}Å")

        # BoltzGen config
        triad_str = ",".join(str(r) for r in placement["assigned_residues"])

        yaml_config = {
            "entities": [
                {
                    "file": {
                        "path": str(output_path),
                        "include": [{"chain": {"id": RHODOPSIN_CHAIN}}],
                        "design": [{"chain": {"id": RHODOPSIN_CHAIN, "res_index": DESIGNABLE_REGIONS}}],
                        "not_design": [{"chain": {"id": RHODOPSIN_CHAIN, "res_index": triad_str}}],
                        "structure_groups": [
                            {"group": {"id": "all", "visibility": 1}},
                            {"group": {"id": RHODOPSIN_CHAIN, "res_index": DESIGNABLE_REGIONS, "visibility": 0}},
                            {"group": {"id": RHODOPSIN_CHAIN, "res_index": triad_str, "visibility": 1}},
                        ],
                        "binding_types": [{"chain": {"id": RHODOPSIN_CHAIN, "binding": triad_str}}],
                    }
                },
                {"ligand": {"id": "SUB", "smiles": SUBSTRATE_SMILES}},
                {"ligand": {"id": "RET", "smiles": RETINAL_SMILES}},
            ],
            "num_designs": 8,  # More designs for refined placements
        }

        yaml_path = yaml_dir / f"{config_name}.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(yaml_config, f, default_flow_style=False)

        config = {
            "rank": rank,
            "name": config_name,
            "rotation_xyz_deg": [float(x) for x in rot_angles],
            "assigned_residues": placement["assigned_residues"],
            "distances": placement["distances"],
            "rmsd": placement["rmsd"],
            "original_rmsd": placement["original_rmsd"],
            "improvement": placement["improvement"],
            "method": placement["method"],
            "structure_path": str(output_path),
            "yaml_path": str(yaml_path),
        }
        all_configs.append(config)

    # ==========================================================================
    # Step 6: Save summary
    # ==========================================================================
    print("\n[5] Saving summary...")

    summary = {
        "workflow": "refined_alignment",
        "timestamp": datetime.now().isoformat(),
        "refinement_params": {
            "search_range_degrees": REFINE_RANGE,
            "step_degrees": REFINE_STEP,
        },
        "placements": all_configs,
    }

    summary_path = OUTPUT_DIR / "refined_alignment_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # ==========================================================================
    # Done
    # ==========================================================================
    print("\n" + "=" * 70)
    print("REFINEMENT COMPLETE")
    print("=" * 70)
    print()
    print(f"Grafted structures: {grafted_dir}")
    print(f"BoltzGen configs: {yaml_dir}")
    print()
    print("Best refined placement:")
    best = all_configs[0]
    print(f"  Rotation: [{best['rotation_xyz_deg'][0]:.1f}, {best['rotation_xyz_deg'][1]:.1f}, {best['rotation_xyz_deg'][2]:.1f}]")
    print(f"  Residues: {best['assigned_residues']}")
    print(f"  RMSD: {best['rmsd']:.3f}Å (improved from {best['original_rmsd']:.3f}Å)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
