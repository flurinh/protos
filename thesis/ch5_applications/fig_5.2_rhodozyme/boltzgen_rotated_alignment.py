#!/usr/bin/env python3
"""BoltzGen Rotated Alignment - Optimize triad orientation in binding pocket.

Problem: Simple centroid alignment leaves catalytic triad pointing outward.

Solution:
1. Translate substrate centroid → G-protein peptide centroid
2. ROTATE the entire assembly (substrate + triad) around that centroid
3. Find optimal rotation minimizing CA-distance to designable residues
4. Generate multiple rotation variants for BoltzGen to evaluate
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

# Rotation grid - sample orientations
ROTATION_ANGLES = list(range(0, 360, 30))  # Every 30 degrees


def parse_designable_regions(region_str: str) -> List[int]:
    """Parse designable region string to list of residue numbers."""
    residues = []
    for part in region_str.split(","):
        if ".." in part:
            start, end = map(int, part.split(".."))
            residues.extend(range(start, end + 1))
        else:
            residues.append(int(part))
    return residues


def get_ca_positions(df: pd.DataFrame, chain: str, residues: List[int] = None) -> Tuple[np.ndarray, List[int]]:
    """Get CA positions for specified chain/residues."""
    mask = (df["auth_chain_id"] == chain) & (df["atom_name"] == "CA")
    if residues is not None:
        mask = mask & (df["auth_seq_id"].isin(residues))

    ca_df = df[mask].sort_values("auth_seq_id")
    coords = ca_df[["x", "y", "z"]].values
    resids = ca_df["auth_seq_id"].values.tolist()
    return coords, resids


def get_residue_atoms(df: pd.DataFrame, chain: str, res_num: int) -> pd.DataFrame:
    """Get all atoms for a specific residue."""
    mask = (df["auth_chain_id"] == chain) & (df["auth_seq_id"] == res_num)
    return df[mask].copy()


def rotate_points(points: np.ndarray, center: np.ndarray, rotation: Rotation) -> np.ndarray:
    """Rotate points around a center."""
    centered = points - center
    rotated = rotation.apply(centered)
    return rotated + center


def compute_placement_score(
    triad_ca: np.ndarray,
    target_ca: np.ndarray,
) -> Tuple[float, List[Tuple[int, float]]]:
    """
    Compute how well the triad fits the target positions.

    Returns: (total_distance, [(nearest_idx, distance), ...])
    """
    assignments = []
    used = set()

    for triad_pos in triad_ca:
        dists = np.linalg.norm(target_ca - triad_pos, axis=1)
        # Find nearest unused target
        for idx in np.argsort(dists):
            if idx not in used:
                assignments.append((idx, dists[idx]))
                used.add(idx)
                break

    total_dist = sum(d for _, d in assignments)
    return total_dist, assignments


def optimize_rotation(
    triad_ca_original: np.ndarray,
    substrate_centroid: np.ndarray,
    peptide_centroid: np.ndarray,
    target_ca: np.ndarray,
    target_resids: List[int],
    angle_step: int = 30,
) -> List[Dict]:
    """
    Find optimal rotations that minimize distance from triad CA to target CA.

    Returns list of top placements sorted by score.
    """
    # First translate triad so substrate centroid aligns with peptide centroid
    translation = peptide_centroid - substrate_centroid
    triad_ca_translated = triad_ca_original + translation

    results = []
    angles = list(range(0, 360, angle_step))

    for rx, ry, rz in product(angles, angles, angles):
        # Create rotation
        rot = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True)

        # Rotate triad around peptide centroid
        triad_ca_rotated = rotate_points(triad_ca_translated, peptide_centroid, rot)

        # Score this placement
        score, assignments = compute_placement_score(triad_ca_rotated, target_ca)

        # Get assigned residues
        assigned_resids = [target_resids[idx] for idx, _ in assignments]
        distances = [d for _, d in assignments]

        results.append({
            "rotation_xyz_deg": [rx, ry, rz],
            "score": score,
            "assigned_residues": assigned_resids,
            "distances": distances,
            "triad_ca_positions": triad_ca_rotated.tolist(),
            "rmsd": np.sqrt(np.mean(np.array(distances) ** 2)),
        })

    # Sort by score (lower is better)
    results.sort(key=lambda x: x["score"])

    return results


def transform_residue_atoms(
    atoms_df: pd.DataFrame,
    original_centroid: np.ndarray,
    target_centroid: np.ndarray,
    rotation: Rotation,
) -> pd.DataFrame:
    """Transform all atoms of a residue."""
    atoms_df = atoms_df.copy()
    coords = atoms_df[["x", "y", "z"]].values

    # Translate so original centroid -> target centroid
    translation = target_centroid - original_centroid
    coords_translated = coords + translation

    # Rotate around target centroid
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
    """Create grafted structure."""
    grafted_df = rhodopsin_df.copy()

    # Remove peptide chain
    grafted_df = grafted_df[grafted_df["auth_chain_id"] != PEPTIDE_CHAIN].copy()

    # Build mapping
    existing_mapping = grafted_df[["auth_seq_id", "label_seq_id"]].drop_duplicates()
    auth_to_label = dict(zip(existing_mapping["auth_seq_id"], existing_mapping["label_seq_id"]))

    for i, (triad_atoms, target_res) in enumerate(zip(transformed_triad, target_residues)):
        # Remove original residue
        mask = ~((grafted_df["auth_seq_id"] == target_res) &
                 (grafted_df["auth_chain_id"] == RHODOPSIN_CHAIN))
        grafted_df = grafted_df[mask].copy()

        # Insert triad residue
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

    # Sort and clean
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
    """Run rotated alignment optimization."""
    print("=" * 70)
    print("ROTATED ALIGNMENT OPTIMIZATION")
    print("Find optimal triad orientation in binding pocket")
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

    print(f"  Trypsin (2AGE): {len(trypsin_df)} atoms")
    print(f"  Rhodopsin (3PQR): {len(rhodopsin_df)} atoms")

    # ==========================================================================
    # Step 2: Extract key positions
    # ==========================================================================
    print("\n[2] Extracting key positions...")

    # Triad CA positions
    triad_ca = []
    triad_atoms_list = []
    for res_info in TRYPSIN_TRIAD:
        atoms = get_residue_atoms(trypsin_df, TRYPSIN_CHAIN, res_info["res_num"])
        triad_atoms_list.append(atoms)
        ca = atoms[atoms["atom_name"] == "CA"][["x", "y", "z"]].values[0]
        triad_ca.append(ca)
    triad_ca = np.array(triad_ca)
    triad_centroid = triad_ca.mean(axis=0)

    print(f"  Triad CA centroid: [{triad_centroid[0]:.1f}, {triad_centroid[1]:.1f}, {triad_centroid[2]:.1f}]")

    # Substrate centroid
    substrate_ca, _ = get_ca_positions(trypsin_df, SUBSTRATE_CHAIN)
    substrate_centroid = substrate_ca.mean(axis=0)
    print(f"  Substrate centroid: [{substrate_centroid[0]:.1f}, {substrate_centroid[1]:.1f}, {substrate_centroid[2]:.1f}]")

    # G-protein peptide centroid
    peptide_ca, _ = get_ca_positions(rhodopsin_df, PEPTIDE_CHAIN)
    peptide_centroid = peptide_ca.mean(axis=0)
    print(f"  Peptide centroid: [{peptide_centroid[0]:.1f}, {peptide_centroid[1]:.1f}, {peptide_centroid[2]:.1f}]")

    # Designable CA positions (targets)
    designable_resids = parse_designable_regions(DESIGNABLE_REGIONS)
    target_ca, target_resids = get_ca_positions(rhodopsin_df, RHODOPSIN_CHAIN, designable_resids)
    print(f"  Designable positions: {len(target_resids)}")

    # ==========================================================================
    # Step 3: Optimize rotation
    # ==========================================================================
    print("\n[3] Optimizing rotation (this may take a moment)...")

    placements = optimize_rotation(
        triad_ca,
        substrate_centroid,
        peptide_centroid,
        target_ca,
        target_resids,
        angle_step=30,  # 30 degree steps = 12^3 = 1728 combinations
    )

    print(f"  Evaluated {len(placements)} orientations")
    print(f"\n  Top 10 placements:")
    print(f"  {'Rank':<5} {'Score':<8} {'RMSD':<8} {'Residues':<20} {'Rotation'}")
    print(f"  {'-'*70}")

    for i, p in enumerate(placements[:10]):
        rot_str = f"({p['rotation_xyz_deg'][0]}, {p['rotation_xyz_deg'][1]}, {p['rotation_xyz_deg'][2]})"
        res_str = ", ".join(str(r) for r in p['assigned_residues'])
        print(f"  {i+1:<5} {p['score']:<8.1f} {p['rmsd']:<8.2f} {res_str:<20} {rot_str}")

    # ==========================================================================
    # Step 4: Create grafted structures for top placements
    # ==========================================================================
    print("\n[4] Creating grafted structures for top 5 placements...")

    grafted_dir = OUTPUT_DIR / "grafted_structures_rotated"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    yaml_dir = OUTPUT_DIR / "boltzgen_configs_rotated"
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

    for rank, placement in enumerate(placements[:5], 1):
        rx, ry, rz = placement["rotation_xyz_deg"]
        rotation = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True)

        # Transform each triad residue
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
        config_name = f"trypsin_rot{rank:02d}_{rx:03d}_{ry:03d}_{rz:03d}"
        output_path = grafted_dir / f"{config_name}.cif"

        grafted_df = create_grafted_structure(
            rhodopsin_chain_a,
            transformed_triad,
            TRYPSIN_TRIAD,
            placement["assigned_residues"],
            output_path,
        )

        print(f"  Rank {rank}: {config_name}")
        print(f"    Residues: {placement['assigned_residues']}, RMSD: {placement['rmsd']:.2f}Å")

        # Generate BoltzGen config
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
            "num_designs": 4,
        }

        yaml_path = yaml_dir / f"{config_name}.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(yaml_config, f, default_flow_style=False)

        config = {
            "rank": rank,
            "name": config_name,
            "rotation_xyz_deg": placement["rotation_xyz_deg"],
            "assigned_residues": placement["assigned_residues"],
            "distances": placement["distances"],
            "rmsd": placement["rmsd"],
            "score": placement["score"],
            "structure_path": str(output_path),
            "yaml_path": str(yaml_path),
        }
        all_configs.append(config)

    # ==========================================================================
    # Step 5: Save summary
    # ==========================================================================
    print("\n[5] Saving summary...")

    summary = {
        "workflow": "rotated_alignment",
        "timestamp": datetime.now().isoformat(),
        "description": "Optimized rotation of triad around substrate centroid",
        "substrate_centroid": substrate_centroid.tolist(),
        "peptide_centroid": peptide_centroid.tolist(),
        "total_orientations_evaluated": len(placements),
        "angle_step_degrees": 30,
        "top_placements": all_configs,
        "all_placements_top20": placements[:20],
    }

    summary_path = OUTPUT_DIR / "rotated_alignment_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Summary: {summary_path}")

    # ==========================================================================
    # Done
    # ==========================================================================
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print()
    print(f"Created {len(all_configs)} grafted structures with optimized rotations")
    print(f"Grafted structures: {grafted_dir}")
    print(f"BoltzGen configs: {yaml_dir}")
    print()
    print("Best placement:")
    best = all_configs[0]
    print(f"  Rotation: {best['rotation_xyz_deg']}")
    print(f"  Residues: {best['assigned_residues']}")
    print(f"  RMSD: {best['rmsd']:.2f}Å")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
