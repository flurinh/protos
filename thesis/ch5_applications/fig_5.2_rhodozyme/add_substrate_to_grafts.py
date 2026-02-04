#!/usr/bin/env python3
"""Add transformed substrate as chain B to grafted structures.

This script:
1. Loads existing grafted structures with catalytic triads
2. Applies the same transformation (translation + rotation) to the substrate from 2AGE
3. Adds the transformed substrate as chain B
4. Updates BoltzGen configs to include chain B with fixed coordinates
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import List

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

import protos
from protos.processing.structure import StructureProcessor
from protos.io.formats.cif_utils import write_cif_file, read_cif_file
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Configuration
# =============================================================================

TRYPSIN_PDB = "2age"
TRYPSIN_CHAIN = "X"
SUBSTRATE_CHAIN = "A"

RHODOPSIN_PDB = "3pqr"
PEPTIDE_CHAIN = "B"

DESIGNABLE_REGIONS = "58..76,134..157,222..256,306..326"
RETINAL_SMILES = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C"


def get_ca_positions(df: pd.DataFrame, chain: str) -> np.ndarray:
    """Get CA positions for a chain."""
    mask = (df["auth_chain_id"] == chain) & (df["atom_name"] == "CA")
    ca_df = df[mask].sort_values("auth_seq_id")
    return ca_df[["x", "y", "z"]].values


def rotate_points(points: np.ndarray, center: np.ndarray, rotation: Rotation) -> np.ndarray:
    """Rotate points around a center."""
    centered = points - center
    rotated = rotation.apply(centered)
    return rotated + center


def transform_substrate_atoms(
    substrate_df: pd.DataFrame,
    substrate_centroid: np.ndarray,
    peptide_centroid: np.ndarray,
    rotation: Rotation,
) -> pd.DataFrame:
    """Apply translation + rotation to substrate atoms."""
    df = substrate_df.copy()
    coords = df[["x", "y", "z"]].values

    # Translation: move substrate centroid to peptide centroid
    translation = peptide_centroid - substrate_centroid
    coords_translated = coords + translation

    # Rotation: rotate around peptide centroid
    coords_rotated = rotate_points(coords_translated, peptide_centroid, rotation)

    df["x"] = coords_rotated[:, 0]
    df["y"] = coords_rotated[:, 1]
    df["z"] = coords_rotated[:, 2]

    return df


def main() -> int:
    print("=" * 70)
    print("ADDING SUBSTRATE TO GRAFTED STRUCTURES")
    print("=" * 70)

    # Initialize
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))

    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)

    # ==========================================================================
    # Step 1: Load source structures
    # ==========================================================================
    print("\n[1] Loading source structures...")

    loader.download_and_register(TRYPSIN_PDB, name="trypsin_2age")
    loader.download_and_register(RHODOPSIN_PDB, name=RHODOPSIN_PDB)

    trypsin_df = sp.load_entity("trypsin_2age").reset_index()
    rhodopsin_df = sp.load_entity(RHODOPSIN_PDB).reset_index()

    # Get substrate atoms from 2AGE chain A (exclude HOH water molecules)
    substrate_df = trypsin_df[
        (trypsin_df["auth_chain_id"] == SUBSTRATE_CHAIN) &
        (trypsin_df["res_name"] != "HOH")
    ].copy()
    print(f"  Substrate atoms: {len(substrate_df)}")
    print(f"  Substrate residues: {substrate_df['auth_seq_id'].unique().tolist()}")
    print(f"  Substrate res names: {substrate_df['res_name'].unique().tolist()}")

    # Get centroids
    substrate_ca = get_ca_positions(trypsin_df, SUBSTRATE_CHAIN)
    substrate_centroid = substrate_ca.mean(axis=0)

    peptide_ca = get_ca_positions(rhodopsin_df, PEPTIDE_CHAIN)
    peptide_centroid = peptide_ca.mean(axis=0)

    print(f"  Substrate centroid: [{substrate_centroid[0]:.1f}, {substrate_centroid[1]:.1f}, {substrate_centroid[2]:.1f}]")
    print(f"  Peptide centroid: [{peptide_centroid[0]:.1f}, {peptide_centroid[1]:.1f}, {peptide_centroid[2]:.1f}]")

    # ==========================================================================
    # Step 2: Load refinement summary
    # ==========================================================================
    print("\n[2] Loading refinement summary...")

    summary_path = OUTPUT_DIR / "refined_alignment_summary.json"
    with open(summary_path) as f:
        summary = json.load(f)

    placements = summary["placements"]
    print(f"  Found {len(placements)} refined placements")

    # ==========================================================================
    # Step 3: Add substrate to each grafted structure
    # ==========================================================================
    print("\n[3] Adding substrate to grafted structures...")

    grafted_dir = OUTPUT_DIR / "grafted_structures_with_substrate"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    yaml_dir = OUTPUT_DIR / "boltzgen_configs_with_substrate"
    yaml_dir.mkdir(parents=True, exist_ok=True)

    updated_configs = []

    for placement in placements:
        rank = placement["rank"]
        rot_angles = placement["rotation_xyz_deg"]
        assigned_residues = placement["assigned_residues"]
        rmsd = placement["rmsd"]

        print(f"\n  Placement {rank}: residues {assigned_residues}, RMSD {rmsd:.3f}Å")
        print(f"    Rotation: [{rot_angles[0]:.1f}, {rot_angles[1]:.1f}, {rot_angles[2]:.1f}]")

        # Load existing grafted structure directly from CIF
        original_path = Path(placement["structure_path"])
        grafted_df = read_cif_file(str(original_path))
        print(f"    Loaded {len(grafted_df)} atoms from original graft")

        # Transform substrate
        rotation = Rotation.from_euler('xyz', rot_angles, degrees=True)
        transformed_substrate = transform_substrate_atoms(
            substrate_df,
            substrate_centroid,
            peptide_centroid,
            rotation,
        )

        # Assign chain B to substrate
        transformed_substrate["auth_chain_id"] = "B"
        transformed_substrate["label_chain_id"] = "B"
        transformed_substrate["entity_id"] = "2"

        # Renumber substrate residues starting from 1
        unique_sub_res = sorted(transformed_substrate["auth_seq_id"].unique())
        res_renumber = {old: i + 1 for i, old in enumerate(unique_sub_res)}
        transformed_substrate["auth_seq_id"] = transformed_substrate["auth_seq_id"].map(res_renumber)
        transformed_substrate["label_seq_id"] = transformed_substrate["auth_seq_id"]
        transformed_substrate["gen_seq_id"] = transformed_substrate["auth_seq_id"]

        print(f"    Transformed substrate: {len(transformed_substrate)} atoms")
        print(f"    Substrate residues (renumbered): {sorted(transformed_substrate['auth_seq_id'].unique())}")

        # Combine grafted structure with substrate
        combined_df = pd.concat([grafted_df, transformed_substrate], ignore_index=True)
        combined_df = combined_df.sort_values(["auth_chain_id", "auth_seq_id", "atom_id"])

        # Write combined structure
        config_name = f"trypsin_refined{rank:02d}_with_substrate"
        output_path = grafted_dir / f"{config_name}.cif"
        write_cif_file(str(output_path), combined_df, force_overwrite=True)
        print(f"    Wrote: {output_path.name}")

        # Create BoltzGen config with chain B
        triad_str = ",".join(str(r) for r in assigned_residues)

        # Get substrate residue range for chain B
        sub_residues = sorted(transformed_substrate["auth_seq_id"].unique())
        sub_res_str = f"1..{max(sub_residues)}"

        yaml_config = {
            "entities": [
                {
                    "file": {
                        "path": str(output_path),
                        "include": [
                            {"chain": {"id": "A"}},
                            {"chain": {"id": "B"}},  # Include substrate chain
                        ],
                        "design": [
                            {"chain": {"id": "A", "res_index": DESIGNABLE_REGIONS}},
                        ],
                        "not_design": [
                            {"chain": {"id": "A", "res_index": triad_str}},
                            {"chain": {"id": "B"}},  # Fix substrate sequence
                        ],
                        "structure_groups": [
                            {"group": {"id": "all", "visibility": 1}},
                            {"group": {"id": "A", "res_index": DESIGNABLE_REGIONS, "visibility": 0}},
                            {"group": {"id": "A", "res_index": triad_str, "visibility": 1}},
                            {"group": {"id": "B", "visibility": 1}},  # Fix substrate coords
                        ],
                        "binding_types": [
                            {"chain": {"id": "A", "binding": triad_str}},
                        ],
                    }
                },
                {"ligand": {"id": "RET", "smiles": RETINAL_SMILES}},
            ],
            "num_designs": 8,
        }

        yaml_path = yaml_dir / f"{config_name}.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(yaml_config, f, default_flow_style=False)
        print(f"    Config: {yaml_path.name}")

        updated_configs.append({
            "rank": rank,
            "name": config_name,
            "rotation_xyz_deg": rot_angles,
            "assigned_residues": assigned_residues,
            "rmsd": rmsd,
            "structure_path": str(output_path),
            "yaml_path": str(yaml_path),
            "substrate_residues": [int(r) for r in sub_residues],
        })

    # ==========================================================================
    # Step 4: Save updated summary
    # ==========================================================================
    print("\n[4] Saving summary...")

    updated_summary = {
        "workflow": "grafts_with_substrate",
        "source_summary": str(summary_path),
        "substrate_source": f"2AGE chain {SUBSTRATE_CHAIN}",
        "placements": updated_configs,
    }

    updated_summary_path = OUTPUT_DIR / "grafts_with_substrate_summary.json"
    with open(updated_summary_path, "w") as f:
        json.dump(updated_summary, f, indent=2)

    # ==========================================================================
    # Done
    # ==========================================================================
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print()
    print(f"Grafted structures: {grafted_dir}")
    print(f"BoltzGen configs: {yaml_dir}")
    print()
    print("Key changes in BoltzGen configs:")
    print("  - Chain B (substrate) added with fixed coordinates (visibility: 1)")
    print("  - Chain B marked as not_design (fixed sequence)")
    print("  - Designable regions will be generated around both triad AND substrate")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
