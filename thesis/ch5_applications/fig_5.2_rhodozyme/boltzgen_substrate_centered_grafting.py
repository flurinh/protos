#!/usr/bin/env python3
"""BoltzGen Substrate-Centered Catalytic Site Grafting.

Approach:
1. Find substrate position from where G-protein peptide (chain B) binds in 3PQR
2. Place catalytic triad around substrate with correct geometry (sidechains pointing inward)
3. Find nearest CA atoms in rhodopsin for each triad residue
4. OVERWRITE those positions with full catalytic residue atoms
5. BoltzGen designs EVERYTHING except the fixed catalytic residues

The resulting structure is NOT a valid protein - but that's fine because
BoltzGen will redesign the fold around the fixed catalytic anchor points.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from datetime import datetime
from dataclasses import dataclass
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize

import numpy as np
import pandas as pd
import yaml

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

RHODOPSIN_PDB = "3pqr"
RHODOPSIN_CHAIN = "A"
PEPTIDE_CHAIN = "B"  # G-protein peptide - this is where substrate should go

# Residues that CAN be designed (cytoplasmic regions)
DESIGNABLE_REGIONS = "58..76,134..157,222..256,306..326"

# Keep retinal
RETINAL_CODES = ["RET", "LYR"]

# Reference enzymes
REFERENCE_ENZYMES = {
    "trypsin": {
        "pdb": "1S81",
        "chain": "A",
        "triad": [
            {"res_num": 195, "res_name": "SER", "role": "nucleophile"},
            {"res_num": 57, "res_name": "HIS", "role": "base"},
            {"res_num": 102, "res_name": "ASP", "role": "electrostatic"},
        ],
        "substrate_smiles": "NC(=N)c1ccccc1",  # Benzamidine
    },
    "subtilisin": {
        "pdb": "1SBC",
        "chain": "A",
        "triad": [
            {"res_num": 221, "res_name": "SER", "role": "nucleophile"},
            {"res_num": 64, "res_name": "HIS", "role": "base"},
            {"res_num": 32, "res_name": "ASP", "role": "electrostatic"},
        ],
        "substrate_smiles": "CC(C)CC(N)C(=O)O",  # Leucine
    },
    "papain": {
        "pdb": "1PPN",
        "chain": "A",
        "triad": [
            {"res_num": 25, "res_name": "CYS", "role": "nucleophile"},
            {"res_num": 159, "res_name": "HIS", "role": "base"},
            {"res_num": 175, "res_name": "ASN", "role": "electrostatic"},
        ],
        "substrate_smiles": "NC(Cc1ccccc1)C(=O)O",  # Phenylalanine
    },
}

RETINAL_SMILES = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C"


@dataclass
class CatalyticTriad:
    """Full catalytic triad with all atom coordinates."""
    enzyme_name: str
    residue_names: List[str]  # e.g., ["SER", "HIS", "ASP"]
    residue_roles: List[str]  # e.g., ["nucleophile", "base", "electrostatic"]
    atoms_per_residue: List[pd.DataFrame]  # Full atom data
    ca_coords: np.ndarray  # (3, 3) CA positions
    centroid: np.ndarray  # Centroid of CA positions
    substrate_smiles: str


def extract_catalytic_triad(
    enzyme_df: pd.DataFrame,
    enzyme_config: Dict,
    chain_id: str,
) -> CatalyticTriad:
    """Extract full catalytic triad from enzyme structure."""
    atoms_per_residue = []
    ca_coords = []
    residue_names = []
    residue_roles = []

    for res_info in enzyme_config["triad"]:
        res_num = res_info["res_num"]
        res_name = res_info["res_name"]
        role = res_info["role"]

        mask = (enzyme_df["auth_seq_id"] == res_num) & (enzyme_df["auth_chain_id"] == chain_id)
        res_atoms = enzyme_df[mask].copy()

        if len(res_atoms) == 0:
            raise ValueError(f"Residue {res_name}{res_num} not found")

        atoms_per_residue.append(res_atoms)
        residue_names.append(res_name)
        residue_roles.append(role)

        ca = res_atoms[res_atoms["atom_name"] == "CA"]
        if len(ca) == 0:
            raise ValueError(f"CA not found for {res_name}{res_num}")
        ca_coords.append(ca[["x", "y", "z"]].values[0])

    return CatalyticTriad(
        enzyme_name=enzyme_config.get("name", "unknown"),
        residue_names=residue_names,
        residue_roles=residue_roles,
        atoms_per_residue=atoms_per_residue,
        ca_coords=np.array(ca_coords),
        centroid=np.array(ca_coords).mean(axis=0),
        substrate_smiles=enzyme_config["substrate_smiles"],
    )


def get_peptide_binding_center(rhodopsin_df: pd.DataFrame, peptide_chain: str) -> np.ndarray:
    """Get the center of where the G-protein peptide binds (substrate position)."""
    peptide_atoms = rhodopsin_df[rhodopsin_df["auth_chain_id"] == peptide_chain]

    if len(peptide_atoms) == 0:
        raise ValueError(f"No atoms found for peptide chain {peptide_chain}")

    # Get CA atoms of the peptide
    peptide_ca = peptide_atoms[peptide_atoms["atom_name"] == "CA"]
    if len(peptide_ca) == 0:
        # Fall back to all atoms
        coords = peptide_atoms[["x", "y", "z"]].values
    else:
        coords = peptide_ca[["x", "y", "z"]].values

    center = coords.mean(axis=0)
    print(f"  Peptide binding center: [{center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}]")
    print(f"  Peptide residues: {len(peptide_ca)} CA atoms")

    return center


def get_designable_ca_positions(
    rhodopsin_df: pd.DataFrame,
    designable_regions: str,
    chain_id: str,
) -> Tuple[np.ndarray, List[int]]:
    """Get CA positions of designable residues."""
    # Parse designable regions
    residues = []
    for part in designable_regions.split(","):
        if ".." in part:
            start, end = map(int, part.split(".."))
            residues.extend(range(start, end + 1))
        else:
            residues.append(int(part))

    mask = (
        (rhodopsin_df["auth_chain_id"] == chain_id) &
        (rhodopsin_df["auth_seq_id"].isin(residues)) &
        (rhodopsin_df["atom_name"] == "CA")
    )

    ca_df = rhodopsin_df[mask][["auth_seq_id", "x", "y", "z"]].copy()
    ca_df = ca_df.sort_values("auth_seq_id")

    coords = ca_df[["x", "y", "z"]].values
    resids = ca_df["auth_seq_id"].values.tolist()

    return coords, resids


def optimize_triad_placement(
    triad: CatalyticTriad,
    substrate_center: np.ndarray,
    target_ca_coords: np.ndarray,
    target_resids: List[int],
) -> Dict:
    """
    Find optimal rotation to place triad around substrate center,
    matching triad CA positions to nearest target CA positions.

    Returns placement info including rotation, target residues, and RMSD.
    """
    # Center triad at origin
    triad_centered = triad.ca_coords - triad.centroid

    def objective(params):
        """Minimize distance from rotated triad CAs to nearest target CAs."""
        rx, ry, rz = params
        R = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True).as_matrix()

        # Rotate and translate triad
        triad_transformed = (R @ triad_centered.T).T + substrate_center

        # For each triad CA, find distance to nearest target CA
        total_dist = 0
        for triad_ca in triad_transformed:
            dists = np.linalg.norm(target_ca_coords - triad_ca, axis=1)
            total_dist += np.min(dists)

        return total_dist

    # Grid search for good starting point
    best_result = None
    best_value = float('inf')

    for rx in range(0, 360, 45):
        for ry in range(0, 360, 45):
            for rz in range(0, 360, 45):
                result = minimize(
                    objective,
                    [rx, ry, rz],
                    method='L-BFGS-B',
                    bounds=[(-360, 360), (-360, 360), (-360, 360)]
                )
                if result.fun < best_value:
                    best_value = result.fun
                    best_result = result

    # Get final placement
    rx, ry, rz = best_result.x
    R = Rotation.from_euler('xyz', [rx, ry, rz], degrees=True).as_matrix()
    triad_final = (R @ triad_centered.T).T + substrate_center

    # Find nearest target CA for each triad residue
    assigned_resids = []
    assigned_distances = []
    used_indices = set()

    for triad_ca in triad_final:
        dists = np.linalg.norm(target_ca_coords - triad_ca, axis=1)
        # Find nearest unused target
        for idx in np.argsort(dists):
            if idx not in used_indices:
                assigned_resids.append(target_resids[idx])
                assigned_distances.append(dists[idx])
                used_indices.add(idx)
                break

    rmsd = np.sqrt(np.mean(np.array(assigned_distances) ** 2))

    return {
        "rotation_xyz_deg": [rx, ry, rz],
        "target_residues": assigned_resids,
        "distances_to_target": assigned_distances,
        "rmsd": rmsd,
        "triad_ca_positions": triad_final.tolist(),
        "substrate_center": substrate_center.tolist(),
    }


def transform_triad_atoms(
    triad: CatalyticTriad,
    substrate_center: np.ndarray,
    rotation_xyz_deg: List[float],
) -> List[pd.DataFrame]:
    """Transform all triad atoms to the target position."""
    R = Rotation.from_euler('xyz', rotation_xyz_deg, degrees=True).as_matrix()

    transformed_residues = []

    for res_atoms in triad.atoms_per_residue:
        res_atoms = res_atoms.copy()
        coords = res_atoms[["x", "y", "z"]].values

        # Center on triad centroid, rotate, translate to substrate center
        coords_centered = coords - triad.centroid
        coords_rotated = (R @ coords_centered.T).T
        coords_final = coords_rotated + substrate_center

        res_atoms["x"] = coords_final[:, 0]
        res_atoms["y"] = coords_final[:, 1]
        res_atoms["z"] = coords_final[:, 2]

        transformed_residues.append(res_atoms)

    return transformed_residues


def create_grafted_structure(
    rhodopsin_df: pd.DataFrame,
    transformed_triad: List[pd.DataFrame],
    triad: CatalyticTriad,
    target_residues: List[int],
    output_path: Path,
) -> pd.DataFrame:
    """
    Create grafted structure by OVERWRITING target residues with triad atoms.

    This creates an "invalid" protein - the backbone won't connect properly,
    but that's fine because BoltzGen will redesign everything.
    """
    grafted_df = rhodopsin_df.copy()

    # Remove peptide chain (we're replacing it with substrate)
    grafted_df = grafted_df[grafted_df["auth_chain_id"] != PEPTIDE_CHAIN].copy()

    # Build mapping
    existing_mapping = grafted_df[["auth_seq_id", "label_seq_id"]].drop_duplicates()
    auth_to_label = dict(zip(existing_mapping["auth_seq_id"], existing_mapping["label_seq_id"]))

    for i, (triad_atoms, target_res) in enumerate(zip(transformed_triad, target_residues)):
        # REMOVE original residue at target position
        mask = ~((grafted_df["auth_seq_id"] == target_res) &
                 (grafted_df["auth_chain_id"] == RHODOPSIN_CHAIN))
        grafted_df = grafted_df[mask].copy()

        # INSERT triad residue at that position
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

        # Update residue name
        triad_atoms["res_name"] = triad.residue_names[i]

        grafted_df = pd.concat([grafted_df, triad_atoms], ignore_index=True)

    # Sort
    grafted_df = grafted_df.sort_values(["auth_chain_id", "auth_seq_id", "atom_id"])

    # Clean: keep only amino acids + retinal
    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    ]
    keep_codes = amino_acids + RETINAL_CODES
    grafted_df = grafted_df[grafted_df["res_name"].isin(keep_codes)].copy()

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_cif_file(str(output_path), grafted_df, force_overwrite=True)

    return grafted_df


def main() -> int:
    """Run substrate-centered catalytic site grafting."""
    print("=" * 70)
    print("SUBSTRATE-CENTERED CATALYTIC SITE GRAFTING")
    print("Place triad around G-protein binding site (substrate position)")
    print("=" * 70)

    # Initialize
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))

    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)

    # ==========================================================================
    # Step 1: Load rhodopsin with peptide
    # ==========================================================================
    print("\n[1] Loading rhodopsin (3PQR) with G-protein peptide...")

    loader.download_and_register(RHODOPSIN_PDB, name=RHODOPSIN_PDB)
    rhodopsin_df = sp.load_entity(RHODOPSIN_PDB)

    if rhodopsin_df is None:
        raise RuntimeError("Failed to load rhodopsin")

    rhodopsin_df = rhodopsin_df.reset_index()
    print(f"  Loaded {len(rhodopsin_df)} atoms")

    # Get peptide binding center (this is where substrate should go)
    substrate_center = get_peptide_binding_center(rhodopsin_df, PEPTIDE_CHAIN)

    # Get designable CA positions
    target_ca_coords, target_resids = get_designable_ca_positions(
        rhodopsin_df, DESIGNABLE_REGIONS, RHODOPSIN_CHAIN
    )
    print(f"  Designable residues: {len(target_resids)} positions")

    # Filter to chain A for grafting (remove peptide later)
    rhodopsin_chain_a = rhodopsin_df[rhodopsin_df["auth_chain_id"].isin([RHODOPSIN_CHAIN])].copy()

    # Renumber
    unique_residues = rhodopsin_chain_a[["auth_seq_id", "res_name"]].drop_duplicates()
    unique_residues = unique_residues.sort_values("auth_seq_id")
    auth_to_label = {auth_seq: i + 1 for i, auth_seq in enumerate(unique_residues["auth_seq_id"].values)}
    rhodopsin_chain_a["label_seq_id"] = rhodopsin_chain_a["auth_seq_id"].map(auth_to_label)
    rhodopsin_chain_a["label_chain_id"] = rhodopsin_chain_a["auth_chain_id"]
    rhodopsin_chain_a["entity_id"] = "1"

    # ==========================================================================
    # Step 2: Extract triads from enzymes
    # ==========================================================================
    print("\n[2] Extracting catalytic triads...")

    triads = {}

    for enzyme_name, config in REFERENCE_ENZYMES.items():
        config["name"] = enzyme_name
        pdb_id = config["pdb"].lower()

        print(f"\n  {enzyme_name} ({pdb_id}):")

        loader.download_and_register(pdb_id, name=f"enzyme_{enzyme_name}")
        enzyme_df = sp.load_entity(f"enzyme_{enzyme_name}")

        if enzyme_df is None:
            print(f"    Failed to load")
            continue

        enzyme_df = enzyme_df.reset_index()

        try:
            triad = extract_catalytic_triad(enzyme_df, config, config["chain"])
            triads[enzyme_name] = triad

            # Show geometry
            d01 = np.linalg.norm(triad.ca_coords[0] - triad.ca_coords[1])
            d12 = np.linalg.norm(triad.ca_coords[1] - triad.ca_coords[2])
            d02 = np.linalg.norm(triad.ca_coords[0] - triad.ca_coords[2])

            print(f"    Triad: {'-'.join(triad.residue_names)}")
            print(f"    CA-CA: {d01:.1f}A, {d12:.1f}A, {d02:.1f}A (max: {max(d01,d12,d02):.1f}A)")

        except Exception as e:
            print(f"    Failed: {e}")

    # ==========================================================================
    # Step 3: Optimize placement for each enzyme
    # ==========================================================================
    print("\n[3] Optimizing triad placements around substrate center...")

    placements = {}

    for enzyme_name, triad in triads.items():
        print(f"\n  {enzyme_name}:")

        placement = optimize_triad_placement(
            triad, substrate_center, target_ca_coords, target_resids
        )
        placements[enzyme_name] = placement

        print(f"    Target residues: {placement['target_residues']}")
        print(f"    Distances to target CAs: {[f'{d:.1f}A' for d in placement['distances_to_target']]}")
        print(f"    RMSD: {placement['rmsd']:.2f}A")

    # ==========================================================================
    # Step 4: Create grafted structures
    # ==========================================================================
    print("\n[4] Creating grafted structures...")

    grafted_dir = OUTPUT_DIR / "grafted_structures_substrate_centered"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    all_configs = []

    for enzyme_name, triad in triads.items():
        placement = placements[enzyme_name]

        # Transform triad atoms
        transformed = transform_triad_atoms(
            triad,
            np.array(placement["substrate_center"]),
            placement["rotation_xyz_deg"]
        )

        # Create structure
        config_name = f"{enzyme_name}_substrate_centered"
        output_path = grafted_dir / f"{config_name}.cif"

        grafted_df = create_grafted_structure(
            rhodopsin_chain_a,
            transformed,
            triad,
            placement["target_residues"],
            output_path,
        )

        print(f"  {enzyme_name}: {output_path.name}")
        print(f"    Catalytic residues at: {placement['target_residues']}")

        config = {
            "name": config_name,
            "enzyme": enzyme_name,
            "target_residues": placement["target_residues"],
            "rmsd": placement["rmsd"],
            "triad_residues": triad.residue_names,
            "structure_path": str(output_path),
            "substrate_smiles": triad.substrate_smiles,
            "substrate_center": placement["substrate_center"],
        }
        all_configs.append(config)

    # ==========================================================================
    # Step 5: Generate BoltzGen configs
    # ==========================================================================
    print("\n[5] Generating BoltzGen configurations...")

    yaml_dir = OUTPUT_DIR / "boltzgen_configs_substrate_centered"
    yaml_dir.mkdir(parents=True, exist_ok=True)

    for config in all_configs:
        # EVERYTHING is designable EXCEPT the catalytic residues
        triad_residues_str = ",".join(str(r) for r in config["target_residues"])

        yaml_config = {
            "entities": [
                {
                    "file": {
                        "path": config["structure_path"],
                        "include": [
                            {"chain": {"id": RHODOPSIN_CHAIN}}
                        ],
                        # Design EVERYTHING in cytoplasmic region
                        "design": [
                            {"chain": {"id": RHODOPSIN_CHAIN, "res_index": DESIGNABLE_REGIONS}}
                        ],
                        # EXCEPT the catalytic triad - these are FIXED
                        "not_design": [
                            {"chain": {"id": RHODOPSIN_CHAIN, "res_index": triad_residues_str}}
                        ],
                        # Hide designable regions for inpainting
                        "structure_groups": [
                            {"group": {"id": "all", "visibility": 1}},
                            {"group": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": DESIGNABLE_REGIONS,
                                "visibility": 0  # Hide - will be designed
                            }},
                            {"group": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": triad_residues_str,
                                "visibility": 1  # Show - fixed anchor points
                            }},
                        ],
                        # Binding site = catalytic triad
                        "binding_types": [
                            {"chain": {"id": RHODOPSIN_CHAIN, "binding": triad_residues_str}}
                        ],
                    }
                },
                # Substrate
                {"ligand": {"id": "SUB", "smiles": config["substrate_smiles"]}},
                # Retinal
                {"ligand": {"id": "RET", "smiles": RETINAL_SMILES}},
            ],
            "num_designs": 4,  # More designs since this is a harder problem
        }

        yaml_path = yaml_dir / f"{config['name']}.yaml"
        with open(yaml_path, "w") as f:
            yaml.dump(yaml_config, f, default_flow_style=False)

        config["yaml_path"] = str(yaml_path)
        print(f"  {config['name']}.yaml")

    # ==========================================================================
    # Step 6: Save summary
    # ==========================================================================
    print("\n[6] Saving summary...")

    summary = {
        "workflow": "substrate_centered_grafting",
        "timestamp": datetime.now().isoformat(),
        "description": "Catalytic triad placed around G-protein peptide binding site",
        "substrate_center": substrate_center.tolist(),
        "scaffold": {"pdb": RHODOPSIN_PDB, "chain": RHODOPSIN_CHAIN},
        "designable_regions": DESIGNABLE_REGIONS,
        "configs": all_configs,
    }

    summary_path = OUTPUT_DIR / "substrate_centered_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    # ==========================================================================
    # Done
    # ==========================================================================
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"Grafted structures: {grafted_dir}")
    print(f"BoltzGen configs: {yaml_dir}")
    print(f"Summary: {summary_path}")
    print()
    print("KEY CONCEPT:")
    print("  - Substrate position = where G-protein peptide binds")
    print("  - Catalytic triad placed AROUND this position")
    print("  - Triad residues OVERWRITE nearest CA positions")
    print("  - BoltzGen designs EVERYTHING else (protein fold)")
    print()
    print("The grafted structures are NOT valid proteins!")
    print("BoltzGen will figure out how to fold around the fixed catalytic anchors.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
