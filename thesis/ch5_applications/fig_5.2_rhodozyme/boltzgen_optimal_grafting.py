#!/usr/bin/env python3
"""BoltzGen Optimal Triad Grafting - Using Geometrically Optimized Placements.

This workflow uses the TM3-TM5 interface placements from the geometric
optimization, which achieve near-perfect fits (RMSD < 0.5A) between the
catalytic triad geometry and the TM helix cytoplasmic segments.

Key insight: The trypsin catalytic triad spans ~10.3A (CA-CA distances),
which closely matches the TM3-TM5 cytoplasmic distance (~10.6A).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
from dataclasses import dataclass

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

# Designable regions (TM cytoplasmic ends and loops)
DESIGNABLE_REGIONS = "58..76,134..157,222..256,306..326"

# Keep retinal - essential for rhodopsin function
RETINAL_CODES = ["RET", "LYR"]

# Reference enzymes with catalytic triads
REFERENCE_ENZYMES = {
    "trypsin": {
        "id": "TRP",
        "pdb": "1S81",
        "chain": "A",
        "description": "Bovine trypsin (serine protease)",
        "triad_residues": [
            {"res_num": 195, "res_name": "SER", "role": "nucleophile"},
            {"res_num": 57, "res_name": "HIS", "role": "base"},
            {"res_num": 102, "res_name": "ASP", "role": "electrostatic"},
        ],
        "substrate_smiles": "NC(=N)c1ccccc1",
        "substrate_name": "benzamidine",
    },
    "subtilisin": {
        "id": "SUB",
        "pdb": "1SBC",
        "chain": "A",
        "description": "Subtilisin (serine protease, different fold)",
        "triad_residues": [
            {"res_num": 221, "res_name": "SER", "role": "nucleophile"},
            {"res_num": 64, "res_name": "HIS", "role": "base"},
            {"res_num": 32, "res_name": "ASP", "role": "electrostatic"},
        ],
        "substrate_smiles": "CC(C)CC(N)C(=O)O",
        "substrate_name": "leucine",
    },
    "papain": {
        "id": "PAP",
        "pdb": "1PPN",  # Higher resolution than 9PAP
        "chain": "A",
        "description": "Papain (cysteine protease)",
        "triad_residues": [
            {"res_num": 25, "res_name": "CYS", "role": "nucleophile"},
            {"res_num": 159, "res_name": "HIS", "role": "base"},
            {"res_num": 175, "res_name": "ASN", "role": "electrostatic"},
        ],
        "substrate_smiles": "NC(Cc1ccccc1)C(=O)O",
        "substrate_name": "phenylalanine",
    },
}

# Retinal SMILES for BoltzGen config
RETINAL_SMILES = "CC1=C(C(CCC1)(C)C)/C=C/C(=C/C=C/C(=C/C=O)/C)/C"


@dataclass
class TriadGeometry:
    """Stores the geometry of a catalytic triad."""
    enzyme_key: str
    enzyme_id: str
    residue_names: List[str]
    residue_roles: List[str]
    atoms_per_residue: List[pd.DataFrame]
    ca_coords: np.ndarray
    centroid: np.ndarray
    substrate_smiles: str
    substrate_name: str


def rotation_matrix_from_euler(rx: float, ry: float, rz: float) -> np.ndarray:
    """Create rotation matrix from Euler angles (in degrees)."""
    rx, ry, rz = np.radians([rx, ry, rz])

    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rx), -np.sin(rx)],
        [0, np.sin(rx), np.cos(rx)]
    ])

    Ry = np.array([
        [np.cos(ry), 0, np.sin(ry)],
        [0, 1, 0],
        [-np.sin(ry), 0, np.cos(ry)]
    ])

    Rz = np.array([
        [np.cos(rz), -np.sin(rz), 0],
        [np.sin(rz), np.cos(rz), 0],
        [0, 0, 1]
    ])

    return Rz @ Ry @ Rx


def extract_triad_geometry(
    enzyme_df: pd.DataFrame,
    enzyme_config: Dict,
    chain_id: str,
) -> TriadGeometry:
    """Extract the catalytic triad geometry from an enzyme structure."""
    triad_residues = enzyme_config["triad_residues"]

    atoms_per_residue = []
    ca_coords = []
    residue_names = []
    residue_roles = []

    for res_info in triad_residues:
        res_num = res_info["res_num"]
        res_name = res_info["res_name"]
        role = res_info["role"]

        mask = (enzyme_df["auth_seq_id"] == res_num) & (enzyme_df["auth_chain_id"] == chain_id)
        res_atoms = enzyme_df[mask].copy()

        if len(res_atoms) == 0:
            raise ValueError(f"Residue {res_name}{res_num} not found in chain {chain_id}")

        atoms_per_residue.append(res_atoms)
        residue_names.append(res_name)
        residue_roles.append(role)

        ca = res_atoms[res_atoms["atom_name"] == "CA"]
        if len(ca) == 0:
            raise ValueError(f"CA not found for {res_name}{res_num}")
        ca_coords.append(ca[["x", "y", "z"]].values[0])

    ca_coords = np.array(ca_coords)
    centroid = ca_coords.mean(axis=0)

    return TriadGeometry(
        enzyme_key=enzyme_config.get("key", "unknown"),
        enzyme_id=enzyme_config["id"],
        residue_names=residue_names,
        residue_roles=residue_roles,
        atoms_per_residue=atoms_per_residue,
        ca_coords=ca_coords,
        centroid=centroid,
        substrate_smiles=enzyme_config["substrate_smiles"],
        substrate_name=enzyme_config["substrate_name"],
    )


def transform_triad_to_placement(
    triad: TriadGeometry,
    target_centroid: np.ndarray,
    rotation_angles: Tuple[float, float, float],
) -> List[pd.DataFrame]:
    """Transform triad atoms to the target centroid with given rotation.

    1. Translate triad so centroid is at origin
    2. Apply rotation
    3. Translate to target centroid
    """
    R = rotation_matrix_from_euler(*rotation_angles)

    transformed_residues = []

    for res_atoms in triad.atoms_per_residue:
        res_atoms = res_atoms.copy()
        coords = res_atoms[["x", "y", "z"]].values

        # Center on triad centroid
        coords_centered = coords - triad.centroid

        # Rotate
        coords_rotated = (R @ coords_centered.T).T

        # Translate to target
        coords_final = coords_rotated + target_centroid

        res_atoms["x"] = coords_final[:, 0]
        res_atoms["y"] = coords_final[:, 1]
        res_atoms["z"] = coords_final[:, 2]

        transformed_residues.append(res_atoms)

    return transformed_residues


def clean_structure(df: pd.DataFrame, keep_codes: List[str]) -> pd.DataFrame:
    """Remove unwanted ligands/molecules, keeping only specified residue types."""
    # Standard amino acids
    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
    ]

    all_keep = amino_acids + keep_codes
    mask = df["res_name"].isin(all_keep)
    return df[mask].copy()


def create_grafted_structure(
    rhodopsin_df: pd.DataFrame,
    triad_residues: List[pd.DataFrame],
    triad_info: TriadGeometry,
    insertion_points: List[int],
    output_path: Path,
) -> pd.DataFrame:
    """Create a grafted structure with triad residues at specified positions."""
    grafted_df = rhodopsin_df.copy()

    # Build auth_seq_id -> label_seq_id mapping
    existing_mapping = grafted_df[["auth_seq_id", "label_seq_id"]].drop_duplicates()
    auth_to_label = dict(zip(existing_mapping["auth_seq_id"], existing_mapping["label_seq_id"]))

    for i, (triad_atoms, insert_pos) in enumerate(zip(triad_residues, insertion_points)):
        # Remove original residue at insertion point
        mask = ~((grafted_df["auth_seq_id"] == insert_pos) &
                 (grafted_df["auth_chain_id"] == RHODOPSIN_CHAIN))
        grafted_df = grafted_df[mask].copy()

        # Prepare triad residue for insertion
        triad_atoms = triad_atoms.copy()
        triad_atoms["auth_seq_id"] = insert_pos
        triad_atoms["auth_chain_id"] = RHODOPSIN_CHAIN
        triad_atoms["label_chain_id"] = RHODOPSIN_CHAIN
        triad_atoms["entity_id"] = "1"
        triad_atoms["gen_seq_id"] = insert_pos

        if insert_pos in auth_to_label:
            triad_atoms["label_seq_id"] = auth_to_label[insert_pos]
        else:
            triad_atoms["label_seq_id"] = insert_pos

        triad_atoms["res_name"] = triad_info.residue_names[i]
        grafted_df = pd.concat([grafted_df, triad_atoms], ignore_index=True)

    # Sort
    grafted_df = grafted_df.sort_values(["auth_chain_id", "auth_seq_id", "atom_id"])

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_cif_file(str(output_path), grafted_df, force_overwrite=True)

    return grafted_df


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


def load_optimal_placements(placements_file: Path, top_n: int = 4) -> List[Dict]:
    """Load the top N optimal placements from the geometric optimization."""
    with open(placements_file) as f:
        data = json.load(f)

    placements = data["placements"][:top_n]
    print(f"  Loaded {len(placements)} optimal placements (top {top_n})")

    for p in placements:
        print(f"    Rank {p['rank']}: RMSD {p['rmsd']:.2f}A - "
              f"{p['helix_assignment']} -> residues {p['target_residues']}")

    return placements


def main() -> int:
    """Run the optimal triad grafting workflow."""
    print("=" * 70)
    print("BOLTZGEN OPTIMAL TRIAD GRAFTING")
    print("Using Geometrically Optimized TM3-TM5 Placements")
    print("=" * 70)
    print()

    # Initialize Protos
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    print(f"Data root: {data_root}")

    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)

    # ==========================================================================
    # Step 1: Load optimal placements
    # ==========================================================================
    print("\n[1] Loading optimal placements...")

    placements_file = OUTPUT_DIR / "correct_triad_placements.json"
    if not placements_file.exists():
        print(f"ERROR: Placements file not found: {placements_file}")
        return 1

    optimal_placements = load_optimal_placements(placements_file, top_n=4)

    # ==========================================================================
    # Step 2: Load rhodopsin scaffold
    # ==========================================================================
    print("\n[2] Loading rhodopsin scaffold...")

    loader.download_and_register(RHODOPSIN_PDB, name=RHODOPSIN_PDB)
    rhodopsin_df = sp.load_entity(RHODOPSIN_PDB)

    if rhodopsin_df is None:
        raise RuntimeError("Failed to load rhodopsin")

    rhodopsin_df = rhodopsin_df.reset_index()
    print(f"  Loaded {len(rhodopsin_df)} atoms (raw)")

    # Filter to chain A only
    rhodopsin_df = rhodopsin_df[rhodopsin_df["auth_chain_id"] == RHODOPSIN_CHAIN].copy()
    print(f"  Filtered to chain A: {len(rhodopsin_df)} atoms")

    # Clean: keep amino acids + retinal only
    rhodopsin_df = clean_structure(rhodopsin_df, RETINAL_CODES)
    print(f"  Cleaned (AA + RET): {len(rhodopsin_df)} atoms")

    # Renumber label_seq_id to be sequential
    unique_residues = rhodopsin_df[["auth_seq_id", "res_name"]].drop_duplicates()
    unique_residues = unique_residues.sort_values("auth_seq_id")
    auth_to_label = {auth_seq: i + 1 for i, auth_seq in enumerate(unique_residues["auth_seq_id"].values)}
    rhodopsin_df["label_seq_id"] = rhodopsin_df["auth_seq_id"].map(auth_to_label)
    rhodopsin_df["label_chain_id"] = rhodopsin_df["auth_chain_id"]
    rhodopsin_df["entity_id"] = "1"
    print(f"  Renumbered label_seq_id: 1-{max(auth_to_label.values())}")

    # ==========================================================================
    # Step 3: Extract triads from enzymes
    # ==========================================================================
    print("\n[3] Extracting catalytic triads from enzymes...")

    triads = {}

    for enzyme_key, enzyme_config in REFERENCE_ENZYMES.items():
        enzyme_config["key"] = enzyme_key
        pdb_id = enzyme_config["pdb"].lower()
        chain_id = enzyme_config["chain"]

        print(f"\n  {enzyme_config['description']}:")

        loader.download_and_register(pdb_id, name=f"enzyme_{enzyme_key}")
        enzyme_df = sp.load_entity(f"enzyme_{enzyme_key}")

        if enzyme_df is None:
            print(f"    Failed to load {pdb_id}")
            continue

        enzyme_df = enzyme_df.reset_index()

        try:
            triad = extract_triad_geometry(enzyme_df, enzyme_config, chain_id)
            triads[enzyme_key] = triad

            # Compute CA-CA distances
            d01 = np.linalg.norm(triad.ca_coords[0] - triad.ca_coords[1])
            d12 = np.linalg.norm(triad.ca_coords[1] - triad.ca_coords[2])
            d02 = np.linalg.norm(triad.ca_coords[0] - triad.ca_coords[2])

            print(f"    Triad: {'-'.join(triad.residue_names)}")
            print(f"    CA distances: {d01:.1f}A, {d12:.1f}A, {d02:.1f}A")
            print(f"    Max span: {max(d01, d12, d02):.1f}A")

        except Exception as e:
            print(f"    Failed: {e}")

    print(f"\n  Total triads extracted: {len(triads)}")

    # ==========================================================================
    # Step 4: Create grafted structures for each placement x enzyme
    # ==========================================================================
    print("\n[4] Creating grafted structures...")

    grafted_dir = OUTPUT_DIR / "grafted_structures_optimal"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    all_configs = []

    for enzyme_key, triad in triads.items():
        print(f"\n  {enzyme_key}:")
        enzyme_dir = grafted_dir / enzyme_key
        enzyme_dir.mkdir(parents=True, exist_ok=True)

        for placement in optimal_placements:
            rank = placement["rank"]
            centroid = np.array(placement["centroid"])
            rotation = placement["rotation_xyz_deg"]
            target_residues = placement["target_residues"]

            # Transform triad to this placement
            transformed_residues = transform_triad_to_placement(
                triad, centroid, tuple(rotation)
            )

            # Create config name
            helix_str = "-".join(placement["helix_assignment"])
            config_name = f"{enzyme_key}_optimal{rank:02d}_{helix_str}"
            output_path = enzyme_dir / f"{config_name}.cif"

            try:
                grafted_df = create_grafted_structure(
                    rhodopsin_df,
                    transformed_residues,
                    triad,
                    target_residues,
                    output_path,
                )

                config = {
                    "name": config_name,
                    "enzyme": enzyme_key,
                    "placement_rank": rank,
                    "helix_assignment": placement["helix_assignment"],
                    "target_residues": target_residues,
                    "rmsd": placement["rmsd"],
                    "rotation": rotation,
                    "centroid": centroid.tolist(),
                    "triad_residues": triad.residue_names,
                    "structure_path": str(output_path),
                    "substrate_smiles": triad.substrate_smiles,
                    "substrate_name": triad.substrate_name,
                }
                all_configs.append(config)
                print(f"    Rank {rank}: {output_path.name}")

            except Exception as e:
                print(f"    Rank {rank}: FAILED - {e}")

    print(f"\n  Total grafted structures: {len(all_configs)}")

    # ==========================================================================
    # Step 5: Generate BoltzGen YAML configs
    # ==========================================================================
    print("\n[5] Generating BoltzGen configurations...")

    yaml_dir = OUTPUT_DIR / "boltzgen_configs_optimal"
    yaml_dir.mkdir(parents=True, exist_ok=True)

    boltzgen_configs = []

    for config in all_configs:
        # Define binding residues - the triad + nearby residues
        binding_residues = config["target_residues"]

        yaml_config = {
            "entities": [
                # Grafted rhodopsin structure
                {
                    "file": {
                        "path": config["structure_path"],
                        "include": [
                            {"chain": {"id": RHODOPSIN_CHAIN}}
                        ],
                        # Design the cytoplasmic loops
                        "design": [
                            {"chain": {"id": RHODOPSIN_CHAIN, "res_index": DESIGNABLE_REGIONS}}
                        ],
                        # FIX the triad positions (don't change these residues)
                        "not_design": [
                            {"chain": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": ",".join(str(p) for p in config["target_residues"])
                            }}
                        ],
                        # Structure visibility for inpainting
                        "structure_groups": [
                            {"group": {"id": "all", "visibility": 1}},
                            # Hide designable regions
                            {"group": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": DESIGNABLE_REGIONS,
                                "visibility": 0
                            }},
                            # But SHOW the triad positions (anchor points)
                            {"group": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": ",".join(str(p) for p in config["target_residues"]),
                                "visibility": 1
                            }},
                        ],
                        # Define binding site around triad
                        "binding_types": [
                            {"chain": {
                                "id": RHODOPSIN_CHAIN,
                                "binding": ",".join(str(p) for p in config["target_residues"])
                            }}
                        ],
                    }
                },
                # Substrate ligand
                {
                    "ligand": {
                        "id": "SUB",
                        "smiles": config["substrate_smiles"],
                    }
                },
                # Retinal ligand (essential for rhodopsin)
                {
                    "ligand": {
                        "id": "RET",
                        "smiles": RETINAL_SMILES,
                    }
                },
            ],
            "num_designs": 2,  # Start with 2 per config, increase for top configs later
        }

        yaml_path = yaml_dir / f"{config['name']}.yaml"

        with open(yaml_path, "w") as f:
            yaml.dump(yaml_config, f, default_flow_style=False)

        boltzgen_configs.append({
            **config,
            "yaml_path": str(yaml_path),
        })

    print(f"  Created {len(boltzgen_configs)} BoltzGen configs")
    print(f"  Saved to: {yaml_dir}")

    # ==========================================================================
    # Step 6: Save workflow summary
    # ==========================================================================
    print("\n[6] Saving workflow summary...")

    summary = {
        "workflow": "boltzgen_optimal_grafting",
        "timestamp": datetime.now().isoformat(),
        "description": "Optimal triad placements at TM3-TM5 interface",
        "scaffold": {
            "pdb": RHODOPSIN_PDB,
            "chain": RHODOPSIN_CHAIN,
        },
        "designable_regions": DESIGNABLE_REGIONS,
        "placements_used": len(optimal_placements),
        "enzymes": {
            key: {
                "id": triad.enzyme_id,
                "triad": triad.residue_names,
                "substrate": triad.substrate_name,
            }
            for key, triad in triads.items()
        },
        "total_configs": len(boltzgen_configs),
        "configs": [
            {
                "name": c["name"],
                "enzyme": c["enzyme"],
                "placement_rank": c["placement_rank"],
                "helix_assignment": c["helix_assignment"],
                "target_residues": c["target_residues"],
                "rmsd": c["rmsd"],
                "structure_path": c["structure_path"],
                "yaml_path": c["yaml_path"],
            }
            for c in boltzgen_configs
        ],
    }

    summary_path = OUTPUT_DIR / "optimal_grafting_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  Summary: {summary_path}")

    # ==========================================================================
    # Final summary
    # ==========================================================================
    print("\n" + "=" * 70)
    print("WORKFLOW COMPLETE")
    print("=" * 70)
    print(f"Scaffold: {RHODOPSIN_PDB}")
    print(f"Enzymes: {list(triads.keys())}")
    print(f"Placements per enzyme: {len(optimal_placements)}")
    print(f"Total configurations: {len(boltzgen_configs)}")
    print()
    print("OUTPUTS:")
    print(f"  Grafted structures: {grafted_dir}")
    print(f"  BoltzGen configs: {yaml_dir}")
    print(f"  Summary: {summary_path}")
    print()
    print("KEY IMPROVEMENTS:")
    print("  - TM3-TM5 interface (RMSD < 0.5A) vs ICL3")
    print("  - Geometrically optimized placement")
    print("  - Includes retinal in BoltzGen configs")
    print()
    print("NEXT STEPS:")
    print("  1. Review grafted structures in PyMOL or HTML viewer")
    print("  2. Run BoltzGen: python run_boltzgen_optimal.py")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
