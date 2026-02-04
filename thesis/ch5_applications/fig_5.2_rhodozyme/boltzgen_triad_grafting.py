#!/usr/bin/env python3
"""BoltzGen Triad Grafting Workflow - Direct Catalytic Site Insertion.

Approach:
1. Extract the catalytic triad (full residues with coordinates) from enzyme CIF
2. Place the triad at the CENTER of the rhodopsin cytoplasmic pocket
3. Try different ORIENTATIONS of the triad (rotation grid search)
4. Let BoltzGen DESIGN the connecting loops around the fixed triad
5. Include SUBSTRATE at the center for binding pocket formation

The triad residues become "anchor points" - their Cα positions are fixed,
and BoltzGen generates the surrounding protein to connect them.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from datetime import datetime
from dataclasses import dataclass
import itertools

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.structure.structure_utils import load_structure
from protos.io.formats.cif_utils import write_cif_file
from protos.io.ingest.structure_loader import StructureLoader
from protos.analysis.structure.alignment import kabsch_alignment


# =============================================================================
# Configuration
# =============================================================================

RHODOPSIN_PDB = "3pqr"
RHODOPSIN_CHAIN = "A"

# Cytoplasmic pocket center: defined by TM5/TM6 cytoplasmic ends
# These residues form the "opening" that moves during activation
POCKET_CENTER_RESIDUES = [230, 231, 232, 233, 234, 235,  # TM5 cyto end
                          244, 245, 246, 247, 248, 249]  # TM6 cyto start

# Designable regions (loops to be generated)
DESIGNABLE_REGIONS = "58..76,134..157,222..256,306..326"

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
    "papain": {
        "id": "PAP",
        "pdb": "9PAP",
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
}

# Rotation grid for orientation search
# Euler angles (degrees) around X, Y, Z axes
ROTATION_ANGLES = [0, 60, 120, 180, 240, 300]  # 6 angles per axis = 216 combinations
# For faster testing, use fewer angles:
ROTATION_ANGLES_FAST = [0, 90, 180, 270]  # 4 angles = 64 combinations


@dataclass
class TriadGeometry:
    """Stores the geometry of a catalytic triad."""
    enzyme_key: str
    enzyme_id: str

    # Residue info
    residue_names: List[str]  # e.g., ["SER", "HIS", "ASP"]
    residue_roles: List[str]  # e.g., ["nucleophile", "base", "electrostatic"]

    # Coordinates (all atoms for each residue)
    atoms_per_residue: List[pd.DataFrame]  # Full atom data for each triad residue

    # Key reference points
    ca_coords: np.ndarray  # (3, 3) - CA coordinates of each residue
    centroid: np.ndarray   # (3,) - centroid of the triad

    # Substrate
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
    """Extract the catalytic triad geometry from an enzyme structure.

    Returns full atom coordinates for each residue, plus CA positions
    and centroid for alignment.
    """
    triad_residues = enzyme_config["triad_residues"]

    atoms_per_residue = []
    ca_coords = []
    residue_names = []
    residue_roles = []

    for res_info in triad_residues:
        res_num = res_info["res_num"]
        res_name = res_info["res_name"]
        role = res_info["role"]

        # Get all atoms for this residue
        mask = (enzyme_df["auth_seq_id"] == res_num) & (enzyme_df["auth_chain_id"] == chain_id)
        res_atoms = enzyme_df[mask].copy()

        if len(res_atoms) == 0:
            raise ValueError(f"Residue {res_name}{res_num} not found in chain {chain_id}")

        atoms_per_residue.append(res_atoms)
        residue_names.append(res_name)
        residue_roles.append(role)

        # Get CA coordinate
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


def compute_pocket_center(
    rhodopsin_df: pd.DataFrame,
    pocket_residues: List[int],
    chain_id: str,
) -> np.ndarray:
    """Compute the center of the cytoplasmic binding pocket."""
    mask = (
        (rhodopsin_df["auth_chain_id"] == chain_id) &
        (rhodopsin_df["auth_seq_id"].isin(pocket_residues)) &
        (rhodopsin_df["atom_name"] == "CA")
    )

    pocket_ca = rhodopsin_df[mask][["x", "y", "z"]].values

    if len(pocket_ca) == 0:
        raise ValueError("No CA atoms found for pocket residues")

    return pocket_ca.mean(axis=0)


def transform_triad_to_pocket(
    triad: TriadGeometry,
    pocket_center: np.ndarray,
    rotation_angles: Tuple[float, float, float],
) -> List[pd.DataFrame]:
    """Transform triad atoms to the pocket center with given rotation.

    1. Translate triad so centroid is at origin
    2. Apply rotation
    3. Translate to pocket center

    Returns transformed atom DataFrames for each residue.
    """
    # Rotation matrix
    R = rotation_matrix_from_euler(*rotation_angles)

    # Translation: centroid -> origin -> pocket_center
    translation = pocket_center - triad.centroid

    transformed_residues = []

    for res_atoms in triad.atoms_per_residue:
        res_atoms = res_atoms.copy()

        # Get coordinates
        coords = res_atoms[["x", "y", "z"]].values

        # Center on triad centroid
        coords_centered = coords - triad.centroid

        # Rotate
        coords_rotated = (R @ coords_centered.T).T

        # Translate to pocket
        coords_final = coords_rotated + pocket_center

        # Update DataFrame
        res_atoms["x"] = coords_final[:, 0]
        res_atoms["y"] = coords_final[:, 1]
        res_atoms["z"] = coords_final[:, 2]

        transformed_residues.append(res_atoms)

    return transformed_residues


def create_grafted_structure(
    rhodopsin_df: pd.DataFrame,
    triad_residues: List[pd.DataFrame],
    triad_info: TriadGeometry,
    insertion_points: List[int],  # Where to insert each triad residue in the sequence
    output_path: Path,
) -> pd.DataFrame:
    """Create a grafted structure with triad residues inserted into rhodopsin.

    The triad residues are inserted at specified sequence positions,
    replacing the original residues at those positions.

    Args:
        rhodopsin_df: Full rhodopsin structure (already filtered to chain A only)
        triad_residues: Transformed triad residue DataFrames
        triad_info: TriadGeometry with residue names
        insertion_points: Sequence positions to insert each triad residue
        output_path: Where to save the grafted CIF

    Returns:
        DataFrame of the grafted structure
    """
    grafted_df = rhodopsin_df.copy()

    # Build auth_seq_id -> label_seq_id mapping from existing structure
    existing_mapping = grafted_df[["auth_seq_id", "label_seq_id"]].drop_duplicates()
    auth_to_label = dict(zip(existing_mapping["auth_seq_id"], existing_mapping["label_seq_id"]))

    # For each triad residue, replace the residue at the insertion point
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

        # Set label_seq_id to match the sequence position in entity_poly_seq
        if insert_pos in auth_to_label:
            triad_atoms["label_seq_id"] = auth_to_label[insert_pos]
        else:
            triad_atoms["label_seq_id"] = insert_pos

        # Update residue name to the triad residue
        triad_atoms["res_name"] = triad_info.residue_names[i]

        # Append triad residue
        grafted_df = pd.concat([grafted_df, triad_atoms], ignore_index=True)

    # Sort by chain and residue number
    grafted_df = grafted_df.sort_values(["auth_chain_id", "auth_seq_id", "atom_id"])

    # Ensure output directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Save structure using protos CIF writer
    write_cif_file(str(output_path), grafted_df, force_overwrite=True)

    return grafted_df


def generate_orientation_grid(
    angles: List[float] = ROTATION_ANGLES_FAST,
    filter_redundant: bool = True,
) -> List[Tuple[float, float, float]]:
    """Generate grid of rotation angles to try.

    Args:
        angles: List of angles (degrees) to try for each axis
        filter_redundant: Remove some redundant rotations

    Returns:
        List of (rx, ry, rz) tuples
    """
    # Full grid
    grid = list(itertools.product(angles, angles, angles))

    if filter_redundant:
        # Remove some obviously redundant cases
        # (This is approximate - true symmetry depends on the triad shape)
        filtered = []
        for rx, ry, rz in grid:
            # Keep if rx <= 180 (half sphere)
            if rx <= 180:
                filtered.append((rx, ry, rz))
        grid = filtered

    return grid


def compute_triad_distances(ca_coords: np.ndarray) -> Tuple[float, float, float]:
    """Compute pairwise CA-CA distances for a triad."""
    d01 = np.linalg.norm(ca_coords[0] - ca_coords[1])
    d12 = np.linalg.norm(ca_coords[1] - ca_coords[2])
    d02 = np.linalg.norm(ca_coords[0] - ca_coords[2])
    return d01, d12, d02


def find_best_insertion_points(
    rhodopsin_df: pd.DataFrame,
    pocket_center: np.ndarray,
    triad_ca_transformed: np.ndarray,
    designable_residues: List[int],
) -> List[int]:
    """Find the best residue positions to insert each triad residue.

    For each triad CA position, find the closest designable residue
    in rhodopsin.
    """
    # Get CA atoms of designable residues
    mask = (
        (rhodopsin_df["auth_chain_id"] == RHODOPSIN_CHAIN) &
        (rhodopsin_df["auth_seq_id"].isin(designable_residues)) &
        (rhodopsin_df["atom_name"] == "CA")
    )

    design_ca = rhodopsin_df[mask][["auth_seq_id", "x", "y", "z"]].copy()
    design_coords = design_ca[["x", "y", "z"]].values
    design_resids = design_ca["auth_seq_id"].values

    insertion_points = []
    used_resids = set()

    for triad_ca in triad_ca_transformed:
        # Find closest designable residue
        distances = np.linalg.norm(design_coords - triad_ca, axis=1)

        # Skip already used positions
        for idx in np.argsort(distances):
            resid = design_resids[idx]
            if resid not in used_resids:
                insertion_points.append(int(resid))
                used_resids.add(resid)
                break

    return insertion_points


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


def main() -> int:
    """Run the triad grafting workflow."""
    print("=" * 70)
    print("BOLTZGEN TRIAD GRAFTING WORKFLOW")
    print("Direct Catalytic Site Insertion with Orientation Grid Search")
    print("=" * 70)
    print()
    print("APPROACH:")
    print("  1. Extract catalytic triad from enzyme structure")
    print("  2. Place triad at CENTER of rhodopsin cytoplasmic pocket")
    print("  3. Try different ORIENTATIONS (rotation grid)")
    print("  4. Insert triad residues at optimal positions")
    print("  5. BoltzGen designs connecting loops")
    print()
    print(f"  Scaffold: {RHODOPSIN_PDB}")
    print(f"  Enzymes: {list(REFERENCE_ENZYMES.keys())}")
    print("=" * 70)

    # Initialize Protos - use project data directory
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    print(f"  Data root: {data_root}")

    sp = StructureProcessor()
    loader = StructureLoader(processor=sp)

    # ==========================================================================
    # Step 1: Download and load rhodopsin
    # ==========================================================================
    print("\n[1] Loading rhodopsin scaffold...")

    loader.download_and_register(RHODOPSIN_PDB, name=RHODOPSIN_PDB)
    rhodopsin_df = sp.load_entity(RHODOPSIN_PDB)

    if rhodopsin_df is None:
        raise RuntimeError("Failed to load rhodopsin")

    rhodopsin_df = rhodopsin_df.reset_index()
    print(f"  ✓ Loaded {len(rhodopsin_df)} atoms (raw)")

    # Filter to only chain A (the main rhodopsin chain)
    # Exclude chain B (H5 peptide) and other ligands except retinal (RET)
    retinal_codes = ["RET", "LYR", "LYS"]  # Retinal and Schiff base lysine
    chain_a_mask = rhodopsin_df["auth_chain_id"] == RHODOPSIN_CHAIN
    retinal_mask = rhodopsin_df["res_name"].isin(retinal_codes)

    # Keep chain A residues (including connected retinal)
    rhodopsin_df = rhodopsin_df[chain_a_mask].copy()
    print(f"  ✓ Filtered to chain A: {len(rhodopsin_df)} atoms")

    # Renumber label_seq_id to be sequential (required for BoltzGen entity_poly_seq)
    # Get unique residues in auth_seq_id order
    unique_residues = rhodopsin_df[["auth_seq_id", "res_name"]].drop_duplicates()
    unique_residues = unique_residues.sort_values("auth_seq_id")

    # Create mapping from auth_seq_id to sequential label_seq_id
    auth_to_label = {auth_seq: i + 1 for i, auth_seq in enumerate(unique_residues["auth_seq_id"].values)}
    rhodopsin_df["label_seq_id"] = rhodopsin_df["auth_seq_id"].map(auth_to_label)

    # Set label_chain_id = auth_chain_id for single chain
    rhodopsin_df["label_chain_id"] = rhodopsin_df["auth_chain_id"]
    rhodopsin_df["entity_id"] = "1"  # Single entity

    print(f"  ✓ Renumbered label_seq_id: 1-{max(auth_to_label.values())}")

    # Compute pocket center
    pocket_center = compute_pocket_center(
        rhodopsin_df, POCKET_CENTER_RESIDUES, RHODOPSIN_CHAIN
    )
    print(f"  ✓ Pocket center: [{pocket_center[0]:.1f}, {pocket_center[1]:.1f}, {pocket_center[2]:.1f}]")

    # Parse designable residues
    designable_residues = parse_designable_regions(DESIGNABLE_REGIONS)
    print(f"  ✓ Designable residues: {len(designable_residues)}")

    # ==========================================================================
    # Step 2: Extract triads from enzymes
    # ==========================================================================
    print("\n[2] Extracting catalytic triads from enzymes...")

    triads = {}

    for enzyme_key, enzyme_config in REFERENCE_ENZYMES.items():
        enzyme_config["key"] = enzyme_key
        pdb_id = enzyme_config["pdb"].lower()
        chain_id = enzyme_config["chain"]

        print(f"\n  {enzyme_config['description']}:")

        # Download enzyme
        loader.download_and_register(pdb_id, name=f"enzyme_{enzyme_key}")
        enzyme_df = sp.load_entity(f"enzyme_{enzyme_key}")

        if enzyme_df is None:
            print(f"    ✗ Failed to load {pdb_id}")
            continue

        enzyme_df = enzyme_df.reset_index()

        # Extract triad
        try:
            triad = extract_triad_geometry(enzyme_df, enzyme_config, chain_id)
            triads[enzyme_key] = triad

            # Show triad geometry
            distances = compute_triad_distances(triad.ca_coords)
            print(f"    ✓ Triad: {'-'.join(triad.residue_names)}")
            print(f"    ✓ CA distances: {distances[0]:.1f}Å, {distances[1]:.1f}Å, {distances[2]:.1f}Å")
            print(f"    ✓ Centroid: [{triad.centroid[0]:.1f}, {triad.centroid[1]:.1f}, {triad.centroid[2]:.1f}]")

        except Exception as e:
            print(f"    ✗ Failed: {e}")

    print(f"\n  Total triads extracted: {len(triads)}")

    # ==========================================================================
    # Step 3: Generate orientation grid
    # ==========================================================================
    print("\n[3] Generating orientation grid...")

    orientations = generate_orientation_grid(ROTATION_ANGLES_FAST, filter_redundant=True)
    print(f"  ✓ {len(orientations)} orientations to try per enzyme")
    print(f"  ✓ Total configurations: {len(orientations) * len(triads)}")

    # ==========================================================================
    # Step 4: Create grafted structures for each orientation
    # ==========================================================================
    print("\n[4] Creating grafted structures...")

    grafted_dir = OUTPUT_DIR / "grafted_structures"
    grafted_dir.mkdir(parents=True, exist_ok=True)

    all_configs = []

    for enzyme_key, triad in triads.items():
        print(f"\n  {enzyme_key}:")

        enzyme_dir = grafted_dir / enzyme_key
        enzyme_dir.mkdir(parents=True, exist_ok=True)

        for i, (rx, ry, rz) in enumerate(orientations):
            # Transform triad to pocket
            transformed_residues = transform_triad_to_pocket(
                triad, pocket_center, (rx, ry, rz)
            )

            # Get transformed CA coordinates for insertion point finding
            transformed_ca = []
            for res_atoms in transformed_residues:
                ca = res_atoms[res_atoms["atom_name"] == "CA"]
                transformed_ca.append(ca[["x", "y", "z"]].values[0])
            transformed_ca = np.array(transformed_ca)

            # Find best insertion points
            insertion_points = find_best_insertion_points(
                rhodopsin_df, pocket_center, transformed_ca, designable_residues
            )

            # Create grafted structure
            config_name = f"{enzyme_key}_rot{rx:03.0f}_{ry:03.0f}_{rz:03.0f}"
            output_path = enzyme_dir / f"{config_name}.cif"

            try:
                grafted_df = create_grafted_structure(
                    rhodopsin_df,
                    transformed_residues,
                    triad,
                    insertion_points,
                    output_path,
                )

                config = {
                    "name": config_name,
                    "enzyme": enzyme_key,
                    "rotation": (rx, ry, rz),
                    "insertion_points": insertion_points,
                    "triad_residues": triad.residue_names,
                    "structure_path": str(output_path),
                    "substrate_smiles": triad.substrate_smiles,
                    "substrate_name": triad.substrate_name,
                }
                all_configs.append(config)

            except Exception as e:
                print(f"    ✗ {config_name}: {e}")

        print(f"    ✓ Created {len(orientations)} grafted structures")

    print(f"\n  Total grafted structures: {len(all_configs)}")

    # ==========================================================================
    # Step 5: Generate BoltzGen YAML configs
    # ==========================================================================
    print("\n[5] Generating BoltzGen configurations...")

    yaml_dir = OUTPUT_DIR / "boltzgen_configs"
    yaml_dir.mkdir(parents=True, exist_ok=True)

    boltzgen_configs = []

    for config in all_configs:
        # BoltzGen config for this grafted structure
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
                        # FIX the triad positions
                        "not_design": [
                            {"chain": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": ",".join(str(p) for p in config["insertion_points"])
                            }}
                        ],
                        # Hide designable regions for inpainting
                        "structure_groups": [
                            {"group": {"id": "all", "visibility": 1}},
                            {"group": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": DESIGNABLE_REGIONS,
                                "visibility": 0
                            }},
                            # But SHOW the triad positions
                            {"group": {
                                "id": RHODOPSIN_CHAIN,
                                "res_index": ",".join(str(p) for p in config["insertion_points"]),
                                "visibility": 1
                            }},
                        ],
                        # Binding site around triad
                        "binding_types": [
                            {"chain": {
                                "id": RHODOPSIN_CHAIN,
                                "binding": ",".join(str(p) for p in config["insertion_points"])
                            }}
                        ],
                    }
                },
                # Substrate
                {
                    "ligand": {
                        "id": "SUB",
                        "smiles": config["substrate_smiles"],
                    }
                },
            ],
            "num_designs": 2,
        }

        yaml_path = yaml_dir / f"{config['name']}.yaml"

        import yaml
        with open(yaml_path, "w") as f:
            yaml.dump(yaml_config, f, default_flow_style=False)

        boltzgen_configs.append({
            **config,
            "yaml_path": str(yaml_path),
        })

    print(f"  ✓ Created {len(boltzgen_configs)} BoltzGen configs")
    print(f"  ✓ Saved to: {yaml_dir}")

    # ==========================================================================
    # Step 6: Save workflow summary
    # ==========================================================================
    print("\n[6] Saving workflow summary...")

    summary = {
        "workflow": "boltzgen_triad_grafting",
        "timestamp": datetime.now().isoformat(),
        "scaffold": {
            "pdb": RHODOPSIN_PDB,
            "chain": RHODOPSIN_CHAIN,
            "pocket_center": pocket_center.tolist(),
            "pocket_residues": POCKET_CENTER_RESIDUES,
        },
        "designable_regions": DESIGNABLE_REGIONS,
        "enzymes": {
            key: {
                "id": triad.enzyme_id,
                "triad": triad.residue_names,
                "ca_distances": compute_triad_distances(triad.ca_coords),
                "substrate": triad.substrate_name,
            }
            for key, triad in triads.items()
        },
        "orientations": {
            "num_angles": len(ROTATION_ANGLES_FAST),
            "total_orientations": len(orientations),
        },
        "total_configs": len(boltzgen_configs),
        "configs": [
            {
                "name": c["name"],
                "enzyme": c["enzyme"],
                "rotation": c["rotation"],
                "insertion_points": c["insertion_points"],
                "structure_path": c["structure_path"],
                "yaml_path": c["yaml_path"],
            }
            for c in boltzgen_configs
        ],
    }

    summary_path = OUTPUT_DIR / "triad_grafting_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"  ✓ Summary: {summary_path}")

    # ==========================================================================
    # Final summary
    # ==========================================================================
    print("\n" + "=" * 70)
    print("WORKFLOW COMPLETE")
    print("=" * 70)
    print(f"Scaffold: {RHODOPSIN_PDB}")
    print(f"Enzymes: {list(triads.keys())}")
    print(f"Orientations per enzyme: {len(orientations)}")
    print(f"Total configurations: {len(boltzgen_configs)}")
    print()
    print("OUTPUTS:")
    print(f"  Grafted structures: {grafted_dir}")
    print(f"  BoltzGen configs: {yaml_dir}")
    print(f"  Summary: {summary_path}")
    print()
    print("NEXT STEPS:")
    print("  1. Review grafted structures in PyMOL")
    print("  2. Run BoltzGen on selected configs")
    print("  3. Compare designed loops by:")
    print("     - Catalytic geometry preservation")
    print("     - Substrate binding pose")
    print("     - Loop stability (pLDDT)")
    print()
    print("To run BoltzGen on a config:")
    print("  mm = ModelManager()")
    print("  config = yaml.safe_load(open('config.yaml'))")
    print("  invocation = mm.prepare('boltzgen', config=config)")
    print("  state = mm.submit_job(invocation)")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
