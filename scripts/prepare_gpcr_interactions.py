#!/usr/bin/env python3
"""Calculate ligand-protein interactions for GPCR structures and export for PyMOL.

This script:
1. Loads 8 GPCR structures
2. Calculates ligand-protein interactions
3. Annotates binding residues with GRN labels
4. Exports PyMOL scripts for visualization
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional

import pandas as pd

# Add src to path
REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor
from protos.analysis.structure_ligand_analysis import (
    calculate_ligand_interactions,
    get_ligand_by_id,
)

DATA_ROOT = REPO_ROOT / "data"
OUTPUT_DIR = DATA_ROOT / "visualizations"

# Structure metadata with ligand codes and resolution
STRUCTURES: Dict[str, Dict[str, Any]] = {
    "2RH1": {"chain": "A", "ligand": "CAU", "ligand_name": "carazolol", "ligand_type": "inverse_agonist", "state": "inactive", "receptor": "ADRB2", "resolution": 2.4},
    "3NY9": {"chain": "A", "ligand": "JSZ", "ligand_name": "ICI-118551", "ligand_type": "inverse_agonist", "state": "inactive", "receptor": "ADRB2", "resolution": 2.8},
    "3SN6": {"chain": "R", "ligand": "P0G", "ligand_name": "BI-167107", "ligand_type": "full_agonist", "state": "active", "receptor": "ADRB2", "resolution": 3.2},
    "4LDO": {"chain": "A", "ligand": "ALE", "ligand_name": "adrenaline", "ligand_type": "full_agonist", "state": "active", "receptor": "ADRB2", "resolution": 2.79},
    "2VT4": {"chain": "A", "ligand": "P32", "ligand_name": "cyanopindolol", "ligand_type": "antagonist", "state": "inactive", "receptor": "ADRB1", "resolution": 2.7},
    "2Y02": {"chain": "A", "ligand": "WHJ", "ligand_name": "isoprenaline", "ligand_type": "full_agonist", "state": "active_like", "receptor": "ADRB1", "resolution": 2.6},
    "2Y04": {"chain": "A", "ligand": "68H", "ligand_name": "salbutamol", "ligand_type": "partial_agonist", "state": "intermediate", "receptor": "ADRB1", "resolution": 2.6},
    "2Y00": {"chain": "A", "ligand": "Y00", "ligand_name": "dobutamine", "ligand_type": "partial_agonist", "state": "intermediate", "receptor": "ADRB1", "resolution": 2.5},
}

# Structure colors (same as main visualization script)
STRUCTURE_COLORS = {
    "2RH1": "marine",
    "3NY9": "slate",
    "3SN6": "forest",
    "4LDO": "lime",
    "2VT4": "deepteal",
    "2Y02": "chartreuse",
    "2Y04": "yellow",
    "2Y00": "orange",
}

# Ligand colors - two separate scales
# Inverse agonist/antagonist: magenta scale
# Agonist: green scale
LIGAND_COLORS = {
    "2RH1": "magenta",     # inverse agonist
    "3NY9": "hotpink",     # inverse agonist
    "2VT4": "salmon",      # antagonist
    "3SN6": "forest",      # full agonist
    "4LDO": "lime",        # full agonist
    "2Y02": "chartreuse",  # full agonist
    "2Y04": "splitpea",    # partial agonist
    "2Y00": "olive",       # partial agonist
}

INTERACTION_CUTOFF = 4.0  # Angstroms


def calculate_all_interactions(sp: StructureProcessor) -> Dict[str, Any]:
    """Calculate interactions for all structures."""
    results = {}

    for pdb_id, meta in STRUCTURES.items():
        chain = meta["chain"]
        ligand_code = meta["ligand"]

        print(f"\n[{pdb_id}] Processing...")

        try:
            # Load structure
            structure_df = sp.load_entity(pdb_id)
            if structure_df is None:
                print(f"  ERROR: Could not load structure")
                continue

            structure_df = structure_df.reset_index()

            # Get ligand atoms
            ligand_atoms = structure_df[
                (structure_df["group"] == "HETATM") &
                (structure_df["res_name3l"] == ligand_code)
            ]

            if ligand_atoms.empty:
                print(f"  ERROR: Ligand {ligand_code} not found")
                continue

            print(f"  Ligand: {ligand_code} ({len(ligand_atoms)} atoms)")

            # Calculate interactions
            interactions = calculate_ligand_interactions(
                sp, pdb_id, ligand_atoms, detailed=True, cutoff=INTERACTION_CUTOFF
            )

            binding_residues = interactions.get("binding_residues", [])
            print(f"  Binding residues: {len(binding_residues)}")

            # Try to get GRN annotations - filter by chain to avoid antibody/fusion protein rows
            grn_map = {}
            try:
                # Check if structure has GRN column
                if "grn" in structure_df.columns:
                    chain_df = structure_df[structure_df["auth_chain_id"] == chain]
                    grn_df = chain_df[chain_df["grn"].notna()][["auth_seq_id", "grn"]].drop_duplicates()
                    for _, row in grn_df.iterrows():
                        grn_map[int(row["auth_seq_id"])] = row["grn"]
                else:
                    # Try to annotate
                    sp.annotate_with_grn(pdb_id, chains=[chain])
                    structure_df = sp.load_entity(pdb_id).reset_index()
                    if "grn" in structure_df.columns:
                        chain_df = structure_df[structure_df["auth_chain_id"] == chain]
                        grn_df = chain_df[chain_df["grn"].notna()][["auth_seq_id", "grn"]].drop_duplicates()
                        for _, row in grn_df.iterrows():
                            grn_map[int(row["auth_seq_id"])] = row["grn"]
            except Exception as e:
                print(f"  Warning: Could not get GRN annotations: {e}")

            # Add GRN and closest ligand atom to binding residues
            import numpy as np
            ligand_coords = ligand_atoms[['x', 'y', 'z']].values
            ligand_atom_names = ligand_atoms['atom_name'].values

            for res in binding_residues:
                res_id = int(res.get("res_id", 0))
                res["grn"] = grn_map.get(res_id, "")

                # Find closest ligand atom to this residue's contact atom
                contact_atoms_list = res.get("contact_atoms", [])
                if contact_atoms_list:
                    # Get the sidechain atoms for this residue
                    res_df = structure_df[
                        (structure_df["auth_seq_id"] == res_id) &
                        (structure_df["auth_chain_id"] == chain)
                    ]
                    # Filter to contact atoms
                    contact_df = res_df[res_df["atom_name"].isin(contact_atoms_list)]

                    if not contact_df.empty:
                        # Find the closest ligand atom to any of the contact atoms
                        min_dist = float('inf')
                        closest_lig_atom = None
                        closest_res_atom = None

                        for _, res_row in contact_df.iterrows():
                            res_coord = np.array([res_row['x'], res_row['y'], res_row['z']])
                            dists = np.sqrt(np.sum((ligand_coords - res_coord)**2, axis=1))
                            min_idx = np.argmin(dists)
                            if dists[min_idx] < min_dist:
                                min_dist = dists[min_idx]
                                closest_lig_atom = ligand_atom_names[min_idx]
                                closest_res_atom = res_row['atom_name']

                        res["closest_ligand_atom"] = closest_lig_atom
                        res["closest_residue_atom"] = closest_res_atom

            print(f"  GRN annotations: {sum(1 for r in binding_residues if r.get('grn'))}")

            # Count crystallographic waters and calculate water network
            waters = structure_df[structure_df["res_name3l"] == "HOH"]
            water_count = len(waters["auth_seq_id"].unique()) if not waters.empty else 0
            print(f"  Crystallographic waters: {water_count}")

            # Calculate water network interactions
            water_network = []
            water_sidechain = []
            WATER_DIST_CUTOFF = 3.5  # Angstroms for H-bond distance

            if not waters.empty:
                # Get water oxygen positions (one per water)
                water_oxygens = waters[waters["atom_name"] == "O"][["auth_seq_id", "x", "y", "z"]].drop_duplicates("auth_seq_id")

                if len(water_oxygens) > 1:
                    water_coords = water_oxygens[["x", "y", "z"]].values
                    water_ids = water_oxygens["auth_seq_id"].values

                    # Water-water interactions
                    for i in range(len(water_coords)):
                        for j in range(i + 1, len(water_coords)):
                            dist = np.sqrt(np.sum((water_coords[i] - water_coords[j])**2))
                            if dist <= WATER_DIST_CUTOFF:
                                water_network.append({
                                    "water1_id": int(water_ids[i]),
                                    "water2_id": int(water_ids[j]),
                                    "distance": float(dist),
                                })

                # Water-sidechain interactions (to binding pocket residues)
                binding_res_ids = [int(r.get("res_id", 0)) for r in binding_residues if r.get("grn")]
                if binding_res_ids:
                    pocket_atoms = structure_df[
                        (structure_df["auth_seq_id"].isin(binding_res_ids)) &
                        (structure_df["auth_chain_id"] == chain) &
                        (~structure_df["atom_name"].isin(["C", "N", "O", "CA"]))  # Sidechain only
                    ]

                    if not pocket_atoms.empty:
                        pocket_coords = pocket_atoms[["x", "y", "z"]].values
                        pocket_res_ids = pocket_atoms["auth_seq_id"].values
                        pocket_atom_names = pocket_atoms["atom_name"].values

                        for i, (wx, wy, wz) in enumerate(water_coords):
                            water_coord = np.array([wx, wy, wz])
                            dists = np.sqrt(np.sum((pocket_coords - water_coord)**2, axis=1))
                            close_idx = np.where(dists <= WATER_DIST_CUTOFF)[0]

                            for idx in close_idx:
                                water_sidechain.append({
                                    "water_id": int(water_ids[i]),
                                    "res_id": int(pocket_res_ids[idx]),
                                    "atom_name": pocket_atom_names[idx],
                                    "distance": float(dists[idx]),
                                })

            print(f"  Water-water interactions: {len(water_network)}")
            print(f"  Water-sidechain interactions: {len(water_sidechain)}")

            # Store results
            results[pdb_id] = {
                "meta": meta,
                "binding_residues": binding_residues,
                "hydrogen_bonds": interactions.get("hydrogen_bonds", []),
                "hydrophobic": interactions.get("hydrophobic", []),
                "pi_stacking": interactions.get("pi_stacking", []),
                "salt_bridges": interactions.get("salt_bridges", []),
                "water_mediated": interactions.get("water_mediated", []),
                "water_network": water_network,
                "water_sidechain": water_sidechain,
                "summary": interactions.get("summary", {}),
                "water_count": water_count,
                "water_note": "No crystallographic waters" if water_count == 0 else None,
            }

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

    return results


def generate_pocket_script(
    pdb_ids: List[str],
    interactions: Dict[str, Any],
    output_path: Path,
    title: str,
    scene_name: str,
    show_grn: bool = True,
    show_water_network: bool = False,
) -> None:
    """Generate PyMOL script for a specific pocket view with interactions.

    Args:
        show_grn: If True, include GRN labels on binding residues.
        show_water_network: If True, include water network interactions (water-water, water-sidechain).
    """

    lines = [
        f"# GPCR Binding Pocket View - {title}",
        f"# Structures: {', '.join(pdb_ids)}",
        "#",
        "# Run after: pymol scripts/gpcr_visualization.pml",
        "",
        f'print("Setting up {title} binding pocket view...")',
        "",
        "# Disable all structures",
        "disable all",
        "",
        "# Enable selected structures",
    ]

    for pdb_id in pdb_ids:
        lines.append(f"enable {pdb_id}_gpcr")

    # Show waters
    sel = " or ".join([f"{pdb_id}_gpcr" for pdb_id in pdb_ids])

    # Build ligand selection for zoom
    ligand_sels = []
    for pdb_id in pdb_ids:
        if pdb_id in interactions:
            ligand = interactions[pdb_id]["meta"]["ligand"]
            ligand_sels.append(f"({pdb_id}_gpcr and resn {ligand})")
    ligand_sel = " or ".join(ligand_sels) if ligand_sels else sel

    lines.extend([
        "",
        "# Show waters as spheres (colored by structure)",
        f"show spheres, ({sel}) and resn HOH",
        "set sphere_scale, 0.3, resn HOH",
    ])

    # Color waters by structure (same as protein)
    for pdb_id in pdb_ids:
        color = STRUCTURE_COLORS.get(pdb_id, "cyan")
        lines.append(f"color {color}, {pdb_id}_gpcr and resn HOH")

    lines.extend([
        "",
        "# Use stored pocket view",
        "view pocket_view, recall",
        "",
    ])

    # Add binding residue selections with GRN labels
    lines.extend([
        "# =============================================================================",
        "# Binding residue selections and labels",
        "# =============================================================================",
        "",
    ])

    all_grn_residues = {}  # grn -> [(pdb_id, chain, res_id, res_name)]

    for pdb_id in pdb_ids:
        if pdb_id not in interactions:
            continue

        data = interactions[pdb_id]
        chain = data["meta"]["chain"]
        binding_residues = data["binding_residues"]

        # Sort by distance
        binding_residues = sorted(binding_residues, key=lambda x: x.get("min_distance", 999))

        # Take top 10 closest residues
        top_residues = binding_residues[:10]

        for res in top_residues:
            res_id = res.get("res_id")
            res_name = res.get("res_name", "UNK")
            grn = res.get("grn", "")

            if grn:
                if grn not in all_grn_residues:
                    all_grn_residues[grn] = []
                all_grn_residues[grn].append((pdb_id, chain, res_id, res_name))

    # Create selections for each GRN position
    for grn in sorted(all_grn_residues.keys()):
        grn_safe = grn.replace(".", "_")
        residues = all_grn_residues[grn]

        sel_parts = []
        for pdb_id, chain, res_id, res_name in residues:
            sel_parts.append(f"({pdb_id}_gpcr and chain {chain} and resi {res_id})")

        if sel_parts:
            selection = " or ".join(sel_parts)
            lines.append(f"select grn_{grn_safe}, {selection}")

    # Show sidechains
    lines.extend([
        "",
        "# Show binding pocket sidechains",
    ])

    grn_sels = [f"grn_{grn.replace('.', '_')}" for grn in sorted(all_grn_residues.keys())]
    if grn_sels:
        all_grn_sel = " or ".join(grn_sels)
        lines.extend([
            f"show sticks, ({all_grn_sel}) and sidechain",
            f"set stick_radius, 0.15, ({all_grn_sel})",
            "",
        ])

    # Color sidechains by structure (same as protein)
    lines.append("# Color sidechains by structure (same as protein)")
    for pdb_id in pdb_ids:
        color = STRUCTURE_COLORS.get(pdb_id, "gray")
        lines.append(f"color {color}, {pdb_id}_gpcr and sidechain")

    # Add labels (one per GRN, using first structure's CA as representative)
    if show_grn:
        lines.extend([
            "",
            "# GRN labels at CA positions (one label per GRN)",
            "set label_size, 36",
            "set label_color, black",
            "set label_font_id, 7",
            "set label_position, (2.0, 2.0, 0)",
            "",
        ])

        for grn in sorted(all_grn_residues.keys()):
            # Use first structure's residue as the label position (all CAs are aligned)
            residues = all_grn_residues[grn]
            pdb_id, chain, res_id, res_name = residues[0]
            lines.append(f'label {pdb_id}_gpcr and chain {chain} and resi {res_id} and name CA, "{grn}"')

    # Add interaction distance lines - ONE per binding residue (closest atom pair)
    lines.extend([
        "",
        "# =============================================================================",
        "# Ligand-sidechain interactions (pink dashes)",
        "# =============================================================================",
        "",
        "set dash_color, lightpink",
        "set dash_width, 1.5",
        "set dash_gap, 0.3",
        "",
    ])

    dist_count = 0
    for pdb_id in pdb_ids:
        if pdb_id not in interactions:
            continue

        data = interactions[pdb_id]
        chain = data["meta"]["chain"]
        ligand = data["meta"]["ligand"]

        # Use binding_residues - one distance line per UNIQUE residue with GRN
        binding_residues = data.get("binding_residues", [])

        # Deduplicate by res_id, keeping the one with smallest min_distance
        seen_residues = {}
        for res in binding_residues:
            if not res.get("grn"):
                continue
            res_id = res.get("res_id")
            if res_id not in seen_residues or res.get("min_distance", 999) < seen_residues[res_id].get("min_distance", 999):
                seen_residues[res_id] = res

        # Sort by distance and take closest 10
        grn_residues = sorted(seen_residues.values(), key=lambda x: x.get("min_distance", 999))

        for res in grn_residues[:10]:  # Top 10 closest unique GRN-annotated residues
            res_id = res.get("res_id")

            # Get the specific closest atom pair
            closest_lig_atom = res.get("closest_ligand_atom")
            closest_res_atom = res.get("closest_residue_atom")

            if not res_id or not closest_lig_atom or not closest_res_atom:
                continue

            dist_name = f"{pdb_id}_int_{dist_count}"
            # Distance from specific ligand atom to specific sidechain atom
            sel1 = f"({pdb_id}_gpcr and resn {ligand} and name {closest_lig_atom})"
            sel2 = f"({pdb_id}_gpcr and chain {chain} and resi {res_id} and name {closest_res_atom})"

            lines.append(f"distance {dist_name}, {sel1}, {sel2}")
            dist_count += 1

    # Hide ligand-sidechain distance labels
    lines.extend([
        "",
        "# Hide distance labels (show just dashes)",
        "hide labels, *_int_*",
        "",
    ])

    # Add water network interactions (water-water and water-sidechain) in red
    if show_water_network:
        lines.extend([
            "# =============================================================================",
            "# Water network interactions (red dashes)",
            "# =============================================================================",
            "",
            "set dash_color, red",
            "set dash_width, 1.0",
            "set dash_gap, 0.2",
            "",
        ])

        water_count = 0
        for pdb_id in pdb_ids:
            if pdb_id not in interactions:
                continue

            data = interactions[pdb_id]
            chain = data["meta"]["chain"]

            # Water-water interactions
            water_network = data.get("water_network", [])
            for ww in water_network:
                w1_id = ww["water1_id"]
                w2_id = ww["water2_id"]
                dist_name = f"{pdb_id}_wat_{water_count}"
                sel1 = f"({pdb_id}_gpcr and resn HOH and resi {w1_id} and name O)"
                sel2 = f"({pdb_id}_gpcr and resn HOH and resi {w2_id} and name O)"
                lines.append(f"distance {dist_name}, {sel1}, {sel2}")
                water_count += 1

            # Water-sidechain interactions
            water_sidechain = data.get("water_sidechain", [])
            # Deduplicate: one interaction per water-residue pair
            seen_ws = set()
            for ws in water_sidechain:
                key = (ws["water_id"], ws["res_id"])
                if key in seen_ws:
                    continue
                seen_ws.add(key)

                w_id = ws["water_id"]
                res_id = ws["res_id"]
                atom_name = ws["atom_name"]
                dist_name = f"{pdb_id}_wat_{water_count}"
                sel1 = f"({pdb_id}_gpcr and resn HOH and resi {w_id} and name O)"
                sel2 = f"({pdb_id}_gpcr and chain {chain} and resi {res_id} and name {atom_name})"
                lines.append(f"distance {dist_name}, {sel1}, {sel2}")
                water_count += 1

        # Hide water distance labels
        lines.extend([
            "",
            "# Hide water distance labels",
            "hide labels, *_wat_*",
            "",
        ])

    # Group all distance objects
    lines.extend([
        "# Group all distance objects",
        f"group {scene_name}_ligand_int, *_int_*",
    ])
    if show_water_network:
        lines.append(f"group {scene_name}_water_int, *_wat_*")
    lines.append("")

    # Store scene
    lines.extend([
        f"# Store as scene",
        f"scene {scene_name}, store",
        "",
        f'print("{title} pocket view ready.")',
        f'print("Structures shown: {", ".join(pdb_ids)}")',
        f'print("GRN-labeled residues: {len(all_grn_residues)}")',
        "",
        "deselect",
    ])

    # Write file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    print(f"  Saved: {output_path}")


def ensure_structures_downloaded(sp: StructureProcessor) -> None:
    """Download structures if they don't exist."""
    loader = StructureLoader()
    pdb_ids = list(STRUCTURES.keys())

    # Check which structures need downloading
    missing = []
    for pdb_id in pdb_ids:
        try:
            df = sp.load_entity(pdb_id)
            if df is None or df.empty:
                missing.append(pdb_id)
        except:
            missing.append(pdb_id)

    if missing:
        print(f"Downloading {len(missing)} structures: {', '.join(missing)}")
        loader.download_batch(missing, dataset_name="gpcr_structures", create_dataset=True, overwrite=False)


def main():
    """Generate interaction data and PyMOL scripts."""
    parser = argparse.ArgumentParser(description="Generate GPCR interaction data and PyMOL scripts")
    parser.add_argument(
        "--grn",
        action="store_true",
        default=True,
        help="Include GRN labels in pocket visualizations (default: True)",
    )
    parser.add_argument(
        "--no-grn",
        action="store_false",
        dest="grn",
        help="Exclude GRN labels from pocket visualizations",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("GPCR Ligand Interaction Analysis")
    print(f"GRN labels: {'enabled' if args.grn else 'disabled'}")
    print("=" * 60)

    protos.set_data_path(str(DATA_ROOT))
    sp = StructureProcessor()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Ensure structures are downloaded
    print("\n[0/3] Checking structures...")
    ensure_structures_downloaded(sp)

    # Calculate all interactions
    print("\n[1/3] Calculating ligand-protein interactions...")
    interactions = calculate_all_interactions(sp)

    # Save interactions to JSON
    json_path = OUTPUT_DIR / "gpcr_interactions.json"

    # Convert to JSON-serializable format
    json_data = {
        "_caveat": {
            "water_analysis": (
                "Crystallographic water counts CANNOT be compared between structures "
                "at different resolutions. Structures with 0 waters (e.g., 3SN6 at 3.2A) "
                "lack waters due to resolution limits, not biological absence. "
                "See Yuan et al. (2014) Nat Commun 5:4733 for MD-based water analysis."
            ),
        },
    }
    for pdb_id, data in interactions.items():
        json_data[pdb_id] = {
            "meta": data["meta"],
            "binding_residues": data["binding_residues"],
            "hydrogen_bonds": data["hydrogen_bonds"],
            "hydrophobic": data["hydrophobic"],
            "pi_stacking": data["pi_stacking"],
            "salt_bridges": data["salt_bridges"],
            "summary": data["summary"],
            "water_count": data.get("water_count", 0),
            "water_note": data.get("water_note"),
        }

    with open(json_path, "w") as f:
        json.dump(json_data, f, indent=2, default=str)
    print(f"\nInteraction data saved to: {json_path}")

    # Analyze water networks by receptor state
    print("\n" + "=" * 60)
    print("WATER NETWORK ANALYSIS BY RECEPTOR STATE")
    print("=" * 60)

    # Define active and inactive structures
    active_structures = ["3SN6", "4LDO", "2Y02"]  # Active states
    inactive_structures = ["2RH1", "3NY9", "2VT4"]  # Inactive states

    # Collect GRNs involved in water-sidechain interactions
    active_water_grns = {}  # grn -> list of (pdb_id, water_id, res_id)
    inactive_water_grns = {}

    for pdb_id, data in interactions.items():
        water_count = data.get("water_count", 0)

        # Skip structures with no waters
        if water_count == 0:
            print(f"\n[{pdb_id}] EXCLUDED - no crystallographic waters (resolution: {data['meta']['resolution']}Å)")
            continue

        # Get GRN map for binding residues
        grn_map = {int(r.get("res_id", 0)): r.get("grn", "") for r in data.get("binding_residues", []) if r.get("grn")}

        water_sidechain = data.get("water_sidechain", [])

        if pdb_id in active_structures:
            target_dict = active_water_grns
            state = "active"
        elif pdb_id in inactive_structures:
            target_dict = inactive_water_grns
            state = "inactive"
        else:
            continue  # Skip intermediate states for this analysis

        print(f"\n[{pdb_id}] {state.upper()} - {water_count} waters, {len(water_sidechain)} water-sidechain contacts")

        # Group by GRN
        for ws in water_sidechain:
            res_id = ws["res_id"]
            grn = grn_map.get(res_id, "")
            if grn:
                if grn not in target_dict:
                    target_dict[grn] = []
                target_dict[grn].append({
                    "pdb_id": pdb_id,
                    "water_id": ws["water_id"],
                    "res_id": res_id,
                    "atom": ws["atom_name"],
                    "distance": ws["distance"],
                })

    # Summary comparison
    print("\n" + "-" * 60)
    print("SUMMARY: GRN positions with water contacts")
    print("-" * 60)

    all_grns = sorted(set(active_water_grns.keys()) | set(inactive_water_grns.keys()))

    print(f"\n{'GRN':<10} {'Active':<20} {'Inactive':<20} {'Difference'}")
    print("-" * 60)

    for grn in all_grns:
        active_count = len(active_water_grns.get(grn, []))
        inactive_count = len(inactive_water_grns.get(grn, []))

        active_structs = set(x["pdb_id"] for x in active_water_grns.get(grn, []))
        inactive_structs = set(x["pdb_id"] for x in inactive_water_grns.get(grn, []))

        diff = ""
        if active_count > 0 and inactive_count == 0:
            diff = "ACTIVE ONLY"
        elif inactive_count > 0 and active_count == 0:
            diff = "INACTIVE ONLY"
        elif active_count > inactive_count:
            diff = f"+{active_count - inactive_count} in active"
        elif inactive_count > active_count:
            diff = f"+{inactive_count - active_count} in inactive"

        active_str = f"{active_count} ({','.join(active_structs)})" if active_count else "-"
        inactive_str = f"{inactive_count} ({','.join(inactive_structs)})" if inactive_count else "-"

        print(f"{grn:<10} {active_str:<20} {inactive_str:<20} {diff}")

    print("\n" + "=" * 60)

    # Generate PyMOL scripts
    print("\n[2/3] Generating PyMOL pocket scripts...")

    # Pocket 1: ADRB1
    generate_pocket_script(
        ["2VT4", "2Y02", "2Y04", "2Y00"],
        interactions,
        OUTPUT_DIR.parent.parent / "scripts" / "gpcr_visualization_pocket1.pml",
        "ADRB1",
        "adrb1_pocket",
        show_grn=args.grn,
    )

    # Pocket 2: ADRB2
    generate_pocket_script(
        ["2RH1", "3NY9", "3SN6", "4LDO"],
        interactions,
        OUTPUT_DIR.parent.parent / "scripts" / "gpcr_visualization_pocket2.pml",
        "ADRB2",
        "adrb2_pocket",
        show_grn=args.grn,
    )

    # Pocket 3: Active (with water network)
    generate_pocket_script(
        ["3SN6", "4LDO", "2Y02"],
        interactions,
        OUTPUT_DIR.parent.parent / "scripts" / "gpcr_visualization_pocket3.pml",
        "Active States",
        "active_pocket",
        show_grn=args.grn,
        show_water_network=True,
    )

    # Pocket 4: Inactive (with water network)
    generate_pocket_script(
        ["2RH1", "3NY9", "2VT4"],
        interactions,
        OUTPUT_DIR.parent.parent / "scripts" / "gpcr_visualization_pocket4.pml",
        "Inactive States",
        "inactive_pocket",
        show_grn=args.grn,
        show_water_network=True,
    )

    # Pocket 5: Agonist
    generate_pocket_script(
        ["3SN6", "4LDO", "2Y02", "2Y04", "2Y00"],
        interactions,
        OUTPUT_DIR.parent.parent / "scripts" / "gpcr_visualization_pocket5.pml",
        "Agonist-Bound",
        "agonist_pocket",
        show_grn=args.grn,
    )

    # Pocket 6: Inverse Agonist
    generate_pocket_script(
        ["2RH1", "3NY9"],
        interactions,
        OUTPUT_DIR.parent.parent / "scripts" / "gpcr_visualization_pocket6.pml",
        "Inverse Agonist-Bound",
        "inverse_agonist_pocket",
        show_grn=args.grn,
    )

    print("\n" + "=" * 60)
    print("DONE")
    print("=" * 60)
    print(f"\nOutput files:")
    print(f"  - {json_path}")
    print(f"  - scripts/gpcr_visualization_pocket[1-6].pml")
    print("\nTo visualize:")
    print("  pymol scripts/gpcr_visualization.pml")
    print("  @scripts/gpcr_visualization_pocket1.pml")


if __name__ == "__main__":
    main()
