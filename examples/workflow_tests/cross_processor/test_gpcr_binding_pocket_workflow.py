#!/usr/bin/env python3
"""
GPCR Binding Pocket Comparison Workflow

This workflow analyzes 8 adrenergic receptor structures to compare:
- Agonist vs antagonist/inverse agonist binding
- Active vs inactive state conformations
- Water network differences
- GRN-annotated ligand-protein interactions

Structures:
- ADRB2: 2RH1, 3NY9 (inverse agonist), 3SN6, 4LDO (agonist)
- ADRB1: 2VT4 (antagonist), 2Y02, 2Y04, 2Y00 (agonist/partial agonist)

Output:
- GRN-annotated structures
- Ligand interaction data with GRN labels
- Water network analysis
- JSON data for PyMOL visualization

IMPORTANT CAVEAT - Water Network Analysis:
    Crystallographic water counts CANNOT be directly compared between structures
    at different resolutions. For example:
    - 2RH1 (inactive): 2.4A resolution, 48 crystallographic waters
    - 3SN6 (active): 3.2A resolution, 0 crystallographic waters

    The absence of waters in lower-resolution structures reflects crystallographic
    limitations, NOT the actual hydration state of the protein.

    For proper water network analysis in GPCRs, see:
    Yuan et al. (2014) Nature Communications 5:4733
    "Activation of G-protein-coupled receptors correlates with the formation
    of a continuous internal water pathway"

    Key findings from Yuan et al.:
    - Used MD simulations (not crystal structure water counts)
    - Inactive states have hydrophobic barriers BLOCKING water penetration
    - Active states form continuous water channels through the receptor
    - The hydrophobic layer near NPxxY motif opens during activation

    Structures analyzed in Yuan et al. (2014):
    - Beta2-AR: 2RH1 (inactive), 3SN6 (active) - same as this workflow
    - A2A-AR: 4EIY (inactive), 2YDV (intermediate)
    - Rhodopsin: 1U19 (inactive), 2X72 (active)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Dict, List, Any, Optional

import pandas as pd

# Setup paths
REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.structure import StructureProcessor
from protos.processing.grn import GRNProcessor
from protos.analysis.structure_ligand_analysis import calculate_ligand_interactions
from protos.analysis.structure_water_networks import summarize_water_networks

# Configuration
DATA_ROOT = REPO_ROOT / "data"
OUTPUT_DIR = DATA_ROOT / "visualizations"

# Structure metadata
STRUCTURES: Dict[str, Dict[str, Any]] = {
    "2RH1": {
        "chain": "A",
        "ligand": "CAU",
        "ligand_name": "carazolol",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "receptor": "ADRB2",
    },
    "3NY9": {
        "chain": "A",
        "ligand": "JSZ",
        "ligand_name": "ICI-118551",
        "ligand_type": "inverse_agonist",
        "state": "inactive",
        "receptor": "ADRB2",
    },
    "3SN6": {
        "chain": "R",
        "ligand": "P0G",
        "ligand_name": "BI-167107",
        "ligand_type": "full_agonist",
        "state": "active",
        "receptor": "ADRB2",
    },
    "4LDO": {
        "chain": "A",
        "ligand": "ALE",
        "ligand_name": "adrenaline",
        "ligand_type": "full_agonist",
        "state": "active",
        "receptor": "ADRB2",
    },
    "2VT4": {
        "chain": "A",
        "ligand": "P32",
        "ligand_name": "cyanopindolol",
        "ligand_type": "antagonist",
        "state": "inactive",
        "receptor": "ADRB1",
    },
    "2Y02": {
        "chain": "A",
        "ligand": "WHJ",
        "ligand_name": "isoprenaline",
        "ligand_type": "full_agonist",
        "state": "active_like",
        "receptor": "ADRB1",
    },
    "2Y04": {
        "chain": "A",
        "ligand": "68H",
        "ligand_name": "salbutamol",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "receptor": "ADRB1",
    },
    "2Y00": {
        "chain": "A",
        "ligand": "Y00",
        "ligand_name": "dobutamine",
        "ligand_type": "partial_agonist",
        "state": "intermediate",
        "receptor": "ADRB1",
    },
}

INTERACTION_CUTOFF = 4.0
GRN_TABLE_NAME = "gpcr_binding_pocket_grn"
DATASET_NAME = "gpcr_binding_pocket"


def ensure_data_root() -> Path:
    """Initialize protos data path."""
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_ROOT))
    return DATA_ROOT


def download_structures(sp: StructureProcessor) -> List[str]:
    """Download structures if not already present."""
    loader = StructureLoader()
    pdb_ids = list(STRUCTURES.keys())

    # Check which are missing
    missing = []
    available = []
    for pdb_id in pdb_ids:
        try:
            df = sp.load_entity(pdb_id)
            if df is None or df.empty:
                missing.append(pdb_id)
            else:
                available.append(pdb_id)
        except Exception:
            missing.append(pdb_id)

    if missing:
        print(f"  Downloading {len(missing)} structures: {', '.join(missing)}")
        loader.download_batch(
            missing,
            dataset_name=DATASET_NAME,
            create_dataset=True,
            overwrite=False,
        )
        available.extend(missing)

    return available


def annotate_with_grn(sp: StructureProcessor, pdb_ids: List[str]) -> Dict[str, Dict[int, str]]:
    """Annotate structures with GRN positions."""
    grn_maps = {}

    print("  Annotating structures with GRN...")

    for pdb_id in pdb_ids:
        chain = STRUCTURES[pdb_id]["chain"]

        try:
            # Use the structure processor's GRN annotation
            sp.annotate_with_grn(pdb_id, chains=[chain])

            # Load the annotated structure
            df = sp.load_entity(pdb_id)
            if df is None:
                continue

            df = df.reset_index()

            # Extract GRN mapping - filter by chain to avoid antibody/fusion protein rows
            grn_map = {}
            if "grn" in df.columns:
                chain_df = df[df["auth_chain_id"] == chain]
                grn_df = chain_df[chain_df["grn"].notna()][["auth_seq_id", "grn"]].drop_duplicates()
                for _, row in grn_df.iterrows():
                    grn_map[int(row["auth_seq_id"])] = row["grn"]

            grn_maps[pdb_id] = grn_map
            print(f"    {pdb_id}: {len(grn_map)} GRN positions")

        except Exception as e:
            print(f"    {pdb_id}: GRN annotation failed - {e}")
            grn_maps[pdb_id] = {}

    return grn_maps


def calculate_interactions(
    sp: StructureProcessor,
    pdb_ids: List[str],
    grn_maps: Dict[str, Dict[int, str]],
) -> Dict[str, Dict[str, Any]]:
    """Calculate ligand-protein interactions for all structures."""
    results = {}

    for pdb_id in pdb_ids:
        meta = STRUCTURES[pdb_id]
        chain = meta["chain"]
        ligand_code = meta["ligand"]
        grn_map = grn_maps.get(pdb_id, {})

        try:
            # Load structure
            df = sp.load_entity(pdb_id)
            if df is None:
                print(f"    {pdb_id}: Could not load structure")
                continue

            df = df.reset_index()

            # Get ligand atoms
            ligand_atoms = df[
                (df["group"] == "HETATM") & (df["res_name3l"] == ligand_code)
            ]

            if ligand_atoms.empty:
                print(f"    {pdb_id}: Ligand {ligand_code} not found")
                continue

            # Calculate interactions
            interactions = calculate_ligand_interactions(
                sp, pdb_id, ligand_atoms, detailed=True, cutoff=INTERACTION_CUTOFF
            )

            binding_residues = interactions.get("binding_residues", [])

            # Add GRN annotations to binding residues
            for res in binding_residues:
                res_id = int(res.get("res_id", 0))
                res["grn"] = grn_map.get(res_id, "")

            # Count GRN-annotated residues
            grn_count = sum(1 for r in binding_residues if r.get("grn"))

            results[pdb_id] = {
                "meta": meta,
                "binding_residues": binding_residues,
                "hydrogen_bonds": interactions.get("hydrogen_bonds", []),
                "hydrophobic": interactions.get("hydrophobic", []),
                "pi_stacking": interactions.get("pi_stacking", []),
                "salt_bridges": interactions.get("salt_bridges", []),
                "water_mediated": interactions.get("water_mediated", []),
                "summary": interactions.get("summary", {}),
            }

            print(
                f"    {pdb_id}: {len(binding_residues)} binding residues, "
                f"{grn_count} with GRN, "
                f"{len(interactions.get('hydrogen_bonds', []))} H-bonds"
            )

        except Exception as e:
            print(f"    {pdb_id}: Interaction calculation failed - {e}")
            import traceback
            traceback.print_exc()

    return results


def analyze_water_networks(
    sp: StructureProcessor, pdb_ids: List[str]
) -> Dict[str, Dict[str, Any]]:
    """Analyze water networks in binding pockets.

    CAVEAT: Crystallographic water counts cannot be compared between structures
    at different resolutions. See module docstring for details.
    """
    print("  Computing water networks...")
    print("  WARNING: Crystallographic water counts vary with resolution.")
    print("           See Yuan et al. (2014) Nat Commun for MD-based analysis.")

    try:
        result = sp.compute_water_networks(
            pdb_ids,
            residue_cutoff=3.4,
            water_water_cutoff=3.4,
            hydrogen_bond_cutoff=3.2,
        )

        structures = result.get("structures", {})
        summary = summarize_water_networks(structures)

        for pdb_id in pdb_ids:
            if pdb_id in summary:
                s = summary[pdb_id].get("summary", {})
                print(
                    f"    {pdb_id}: {s.get('network_count', 0)} networks, "
                    f"{s.get('water_count', 0)} waters"
                )

        return summary

    except Exception as e:
        print(f"    Water network analysis failed: {e}")
        return {}


def save_results(
    interactions: Dict[str, Dict[str, Any]],
    water_networks: Dict[str, Dict[str, Any]],
    output_path: Path,
) -> None:
    """Save all results to JSON."""
    output = {
        "metadata": {
            "structures": STRUCTURES,
            "interaction_cutoff": INTERACTION_CUTOFF,
        },
        "interactions": {},
        "water_networks": {},
    }

    # Add interaction data
    for pdb_id, data in interactions.items():
        output["interactions"][pdb_id] = {
            "meta": data["meta"],
            "binding_residues": data["binding_residues"],
            "hydrogen_bonds": data["hydrogen_bonds"],
            "hydrophobic": data["hydrophobic"],
            "pi_stacking": data["pi_stacking"],
            "salt_bridges": data["salt_bridges"],
            "summary": data["summary"],
        }

    # Add water network data
    for pdb_id, data in water_networks.items():
        output["water_networks"][pdb_id] = {
            "summary": data.get("summary", {}),
            "network_count": data.get("summary", {}).get("network_count", 0),
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\nResults saved to: {output_path}")


def print_summary(interactions: Dict[str, Dict[str, Any]]) -> None:
    """Print summary of binding pocket analysis."""
    print("\n" + "=" * 70)
    print("BINDING POCKET COMPARISON SUMMARY")
    print("=" * 70)

    # Group by category
    categories = {
        "Active States": ["3SN6", "4LDO", "2Y02"],
        "Inactive States": ["2RH1", "3NY9", "2VT4"],
        "Intermediate States": ["2Y04", "2Y00"],
    }

    for category, pdb_ids in categories.items():
        print(f"\n{category}:")
        for pdb_id in pdb_ids:
            if pdb_id not in interactions:
                continue
            data = interactions[pdb_id]
            meta = data["meta"]
            binding = data["binding_residues"]
            grn_residues = [r for r in binding if r.get("grn")]

            # Get key GRN positions
            key_grns = set()
            for r in grn_residues[:10]:
                grn = r.get("grn", "")
                if grn:
                    key_grns.add(grn)

            print(
                f"  {pdb_id} ({meta['receptor']}, {meta['ligand_name']}): "
                f"{len(binding)} contacts, GRN: {', '.join(sorted(key_grns)[:5])}"
            )


def main() -> int:
    """Run the GPCR binding pocket workflow."""
    print("=" * 70)
    print("GPCR BINDING POCKET COMPARISON WORKFLOW")
    print("=" * 70)

    ensure_data_root()
    sp = StructureProcessor()

    # Step 1: Download structures
    print("\n[1/5] Downloading structures...")
    pdb_ids = download_structures(sp)
    print(f"  Available: {len(pdb_ids)} structures")

    # Step 2: Annotate with GRN
    print("\n[2/5] Annotating with GRN...")
    grn_maps = annotate_with_grn(sp, pdb_ids)

    # Step 3: Calculate ligand interactions
    print("\n[3/5] Calculating ligand-protein interactions...")
    interactions = calculate_interactions(sp, pdb_ids, grn_maps)

    # Step 4: Analyze water networks
    print("\n[4/5] Analyzing water networks...")
    water_networks = analyze_water_networks(sp, pdb_ids)

    # Step 5: Save results
    print("\n[5/5] Saving results...")
    output_path = OUTPUT_DIR / "gpcr_binding_pocket_data.json"
    save_results(interactions, water_networks, output_path)

    # Print summary
    print_summary(interactions)

    print("\n" + "=" * 70)
    print("WORKFLOW COMPLETE")
    print("=" * 70)
    print(f"\nOutput: {output_path}")
    print("\nTo generate PyMOL visualizations, run:")
    print("  python scripts/prepare_gpcr_interactions.py")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
