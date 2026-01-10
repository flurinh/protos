#!/usr/bin/env python3
"""GPCR agonist/inverse agonist binding pocket comparison workflow.

This workflow uses the consolidated StructureProcessor API to:
1. Download adrenergic receptor structures (beta1 and beta2)
2. Annotate structures with GRN directly using annotate_with_grn()
3. Extract ligand-protein interactions using get_ligand_interactions()
4. Store all results in a property table with GRN columns

The output is a property table where:
- Rows = receptor structures (proteins)
- Columns = metadata + GRN positions
- Values = interaction distances (or NaN if no interaction)
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd

import protos
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.property import PropertyProcessor
from protos.processing.structure import StructureProcessor

DATA_ROOT = Path(__file__).resolve().parents[2] / "data"
OUTPUT_DIR = DATA_ROOT / "reports"

# Dataset and table names
STRUCTURE_DATASET = "adrenergic_agonist_antagonist_structures"
PROPERTY_TABLE = "adrenergic_binding_pocket_comparison"

# Interaction distance cutoff (Angstroms)
INTERACTION_CUTOFF = 4.5

# Structure metadata: PDB ID -> metadata
ADRENERGIC_STRUCTURES: Dict[str, Dict[str, Any]] = {
    # Beta2 Adrenergic Receptor (Human)
    "2RH1": {
        "receptor": "beta2_AR",
        "species": "human",
        "ligand": "carazolol",
        "ligand_code": "CAU",
        "ligand_type": "inverse_agonist",
        "receptor_state": "inactive",
        "chain": "A",
        "resolution": 2.4,
        "notes": "First human GPCR structure, T4L fusion",
    },
    "3NY9": {
        "receptor": "beta2_AR",
        "species": "human",
        "ligand": "ICI-118551",
        "ligand_code": "268",
        "ligand_type": "inverse_agonist",
        "receptor_state": "inactive",
        "chain": "A",
        "resolution": 2.8,
        "notes": "High-affinity inverse agonist",
    },
    "3SN6": {
        "receptor": "beta2_AR",
        "species": "human",
        "ligand": "BI-167107",
        "ligand_code": "P0G",
        "ligand_type": "full_agonist",
        "receptor_state": "active",
        "chain": "R",
        "resolution": 3.2,
        "notes": "Active state with Gs protein, 14A TM6 movement",
    },
    "4LDO": {
        "receptor": "beta2_AR",
        "species": "human",
        "ligand": "adrenaline",
        "ligand_code": "ALE",
        "ligand_type": "full_agonist",
        "receptor_state": "active",
        "chain": "A",
        "resolution": 3.2,
        "notes": "Endogenous agonist, nanobody stabilized",
    },
    # Beta1 Adrenergic Receptor (Turkey - high homology to human)
    "2VT4": {
        "receptor": "beta1_AR",
        "species": "turkey",
        "ligand": "cyanopindolol",
        "ligand_code": "CYP",
        "ligand_type": "antagonist",
        "receptor_state": "inactive",
        "chain": "A",
        "resolution": 2.7,
        "notes": "Reference antagonist structure",
    },
    "2Y03": {
        "receptor": "beta1_AR",
        "species": "turkey",
        "ligand": "isoprenaline",
        "ligand_code": "68H",
        "ligand_type": "full_agonist",
        "receptor_state": "active_like",
        "chain": "A",
        "resolution": 2.85,
        "notes": "1A pocket contraction, Ser5.42+5.46 bonds",
    },
    "2Y04": {
        "receptor": "beta1_AR",
        "species": "turkey",
        "ligand": "salbutamol",
        "ligand_code": "SAL",
        "ligand_type": "partial_agonist",
        "receptor_state": "intermediate",
        "chain": "A",
        "resolution": 3.05,
        "notes": "Only Ser5.42 interaction",
    },
    "2Y00": {
        "receptor": "beta1_AR",
        "species": "turkey",
        "ligand": "dobutamine",
        "ligand_code": "0GR",
        "ligand_type": "partial_agonist",
        "receptor_state": "intermediate",
        "chain": "A",
        "resolution": 2.6,
        "notes": "Lacks beta-hydroxyl group",
    },
}


def main() -> None:
    """Run the binding pocket comparison workflow."""
    print("=" * 60)
    print("GPCR Binding Pocket Comparison Workflow")
    print("=" * 60)

    # Initialize
    DATA_ROOT.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_ROOT))

    loader = StructureLoader()
    sp = StructureProcessor()
    prop = PropertyProcessor()

    pdb_ids = list(ADRENERGIC_STRUCTURES.keys())

    # Step 1: Download structures
    print("\n[1/4] Downloading structures...")
    print(f"  Downloading {len(pdb_ids)} structures: {', '.join(pdb_ids)}")
    loader.download_batch(pdb_ids, dataset_name=STRUCTURE_DATASET)

    # Step 2: Annotate with GRN and extract interactions
    print("\n[2/4] Annotating with GRN and extracting ligand interactions...")

    all_results: List[Dict[str, Any]] = []
    all_grns: Set[str] = set()

    for pdb_id in pdb_ids:
        metadata = ADRENERGIC_STRUCTURES[pdb_id]
        chain_id = metadata["chain"]
        ligand_code = metadata.get("ligand_code")

        print(f"\n  {pdb_id} ({metadata['ligand_type']}):")

        # Annotate structure with GRN (save=True so water networks get GRN info)
        try:
            sp.annotate_with_grn(pdb_id, chains=[chain_id], save=True)
            print(f"    GRN annotation: OK")
        except Exception as e:
            print(f"    GRN annotation failed: {e}")

        # List ligands
        ligands = sp.list_ligands(pdb_id)
        print(f"    Ligands found: {[l['res_name'] for l in ligands]}")

        # Get ligand interactions (use specific ligand code or largest)
        result = sp.get_ligand_interactions(
            pdb_id,
            ligand_id=ligand_code,
            cutoff=INTERACTION_CUTOFF,
        )

        if not result["ligand"]:
            # Try without specific ligand code
            result = sp.get_ligand_interactions(pdb_id, cutoff=INTERACTION_CUTOFF)

        if result["ligand"]:
            detected_ligand = result["ligand"]["res_name"]
            binding_residues = result["binding_residues"]

            # Filter to target chain
            binding_residues = [r for r in binding_residues if r["chain_id"] == chain_id]

            print(f"    Detected ligand: {detected_ligand}")
            print(f"    Binding residues: {len(binding_residues)}")

            # Build GRN -> distance mapping
            grn_distances: Dict[str, float] = {}
            for res in binding_residues:
                grn = res.get("grn")
                if grn:
                    dist = res["min_distance"]
                    if grn not in grn_distances or dist < grn_distances[grn]:
                        grn_distances[grn] = dist
                    all_grns.add(grn)

            grn_count = len(grn_distances)
            print(f"    GRN-mapped contacts: {grn_count}")

            # Get water contacts for this structure
            water_contacts = sp.get_water_contacts(pdb_id, protein_chain=chain_id)
            water_grns = set()
            for contact in water_contacts:
                if contact.get("grn"):
                    water_grns.add(contact["grn"])
            print(f"    Water contacts: {len(water_contacts)} residues ({len(water_grns)} with GRN)")

            # Store result
            all_results.append({
                "pdb_id": pdb_id,
                "detected_ligand": detected_ligand,
                "n_binding_residues": len(binding_residues),
                "n_grn_contacts": grn_count,
                "n_water_contacts": len(water_contacts),
                "water_contact_grns": water_grns,
                "grn_distances": grn_distances,
                **metadata,
            })
        else:
            print(f"    No ligand interactions found")
            # Still get water contacts
            water_contacts = sp.get_water_contacts(pdb_id, protein_chain=chain_id)
            water_grns = set()
            for contact in water_contacts:
                if contact.get("grn"):
                    water_grns.add(contact["grn"])

            all_results.append({
                "pdb_id": pdb_id,
                "detected_ligand": None,
                "n_binding_residues": 0,
                "n_grn_contacts": 0,
                "n_water_contacts": len(water_contacts),
                "water_contact_grns": water_grns,
                "grn_distances": {},
                **metadata,
            })

    # Compute water networks for all structures
    print("\n  Computing water networks...")
    water_network_output = sp.compute_water_networks(pdb_ids)
    water_network_results = water_network_output.get("structures", {})

    # Store water-mediated GRN pairs for each structure
    water_grn_pairs: Dict[str, Set[tuple]] = {}

    for pdb_id in pdb_ids:
        if pdb_id in water_network_results:
            networks = water_network_results[pdb_id].get("networks", [])
            summary = water_network_results[pdb_id].get("summary", {})

            # Extract GRN pairs from residue_pair_paths
            grn_pairs = set()
            all_network_grns = set()

            for network in networks:
                # Get GRNs from residue_pair_paths
                paths = network.get("residue_pair_paths", [])
                for path in paths:
                    source_grns = path.get("source_grn", [])
                    target_grns = path.get("target_grn", [])
                    for sg in source_grns:
                        all_network_grns.add(sg)
                    for tg in target_grns:
                        all_network_grns.add(tg)
                    # Create pairs
                    for sg in source_grns:
                        for tg in target_grns:
                            if sg != tg:
                                pair = tuple(sorted([sg, tg]))
                                grn_pairs.add(pair)

                # Also get GRNs from residues directly
                for res in network.get("residues", []):
                    for grn in res.get("grn_labels", []):
                        all_network_grns.add(grn)

            water_grn_pairs[pdb_id] = grn_pairs

            # Find the result and update it
            for result in all_results:
                if result["pdb_id"] == pdb_id:
                    result["n_water_networks"] = summary.get("network_count", 0)
                    result["n_bridging_waters"] = summary.get("bridging_water_count", 0)
                    result["water_network_grns"] = all_network_grns
                    result["water_grn_pairs"] = grn_pairs
                    break

            print(f"    {pdb_id}: {summary.get('network_count', 0)} networks, "
                  f"{summary.get('bridging_water_count', 0)} bridging waters, "
                  f"{len(all_network_grns)} GRN positions")

    # Step 3: Build property table
    print("\n[3/4] Building property table...")

    # Sort GRN positions
    def grn_sort_key(grn: str) -> tuple:
        parts = grn.split(".")
        try:
            return (int(parts[0]), int(parts[1]) if len(parts) > 1 else 0)
        except ValueError:
            return (999, 0)

    sorted_grns = sorted(all_grns, key=grn_sort_key)
    print(f"  Unique GRN positions: {len(sorted_grns)}")

    # Build DataFrame rows
    rows = []
    for result in all_results:
        row = {
            "pdb_id": result["pdb_id"],
            "receptor": result["receptor"],
            "species": result["species"],
            "ligand": result["ligand"],
            "ligand_type": result["ligand_type"],
            "receptor_state": result["receptor_state"],
            "resolution": result["resolution"],
            "detected_ligand": result["detected_ligand"],
            "n_binding_residues": result["n_binding_residues"],
            "n_grn_contacts": result["n_grn_contacts"],
            "n_water_contacts": result.get("n_water_contacts", 0),
            "n_water_networks": result.get("n_water_networks", 0),
            "n_bridging_waters": result.get("n_bridging_waters", 0),
        }

        # Add GRN columns for ligand contacts
        grn_distances = result["grn_distances"]
        for grn in sorted_grns:
            row[f"grn_{grn}"] = grn_distances.get(grn)

        rows.append(row)

    df = pd.DataFrame(rows)

    # Reorder columns
    meta_cols = [
        "pdb_id", "receptor", "species", "ligand", "ligand_type",
        "receptor_state", "resolution", "detected_ligand",
        "n_binding_residues", "n_grn_contacts",
        "n_water_contacts", "n_water_networks", "n_bridging_waters",
    ]
    grn_cols = [f"grn_{grn}" for grn in sorted_grns]
    df = df[meta_cols + grn_cols]

    # Save to CSV
    csv_path = DATA_ROOT / "property" / "tables" / f"{PROPERTY_TABLE}.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)
    print(f"  Saved: {csv_path}")

    # Step 4: Register with PropertyProcessor (optional)
    print("\n[4/4] Registering property table...")

    try:
        import json
        prop_rows = []
        for _, row_data in df.iterrows():
            prop_row = row_data.to_dict()
            prop_row["entity_name"] = prop_row["pdb_id"]
            # Scope must be JSON-encoded for PropertyProcessor
            prop_row["scope"] = json.dumps([{"format": "structure", "name": prop_row["pdb_id"]}])
            prop_rows.append(prop_row)

        prop.record_properties(PROPERTY_TABLE, prop_rows, allow_create=True)
        print(f"  Registered: {PROPERTY_TABLE}")
    except Exception as e:
        print(f"  PropertyProcessor registration skipped: {e}")
        print(f"  (Table saved to CSV: {csv_path})")

    # Display results
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    display_cols = ["pdb_id", "ligand_type", "receptor_state", "n_grn_contacts", "n_water_contacts", "n_bridging_waters"]
    print(df[display_cols].to_string(index=False))

    # Analyze binding pocket differences
    print("\n" + "-" * 60)
    print("BINDING POCKET ANALYSIS")
    print("-" * 60)

    grn_cols = [c for c in df.columns if c.startswith("grn_")]

    # Group by ligand type
    for ligand_type in df["ligand_type"].unique():
        subset = df[df["ligand_type"] == ligand_type]
        print(f"\n{ligand_type.upper()} ({len(subset)} structures):")

        # Common binding GRNs
        common = []
        for col in grn_cols:
            if subset[col].notna().all():
                grn = col.replace("grn_", "")
                avg_dist = subset[col].mean()
                common.append((grn, avg_dist))

        if common:
            print(f"  Common GRN contacts: {len(common)}")
            for grn, dist in sorted(common, key=lambda x: x[1])[:5]:
                print(f"    {grn}: {dist:.2f} A")

    # Compare agonist vs antagonist
    print("\n" + "-" * 60)
    print("AGONIST vs ANTAGONIST/INVERSE AGONIST")
    print("-" * 60)

    active_types = ["full_agonist", "partial_agonist"]
    inactive_types = ["inverse_agonist", "antagonist"]

    active_df = df[df["ligand_type"].isin(active_types)]
    inactive_df = df[df["ligand_type"].isin(inactive_types)]

    if len(active_df) > 0 and len(inactive_df) > 0:
        # Find differentially engaged positions
        agonist_only = []
        antagonist_only = []
        shared = []

        for col in grn_cols:
            grn = col.replace("grn_", "")
            active_has = active_df[col].notna().any()
            inactive_has = inactive_df[col].notna().any()

            if active_has and not inactive_has:
                agonist_only.append(grn)
            elif inactive_has and not active_has:
                antagonist_only.append(grn)
            elif active_has and inactive_has:
                active_dist = active_df[col].mean()
                inactive_dist = inactive_df[col].mean()
                shared.append((grn, active_dist, inactive_dist))

        print(f"\nAgonist-specific GRNs: {agonist_only if agonist_only else 'None'}")
        print(f"Antagonist-specific GRNs: {antagonist_only if antagonist_only else 'None'}")

        print(f"\nShared positions with distance differences:")
        for grn, act_d, inact_d in sorted(shared, key=lambda x: abs(x[1] - x[2]), reverse=True)[:10]:
            diff = act_d - inact_d
            direction = "closer in agonist" if diff < 0 else "closer in antagonist"
            print(f"  {grn}: agonist={act_d:.2f}A, antagonist={inact_d:.2f}A ({direction})")

    # Water Network GRN Comparison
    print("\n" + "-" * 60)
    print("WATER NETWORK GRN COMPARISON")
    print("-" * 60)

    # Group results by state
    active_results = [r for r in all_results if r["receptor_state"] in ["active", "active_like"]]
    inactive_results = [r for r in all_results if r["receptor_state"] in ["inactive", "intermediate"]]

    agonist_results = [r for r in all_results if r["ligand_type"] in ["full_agonist", "partial_agonist"]]
    antagonist_results = [r for r in all_results if r["ligand_type"] in ["inverse_agonist", "antagonist"]]

    def get_combined_grns(results: list) -> Set[str]:
        """Get union of all water network GRNs from a list of results."""
        combined = set()
        for r in results:
            combined.update(r.get("water_network_grns", set()))
        return combined

    def get_common_grns(results: list) -> Set[str]:
        """Get intersection of water network GRNs from results."""
        if not results:
            return set()
        grn_sets = [r.get("water_network_grns", set()) for r in results if r.get("water_network_grns")]
        if not grn_sets:
            return set()
        return set.intersection(*grn_sets) if len(grn_sets) > 1 else grn_sets[0]

    # Compare Active vs Inactive
    print("\n--- Active vs Inactive State ---")
    active_grns = get_combined_grns(active_results)
    inactive_grns = get_combined_grns(inactive_results)

    active_only = active_grns - inactive_grns
    inactive_only = inactive_grns - active_grns
    shared_grns = active_grns & inactive_grns

    print(f"Active state water network GRNs: {len(active_grns)}")
    print(f"Inactive state water network GRNs: {len(inactive_grns)}")
    print(f"Shared: {len(shared_grns)}")

    if active_only:
        print(f"\nActive-state specific GRNs (water networks): {sorted(active_only, key=grn_sort_key)}")
    else:
        print("\nActive-state specific GRNs: None (likely fewer waters in active state)")

    if inactive_only:
        print(f"Inactive-state specific GRNs (water networks): {sorted(inactive_only, key=grn_sort_key)[:20]}")
        if len(inactive_only) > 20:
            print(f"  ... and {len(inactive_only) - 20} more")

    # Compare Agonist vs Antagonist/Inverse Agonist
    print("\n--- Agonist vs Antagonist/Inverse Agonist ---")
    agonist_grns = get_combined_grns(agonist_results)
    antagonist_grns = get_combined_grns(antagonist_results)

    agonist_only = agonist_grns - antagonist_grns
    antagonist_only = antagonist_grns - agonist_grns
    shared_ligand_grns = agonist_grns & antagonist_grns

    print(f"Agonist water network GRNs: {len(agonist_grns)}")
    print(f"Antagonist/Inverse agonist water network GRNs: {len(antagonist_grns)}")
    print(f"Shared: {len(shared_ligand_grns)}")

    if agonist_only:
        print(f"\nAgonist-specific GRNs (water networks): {sorted(agonist_only, key=grn_sort_key)}")
    if antagonist_only:
        print(f"Antagonist-specific GRNs (water networks): {sorted(antagonist_only, key=grn_sort_key)[:20]}")
        if len(antagonist_only) > 20:
            print(f"  ... and {len(antagonist_only) - 20} more")

    # Key functional water positions
    print("\n--- Key Functional Water Positions ---")
    print("(GRNs involved in water networks that differ between states)")

    # TM3-TM6-TM7 interface (important for activation)
    key_tm_regions = {
        "TM3": ["3.32", "3.36", "3.40", "3.43", "3.46", "3.50"],
        "TM5": ["5.42", "5.43", "5.46", "5.50"],
        "TM6": ["6.44", "6.47", "6.48", "6.51", "6.52"],
        "TM7": ["7.45", "7.49", "7.50", "7.53"],
    }

    for tm, positions in key_tm_regions.items():
        in_inactive = [p for p in positions if p in inactive_grns]
        in_active = [p for p in positions if p in active_grns]
        if in_inactive or in_active:
            print(f"  {tm}: Inactive={in_inactive}, Active={in_active}")

    # Save analysis report
    report_path = OUTPUT_DIR / "binding_pocket_analysis.txt"
    with open(report_path, "w") as f:
        f.write("GPCR Binding Pocket Comparison\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Structures analyzed: {len(df)}\n")
        f.write(f"GRN positions detected: {len(sorted_grns)}\n\n")
        f.write("Full table:\n")
        f.write(df.to_string())
        f.write("\n")
    print(f"\n  Report saved: {report_path}")

    print("\n" + "=" * 60)
    print("WORKFLOW COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
