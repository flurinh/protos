#!/usr/bin/env python3
"""
GPCR Mechanism Analysis: Water Networks and Activation States

Based on Yuan et al. (2014) Nat Commun 5:4733 - "Activation of G-protein-coupled
receptors correlates with the formation of a continuous internal water pathway"

Key hypotheses to test:
1. Active states have more water contacts at NPxxY region (7.50, 7.53)
2. Active states show water contacts at TM6 positions (6.48-6.55) due to outward movement
3. Ionic lock region (3.50, 6.30) shows state-dependent water accessibility
4. Agonists vs inverse agonists show different water-mediated interaction patterns
"""
from __future__ import annotations

import sys
from pathlib import Path
from collections import defaultdict
import numpy as np

# Add src to path
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.structure import StructureProcessor

DATA_ROOT = REPO_ROOT / "data"

# Structure classifications with DOI references for verification
# Format: PDB: {chain, ligand, ligand_type, state, receptor, doi, notes}
STRUCTURES = {
    "2RH1": {
        "chain": "A", "ligand": "CAU", "ligand_type": "inverse_agonist",
        "state": "inactive", "receptor": "ADRB2",
        "doi": "10.1126/science.1150577",
        "notes": "Carazolol-bound, first high-res GPCR structure (Cherezov 2007)"
    },
    "3NY9": {
        "chain": "A", "ligand": "JSZ", "ligand_type": "inverse_agonist",
        "state": "inactive", "receptor": "ADRB2",
        "doi": "10.1021/ja105108q",
        "notes": "ICI 118,551-bound inverse agonist (Wacker 2010)"
    },
    "3SN6": {
        "chain": "R", "ligand": "P0G", "ligand_type": "full_agonist",
        "state": "active", "receptor": "ADRB2",
        "doi": "10.1038/nature10361",
        "notes": "BI-167107 + Gs protein, FULLY ACTIVE (Rasmussen 2011)"
    },
    "4LDO": {
        "chain": "A", "ligand": "ALE", "ligand_type": "full_agonist",
        "state": "active", "receptor": "ADRB2",
        "doi": "10.1038/nature12572",
        "notes": "Adrenaline + Nb6B9 nanobody, active state (Ring 2013)"
    },
    "2VT4": {
        "chain": "A", "ligand": "P32", "ligand_type": "antagonist",
        "state": "inactive", "receptor": "ADRB1",
        "doi": "10.1038/nature07101",
        "notes": "Cyanopindolol-bound, inactive (Warne 2008)"
    },
    "2Y02": {
        "chain": "A", "ligand": "WHJ", "ligand_type": "full_agonist",
        "state": "active_like", "receptor": "ADRB1",
        "doi": "10.1038/nature09746",
        "notes": "Isoprenaline-bound, agonist-induced changes but NO G protein (Warne 2011)"
    },
    "2Y04": {
        "chain": "A", "ligand": "68H", "ligand_type": "partial_agonist",
        "state": "intermediate", "receptor": "ADRB1",
        "doi": "10.1038/nature09746",
        "notes": "Salbutamol-bound, intermediate state (Warne 2011)"
    },
    "2Y00": {
        "chain": "A", "ligand": "Y00", "ligand_type": "partial_agonist",
        "state": "intermediate", "receptor": "ADRB1",
        "doi": "10.1038/nature09746",
        "notes": "Dobutamine-bound, intermediate state (Warne 2011)"
    },
}

# Helper to maintain backward compatibility with tuple access
def get_structure_tuple(pdb_id):
    """Get (chain, ligand, ligand_type, state, receptor) tuple for a structure."""
    s = STRUCTURES[pdb_id]
    return (s["chain"], s["ligand"], s["ligand_type"], s["state"], s["receptor"])

# Key GRN positions from Yuan et al. for GPCR activation
KEY_GRNS = {
    # NPxxY motif - the "gate"
    "7.49": "NPxxY (N)",
    "7.50": "NPxxY (P) - conserved proline",
    "7.53": "NPxxY (Y) - key rotamer switch",

    # DRY motif / Ionic lock
    "3.49": "DRY (D)",
    "3.50": "DRY (R) - ionic lock partner",
    "3.51": "DRY (Y)",

    # TM6 - moves outward in active state
    "6.30": "Ionic lock partner (E/K)",
    "6.34": "TM6 activation",
    "6.37": "TM6 activation",
    "6.40": "TM6 hydrophobic layer",
    "6.44": "TM6 connector switch",
    "6.48": "Toggle switch (W/F)",
    "6.51": "TM6 binding pocket",
    "6.52": "TM6 binding pocket",
    "6.55": "TM6 binding pocket (N)",

    # Allosteric sodium site
    "2.50": "Na+ allosteric site (D)",

    # TM5 connection to NPxxY
    "5.58": "Connects to Y7.53",

    # Binding pocket positions
    "3.32": "Orthosteric binding (D)",
    "3.33": "Orthosteric binding",
    "5.42": "Binding pocket",
    "5.43": "Binding pocket (S)",
    "5.46": "Binding pocket",
}


def analyze_water_networks(sp: StructureProcessor):
    """Analyze water network contacts at key GRN positions."""

    print("=" * 80)
    print("GPCR ACTIVATION MECHANISM ANALYSIS")
    print("Based on Yuan et al. (2014) Nat Commun 5:4733")
    print("=" * 80)

    WATER_CUTOFF = 3.5  # Angstroms for H-bond

    # Store results by state and ligand type
    results = {}

    for pdb_id, struct_info in STRUCTURES.items():
        chain = struct_info["chain"]
        ligand = struct_info["ligand"]
        lig_type = struct_info["ligand_type"]
        state = struct_info["state"]
        receptor = struct_info["receptor"]
        print(f"\n{'='*60}")
        print(f"[{pdb_id}] {receptor} - {state.upper()} - {lig_type}")
        print("=" * 60)

        df = sp.load_entity(pdb_id)
        if df is None:
            print("  ERROR: Could not load structure")
            continue
        df = df.reset_index()

        # Check for waters
        waters = df[df["res_name3l"] == "HOH"]
        water_count = len(waters["auth_seq_id"].unique()) if not waters.empty else 0

        if water_count == 0:
            print(f"  EXCLUDED: No crystallographic waters")
            continue

        print(f"  Waters: {water_count}")

        # Get water oxygen positions
        water_oxygens = waters[waters["atom_name"] == "O"][["auth_seq_id", "x", "y", "z"]]
        water_coords = water_oxygens[["x", "y", "z"]].values
        water_ids = water_oxygens["auth_seq_id"].values

        # Get GRN-annotated residues on receptor chain
        chain_df = df[(df["auth_chain_id"] == chain) & (df["grn"].notna())]

        # Build GRN to residue mapping
        grn_residues = {}
        for _, row in chain_df[["auth_seq_id", "grn"]].drop_duplicates().iterrows():
            grn = row["grn"]
            res_id = int(row["auth_seq_id"])
            if grn not in grn_residues:
                grn_residues[grn] = []
            grn_residues[grn].append(res_id)

        # Analyze water contacts at key positions
        water_contacts = {}

        for grn, description in KEY_GRNS.items():
            if grn not in grn_residues:
                continue

            res_ids = grn_residues[grn]

            # Get heavy atoms for this residue (exclude H)
            res_atoms = chain_df[
                (chain_df["auth_seq_id"].isin(res_ids)) &
                (~chain_df["element"].isin(["H"]))
            ]

            if res_atoms.empty:
                continue

            res_coords = res_atoms[["x", "y", "z"]].values
            res_atom_names = res_atoms["atom_name"].values

            # Find closest water to any atom of this residue
            min_dist = float('inf')
            closest_water = None
            closest_atom = None
            contacts = []

            for i, water_coord in enumerate(water_coords):
                dists = np.sqrt(np.sum((res_coords - water_coord)**2, axis=1))
                min_idx = np.argmin(dists)

                if dists[min_idx] < WATER_CUTOFF:
                    contacts.append({
                        "water_id": int(water_ids[i]),
                        "atom": res_atom_names[min_idx],
                        "distance": float(dists[min_idx])
                    })

                if dists[min_idx] < min_dist:
                    min_dist = dists[min_idx]
                    closest_water = water_ids[i]
                    closest_atom = res_atom_names[min_idx]

            water_contacts[grn] = {
                "description": description,
                "min_distance": min_dist,
                "closest_water": closest_water,
                "closest_atom": closest_atom,
                "n_contacts": len(contacts),
                "has_contact": min_dist <= WATER_CUTOFF,
            }

        results[pdb_id] = {
            "state": state,
            "ligand_type": lig_type,
            "receptor": receptor,
            "water_count": water_count,
            "contacts": water_contacts,
        }

        # Print contacts for this structure
        print(f"\n  Key GRN water contacts (<{WATER_CUTOFF}Å):")
        for grn in sorted(water_contacts.keys()):
            c = water_contacts[grn]
            if c["has_contact"]:
                print(f"    {grn:6s} {c['description'][:30]:32s} {c['min_distance']:.2f}Å ({c['n_contacts']} contacts)")

    # Summary comparison
    print("\n" + "=" * 80)
    print("COMPARATIVE ANALYSIS BY RECEPTOR STATE")
    print("=" * 80)

    # Group by state
    active_structs = {k: v for k, v in results.items() if v["state"] in ["active", "active_like"]}
    inactive_structs = {k: v for k, v in results.items() if v["state"] == "inactive"}

    print(f"\nActive structures: {list(active_structs.keys())}")
    print(f"Inactive structures: {list(inactive_structs.keys())}")

    # Compare water contacts at each key position
    print("\n" + "-" * 80)
    print(f"{'GRN':<8} {'Description':<35} {'Active':<15} {'Inactive':<15} {'Interpretation'}")
    print("-" * 80)

    all_grns = set()
    for data in results.values():
        all_grns.update(data["contacts"].keys())

    for grn in sorted(all_grns):
        desc = KEY_GRNS.get(grn, "")[:35]

        # Count contacts in each state
        active_contacts = sum(1 for k, v in active_structs.items()
                             if grn in v["contacts"] and v["contacts"][grn]["has_contact"])
        inactive_contacts = sum(1 for k, v in inactive_structs.items()
                               if grn in v["contacts"] and v["contacts"][grn]["has_contact"])

        active_total = sum(1 for k, v in active_structs.items() if grn in v["contacts"])
        inactive_total = sum(1 for k, v in inactive_structs.items() if grn in v["contacts"])

        active_str = f"{active_contacts}/{active_total}" if active_total > 0 else "-"
        inactive_str = f"{inactive_contacts}/{inactive_total}" if inactive_total > 0 else "-"

        # Interpretation
        interp = ""
        if active_contacts > 0 and inactive_contacts == 0:
            interp = "** ACTIVE ONLY **"
        elif inactive_contacts > 0 and active_contacts == 0:
            interp = "** INACTIVE ONLY **"
        elif active_contacts > inactive_contacts:
            interp = f"+{active_contacts - inactive_contacts} active"
        elif inactive_contacts > active_contacts:
            interp = f"+{inactive_contacts - active_contacts} inactive"

        if active_contacts > 0 or inactive_contacts > 0:
            print(f"{grn:<8} {desc:<35} {active_str:<15} {inactive_str:<15} {interp}")

    # NPxxY analysis
    print("\n" + "=" * 80)
    print("HYPOTHESIS 1: NPxxY REGION WATER CONTACTS")
    print("Yuan et al.: Y7.53 rotamer controls water pathway gate")
    print("=" * 80)

    npxxy_grns = ["7.49", "7.50", "7.53"]
    for grn in npxxy_grns:
        print(f"\n{grn} ({KEY_GRNS.get(grn, '')}):")
        for pdb_id, data in results.items():
            if grn in data["contacts"]:
                c = data["contacts"][grn]
                status = "CONTACT" if c["has_contact"] else "no contact"
                print(f"  {pdb_id} ({data['state']:12s}): {c['min_distance']:.2f}Å - {status}")

    # TM6 analysis
    print("\n" + "=" * 80)
    print("HYPOTHESIS 2: TM6 WATER ACCESS (outward movement in active)")
    print("Yuan et al.: TM6 moves outward, exposing cytoplasmic region")
    print("=" * 80)

    tm6_grns = ["6.30", "6.34", "6.40", "6.44", "6.48", "6.51", "6.52", "6.55"]
    for grn in tm6_grns:
        if grn not in KEY_GRNS:
            continue

        active_with = []
        inactive_with = []

        for pdb_id, data in results.items():
            if grn in data["contacts"] and data["contacts"][grn]["has_contact"]:
                if data["state"] in ["active", "active_like"]:
                    active_with.append(pdb_id)
                elif data["state"] == "inactive":
                    inactive_with.append(pdb_id)

        if active_with or inactive_with:
            print(f"\n{grn} ({KEY_GRNS.get(grn, '')}):")
            print(f"  Active with contact: {active_with if active_with else 'none'}")
            print(f"  Inactive with contact: {inactive_with if inactive_with else 'none'}")

    # Ligand type analysis
    print("\n" + "=" * 80)
    print("HYPOTHESIS 3: AGONIST vs INVERSE AGONIST WATER PATTERNS")
    print("=" * 80)

    agonist_structs = {k: v for k, v in results.items()
                       if v["ligand_type"] in ["full_agonist", "partial_agonist"]}
    inverse_structs = {k: v for k, v in results.items()
                       if v["ligand_type"] == "inverse_agonist"}

    print(f"\nAgonist structures: {list(agonist_structs.keys())}")
    print(f"Inverse agonist structures: {list(inverse_structs.keys())}")

    # Find positions unique to each
    agonist_only_positions = set()
    inverse_only_positions = set()

    for grn in sorted(all_grns):
        agonist_contacts = sum(1 for k, v in agonist_structs.items()
                              if grn in v["contacts"] and v["contacts"][grn]["has_contact"])
        inverse_contacts = sum(1 for k, v in inverse_structs.items()
                              if grn in v["contacts"] and v["contacts"][grn]["has_contact"])

        if agonist_contacts > 0 and inverse_contacts == 0:
            agonist_only_positions.add(grn)
        elif inverse_contacts > 0 and agonist_contacts == 0:
            inverse_only_positions.add(grn)

    print(f"\nPositions with water contacts ONLY in agonist-bound:")
    for grn in sorted(agonist_only_positions):
        print(f"  {grn}: {KEY_GRNS.get(grn, '')}")

    print(f"\nPositions with water contacts ONLY in inverse agonist-bound:")
    for grn in sorted(inverse_only_positions):
        print(f"  {grn}: {KEY_GRNS.get(grn, '')}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY: CHARACTERISTIC POSITIONS FOR ACTIVATION STATE")
    print("=" * 80)

    print("""
Based on our crystallographic water analysis and Yuan et al. (2014):

1. ACTIVE STATE MARKERS (water contacts present):
   - These positions show water contacts preferentially in active structures
   - Consistent with opening of the internal water pathway upon activation

2. INACTIVE STATE MARKERS (water contacts present):
   - These positions show water contacts preferentially in inactive structures
   - May indicate stable water molecules in the closed/blocked state

3. CAVEATS:
   - 3SN6 excluded (no crystallographic waters at 3.2Å resolution)
   - Crystal waters represent static snapshot, not full dynamics
   - MD simulations (as in Yuan et al.) needed for complete picture

4. VERIFIABLE PREDICTIONS:
   - Position 6.55 should show state-dependent water accessibility
   - NPxxY region (7.50, 7.53) should correlate with activation
   - TM6 positions (6.48-6.55) should become more water-accessible in active state
""")

    return results


def main():
    protos.set_data_path(str(DATA_ROOT))
    sp = StructureProcessor()

    results = analyze_water_networks(sp)


if __name__ == "__main__":
    main()
