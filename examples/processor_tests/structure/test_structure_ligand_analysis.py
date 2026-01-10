#!/usr/bin/env python3
"""StructureProcessor ligand analysis demonstration.

This script demonstrates ligand analysis methods:
- Listing ligands in structures
- Summarizing ligand information
- Analyzing ligand-protein interactions
- Ion contact analysis
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def ensure_data_root() -> Path:
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


# Structure with interesting ligands (adenosine receptor with antagonist)
TARGET_STRUCTURE = "5d5a"


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.structure import StructureProcessor
    from protos.io.ingest.structure_loader import StructureLoader

    print("=== StructureProcessor Ligand Analysis Demo ===\n")

    processor = StructureProcessor()
    loader = StructureLoader(processor=processor)

    # Ensure structure is available
    print(f"1. Loading structure '{TARGET_STRUCTURE}'...")
    structure_df = processor.load_entity(TARGET_STRUCTURE)

    if structure_df is None:
        print(f"   Downloading {TARGET_STRUCTURE}...")
        success, failed = loader.download_batch(
            [TARGET_STRUCTURE],
            dataset_name="ligand_analysis_structures",
            create_dataset=True,
        )
        if failed:
            print(f"   Failed to download: {failed}")
            return
        structure_df = processor.load_entity(TARGET_STRUCTURE)

    if structure_df is None:
        print("   Could not load structure")
        return

    print(f"   Loaded structure with {len(structure_df)} atoms\n")

    # 2. List all ligands
    print("2. Listing Ligands")
    print("-" * 40)

    ligands = processor.list_ligands(
        TARGET_STRUCTURE,
        exclude_common=True,  # Exclude water, ions, common buffers
        min_atoms=5,          # Only ligands with 5+ atoms
    )

    if ligands:
        print(f"   Found {len(ligands)} ligand(s):\n")
        for lig in ligands:
            print(f"   - {lig['res_name']} (chain {lig['chain_id']}, res {lig['res_id']})")
            print(f"     Atoms: {lig['num_atoms']}, Centroid: {lig.get('centroid', 'N/A')}")
    else:
        print("   No ligands found (with exclude_common=True, min_atoms=5)")

    print()

    # 3. Summarize ligands
    print("3. Ligand Summary")
    print("-" * 40)

    summary = processor.summarize_ligands(
        TARGET_STRUCTURE,
        min_atoms=1,  # Include smaller molecules for summary
    )

    if summary:
        chains_info = summary.get("chains", {})
        total_ligands = sum(len(ligs) for ligs in chains_info.values())
        print(f"   Total ligands across all chains: {total_ligands}")

        for chain_id, ligands_dict in chains_info.items():
            print(f"\n   Chain {chain_id}:")
            for res_key, info in ligands_dict.items():
                comp_id = info.get("comp_id", "UNK")
                atom_count = info.get("atom_count", 0)
                print(f"     - {comp_id} ({res_key}): {atom_count} atoms")
    else:
        print("   No ligand summary available")

    print()

    # 4. Ligand-protein interactions
    print("4. Ligand-Protein Interactions")
    print("-" * 40)

    if ligands:
        # Analyze interactions for the first (largest) ligand
        target_ligand = ligands[0]
        ligand_id = target_ligand["res_name"]
        chain_id = target_ligand["chain_id"]

        print(f"   Analyzing interactions for {ligand_id} (chain {chain_id})...")

        interactions = processor.get_ligand_interactions(
            TARGET_STRUCTURE,
            ligand_id=ligand_id,
            chain_id=chain_id,
            cutoff=4.5,  # Distance cutoff in Angstroms
            include_water_bridges=False,
        )

        if interactions:
            binding_residues = interactions.get("binding_residues", [])
            summary_info = interactions.get("summary", {})

            print(f"\n   Interaction summary:")
            print(f"     Total binding residues: {summary_info.get('total_residues', len(binding_residues))}")
            print(f"     Unique chains involved: {summary_info.get('chains', [])}")

            if binding_residues:
                print(f"\n   Binding residues (within {4.5}A):")
                for i, res in enumerate(binding_residues[:10]):  # Show first 10
                    res_name = res.get("res_name", "UNK")
                    res_id = res.get("res_id", "?")
                    res_chain = res.get("chain_id", "?")
                    min_dist = res.get("min_distance", 0)
                    grn = res.get("grn", "")
                    grn_str = f" [GRN: {grn}]" if grn else ""
                    print(f"     {i+1}. {res_name}{res_id} (chain {res_chain}) - {min_dist:.2f}A{grn_str}")

                if len(binding_residues) > 10:
                    print(f"     ... and {len(binding_residues) - 10} more")
        else:
            print("   No interactions found")
    else:
        print("   Skipping (no ligands to analyze)")

    print()

    # 5. Ion contact analysis
    print("5. Ion Contact Analysis")
    print("-" * 40)

    ion_contacts = processor.get_ion_contacts(
        TARGET_STRUCTURE,
        cutoff=3.5,
        ion_types=["NA", "K", "CA", "MG", "ZN", "FE", "CL"],
    )

    if ion_contacts:
        print(f"   Found {len(ion_contacts)} ion-protein contacts:\n")
        for i, contact in enumerate(ion_contacts[:10]):
            ion_name = contact.get("ion_name", "UNK")
            ion_chain = contact.get("ion_chain", "?")
            res_name = contact.get("res_name", "UNK")
            res_id = contact.get("res_id", "?")
            res_chain = contact.get("res_chain", "?")
            distance = contact.get("distance", 0)
            print(f"   {i+1}. {ion_name} (chain {ion_chain}) - {res_name}{res_id} (chain {res_chain}): {distance:.2f}A")
    else:
        print("   No ion contacts found in structure")

    print()

    # 6. Water contacts (brief)
    print("6. Water Contact Analysis")
    print("-" * 40)

    water_contacts = processor.get_water_contacts(
        TARGET_STRUCTURE,
        cutoff=3.5,
    )

    if water_contacts:
        print(f"   Found {len(water_contacts)} water-protein contacts")
        # Group by residue
        residues_with_water = set()
        for contact in water_contacts:
            res_key = (contact.get("res_chain"), contact.get("res_name"), contact.get("res_id"))
            residues_with_water.add(res_key)
        print(f"   Unique residues in contact with water: {len(residues_with_water)}")
    else:
        print("   No water contacts found")

    print()
    print("=== StructureProcessor Ligand Analysis Demo Complete ===")


if __name__ == "__main__":
    main()
