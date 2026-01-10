#!/usr/bin/env python3
"""StructureProcessor quality analysis demonstration.

This script demonstrates structure quality assessment methods:
- Missing atoms detection
- Bond length validation
- Clash detection
- Chain break detection
- B-factor statistics
- Chirality validation
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


TARGET_STRUCTURE = "3sn6"  # Beta2-adrenergic receptor structure


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.structure import StructureProcessor
    from protos.io.ingest.structure_loader import StructureLoader

    print("=== StructureProcessor Quality Analysis Demo ===\n")

    processor = StructureProcessor()
    loader = StructureLoader(processor=processor)

    # Ensure structure is available
    print(f"1. Loading structure '{TARGET_STRUCTURE}'...")
    structure_df = processor.load_entity(TARGET_STRUCTURE)

    if structure_df is None:
        print(f"   Downloading {TARGET_STRUCTURE}...")
        success, failed = loader.download_batch(
            [TARGET_STRUCTURE],
            dataset_name="quality_analysis_structures",
            create_dataset=True,
        )
        if failed:
            print(f"   Failed to download: {failed}")
            return
        structure_df = processor.load_entity(TARGET_STRUCTURE)

    if structure_df is None:
        print("   Could not load structure")
        return

    print(f"   Loaded structure with {len(structure_df)} atoms")

    # Get chain information
    chains = structure_df["auth_chain_id"].unique().tolist()
    print(f"   Chains: {chains}\n")

    # Import analysis functions
    from protos.analysis.structure.quality import (
        check_missing_atoms,
        validate_bond_lengths,
        check_chain_breaks,
        calculate_b_factor_statistics,
        validate_chirality,
    )
    # Note: check_clashes and validate_structure_integrity are available but slow for large structures

    # Work with a single chain for efficiency
    protein_chains = structure_df[structure_df["group"].str.upper() == "ATOM"]["auth_chain_id"].unique()
    if len(protein_chains) > 0:
        main_chain = protein_chains[0]
        chain_df = structure_df[structure_df["auth_chain_id"] == main_chain]
        print(f"   Analyzing chain {main_chain} ({len(chain_df)} atoms) for detailed quality checks\n")
    else:
        chain_df = structure_df
        main_chain = "all"
        print(f"   Analyzing all atoms for quality checks\n")

    # 2. Missing atoms analysis
    print("2. Missing Atoms Analysis")
    print("-" * 40)

    missing = check_missing_atoms(chain_df)
    if missing is not None and not missing.empty:
        print(f"   Found {len(missing)} residues with missing atoms:")
        # Group by residue type
        if "residue_name" in missing.columns:
            by_type = missing.groupby("residue_name").size()
            for res_type, count in list(by_type.items())[:5]:
                print(f"     {res_type}: {count} residues")
            if len(by_type) > 5:
                print(f"     ... and {len(by_type) - 5} more residue types")
    else:
        print("   No missing atoms detected (or check not applicable)")

    print()

    # 3. Bond length validation
    print("3. Bond Length Validation")
    print("-" * 40)

    bond_issues = validate_bond_lengths(chain_df)
    if bond_issues is not None and not bond_issues.empty:
        print(f"   Found {len(bond_issues)} unusual bond lengths:")
        cols_to_show = [c for c in ["chain_id", "residue_id", "bond_type", "length", "expected", "deviation_sigma"]
                        if c in bond_issues.columns]
        print(bond_issues.head(5)[cols_to_show].to_string(index=False))
        if len(bond_issues) > 5:
            print(f"   ... and {len(bond_issues) - 5} more")
    else:
        print("   All bond lengths within expected ranges")

    print()

    # 4. Chain break detection
    print("4. Chain Break Detection")
    print("-" * 40)

    chain_breaks = check_chain_breaks(chain_df, max_ca_distance=4.2)
    if chain_breaks is not None and not chain_breaks.empty:
        print(f"   Found {len(chain_breaks)} chain breaks:")
        for _, row in chain_breaks.head(10).iterrows():
            chain = row.get("chain_id", "?")
            pos1 = row.get("residue1_id", "?")
            pos2 = row.get("residue2_id", "?")
            dist = row.get("ca_distance", 0)
            gap = row.get("gap_size", 0)
            print(f"     Chain {chain}: {pos1} -> {pos2} (CA dist: {dist:.2f}A, gap: {gap} residues)")
        if len(chain_breaks) > 10:
            print(f"   ... and {len(chain_breaks) - 10} more")
    else:
        print("   No chain breaks detected (CA distance < 4.2A between consecutive residues)")

    print()

    # 5. B-factor statistics
    print("5. B-factor Statistics")
    print("-" * 40)

    bfactor_stats = calculate_b_factor_statistics(chain_df)
    if bfactor_stats and "error" not in bfactor_stats:
        print(f"   Mean B-factor: {bfactor_stats.get('mean', 0):.2f}")
        print(f"   Std B-factor:  {bfactor_stats.get('std', 0):.2f}")
        print(f"   Min B-factor:  {bfactor_stats.get('min', 0):.2f}")
        print(f"   Max B-factor:  {bfactor_stats.get('max', 0):.2f}")

        # High B-factor residues (flexible regions)
        high_b_residues = bfactor_stats.get("high_b_residues", [])
        if high_b_residues:
            print(f"\n   High B-factor residues (top flexible regions):")
            for res in high_b_residues[:5]:
                print(f"     {res['residue_name']}{res['residue_id']} (chain {res['chain']}): B={res['mean_b_factor']:.1f}")
    else:
        print(f"   B-factor statistics: {bfactor_stats.get('error', 'not available')}")

    print()

    # 6. Chirality validation
    print("6. Chirality Validation")
    print("-" * 40)

    chirality_issues = validate_chirality(chain_df)
    if chirality_issues is not None and not chirality_issues.empty:
        print(f"   Found {len(chirality_issues)} chirality issues (D-amino acids):")
        cols_to_show = [c for c in ["chain_id", "residue_id", "residue_name", "chirality"]
                        if c in chirality_issues.columns]
        print(chirality_issues.head(5)[cols_to_show].to_string(index=False))
    else:
        print("   All chiral centers have correct stereochemistry (L-amino acids)")

    print()

    # 7. Summary by chain
    print("7. Per-Chain Summary")
    print("-" * 40)

    for chain in chains[:3]:  # Limit to first 3 chains
        this_chain_df = structure_df[structure_df["auth_chain_id"] == chain]
        n_atoms = len(this_chain_df)
        n_residues = this_chain_df["auth_seq_id"].nunique()

        # Check for protein vs ligand
        is_protein = (this_chain_df["group"].str.upper() == "ATOM").any()
        chain_type = "protein" if is_protein else "ligand/other"

        print(f"\n   Chain {chain} ({chain_type}):")
        print(f"     Atoms: {n_atoms}")
        print(f"     Residues: {n_residues}")

        if is_protein:
            chain_breaks_this = check_chain_breaks(this_chain_df)
            n_breaks = len(chain_breaks_this) if chain_breaks_this is not None else 0
            print(f"     Chain breaks: {n_breaks}")

    print()
    print("=== StructureProcessor Quality Analysis Demo Complete ===")


if __name__ == "__main__":
    main()
