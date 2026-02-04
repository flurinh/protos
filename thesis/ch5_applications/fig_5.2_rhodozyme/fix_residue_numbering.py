#!/usr/bin/env python3
"""Fix residue numbering in grafted structures for BoltzGen compatibility.

BoltzGen expects contiguous residue numbering (1, 2, 3, ..., N).
The current structures have gaps in auth_seq_id. This script:
1. Renumbers all residues contiguously
2. Creates a mapping file so we know the correspondence
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

from protos.io.formats.cif_utils import read_cif_file, write_cif_file

OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"


def renumber_structure(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Renumber residues contiguously per chain."""
    df = df.copy()
    mapping = {}

    for chain in sorted(df["auth_chain_id"].unique()):
        chain_df = df[df["auth_chain_id"] == chain]
        old_residues = sorted(chain_df["auth_seq_id"].unique())

        # Create mapping: old_res -> new_res (1-indexed)
        chain_mapping = {old: i + 1 for i, old in enumerate(old_residues)}
        mapping[chain] = chain_mapping

        # Apply mapping
        mask = df["auth_chain_id"] == chain
        df.loc[mask, "auth_seq_id"] = df.loc[mask, "auth_seq_id"].map(chain_mapping)
        df.loc[mask, "label_seq_id"] = df.loc[mask, "auth_seq_id"]
        df.loc[mask, "gen_seq_id"] = df.loc[mask, "auth_seq_id"]

    return df, mapping


def main() -> int:
    print("=" * 70)
    print("FIXING RESIDUE NUMBERING FOR BOLTZGEN")
    print("=" * 70)

    # Load summary
    with open(OUTPUT_DIR / "grafts_with_substrate_summary.json") as f:
        summary = json.load(f)

    output_dir = OUTPUT_DIR / "grafted_structures_renumbered"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_mappings = {}
    updated_placements = []

    for placement in summary["placements"]:
        rank = placement["rank"]
        old_path = Path(placement["structure_path"])

        print(f"\n[{rank}] Processing: {old_path.name}")

        # Load and renumber
        df = read_cif_file(str(old_path))
        df_renumbered, mapping = renumber_structure(df)

        # Save renumbered structure
        new_name = f"trypsin_refined{rank:02d}_renumbered.cif"
        new_path = output_dir / new_name
        write_cif_file(str(new_path), df_renumbered, force_overwrite=True)

        # Map old residue numbers to new
        old_triad = placement["assigned_residues"]
        new_triad = [mapping["A"][r] for r in old_triad]

        # Find retinal (last residue with res_name RET)
        chain_a = df[df["auth_chain_id"] == "A"]
        ret_old = int(chain_a[chain_a["res_name"] == "RET"]["auth_seq_id"].iloc[0])
        ret_new = mapping["A"][ret_old]

        print(f"  Chain A: {len(mapping['A'])} residues (1 to {len(mapping['A'])})")
        print(f"  Chain B: {len(mapping['B'])} residues (1 to {len(mapping['B'])})")
        print(f"  Triad mapping: {old_triad} -> {new_triad}")
        print(f"  Retinal mapping: {ret_old} -> {ret_new}")

        all_mappings[rank] = {
            "chain_A": {str(k): v for k, v in mapping["A"].items()},
            "chain_B": {str(k): v for k, v in mapping["B"].items()},
            "triad_old": old_triad,
            "triad_new": new_triad,
            "retinal_old": ret_old,
            "retinal_new": ret_new,
        }

        updated_placements.append({
            **placement,
            "structure_path_renumbered": str(new_path),
            "triad_residues_new": new_triad,
            "retinal_residue_new": ret_new,
        })

    # Save mapping
    mapping_path = output_dir / "residue_mapping.json"
    with open(mapping_path, "w") as f:
        json.dump(all_mappings, f, indent=2)

    # Save updated summary
    updated_summary = {
        **summary,
        "placements": updated_placements,
    }
    summary_path = output_dir / "placements_renumbered.json"
    with open(summary_path, "w") as f:
        json.dump(updated_summary, f, indent=2)

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\nRenumbered structures: {output_dir}")
    print(f"Mapping file: {mapping_path}")

    # Print designable regions mapping for first placement
    print("\nDesignable regions mapping (for BoltzGen config):")
    m = all_mappings[1]["chain_A"]
    regions = [(58, 76), (134, 157), (222, 256), (306, 326)]
    for start, end in regions:
        new_start = m.get(str(start), "?")
        new_end = m.get(str(end), "?")
        print(f"  {start}..{end} -> {new_start}..{new_end}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
