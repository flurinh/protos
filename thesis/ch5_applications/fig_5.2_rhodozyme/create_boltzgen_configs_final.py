#!/usr/bin/env python3
"""Create final BoltzGen configs with correct residue numbering."""
from __future__ import annotations

import json
import sys
import yaml
from pathlib import Path

THESIS_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"


def main() -> int:
    print("=" * 70)
    print("CREATING FINAL BOLTZGEN CONFIGS")
    print("=" * 70)

    # Load renumbered placements
    with open(OUTPUT_DIR / "grafted_structures_renumbered" / "placements_renumbered.json") as f:
        summary = json.load(f)

    config_dir = OUTPUT_DIR / "boltzgen_configs_final"
    config_dir.mkdir(parents=True, exist_ok=True)

    # Designable regions (same in renumbered structure)
    designable_ranges = [(58, 76), (134, 157), (222, 256), (306, 326)]

    for placement in summary["placements"]:
        rank = placement["rank"]
        # Use absolute path that works on host (protos wrapper handles Docker mapping)
        structure_path = placement["structure_path_renumbered"]
        triad = placement["triad_residues_new"]
        retinal = placement["retinal_residue_new"]

        print(f"\n[{rank}] Creating config")
        print(f"  Structure: {Path(structure_path).name}")
        print(f"  Triad: {triad}")
        print(f"  Retinal: {retinal}")

        # Build config
        config = {
            "entities": [
                {
                    "file": {
                        "path": structure_path,
                        # Include both chains
                        "include": [
                            {"chain": {"id": "A"}},
                            {"chain": {"id": "B"}},
                        ],
                        # Design only the binding domain loops
                        "design": [
                            {"chain": {"id": "A", "res_index": f"{start}..{end}"}}
                            for start, end in designable_ranges
                        ],
                        # Fix: catalytic triad, retinal, and entire substrate chain
                        "not_design": [
                            {"chain": {"id": "A", "res_index": str(r)}}
                            for r in triad + [retinal]
                        ] + [
                            {"chain": {"id": "B"}}
                        ],
                        # Structure groups for coordinate constraints
                        "structure_groups": [
                            # Default: everything fixed
                            {"group": {"id": "all", "visibility": 1}},
                            # Designable regions: regenerate structure
                        ] + [
                            {"group": {"id": "A", "res_index": f"{start}..{end}", "visibility": 0}}
                            for start, end in designable_ranges
                        ],
                        # Binding site specification
                        "binding_types": [
                            {"chain": {"id": "A", "binding": str(r)}}
                            for r in triad
                        ],
                    }
                }
            ],
            # Note: num_designs is handled by protos wrapper, not in yaml
        }

        # Save config
        config_path = config_dir / f"trypsin_refined{rank:02d}_final.yaml"
        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

        print(f"  Config: {config_path.name}")

    # Also save a metadata file for the protos wrapper
    metadata = {
        "num_designs": 8,
        "configs": [
            str(config_dir / f"trypsin_refined{p['rank']:02d}_final.yaml")
            for p in summary["placements"]
        ]
    }
    with open(config_dir / "run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\nConfigs saved to: {config_dir}")
    print(f"Total configs: {len(summary['placements'])}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
