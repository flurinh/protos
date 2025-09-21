#!/usr/bin/env python3
"""Showcase LigandProcessor ingestion and interaction recording."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


SMILES_SAMPLE = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
}


def main() -> None:
    ensure_data_root()

    from protos.processing.ligand import LigandProcessor
    from protos.processing.structure import StructureProcessor

    ligand_proc = LigandProcessor()

    dataset, entities = ligand_proc.register_smiles_dataset(SMILES_SAMPLE)
    print(f"Registered {len(entities)} ligands in dataset '{dataset}':")
    for name in entities:
        record = ligand_proc.load_entity(name)
        props = record.get("metadata", {}).get("properties", {})
        mw = props.get("mw")
        print(f"  • {name[:40]}... MW={mw:.1f}" if mw else f"  • {name}")

    struct_proc = StructureProcessor()
    structure_entities = struct_proc.list_entities()
    if not structure_entities:
        print("No structures available for ligand extraction.")
        return

    structure_id = structure_entities[0]
    print(f"Attempting ligand extraction from structure '{structure_id}'...")

    try:
        ligand_dataset, lig_entities = ligand_proc.extract_ligands_from_structure(structure_id)
    except ValueError as exc:
        print(f"  ! Extraction skipped: {exc}")
        lig_entities = []
    else:
        if lig_entities:
            print(f"  • Stored {len(lig_entities)} structure-derived ligands in '{ligand_dataset}'")
        else:
            print("  • Structure contained no hetero ligands")

    try:
        interactions = ligand_proc.compute_interactions(structure_id, ligands=lig_entities or None)
    except ValueError as exc:
        print(f"  ! Interaction computation skipped: {exc}")
        interactions = pd.DataFrame()

    if not interactions.empty:
        print(f"Found {len(interactions)} contacts within 4.0 Å")
        ligand_proc.record_interactions(
            table_name="ligand_interactions",
            interactions=interactions,
            metadata={"structure_id": structure_id},
        )
        print("Recorded interactions under property table 'ligand_interactions'.")
    else:
        print("No ligand-protein contacts detected (or structure lacked ligands).")


if __name__ == "__main__":
    import pandas as pd  # Local import for optional data frame usage

    main()

