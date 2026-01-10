#!/usr/bin/env python3
"""MoleculeProcessor basic data management demonstration.

This script demonstrates core MoleculeProcessor capabilities:
- Registering molecules from SMILES
- Loading and saving molecule entities
- Calculating molecular properties
- Drug-likeness filtering
- Dataset management
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


# Sample molecules for testing
SAMPLE_MOLECULES = {
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "acetaminophen": "CC(=O)NC1=CC=C(C=C1)O",
    "glucose": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "atp": "NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O",
}


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.molecule import MoleculeProcessor

    # Initialize processor
    print("=== MoleculeProcessor Basic Demo ===\n")
    mol_proc = MoleculeProcessor()

    # 1. List existing entities
    print("1. Listing existing molecule entities...")
    entities = mol_proc.list_entities()
    print(f"   Found {len(entities)} entities\n")

    # 2. Register molecules from SMILES map
    print("2. Registering molecules from SMILES...")
    mol_proc.register_smiles_map(
        SAMPLE_MOLECULES,
        dataset_name="demo_molecules",
        metadata={"source": "test_molecule_processor.py"},
    )
    print(f"   Registered {len(SAMPLE_MOLECULES)} molecules in dataset 'demo_molecules'\n")

    # 3. List entities again
    print("3. Verifying registration...")
    entities = mol_proc.list_entities()
    print(f"   Now have {len(entities)} entities")
    for mol_name in list(SAMPLE_MOLECULES.keys())[:3]:
        exists = mol_proc.entity_exists(mol_name)
        print(f"   '{mol_name}' exists: {exists}")
    print()

    # 4. Load a specific molecule
    print("4. Loading molecule 'aspirin'...")
    aspirin = mol_proc.load_entity("aspirin")
    if aspirin:
        print(f"   Loaded: {type(aspirin)}")
        if isinstance(aspirin, dict):
            print(f"   Keys: {list(aspirin.keys())}")
            if "smiles" in aspirin:
                print(f"   SMILES: {aspirin['smiles']}")
    print()

    # 5. Calculate molecular properties
    print("5. Calculating molecular properties...")
    props = mol_proc.calculate_properties(list(SAMPLE_MOLECULES.keys()))
    if props is not None and not props.empty:
        print(f"   Calculated properties for {len(props)} molecules")
        print(f"   Properties: {list(props.columns)[:6]}...")
        print(f"   Sample (aspirin):")
        if "aspirin" in props.index:
            aspirin_props = props.loc["aspirin"]
            for col in ["molecular_weight", "logp", "hbd", "hba"][:4]:
                if col in aspirin_props:
                    print(f"      {col}: {aspirin_props[col]}")
    else:
        print("   Property calculation requires RDKit")
    print()

    # 6. Filter drug-like molecules
    print("6. Filtering drug-like molecules (Lipinski's Rule of Five)...")
    drug_like = mol_proc.filter_drug_like(list(SAMPLE_MOLECULES.keys()))
    if drug_like is not None:
        print(f"   Drug-like molecules: {drug_like}")
    else:
        print("   Drug-likeness filtering requires RDKit")
    print()

    # 7. List datasets
    print("7. Listing datasets...")
    datasets = mol_proc.list_datasets()
    print(f"   Found {len(datasets)} datasets")
    if "demo_molecules" in datasets:
        info = mol_proc.get_dataset_info("demo_molecules")
        print(f"   demo_molecules info: {info}\n")

    # 8. Load dataset
    print("8. Loading dataset 'demo_molecules'...")
    dataset = mol_proc.load_dataset("demo_molecules")
    print(f"   Loaded {len(dataset)} molecules\n")

    # 9. Save a new molecule entity
    print("9. Saving new molecule 'test_molecule'...")
    mol_proc.save_entity(
        "test_molecule",
        {"smiles": "CCCC", "name": "butane", "kind": "smiles_record"},
        metadata={"source": "manual"},
    )
    print("   Saved 'test_molecule'\n")

    # 10. Create a custom dataset
    print("10. Creating custom dataset 'small_drugs'...")
    small_drugs = ["aspirin", "caffeine", "acetaminophen"]
    mol_proc.create_dataset("small_drugs", small_drugs, {"description": "Small drug molecules"})
    print(f"   Created dataset with {len(small_drugs)} molecules\n")

    # 11. Export dataset
    print("11. Exporting dataset...")
    try:
        export_dir = data_root / "exports" / "molecule"
        export_dir.mkdir(parents=True, exist_ok=True)
        export_info = mol_proc.export_dataset(
            "demo_molecules",
            output_dir=export_dir,
            overwrite=True,
        )
        print(f"   Exported {len(export_info)} entities to: {export_dir}\n")
    except Exception as e:
        print(f"   Export not available: {e}\n")

    # 12. List ligands (if structure data available)
    print("12. Listing available ligands...")
    ligands = mol_proc.list_ligands()
    print(f"   Found {len(ligands)} ligands\n")

    print("=== MoleculeProcessor Demo Complete ===")


if __name__ == "__main__":
    main()
