#!/usr/bin/env python3
"""PropertyProcessor basic data management demonstration.

This script demonstrates core PropertyProcessor capabilities:
- Creating and listing property tables
- Recording properties with entity scopes
- Loading and querying properties
- Filtering by property values
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


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.property import PropertyProcessor

    # Initialize processor
    print("=== PropertyProcessor Basic Demo ===\n")
    prop_proc = PropertyProcessor()

    # 1. List existing tables
    print("1. Listing existing property tables...")
    tables = prop_proc.list_tables()
    print(f"   Found {len(tables)} tables: {tables[:5]}{'...' if len(tables) > 5 else ''}\n")

    # 2. Create a property table with sample data
    print("2. Creating property table 'demo_protein_properties'...")
    # Clean up any existing table to avoid duplicate data from previous runs
    table_path = data_root / "property" / "tables" / "demo_protein_properties.csv"
    if table_path.exists():
        table_path.unlink()
    sample_properties = [
        {
            "entity_name": "protein_A",
            "scope": [{"format": "sequence", "name": "protein_A"}],
            "molecular_weight": 45000.0,
            "pi": 6.8,
            "solubility": "high",
        },
        {
            "entity_name": "protein_B",
            "scope": [{"format": "sequence", "name": "protein_B"}],
            "molecular_weight": 32000.0,
            "pi": 5.2,
            "solubility": "medium",
        },
        {
            "entity_name": "protein_C",
            "scope": [{"format": "sequence", "name": "protein_C"}],
            "molecular_weight": 78000.0,
            "pi": 8.1,
            "solubility": "low",
        },
        {
            "entity_name": "structure_1_chain_A",
            "scope": [
                {"format": "structure", "name": "structure_1"},
                {"format": "sequence", "name": "structure_1_chain_A"},
            ],
            "resolution": 2.1,
            "method": "X-ray",
            "r_factor": 0.19,
        },
    ]

    prop_proc.record_properties(
        "demo_protein_properties",
        sample_properties,
        metadata={
            "description": "Demo property table for testing",
            "source": "test_property_processor.py",
        },
        allow_create=True,
    )
    print("   Table created with 4 property records\n")

    # 3. List tables again
    print("3. Verifying table creation...")
    tables = prop_proc.list_tables()
    assert "demo_protein_properties" in tables
    print(f"   'demo_protein_properties' in tables: True\n")

    # 4. Load the property table
    print("4. Loading property table...")
    table_df = prop_proc.load_table("demo_protein_properties")
    print(f"   Loaded table with {len(table_df)} rows")
    print(f"   Columns: {list(table_df.columns)}\n")

    # 5. Get properties for a specific entity
    print("5. Getting properties for 'protein_A'...")
    props = prop_proc.get_properties("protein_A")
    print(f"   Properties:\n{props.to_string(index=False)}\n")

    # 6. Get properties for structure (multi-scope)
    print("6. Getting properties for 'structure_1'...")
    struct_props = prop_proc.get_properties("structure_1")
    print(f"   Properties:\n{struct_props.to_string(index=False)}\n")

    # 7. Filter by property value
    print("7. Filtering proteins by molecular_weight > 40000...")
    filtered = prop_proc.filter_by_property(
        "demo_protein_properties",
        property_name="molecular_weight",
        condition=lambda x: x > 40000,
    )
    print(f"   Found {len(filtered)} matching records")
    if not filtered.empty:
        print(f"   Entities: {filtered['entity_name'].tolist()}\n")

    # 8. Add a new property column
    print("8. Adding new property column 'stability_score'...")
    # add_property_column takes: table_name, property_name, values (dict entity->value)
    stability_values = {
        "protein_A": 0.85,
        "protein_B": 0.72,
        "protein_C": 0.91,
        "structure_1_chain_A": 0.78,
    }
    prop_proc.add_property_column(
        "demo_protein_properties",
        "stability_score",
        stability_values,
    )
    table_df = prop_proc.load_table("demo_protein_properties")
    print(f"   Column added. Columns: {list(table_df.columns)}\n")

    # 9. Create a dataset of entities
    print("9. Creating dataset 'high_mw_proteins'...")
    high_mw = filtered["entity_name"].tolist() if not filtered.empty else ["protein_A"]
    prop_proc.create_dataset("high_mw_proteins", high_mw, {"filter": "mw > 40000"})
    print(f"   Dataset created with {len(high_mw)} entities\n")

    # 10. List datasets
    print("10. Listing datasets...")
    datasets = prop_proc.list_datasets()
    print(f"   Found {len(datasets)} datasets")
    if "high_mw_proteins" in datasets:
        info = prop_proc.get_dataset_info("high_mw_proteins")
        print(f"   high_mw_proteins info: {info}\n")

    # 11. Load dataset
    print("11. Loading dataset 'high_mw_proteins'...")
    dataset_entities = prop_proc.load_dataset("high_mw_proteins")
    print(f"   Loaded {len(dataset_entities)} entities\n")

    # 12. Export table to CSV
    print("12. Exporting table to CSV...")
    export_dir = data_root / "exports" / "property"
    export_dir.mkdir(parents=True, exist_ok=True)
    export_info = prop_proc.export_dataset(
        "demo_protein_properties",
        output_dir=export_dir,
        overwrite=True,
    )
    print(f"   Exported {len(export_info)} entities to: {export_dir}\n")

    print("=== PropertyProcessor Demo Complete ===")


if __name__ == "__main__":
    main()
