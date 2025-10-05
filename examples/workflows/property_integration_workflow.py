#!/usr/bin/env python3
"""Property integration workflow for GPCR ligand binding analysis."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict

import pandas as pd
import protos
from protos.processing.property import PropertyProcessor

DATA_RELATIVE_ROOT = Path(__file__).resolve().parents[2] / "data"
PROPERTY_DATASET = "gpcr_ligand_binding_analysis"
OUTPUT_DIR = DATA_RELATIVE_ROOT / "reports"
SUMMARY_PATH = OUTPUT_DIR / "property_integration_summary.json"


@dataclass
class PropertySummary:
    dataset: str
    entity_count: int
    property_columns: int
    unique_receptors: int
    unique_ligand_types: int
    unique_ligands: int


def ensure_data_root() -> Path:
    DATA_RELATIVE_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_RELATIVE_ROOT))
    return DATA_RELATIVE_ROOT


def load_dataset_metadata() -> Dict[str, object]:
    dataset_path = DATA_RELATIVE_ROOT / "property" / "datasets" / f"{PROPERTY_DATASET}.json"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Property dataset metadata not found: {dataset_path}")
    return json.loads(dataset_path.read_text(encoding="utf-8"))


def summarize_properties(processor: PropertyProcessor, metadata: Dict[str, object]) -> PropertySummary:
    entities = metadata.get("entities", [])
    table_file = metadata.get("table_file")
    table_path = DATA_RELATIVE_ROOT / "property" / "tables" / Path(table_file).name
    table = pd.read_csv(table_path)

    return PropertySummary(
        dataset=metadata.get("name", PROPERTY_DATASET),
        entity_count=len(entities),
        property_columns=len(table.columns),
        unique_receptors=table["receptor"].nunique() if "receptor" in table else 0,
        unique_ligand_types=table["ligand_type"].nunique() if "ligand_type" in table else 0,
        unique_ligands=table["ligand_name"].nunique() if "ligand_name" in table else 0,
    )


def main() -> None:
    ensure_data_root()
    processor = PropertyProcessor()
    metadata = load_dataset_metadata()
    summary = summarize_properties(processor, metadata)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, Dict[str, object]] = {
        "properties": asdict(summary)
    }
    SUMMARY_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
