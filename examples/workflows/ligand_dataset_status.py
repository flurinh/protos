#!/usr/bin/env python3
"""Ligand dataset status and cache inspection workflow."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict

import pandas as pd
import protos
from protos.processing.ligand import LigandProcessor

DATA_RELATIVE_ROOT = Path(__file__).resolve().parents[2] / "data"
LIGAND_DATASET = "P24941_chembl_ligands"
OUTPUT_DIR = DATA_RELATIVE_ROOT / "reports"
SUMMARY_PATH = OUTPUT_DIR / "ligand_dataset_status.json"


@dataclass
class DatasetSummary:
    dataset: str
    entity_count: int
    activity_records: int
    unique_ligands: int
    activity_types: int
    median_activity_nm: float


@dataclass
class CacheSummary:
    ccd_index_present: bool
    qm9_archive_present: bool
    enamine_manifest_present: bool


def ensure_data_root() -> Path:
    DATA_RELATIVE_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_RELATIVE_ROOT))
    return DATA_RELATIVE_ROOT


def load_dataset_metadata() -> Dict[str, object]:
    dataset_path = DATA_RELATIVE_ROOT / "ligand" / "datasets" / f"{LIGAND_DATASET}.json"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Ligand dataset metadata not found: {dataset_path}")
    return json.loads(dataset_path.read_text(encoding="utf-8"))


def summarize_dataset(metadata: Dict[str, object]) -> DatasetSummary:
    entities = metadata.get("entities", [])
    data_file = metadata.get("data_file")
    table_path = DATA_RELATIVE_ROOT / "ligand" / "tables" / Path(data_file).name
    activities = pd.read_csv(table_path)

    median_value = float(activities["value_nm"].median()) if "value_nm" in activities else 0.0

    return DatasetSummary(
        dataset=metadata.get("name", LIGAND_DATASET),
        entity_count=len(entities),
        activity_records=len(activities),
        unique_ligands=activities["chembl_id"].nunique(),
        activity_types=activities["activity_type"].nunique(),
        median_activity_nm=median_value,
    )


def summarize_caches(processor: LigandProcessor) -> CacheSummary:
    databases_dir = Path(processor.get_subdirectory_path("databases"))
    ccd_index = databases_dir / "ccd" / "ccd_index.json"
    qm9_archive = databases_dir / "qm9" / "qm9.tar.bz2"
    enamine_manifest = databases_dir / "enamine" / "datasets.json"

    return CacheSummary(
        ccd_index_present=ccd_index.exists(),
        qm9_archive_present=qm9_archive.exists(),
        enamine_manifest_present=enamine_manifest.exists(),
    )


def main() -> None:
    ensure_data_root()
    metadata = load_dataset_metadata()
    dataset_summary = summarize_dataset(metadata)
    processor = LigandProcessor()
    cache_summary = summarize_caches(processor)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, Dict[str, object]] = {
        "dataset": asdict(dataset_summary),
        "cache": asdict(cache_summary),
    }
    SUMMARY_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
