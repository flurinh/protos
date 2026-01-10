#!/usr/bin/env python3
"""Structure alignment and GRN annotation workflow using bundled datasets."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import mean
from typing import Iterable, List

import protos
from protos.processing.structure import StructureProcessor

DATA_RELATIVE_ROOT = Path(__file__).resolve().parents[2] / "data"
STRUCTURE_DATASET = "rhodopsin_states"
REFERENCE_TABLE = "gpcrdb_ref"
PROTEIN_FAMILY = "gpcr_a"
OUTPUT_DIR = DATA_RELATIVE_ROOT / "reports"
SUMMARY_PATH = OUTPUT_DIR / "structure_alignment_annotation.json"


def ensure_data_root() -> Path:
    DATA_RELATIVE_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_RELATIVE_ROOT))
    return DATA_RELATIVE_ROOT


def ensure_structures(processor: StructureProcessor, dataset_name: str) -> List[str]:
    entity_ids = processor.dataset_manager.get_dataset_entities(dataset_name)
    missing = []
    for structure_id in entity_ids:
        if processor.load_entity(structure_id) is None:
            missing.append(structure_id)
    if missing:
        raise FileNotFoundError(
            "The following structures are missing from the data root: "
            + ", ".join(sorted(missing))
        )
    return list(entity_ids)


@dataclass
class AlignmentMetrics:
    reference: str
    aligned_structures: int
    rmsd_min: float
    rmsd_mean: float
    rmsd_max: float


@dataclass
class StructureAnnotationMetrics:
    reference_table: str
    protein_family: str
    chains_annotated: int
    total_chains: int
    mean_coverage: float
    min_coverage: float
    max_coverage: float


def align_structures(processor: StructureProcessor, structure_ids: List[str]) -> AlignmentMetrics:
    reference_id, *mobile_ids = structure_ids
    summary, _ = processor.align_and_record(
        structure_ids=mobile_ids,
        reference_id=reference_id,
        method="cealign",
        atom_selection="CA",
        apply_transform=True,
        save_aligned=True,
        aligned_dataset_name=f"{STRUCTURE_DATASET}_aligned",
        aligned_dataset_include_reference=True,
        summary_name=f"{STRUCTURE_DATASET}_alignment_summary",
    )

    rmsd_stats = summary.get("rmsd", {}).get("global", {})
    return AlignmentMetrics(
        reference=reference_id,
        aligned_structures=len(mobile_ids) + 1,
        rmsd_min=float(rmsd_stats.get("min", 0.0)),
        rmsd_mean=float(rmsd_stats.get("mean", 0.0)),
        rmsd_max=float(rmsd_stats.get("max", 0.0)),
    )


def annotate_structures(processor: StructureProcessor, structure_ids: Iterable[str]) -> StructureAnnotationMetrics:
    dataframe, summary = processor.annotate_structures_with_grn(
        structure_ids=structure_ids,
        reference_table=REFERENCE_TABLE,
        protein_family=PROTEIN_FAMILY,
        output_table=f"{STRUCTURE_DATASET}_grn",
        return_summary=True,
        allow_create=True,
    )

    coverages = [
        info.get("coverage", 0.0)
        for info in summary.get("per_sequence", {}).values()
    ]

    return StructureAnnotationMetrics(
        reference_table=REFERENCE_TABLE,
        protein_family=PROTEIN_FAMILY,
        chains_annotated=sum(1 for value in coverages if value > 0.0),
        total_chains=len(coverages),
        mean_coverage=mean(coverages) if coverages else 0.0,
        min_coverage=min(coverages) if coverages else 0.0,
        max_coverage=max(coverages) if coverages else 0.0,
    )


def main() -> None:
    ensure_data_root()
    processor = StructureProcessor()

    structure_ids = ensure_structures(processor, STRUCTURE_DATASET)
    alignment_metrics = align_structures(processor, structure_ids)
    annotation_metrics = annotate_structures(processor, structure_ids)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    report = {
        "alignment": asdict(alignment_metrics),
        "annotation": asdict(annotation_metrics),
    }
    SUMMARY_PATH.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
