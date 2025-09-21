#!/usr/bin/env python3
"""Cross-processor demo: classify GPCR chains via sequence alignment and annotate structures."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Tuple

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


def ensure_gpcr_structures() -> None:
    from protos.io.ingest.structure_loader import StructureLoader

    loader = StructureLoader()
    gpcr_ids = ["3sn6", "5d5a", "6b73"]
    try:
        loader.download_batch(gpcr_ids, dataset_name="gpcr_structures")
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: could not refresh GPCR structures: {exc}")


TARGET_REFERENCE = "5d5a_chain_A"


def register_chain_sequences(limit: int | None = None) -> Dict[str, str]:
    from protos.processing.structure import StructureProcessor
    from protos.processing.sequence import SequenceProcessor

    structure_processor = StructureProcessor()
    sequence_processor = SequenceProcessor()

    try:
        structure_ids = structure_processor.get_dataset_entities("gpcr_structures")
    except FileNotFoundError:
        structure_ids = ["3sn6", "5d5a", "6b73"]

    # Ensure canonical PKLs exist by loading the dataset once
    structure_processor.load_dataset("gpcr_structures", return_format="dict")

    structure_processor.register_chain_sequences(
        structure_ids,
        dataset_prefix="gpcr_chain_dataset",
        create_dataset=True,
        overwrite=False,
    )

    related = structure_processor.list_dataset_related_sequences("gpcr_structures")
    chain_ids = []
    for relations in related.values():
        for rel in relations:
            chain_ids.append(rel["name"])

    unique_chain_ids = sorted(set(chain_ids))
    if limit is not None:
        unique_chain_ids = unique_chain_ids[:limit]
    if not unique_chain_ids:
        raise RuntimeError("No chain sequences available.")

    seq_map: Dict[str, str] = {}
    for chain_id in unique_chain_ids:
        data = sequence_processor.load_entity(chain_id)
        if isinstance(data, str):
            seq_map[chain_id] = data
        elif isinstance(data, dict):
            seq_map.update({f"{chain_id}_{k}": v for k, v in data.items()})

    if not seq_map:
        raise RuntimeError("Failed to retrieve chain sequences.")

    return seq_map


def align_and_classify(
    reference_sequences: Dict[str, str],
    all_sequences: Dict[str, str],
) -> Dict[str, Tuple[str, float]]:
    from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine

    engine = SequenceAlignmentEngine()

    results: Dict[str, Tuple[str, float]] = {}

    for seq_id, seq in all_sequences.items():
        best_ref = None
        best_score = -1.0
        for ref_id, ref_seq in reference_sequences.items():
            result = engine.align_pairwise(seq_id, seq, ref_id, ref_seq)
            score = result.score / max(len(seq), 1)
            if score > best_score:
                best_score = score
                best_ref = ref_id
        if best_ref is None:
            continue
        results[seq_id] = (best_ref, best_score)

    return results


def annotate_structures(
    classifications: Dict[str, Tuple[str, float]],
    *,
    threshold: float = 0.3,
    reference_key: str,
) -> None:
    from protos.processing.structure import StructureProcessor

    structure_processor = StructureProcessor()

    per_structure: Dict[str, Dict[str, Dict[str, str]]] = {}
    for seq_id, (ref_id, score) in classifications.items():
        parts = seq_id.split("_chain_")
        if len(parts) != 2:
            continue
        structure_id, chain_id = parts
        is_reference_chain = seq_id == reference_key
        if is_reference_chain:
            classification = "reference"
        else:
            classification = "gpcr_like" if score >= threshold else "low_similarity"
        per_structure.setdefault(structure_id, {}).setdefault(chain_id, {})[
            "classification"
        ] = classification

    for structure_id, chain_map in per_structure.items():
        annotations = {"chains": chain_map}
        structure_processor.annotate_structure(structure_id, annotations)
        df = structure_processor.frames[structure_id].reset_index()
        summary = (
            df.groupby("auth_chain_id")["classification"].first().to_dict()
        )
        print(f"Structure {structure_id} classifications: {summary}")


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    ensure_gpcr_structures()
    sequences = register_chain_sequences(limit=6)

    reference_sequences = {k: v for k, v in sequences.items() if k == TARGET_REFERENCE}
    if not reference_sequences:
        raise RuntimeError(f"Reference chain {TARGET_REFERENCE} not available")

    classifications = align_and_classify(reference_sequences, sequences)
    annotate_structures(
        classifications,
        threshold=0.35,
        reference_key=TARGET_REFERENCE,
    )


if __name__ == "__main__":
    main()
