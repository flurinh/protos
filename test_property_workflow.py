#!/usr/bin/env python3
"""PropertyProcessor workflow demo for sequences and structures."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos

TARGET_REFERENCE = "5d5a_chain_A"
GPCR_IDS = ["3sn6", "5d5a", "6b73"]
LOCAL_GPCR_SEQUENCES = {
    "3sn6_chain_A": "MKTIIALSYIFCLVFADYKDDDDAAAFVVVLG",
    "5d5a_chain_A": "MNTSVYIFCLVFADVTDKDNRTLLGFFVASLL",
    "6b73_chain_A": "MKSVLIFCLVFADYKDDDAAGGMVLLVFVVIL",
}


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def ensure_gpcr_structures() -> None:
    from protos.io.ingest.structure_loader import StructureLoader

    loader = StructureLoader()
    try:
        loader.download_batch(GPCR_IDS, dataset_name="gpcr_structures")
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: could not refresh GPCR structures: {exc}")


def register_chain_sequences() -> Dict[str, str]:
    from protos.processing.sequence import SequenceProcessor
    from protos.io.ingest.sequence_loader import SequenceLoader
    from protos.io.paths import get_protos_paths

    sequence_proc = SequenceProcessor()

    try:
        return sequence_proc.load_dataset("gpcr_sequences")
    except Exception:
        pass

    paths = get_protos_paths()
    input_dir = Path(paths.get_processor_path('input'))
    input_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = input_dir / "gpcr_sequences.fasta"
    with open(fasta_path, "w") as handle:
        for seq_id, seq in LOCAL_GPCR_SEQUENCES.items():
            handle.write(f">{seq_id}\n{seq}\n")

    loader = SequenceLoader(processor=sequence_proc)
    dataset_name = loader.download_and_register(
        str(fasta_path),
        name="gpcr_sequences",
        materialize_entities=False,
    )

    sequences = sequence_proc.load_dataset(dataset_name)
    return sequences


def align_against_reference(sequences: Dict[str, str]) -> Dict[str, Tuple[str, float]]:
    from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine

    engine = SequenceAlignmentEngine()

    if TARGET_REFERENCE not in sequences:
        raise RuntimeError(f"Reference sequence {TARGET_REFERENCE} not found")

    reference_sequences = {TARGET_REFERENCE: sequences[TARGET_REFERENCE]}
    results: Dict[str, Tuple[str, float]] = {}

    for seq_id, seq in sequences.items():
        best_ref = None
        best_score = -1.0
        for ref_id, ref_seq in reference_sequences.items():
            alignment = engine.align_pairwise(seq_id, seq, ref_id, ref_seq)
            score = alignment.score / max(len(seq), 1)
            if score > best_score:
                best_score = score
                best_ref = ref_id
        if best_ref is None:
            continue
        results[seq_id] = (best_ref, best_score)

    return results


def store_sequence_properties(classifications: Dict[str, Tuple[str, float]]) -> None:
    from protos.processing.property import PropertyProcessor

    prop_proc = PropertyProcessor()
    rows = []
    for seq_id, (ref_id, score) in classifications.items():
        rows.append(
            {
                "scope": [{"format": "sequence", "name": seq_id}],
                "reference": ref_id,
                "score": score,
                "entity_name": seq_id,
            }
        )
    prop_proc.record_properties(
        "gpcr_sequence_alignment",
        rows,
        allow_create=True,
    )

    print("Stored sequence properties in 'gpcr_sequence_alignment'.")

    prop_proc = PropertyProcessor()
    sample = next(iter(classifications))
    seq_props = prop_proc.get_properties(sample)
    print(f"Properties for {sample}:\n{seq_props.to_string(index=False)}\n")


def store_structure_properties(classifications: Dict[str, Tuple[str, float]]) -> None:
    from protos.processing.property import PropertyProcessor

    prop_proc = PropertyProcessor()
    rows = []
    threshold = 0.35
    for seq_id, (ref_id, score) in classifications.items():
        structure_id = seq_id.split("_chain_")[0]
        classification = (
            "reference" if seq_id == TARGET_REFERENCE else
            ("gpcr_like" if score >= threshold else "low_similarity")
        )
        rows.append(
            {
                "scope": [
                    {"format": "structure", "name": structure_id},
                    {"format": "sequence", "name": seq_id},
                ],
                "classification": classification,
                "score": score,
                "entity_name": seq_id,
            }
        )
    prop_proc.record_properties(
        "gpcr_structure_chain_scores",
        rows,
        allow_create=True,
    )

    print("Stored structure-chain properties in 'gpcr_structure_chain_scores'.")

    prop_proc = PropertyProcessor()
    sample_structure = TARGET_REFERENCE.split("_chain_")[0]
    struct_props = prop_proc.get_properties(sample_structure)
    print(f"Properties for structure {sample_structure}:\n{struct_props.to_string(index=False)}\n")


def main() -> None:
    ensure_data_root()
    ensure_gpcr_structures()
    sequences = register_chain_sequences()
    classifications = align_against_reference(sequences)
    store_sequence_properties(classifications)
    store_structure_properties(classifications)


if __name__ == "__main__":
    main()
