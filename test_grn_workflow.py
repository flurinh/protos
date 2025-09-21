#!/usr/bin/env python3
"""Demonstrate GRNProcessor recording and retrieval."""

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


REFERENCE_TABLE = "gpcrdb_ref"


def register_sequences() -> None:
    from protos.processing.sequence import SequenceProcessor
    from protos.io.paths import get_protos_paths
    from protos.io.ingest.sequence_loader import SequenceLoader

    processor = SequenceProcessor()

    try:
        processor.load_dataset("gpcr_sequences")
        return
    except Exception:
        pass

    sequences = {
        "3sn6_chain_A": "MKTIIALSYIFCLVFADYKDDDDAAAFVVVLG",
        "5d5a_chain_A": "MNTSVYIFCLVFADVTDKDNRTLLGFFVASLL",
        "6b73_chain_A": "MKSVLIFCLVFADYKDDDAAGGMVLLVFVVIL",
    }

    paths = get_protos_paths()
    input_dir = Path(paths.get_processor_path('input'))
    input_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = input_dir / "gpcr_sequences.fasta"
    with open(fasta_path, "w") as handle:
        for seq_id, seq in sequences.items():
            handle.write(f">{seq_id}\n{seq}\n")

    loader = SequenceLoader(processor=processor)
    loader.download_and_register(
        str(fasta_path),
        name="gpcr_sequences",
        materialize_entities=True,
    )


def record_grn_table() -> None:
    from protos.processing.sequence import SequenceProcessor
    from protos.processing.grn import GRNProcessor

    seq_proc = SequenceProcessor()

    dataset_sequences = seq_proc.load_dataset("gpcr_sequences")
    reference_id = "5d5a_chain_A"
    reference_sequence = dataset_sequences.get(reference_id)
    if reference_sequence is None:
        loaded_ref = seq_proc.load_entity(reference_id)
        if isinstance(loaded_ref, dict):
            reference_sequence = loaded_ref.get(reference_id)
        else:
            reference_sequence = loaded_ref
    if not reference_sequence:
        raise RuntimeError(f"Reference sequence '{reference_id}' could not be resolved")

    similarity_scores = {}
    reference_map = {reference_id: reference_sequence}
    for seq_id, sequence in dataset_sequences.items():
        _, score, _ = seq_proc.find_best_match(
            sequence,
            reference_map,
            use_mmseqs=True,
        )
        if score is None or score == float("-inf"):
            _, score, _ = seq_proc.find_best_match(
                sequence,
                reference_map,
                use_mmseqs=False,
            )
        similarity_scores[seq_id] = score / max(len(sequence), 1)

    threshold = 0.35
    classified = {
        seq_id: (score >= threshold)
        for seq_id, score in similarity_scores.items()
    }

    print("Similarity to 5d5a_chain_A (normalized score):")
    for seq_id, score in similarity_scores.items():
        status = "gpcr_like" if classified[seq_id] else "low_similarity"
        print(f"  {seq_id}: {score:.3f} ({status})")

    grn_table = seq_proc.annotate_with_grn_reference(
        dataset_name="gpcr_sequences",
        reference_table=REFERENCE_TABLE,
        protein_family="gpcr_a",
        output_table_name="gpcr_grn_demo",
        allow_create=True,
    )

    grn_proc = GRNProcessor()
    reference_table = grn_proc.load_reference_table(REFERENCE_TABLE)

    # Ensure column coverage matches the reference definition
    assert list(grn_table.columns) == list(reference_table.columns), (
        "GRN annotation table missing reference columns"
    )

    # Every sequence must have at least one assigned GRN position
    assert all((row != '-').any() for _, row in grn_table.iterrows()), (
        "Expected each sequence to have at least one annotated GRN"
    )

    dataset_info = grn_proc.get_dataset_info("gpcr_grn_demo")
    metadata = dataset_info.get("metadata", {}) if dataset_info else {}
    match_summary = metadata.get("match_summary", {})
    assert match_summary, "GRN metadata missing alignment summary"
    assert all(details.get("assigned_positions", 0) > 0 for details in match_summary.values()), (
        "Alignment summary reports zero assigned positions"
    )
    assert all(
        details.get("coverage_fraction", 0) > 0
        for details in match_summary.values()
    ), "Coverage metadata missing for GRN annotations"

    print("Saved GRN table 'gpcr_grn_demo':")
    print(grn_table)

    seq_annotations = grn_proc.get_annotations("3sn6_chain_A")
    print("\nResolved annotations for 3sn6_chain_A:")
    print(seq_annotations)

    related = seq_proc.resolve_related_entities("3sn6_chain_A", format_type="grn")
    print("\nRegistry relationships for 3sn6_chain_A:")
    print(related)


def main() -> None:
    ensure_data_root()
    register_sequences()
    record_grn_table()


if __name__ == "__main__":
    main()
