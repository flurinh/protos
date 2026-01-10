#!/usr/bin/env python3
"""Sequence alignment workflow demo."""

from __future__ import annotations

from pathlib import Path

import protos

DATASET_FASTA = ">SEQ1\nMKTIIALSYIFCLVFADYKDDDDA\n>SEQ2\nGSHSMRYFYTAMSRPGRGEPRFIAVGYVDDMRFYQRS\n"


def ensure_data_root() -> Path:
    root = Path(__file__).resolve().parents[3] / "data"  # -> repo root/data
    root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(root))
    return root


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.sequence import SequenceProcessor
    from protos.io.ingest.sequence_loader import SequenceLoader
    from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine

    processor = SequenceProcessor()
    loader = SequenceLoader(processor=processor)

    # Local dataset
    from protos.io.paths import get_protos_paths

    input_dir = Path(get_protos_paths().get_processor_path('input'))
    input_dir.mkdir(parents=True, exist_ok=True)
    dataset_path = input_dir / "gpcr_demo.fasta"
    dataset_path.write_text(DATASET_FASTA)

    dataset_name = loader.download_and_register(
        str(dataset_path),
        name="gpcr_demo",
        materialize_entities=False,
    )
    print(f"Registered dataset: {dataset_name}")

    sequences = processor.load_dataset(dataset_name)
    ids = list(sequences.keys())

    engine = SequenceAlignmentEngine()

    if len(ids) >= 2:
        result = engine.align_pairwise(ids[0], sequences[ids[0]], ids[1], sequences[ids[1]])
        print(f"Biopython alignment score: {result.score:.2f}")
        print("Alignment preview:\n" + "\n".join(result.alignment[:3]) + "\n...")

    try:
        mmseqs_output = engine.align_pairwise_mmseqs(sequences)
        print("MMseqs output lines:", len(mmseqs_output))
    except Exception as exc:  # noqa: BLE001
        print(f"MMseqs alignment unavailable: {exc}")

    # Export aligned sequences (biopython)
    export_info = processor.export_dataset(
        dataset_name,
        export_name="gpcr_demo_alignment",
        overwrite=True,
    )
    print(
        "Exported alignment dataset to",
        export_info.get("artifact_path") or export_info.get("file_path"),
    )


if __name__ == "__main__":
    main()

