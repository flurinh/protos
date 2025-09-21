#!/usr/bin/env python3
"""Sequence ingestion and export workflow demo."""

from __future__ import annotations

from pathlib import Path

import protos


DATASET_FASTA = ">SEQ_ALPHA\nMKTIIALSYIFCLVFADYKDDDDA\n>SEQ_BETA\nGSHSMRYFYTAMSRPGRGEPRFIAVGYVDDMRFYQRS\n"
SINGLE_FASTA = ">SINGLE_SEQ\nMPLNVSFTDLEK\n"
GPCR_ACCESSIONS = ["P30542", "P07550", "Q9Y5N6"]


def ensure_data_root() -> Path:
    root = Path(__file__).resolve().parent / "data"
    root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(root))
    return root


def create_fasta_files(base_dir: Path) -> tuple[Path, Path]:
    dataset_path = base_dir / "demo_sequences.fasta"
    single_path = base_dir / "single_sequence.fasta"
    dataset_path.write_text(DATASET_FASTA)
    single_path.write_text(SINGLE_FASTA)
    return dataset_path, single_path


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.sequence import SequenceProcessor
    from protos.io.ingest.sequence_loader import SequenceLoader

    processor = SequenceProcessor()
    loader = SequenceLoader(processor=processor)

    from protos.io.paths import get_protos_paths

    input_dir = Path(get_protos_paths().get_processor_path('input'))
    input_dir.mkdir(parents=True, exist_ok=True)

    dataset_fasta, single_fasta = create_fasta_files(input_dir)

    entity_name = loader.download_and_register(
        str(single_fasta),
        name="single_sequence",
        materialize_entities=True,
    )
    print(f"Registered single sequence entity: {entity_name}")

    dataset_name = loader.download_and_register(
        str(dataset_fasta),
        name="demo_sequences",
        materialize_entities=False,
    )
    print(f"Registered dataset: {dataset_name}")

    gpcr_dataset = None
    try:
        gpcr_dataset = loader.download_and_register(
            "uniprot:" + ",".join(GPCR_ACCESSIONS),
            name="gpcr_sequences",
            materialize_entities=False,
        )
        print(f"Downloaded GPCR dataset: {gpcr_dataset}")
    except Exception as exc:  # noqa: BLE001
        print(f"GPCR dataset fetch failed: {exc}")

    sequences = processor.load_dataset(dataset_name)
    print(f"Loaded {len(sequences)} sequences from dataset")

    if gpcr_dataset:
        gpcr_sequences = processor.load_dataset(gpcr_dataset)
        print(f"Loaded {len(gpcr_sequences)} GPCR sequences")

    ids = list(sequences.keys())
    if len(ids) >= 2:
        score, alignment = processor.align_sequences(
            sequences[ids[0]], sequences[ids[1]], ids[0], ids[1]
        )
        print(f"Alignment score ({ids[0]} vs {ids[1]}): {score:.2f}")
        print("Alignment preview:\n" + "\n".join(alignment[:3]) + "\n...")

    output_dir = data_root / "sequence_output"
    output_dir.mkdir(exist_ok=True)

    processor.export_dataset(dataset_name, output_dir, overwrite=True)
    print(f"Exported full dataset to {output_dir}")

    if gpcr_dataset:
        processor.export_dataset(gpcr_dataset, output_dir, overwrite=True)
        print(f"Exported GPCR dataset to {output_dir}")

    subset_path = processor.export_dataset(
        dataset_name,
        output_dir,
        name_pattern="{name}_subset.fasta",
        sequence_ids=[ids[-1]],
        overwrite=True,
    )[dataset_name]
    print(f"Exported subset to {subset_path}")

    target_seq_id = ids[0]
    processor.save_entity(target_seq_id, sequences[target_seq_id])
    single_path = processor.export_entity(
        target_seq_id,
        output_dir / f"{target_seq_id}.fasta",
        overwrite=True,
    )
    print(f"Exported single sequence to {single_path}")

    # Mutant library and downstream analytics
    library, library_meta = processor.create_mutant_library(
        base_sequence_id=target_seq_id,
        mutation_map={5: ['S', 'T'], 10: ['K']},
        limit=3,
        return_metadata=True,
        include_wildtype=True,
    )
    print(f"Generated mutant library with {len(library)} variants")
    print(library_meta.head().to_string(index=False))

    conservation = processor.compute_conservation(sequences=library)
    print("Conservation (first 5 positions):")
    print(conservation.head()[['position', 'consensus', 'consensus_frequency']].to_string(index=False))

    linkage = processor.compute_linkage(sequences=library, top_k=3)
    if not linkage.empty:
        print("Top linkage pairs:")
        print(linkage[['pos_i', 'pos_j', 'mutual_information']].to_string(index=False))
    else:
        print("Linkage analysis produced no pairs (insufficient diversity).")

    # UniProt single sequence
    try:
        uniprot_entity = loader.download_and_register(
            "uniprot:P69905",
            materialize_entities=True,
        )
        print(f"Downloaded UniProt entity: {uniprot_entity}")
    except Exception as exc:  # noqa: BLE001
        print(f"UniProt single fetch failed: {exc}")

    if gpcr_dataset:
        try:
            saved_path = processor.save_sequences(
                gpcr_sequences,
                output_file="gpcr_sequences_export",
                dataset_name="gpcr_sequences_export",
                materialize_entities=False,
            )
            print(f"Saved GPCR sequences to {saved_path}")
        except Exception as exc:  # noqa: BLE001
            print(f"Saving GPCR dataset failed: {exc}")


if __name__ == "__main__":
    main()
