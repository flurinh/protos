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
    grn_table = seq_proc.annotate_with_grn_reference(
        dataset_name="gpcr_sequences",
        reference_table=REFERENCE_TABLE,
        output_table_name="gpcr_grn_demo",
        allow_create=True,
    )

    print("Saved GRN table 'gpcr_grn_demo':")
    print(grn_table)

    grn_proc = GRNProcessor()
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
