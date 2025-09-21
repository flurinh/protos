#!/usr/bin/env python3
"""Real-data sequence workflow using GPCR chains."""

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


def ensure_gpcr_structures() -> None:
    from protos.io.ingest.structure_loader import StructureLoader

    loader = StructureLoader()
    gpcr_ids = ["3sn6", "5d5a", "6b73"]
    try:
        loader.download_batch(gpcr_ids, dataset_name="gpcr_structures")
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: could not refresh GPCR structures: {exc}")


def collect_gpcr_chain_sequences(limit: int = 3) -> dict[str, str]:
    from protos.processing.structure import StructureProcessor
    from protos.processing.sequence import SequenceProcessor

    structure_processor = StructureProcessor()
    sequence_processor = SequenceProcessor()

    try:
        structure_ids = structure_processor.get_dataset_entities("gpcr_structures")
    except FileNotFoundError:
        structure_ids = ["3sn6", "5d5a", "6b73"]

    structure_processor.register_chain_sequences(
        structure_ids,
        dataset_prefix="gpcr_chain_dataset",
        create_dataset=True,
        overwrite=False,
    )

    related = structure_processor.list_dataset_related_sequences("gpcr_structures")
    chain_ids: list[str] = []
    for relations in related.values():
        for rel in relations:
            chain_ids.append(rel["name"])

    unique_chain_ids = sorted(set(chain_ids))[:limit]
    if not unique_chain_ids:
        raise RuntimeError("No chain sequences registered; ensure structure processing succeeded.")

    seq_map: dict[str, str] = {}
    for chain_id in unique_chain_ids:
        data = sequence_processor.load_entity(chain_id)
        if isinstance(data, str):
            seq_map[chain_id] = data
        elif isinstance(data, dict):
            seq_map.update({f"{chain_id}_{k}": v for k, v in data.items()})

    if not seq_map:
        raise RuntimeError("Failed to collect chain sequences.")

    return seq_map


def register_sequences_via_loader(sequences: dict[str, str]) -> str:
    from protos.io.ingest.sequence_loader import SequenceLoader
    from protos.io.formats.fasta_utils import write_fasta
    from protos.io.paths import get_protos_paths

    paths = get_protos_paths()
    input_dir = Path(paths.get_processor_path("input"))
    input_dir.mkdir(parents=True, exist_ok=True)

    input_file = input_dir / "gpcr_chain_input.fasta"
    write_fasta(sequences, str(input_file))

    loader = SequenceLoader()
    dataset_name = loader.download_and_register(
        str(input_file),
        name="gpcr_chains_real",
        materialize_entities=False,
    )
    return dataset_name


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    ensure_gpcr_structures()

    sequences = collect_gpcr_chain_sequences(limit=4)
    dataset_name = register_sequences_via_loader(sequences)

    from protos.processing.sequence import SequenceProcessor
    from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine

    processor = SequenceProcessor()
    alignment_engine = SequenceAlignmentEngine()

    gpcr_sequences = processor.load_dataset(dataset_name)
    print(f"Registered GPCR dataset '{dataset_name}' with {len(gpcr_sequences)} sequences")

    seq_ids = list(gpcr_sequences.keys())
    if len(seq_ids) >= 2:
        ref_id, mob_id = seq_ids[0], seq_ids[1]
        pair_result = alignment_engine.align_pairwise(
            ref_id,
            gpcr_sequences[ref_id],
            mob_id,
            gpcr_sequences[mob_id],
        )
        print(
            f"Pairwise alignment {pair_result.seq1_id} vs {pair_result.seq2_id}: "
            f"score={pair_result.score:.2f}"
        )
        print("Alignment preview:\n" + "\n".join(pair_result.alignment[:3]) + "\n...")

    try:
        mmseqs_lines = alignment_engine.align_pairwise_mmseqs(gpcr_sequences)
        print(f"MMseqs alignment lines: {len(mmseqs_lines)}")
    except Exception as exc:  # noqa: BLE001
        print(f"MMseqs alignment unavailable: {exc}")

    length_groups: dict[int, list[str]] = {}
    for seq_id, seq in gpcr_sequences.items():
        length_groups.setdefault(len(seq), []).append(seq_id)

    primary_length, primary_ids = max(length_groups.items(), key=lambda item: len(item[1]))

    if len(primary_ids) >= 2:
        aligned_subset = {sid: gpcr_sequences[sid] for sid in primary_ids}
        conservation = processor.compute_conservation(sequences=aligned_subset)
        print(f"Conservation summary (first 5 positions) for length {primary_length}:")
        print(
            conservation.head()[
                ["position", "consensus", "consensus_frequency", "normalized_entropy"]
            ].to_string(index=False)
        )

        linkage = processor.compute_linkage(sequences=aligned_subset, top_k=5)
        if not linkage.empty:
            print("Top linkage pairs:")
            print(linkage.to_string(index=False))
        else:
            print("Linkage analysis produced no significant pairs (aligned subset).")
    else:
        print("Not enough sequences with identical length for conservation/linkage analysis.")

    base_sequence_id = seq_ids[0]
    base_sequence = gpcr_sequences[base_sequence_id]
    seq_len = len(base_sequence)

    pos_candidates = [5, 100, 150]
    normalized_positions: dict[int, list[str]] = {}
    for pos, residues in zip(pos_candidates, (["S", "T"], ["K"], ["E"])):
        upper_bound = max(3, seq_len - 1)
        normalized = min(max(3, pos), upper_bound)
        normalized_positions[normalized] = residues

    library_result = processor.create_mutant_library(
        base_sequence_id=base_sequence_id,
        mutation_map=normalized_positions,
        base_name=f"{base_sequence_id}_mut",
        include_wildtype=True,
        limit=10,
        register=True,
        dataset_name=f"{base_sequence_id}_mutants",
        materialize_entities=False,
        return_metadata=True,
    )

    variants, variant_meta, library_path = library_result
    print(f"Generated {len(variants)} variants (including wildtype)")
    print(variant_meta.head().to_string(index=False))
    if library_path:
        print(f"Mutant library FASTA: {library_path}")

    mutant_conservation = processor.compute_conservation(sequences=variants)
    print("Mutant conservation (first 5 positions):")
    print(
        mutant_conservation.head()[
            ["position", "consensus", "consensus_frequency", "normalized_entropy"]
        ].to_string(index=False)
    )

    mutant_linkage = processor.compute_linkage(sequences=variants, top_k=5)
    if not mutant_linkage.empty:
        print("Mutant linkage pairs:")
        print(mutant_linkage.to_string(index=False))


if __name__ == "__main__":
    main()
