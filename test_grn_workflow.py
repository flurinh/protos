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
        processor.load_dataset("gpcr_seqs")
        return
    except Exception:
        pass

    # Real GPCR sequences - full length examples
    # These are human GPCRs with known structures
    sequences = {
        # Beta-2 adrenergic receptor (PDB: 3SN6) - full sequence
        "ADRB2_HUMAN": "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAKFERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQEAINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRFHVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQDNLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNTGEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL",
        
        # Adenosine A2a receptor (PDB: 5G53) - full sequence  
        "AA2AR_HUMAN": "MPIMGSSVYITVELAIAVLAILGNVLVCWAVWLNSNLQNVTNYFVVSLAAADIAVGVLAIPFAITISTGFCAACHGCLFIACFVLVLTQSSIFSLLAIAIDRYIAIRIPLRYNGLVTGTRAKGIIAICWVLSFAIGLTPMLGWNNCGQPKEGKNHSQGCGEGQVACLFEDVVPMNYMVYFNFFACVLVPLLLMLGVYLRIFLAARRQLKQMESQPLPGERARSTLQKEVHAAKSLAIIVGLFALCWLPLHIINCFTFFCPDCSHAPLWLMYLAIVLSHTNSVVNPFIYAYRIREFRQTFRKIIRSHVLRQQEPFKAAGTSARVLAAHGSDGEQVSLRLNGHPPGVWANGSAPHPERRPNGYALGLVSGGSAQESQGNTGLPDVELLSHELKGVCPEPPGLDDPLAQDGAGVS",
        
        # Mu-opioid receptor (PDB: 5C1M) - full sequence
        "OPRM_HUMAN": "MDSSAAPTNASNCTDALAYSSCSPAPSPGSWVNLSHLDGNLSDPCGPNRTDLGGRDSLCPPTGSPSMITAITIMALYSIVCVVGLFGNFLVMYVIVRYTKMKTATNIYIFNLALADALATSTLPFQSVNYLMGTWPFGTILCKIVISIDYYNMFTSIFTLCTMSVDRYIAVCHPVKALDFRTPRNAKIINVCNWILSSAIGLPVMFMATTKYRQGSIDCTLTFSHPTWYWENLLKICVFIFAFIMPVLIITVCYGLMILRLKSVRMLSGSKEKDRNLRRITRMVLVVVAVFIVCWTPIHIYVIIKALVTIPETTFQTVSWHFCIALGYTNSCLNPVLYAFLDENFKRCFRDFCFPLKMRMERQSTSRVRNTVQDPAYLRDIDGMNKPV"
    }

    paths = get_protos_paths()
    input_dir = Path(paths.get_processor_path('input'))
    input_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = input_dir / "gpcr_seqs.fasta"
    with open(fasta_path, "w") as handle:
        for seq_id, seq in sequences.items():
            handle.write(f">{seq_id}\n{seq}\n")

    loader = SequenceLoader(processor=processor)
    loader.download_and_register(
        str(fasta_path),
        name="gpcr_seqs",
        materialize_entities=True,
    )


def record_grn_table() -> None:
    from protos.processing.sequence import SequenceProcessor
    from protos.processing.grn import GRNProcessor

    seq_proc = SequenceProcessor()

    dataset_sequences = seq_proc.load_dataset("gpcr_seqs")
    reference_id = "ADRB2_HUMAN"  # Use beta-2 adrenergic receptor as reference
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

    print(f"Similarity to {reference_id} (normalized score):")
    for seq_id, score in similarity_scores.items():
        status = "gpcr_like" if classified[seq_id] else "low_similarity"
        print(f"  {seq_id}: {score:.3f} ({status})")

    grn_table, summary = seq_proc.annotate_with_grn(
        dataset_name="gpcr_seqs",
        reference_table=REFERENCE_TABLE,
        protein_family="gpcr_a",
        output_table="gpcr_grn_demo",
        return_summary=True,
        allow_create=True,
    )

    grn_proc = GRNProcessor()
    reference_table = grn_proc.load_reference_table(REFERENCE_TABLE)

    # Ensure column coverage matches the reference definition
    assert list(grn_table.columns) == list(reference_table.columns), (
        "GRN annotation table missing reference columns"
    )

    # Check if we got any annotations (may fail if reference table not available)
    has_annotations = any((row != '-').any() for _, row in grn_table.iterrows())
    if has_annotations:
        print("✓ Successfully annotated at least one sequence")
    else:
        print("⚠ Warning: No GRN annotations were assigned (check reference table)")

    dataset_info = grn_proc.get_dataset_info("gpcr_grn_demo")
    metadata = dataset_info.get("metadata", {}) if dataset_info else {}
    
    # Basic sanity checks
    assert summary["global"]["total"] == len(grn_table), "Summary total mismatch"
    print(f"✓ Processed {summary['global']['total']} sequences")
    
    # Print annotation summary
    print("\nAnnotation Summary:")
    for seq_id, info in summary["per_sequence"].items():
        coverage = info.get("coverage", 0)
        status = info.get("status", "unknown")
        print(f"  {seq_id}: coverage={coverage:.1%}, status={status}")

    print("Saved GRN table 'gpcr_grn_demo':")
    print(grn_table)

    # Load the GRN annotations for a specific entity
    test_entity = "ADRB2_HUMAN"
    seq_annotations = grn_proc.load_entity(test_entity)
    print(f"\nGRN annotations for {test_entity}:")
    print(seq_annotations if seq_annotations else "Not found")
    
    # Alternative: Load the entire table and access specific row
    loaded_table = grn_proc.load_table("gpcr_grn_demo")
    print(f"\nGRN table row for {test_entity}:")
    if test_entity in loaded_table.index:
        print(loaded_table.loc[test_entity])
    else:
        print("Entity not found in table")

    # Check if SequenceProcessor has resolve_related_entities method
    # If not, we'll skip this part
    try:
        related = seq_proc.resolve_related_entities(test_entity, format_type="grn")
        print(f"\nRegistry relationships for {test_entity}:")
        print(related)
    except AttributeError:
        print("\nNote: resolve_related_entities not available in SequenceProcessor")


def main() -> None:
    ensure_data_root()
    register_sequences()
    record_grn_table()


if __name__ == "__main__":
    main()
