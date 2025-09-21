#!/usr/bin/env python3
"""Structure alignment workflow using only StructureProcessor and StructureLoader."""

from pathlib import Path

import protos

GPCR_IDS = ["3sn6", "5d5a", "6b73"]


def ensure_data_root() -> Path:
    root = Path(__file__).resolve().parent / "data"
    root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(root))
    return root


def ensure_gpcr_structures(processor, loader) -> list[str]:
    dataset_name = "gpcr_structures"
    try:
        members = processor.get_dataset_entities(dataset_name)
    except Exception:
        members = []

    if not members:
        success, failed = loader.download_batch(
            GPCR_IDS,
            dataset_name=dataset_name,
            create_dataset=True,
        )
        if failed:
            print(f"✗ Failed to download: {failed}")
        members = success

    return members


def load_structures(processor, identifiers) -> list[str]:
    available: list[str] = []
    for pdb_id in identifiers:
        df = processor.load_entity(pdb_id)
        if df is not None:
            print(f"✓ Loaded {pdb_id}")
            available.append(pdb_id)
        else:
            print(f"✗ Missing structure {pdb_id}")
    return available


def align_and_report(processor, reference, mobiles):
    rmsd_map, results = processor.align_structures(
        structure_ids=mobiles,
        reference_id=reference,
        apply_transform=True,
    )

    print("\nRMSD summary (target -> reference = value):")
    for target, refs in rmsd_map.items():
        for ref, rmsd in refs.items():
            suffix = f"{rmsd:.3f} Å" if rmsd == rmsd else "NaN"
            print(f"  {target} -> {ref}: {suffix}")

    return results


def export_aligned(processor, aligned_results, reference_id: str):
    from protos.io.paths import get_protos_paths

    output_dir = Path(get_protos_paths().data_root) / "alignment_output"
    output_dir.mkdir(exist_ok=True)

    for struct_id, result in aligned_results.items():
        if result is None or result.error:
            continue

        aligned_id = result.aligned_id
        df = processor.frames.get(aligned_id)
        if df is not None:
            processor.save_entity(
                aligned_id,
                df,
                metadata={"aligned_to": reference_id, "rmsd": result.rmsd},
            )

        out_path = output_dir / f"{aligned_id}.cif"
        try:
            processor.export_entity(aligned_id, out_path, format="cif", overwrite=True)
            print(f"✓ Exported {aligned_id} -> {out_path}")
        except Exception as exc:  # noqa: BLE001
            print(f"✗ Export failed for {aligned_id}: {exc}")


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.structure import StructureProcessor
    from protos.io.ingest.structure_loader import StructureLoader

    processor = StructureProcessor()
    loader = StructureLoader(processor=processor)

    members = ensure_gpcr_structures(processor, loader)
    if len(members) < 2:
        print("Need at least two GPCR structures to align.")
        return

    available = load_structures(processor, members)
    if len(available) < 2:
        print("Loaded fewer than two structures; aborting alignment.")
        return

    reference, *mobile_ids = available
    alignment_results = align_and_report(processor, reference, mobile_ids)
    export_aligned(processor, alignment_results, reference)


if __name__ == "__main__":
    main()
