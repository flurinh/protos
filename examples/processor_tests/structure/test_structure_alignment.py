#!/usr/bin/env python3
"""Structure alignment workflow using only StructureProcessor and StructureLoader."""

from __future__ import annotations

from pathlib import Path
from datetime import datetime

import protos

GPCR_IDS = ["3sn6", "5d5a", "6b73"]


def ensure_data_root() -> Path:
    root = Path(__file__).resolve().parents[3] / "data"  # -> repo root/data
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


def align_structures(processor, reference_id: str, mobile_ids: list[str]):
    summary, _ = processor.align_and_record(
        structure_ids=mobile_ids,
        reference_id=reference_id,
        method="cealign",
        atom_selection="CA",
        apply_transform=True,
        save_aligned=True,
        summary_name=f"{reference_id}_gpcr_alignment",
        aligned_dataset_name="gpcr_alignment_aligned",
        aligned_dataset_include_reference=True,
        property_table_name="gpcr_structure_alignment_properties",
    )

    global_stats = summary.get("rmsd", {}).get("global", {})
    print("\nRMSD summary (aligned to reference):")
    for key in ("min", "mean", "max"):
        value = global_stats.get(key)
        if value is not None:
            print(f"  {key}: {value:.3f} Å")

    for row in summary.get("rmsd", {}).get("pairwise", []):
        rmsd = row.get("rmsd")
        suffix = f"{rmsd:.3f} Å" if rmsd is not None else "NaN"
        print(f"  {row['target_id']} -> {row['reference_id']}: {suffix}")

    return summary


OVERWRITE_ORIGINAL = True  # Set to False to create versioned copies instead of overwriting
VERSION_SUFFIX = "_aligned"  # Used when OVERWRITE_ORIGINAL is False


def export_aligned(processor, aligned_entities: list[str], reference_id: str):
    from protos.io.paths import get_protos_paths

    if not aligned_entities:
        print("No aligned entities to export.")
        return

    paths = get_protos_paths()
    mmcif_dir = Path(processor.path_cif_dir)
    mmcif_dir.mkdir(parents=True, exist_ok=True)

    export_map: dict[str, Path] = {}
    for struct_id in aligned_entities:
        df = processor.load_entity(struct_id)
        if df is None:
            print(f"✗ Skipping {struct_id}: aligned frame not available")
            continue

        target_id = struct_id
        metadata = {
            "aligned_to": reference_id,
            "updated_at": datetime.utcnow().isoformat(),
        }

        if OVERWRITE_ORIGINAL:
            processor.save_entity(struct_id, df, metadata=metadata)
        else:
            versioned_id = f"{struct_id}{VERSION_SUFFIX}"
            processor.save_entity(versioned_id, df, metadata={**metadata, "source_entity": struct_id})
            target_id = versioned_id

        out_path = mmcif_dir / f"{target_id}.cif"
        processor.export_entity(target_id, out_path, format="cif", overwrite=True)
        export_map[target_id] = out_path

    for target_id, out_path in export_map.items():
        print(f"✓ Exported {target_id} -> {out_path}")


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
    summary = align_structures(processor, reference, mobile_ids)
    export_aligned(processor, summary.get("aligned_entities", []), reference)

    if summary.get("summary_file"):
        print(f"\nSummary written to: {summary['summary_file']}")
    if summary.get("summary_dataset"):
        print(f"Summary dataset: {summary['summary_dataset']}")


if __name__ == "__main__":
    main()
