#!/usr/bin/env python3
"""Quick structure-loading workflow demo."""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def download_structures() -> None:
    from protos.io.ingest.structure_loader import StructureLoader

    loader = StructureLoader()

    pdb_ids = ["1ubq", "7vvl"]
    print("Downloading individual PDB structures...")
    for pdb_id in pdb_ids:
        try:
            result = loader.download_and_register(pdb_id)
            if result:
                print(f"  ✓ {pdb_id} registered as {result}")
            else:
                print(f"  ✗ No registration for {pdb_id}")
        except Exception as exc:  # noqa: BLE001
            print(f"  ✗ {pdb_id} download error: {exc}")

    batch_ids = ["3sn6", "5d5a", "6b73"]
    print("\nDownloading GPCR batch and creating dataset...")
    try:
        success, failed = loader.download_batch(batch_ids, dataset_name="gpcr_structures")
        print(f"  ✓ Batch complete ({len(success)} succeeded, {len(failed)} failed)")
        if failed:
            print(f"    Failed IDs: {failed}")
    except Exception as exc:  # noqa: BLE001
        print(f"  ✗ Batch download error: {exc}")


def inspect_processor_state() -> None:
    from protos.processing.structure import StructureProcessor

    processor = StructureProcessor()
    entities = processor.list_entities()
    datasets = processor.list_datasets()

    print("\nRegistered structures:")
    for name in sorted(entities)[:10]:
        print(f"  - {name}")
    if len(entities) > 10:
        print(f"  ... {len(entities) - 10} more")

    print("\nDatasets:")
    for dataset in datasets:
        print(f"  - {dataset}")

    if "gpcr_structures" in datasets:
        data = processor.load_dataset("gpcr_structures", return_format="dict")
        print("\nLoaded 'gpcr_structures' dataset members:")
        for name, df in data.items():
            atom_count = len(df) if df is not None else 0
            print(f"  - {name}: {atom_count} atoms")


def register_chain_sequences_demo(dataset_name: str = "gpcr_structures") -> None:
    from protos.processing.structure import StructureProcessor

    processor = StructureProcessor()

    try:
        structure_ids = processor.get_dataset_entities(dataset_name)
    except FileNotFoundError:
        print(f"\nDataset '{dataset_name}' not available; skipping chain registration demo.")
        return

    print(f"\nRegistering chain sequences for dataset '{dataset_name}'...")
    summary = processor.register_chain_sequences(
        structure_ids,
        dataset_prefix=f"{dataset_name}_chains",
        create_dataset=True,
    )

    for structure_id, info in summary.items():
        registered = info.get('registered_entities', [])
        dataset = info.get('dataset')
        print(
            f"  - {structure_id}: {len(registered)} chains registered"
            + (f" (dataset: {dataset})" if dataset else "")
        )

    related = processor.list_dataset_related_sequences(dataset_name)
    if related:
        print("\nChain sequence relationships:")
        for structure_id, relations in related.items():
            if not relations:
                continue
            chain_labels = ", ".join(rel['name'] for rel in relations)
            print(f"  - {structure_id}: {chain_labels}")


def main() -> None:
    from protos.io.paths import get_protos_paths

    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))

    paths = get_protos_paths()
    print("=== Structure Loading Workflow ===\n")
    print(f"Data root: {paths.data_root}\n")

    download_structures()
    inspect_processor_state()
    register_chain_sequences_demo()


if __name__ == "__main__":
    main()
