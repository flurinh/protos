#!/usr/bin/env python3
"""Runtime sanity check for StructureProcessor PKL workflow and CIF export."""

import os
from pathlib import Path

import pandas as pd

from protos.processing.structure import StructureProcessor
from protos.io.export.structure_exporter import StructureExporter


def main() -> None:
    repo_root = Path(__file__).resolve().parent
    data_root = repo_root / "data"

    print("=== StructureProcessor CIF Ingest Metadata Check ===")

    data_root.mkdir(parents=True, exist_ok=True)

    os.environ["PROTOS_DATA_ROOT"] = str(data_root)

    processor = StructureProcessor(name="ingest_check")
    print(f"✓ Initialized StructureProcessor with data root {data_root}")

    # Build a minimal, valid CIF using write_cif_file so the ingest path is deterministic
    structure_id = "ingest_demo"
    structure_dir = Path(processor.paths.get_subdir_path("structure", "structure_dir"))
    structure_dir.mkdir(parents=True, exist_ok=True)

    raw_df = pd.DataFrame(
        {
            "group": ["ATOM", "ATOM"],
            "atom_id": [1, 2],
            "atom_name": ["N", "CA"],
            "res_name": ["MET", "MET"],
            "auth_chain_id": ["A", "A"],
            "auth_seq_id": [1, 1],
            "x": [1.0, 1.5],
            "y": [2.0, 2.5],
            "z": [3.0, 3.5],
        }
    )

    cache_dir = Path(processor.paths.get_subdir_path("structure", "cache_dir"))
    pkl_path = cache_dir / f"{structure_id}.pkl"
    if pkl_path.exists():
        pkl_path.unlink()

    processor.save_entity(
        structure_id,
        raw_df,
        metadata={"source": "unit_test", "note": "saved_via_script"},
    )
    print("✓ Saved structure via save_entity")

    df = processor.load_entity(structure_id)
    if df is None:
        raise RuntimeError("Failed to reload saved structure")
    if df.index.names != ["structure_id", "atom_id"]:
        raise AssertionError(f"Unexpected index names: {df.index.names}")
    print(f"✓ Reloaded structure; canonical shape {df.shape}")

    entity = processor.entity_registry.find_entity(structure_id, format_type="structure")
    if entity is None:
        raise RuntimeError(f"Entity '{structure_id}' was not registered")

    if entity.metadata.get("note") != "saved_via_script":
        raise AssertionError("Entity metadata not persisted as expected")

    exporter = StructureExporter(processor)
    export_path = exporter.export_entity(
        structure_id,
        structure_dir / f"{structure_id}.cif",
        format='cif',
        overwrite=True,
    )
    if not Path(export_path).exists():
        raise AssertionError("CIF export was not created")
    print(f"✓ Exported CIF written to {export_path}")

    print("All checks passed. StructureProcessor PKL workflow is healthy.")


if __name__ == "__main__":
    main()
