#!/usr/bin/env python3
"""Demonstrate graph generation from structure data."""

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


def main() -> None:
    ensure_data_root()

    from protos.processing.structure import StructureProcessor
    from protos.processing.graph import GraphProcessor

    struct_proc = StructureProcessor()
    available_structures = struct_proc.list_entities()

    if not available_structures:
        print("No registered structures found. Run a structure workflow first.")
        return

    structure_ids = available_structures[:1]
    print("Generating atom-level graphs for:", ", ".join(structure_ids))

    structure_dataset = "graph_demo_structures"
    existing_datasets = struct_proc.list_datasets()
    if structure_dataset not in existing_datasets:
        struct_proc.create_dataset(structure_dataset, structure_ids, {})

    graph_proc = GraphProcessor(default_cutoff=5.0)

    graph_dataset, graph_entities = graph_proc.generate_graphs_from_dataset(
        structure_dataset=structure_dataset,
        structure_processor=struct_proc,
        level="atom",
        cutoff=5.0,
        include_hydrogens=False,
        dataset_name=f"{structure_dataset}__atom_graphs",
    )

    print(f"Stored {len(graph_entities)} graph entities in dataset '{graph_dataset}'")

    for name in graph_entities:
        graph = graph_proc.load_graph(name)
        meta = graph["graph_metadata"]
        print(
            f"  • {name}: nodes={meta['node_count']} edges={meta['edge_count']} cutoff={meta['cutoff']}Å"
        )


if __name__ == "__main__":
    main()
