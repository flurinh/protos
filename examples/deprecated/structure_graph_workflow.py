#!/usr/bin/env python3
"""Generate residue-level interaction graphs for the rhodopsin state dataset."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List

import protos
from protos.processing.graph import GraphProcessor
from protos.processing.structure import StructureProcessor

DATA_RELATIVE_ROOT = Path(__file__).resolve().parents[2] / "data"
STRUCTURE_DATASET = "rhodopsin_states"
GRAPH_DATASET_NAME = "rhodopsin_states_residue_graphs"
OUTPUT_DIR = DATA_RELATIVE_ROOT / "reports"
SUMMARY_PATH = OUTPUT_DIR / "structure_graph_summary.json"
GRAPH_LEVEL = "residue"
CUTOFF = 8.0


def ensure_data_root() -> Path:
    DATA_RELATIVE_ROOT.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(DATA_RELATIVE_ROOT))
    return DATA_RELATIVE_ROOT


@dataclass
class GraphStats:
    dataset: str
    generated_graphs: int
    level: str
    cutoff: float
    node_count: int
    edge_count: int


def coerce_count(value) -> int:
    if value is None:
        return 0
    if hasattr(value, "shape"):
        return int(value.shape[0])
    return len(value)


def coerce_edge_count(edge_index) -> int:
    if edge_index is None:
        return 0
    if hasattr(edge_index, "shape"):
        return int(edge_index.shape[-1])
    if isinstance(edge_index, (list, tuple)) and edge_index:
        return len(edge_index[0])
    return 0


def summarize_graphs(graph_processor: GraphProcessor, graph_ids: List[str]) -> GraphStats:
    first_graph = graph_processor.load_graph(graph_ids[0])
    node_positions = first_graph.get("node_positions")
    if node_positions is None:
        node_positions = first_graph.get("node_features")
    node_count = coerce_count(node_positions)
    edge_count = coerce_edge_count(first_graph.get("edge_index"))

    return GraphStats(
        dataset=GRAPH_DATASET_NAME,
        generated_graphs=len(graph_ids),
        level=GRAPH_LEVEL,
        cutoff=CUTOFF,
        node_count=node_count,
        edge_count=edge_count,
    )


def main() -> None:
    ensure_data_root()
    structure_processor = StructureProcessor()
    graph_processor = GraphProcessor()

    dataset_name, graph_ids = graph_processor.generate_graphs_from_dataset(
        STRUCTURE_DATASET,
        structure_processor=structure_processor,
        dataset_name=GRAPH_DATASET_NAME,
        level=GRAPH_LEVEL,
        cutoff=CUTOFF,
    )

    stats = summarize_graphs(graph_processor, graph_ids)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    payload: Dict[str, Dict[str, object]] = {
        "graphs": asdict(stats)
    }
    SUMMARY_PATH.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
