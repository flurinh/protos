#!/usr/bin/env python3
"""Minimal GraphProcessor test — matches draft code example."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
protos.set_data_path(str(REPO_ROOT / "data"))

from protos.processing.graph import GraphProcessor

graph_proc = GraphProcessor(default_cutoff=6.0, default_level="residue")

# Generate binding pocket graph: residues within 7A of retinal, edges at 6A
graph_name = graph_proc.generate_graph(
    "1U19", chain="A", near_hetatm=("RET", 7.0), cutoff=6.0
)
graph_data = graph_proc.load_entity(graph_name)

pos = graph_data.get("node_positions", [])
edge_index = graph_data.get("edge_index", [])
print(f"Graph: {len(pos)} nodes, {edge_index.shape[1]} edges")
