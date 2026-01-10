#!/usr/bin/env python3
"""GraphProcessor basic data management demonstration.

This script demonstrates core GraphProcessor capabilities:
- Generating graphs from structure data
- Loading and saving graph entities
- Graph metadata and properties
- Dataset management
- PyG conversion (if available)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def ensure_data_root() -> Path:
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def create_sample_graph(name: str, num_nodes: int = 10) -> dict:
    """Create a sample graph dictionary for testing."""
    np.random.seed(42)

    # Create random node features (e.g., one-hot encoding of atom types)
    node_features = np.random.randn(num_nodes, 8).astype(np.float32)

    # Create random edges (sparse connectivity)
    num_edges = num_nodes * 2
    edge_index = np.random.randint(0, num_nodes, size=(2, num_edges)).astype(np.int64)

    # Edge features (e.g., bond type encoding)
    edge_features = np.random.randn(num_edges, 4).astype(np.float32)

    # Node positions (3D coordinates)
    positions = np.random.randn(num_nodes, 3).astype(np.float32) * 10

    return {
        "node_features": node_features,
        "edge_index": edge_index,
        "edge_features": edge_features,
        "positions": positions,
        "graph_metadata": {
            "name": name,
            "node_count": num_nodes,
            "edge_count": num_edges,
            "cutoff": 5.0,
            "level": "atom",
        },
    }


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")

    from protos.processing.graph import GraphProcessor

    # Initialize processor
    print("=== GraphProcessor Basic Demo ===\n")
    graph_proc = GraphProcessor(default_cutoff=5.0)

    # 1. List existing graphs
    print("1. Listing existing graphs...")
    graphs = graph_proc.list_graphs()
    print(f"   Found {len(graphs)} graphs\n")

    # 2. Create and save sample graphs
    print("2. Creating sample graphs...")
    sample_graphs = {}
    for i, name in enumerate(["graph_A", "graph_B", "graph_C"]):
        graph_data = create_sample_graph(name, num_nodes=10 + i * 5)
        sample_graphs[name] = graph_data
        graph_proc.save_entity(
            name,
            graph_data,
            metadata={"source": "test_graph_processor.py", "index": i},
        )
        meta = graph_data["graph_metadata"]
        print(f"   Saved {name}: {meta['node_count']} nodes, {meta['edge_count']} edges")
    print()

    # 3. List graphs again
    print("3. Verifying registration...")
    graphs = graph_proc.list_graphs()
    print(f"   Now have {len(graphs)} graphs")
    for name in sample_graphs.keys():
        exists = graph_proc.entity_exists(name)
        print(f"   '{name}' exists: {exists}")
    print()

    # 4. Load a specific graph
    print("4. Loading graph 'graph_A'...")
    loaded_graph = graph_proc.load_graph("graph_A")
    if loaded_graph:
        print(f"   Keys: {list(loaded_graph.keys())}")
        meta = loaded_graph.get("graph_metadata", {})
        print(f"   Nodes: {meta.get('node_count')}")
        print(f"   Edges: {meta.get('edge_count')}")
        print(f"   Node features shape: {loaded_graph.get('node_features', np.array([])).shape}")
    print()

    # 5. Create a dataset
    print("5. Creating dataset 'demo_graphs'...")
    graph_proc.create_dataset(
        "demo_graphs",
        list(sample_graphs.keys()),
        {"description": "Demo graph dataset"},
    )
    print(f"   Created dataset with {len(sample_graphs)} graphs\n")

    # 6. List datasets
    print("6. Listing datasets...")
    datasets = graph_proc.list_datasets()
    print(f"   Found {len(datasets)} datasets")
    if "demo_graphs" in datasets:
        info = graph_proc.get_dataset_info("demo_graphs")
        print(f"   demo_graphs info: {info}\n")

    # 7. Load dataset
    print("7. Loading dataset 'demo_graphs'...")
    dataset = graph_proc.load_dataset("demo_graphs")
    print(f"   Loaded {len(dataset)} graphs\n")

    # 8. Check PyG conversion
    print("8. Checking PyTorch Geometric conversion...")
    try:
        pyg_graph = graph_proc.to_pyg(loaded_graph)
        print(f"   PyG conversion successful: {type(pyg_graph)}")
        if hasattr(pyg_graph, 'x'):
            print(f"   PyG x shape: {pyg_graph.x.shape}")
        if hasattr(pyg_graph, 'edge_index'):
            print(f"   PyG edge_index shape: {pyg_graph.edge_index.shape}")
    except Exception as e:
        print(f"   PyG not available: {e}")
    print()

    # 9. Generate graph from structure (if structures available)
    print("9. Checking structure-based graph generation...")
    from protos.processing.structure import StructureProcessor
    struct_proc = StructureProcessor()
    structures = struct_proc.list_entities()

    if structures:
        struct_id = structures[0]
        print(f"   Generating graph from structure '{struct_id}'...")
        try:
            graph = graph_proc.generate_graph(
                structure_id=struct_id,
                structure_processor=struct_proc,
                level="residue",
                cutoff=8.0,
            )
            if graph:
                meta = graph.get("graph_metadata", {})
                print(f"   Generated: {meta.get('node_count')} nodes, {meta.get('edge_count')} edges")
        except Exception as e:
            print(f"   Generation failed: {e}")
    else:
        print("   No structures available. Run a structure workflow first.")
    print()

    # 10. Export dataset
    print("10. Exporting dataset...")
    try:
        export_dir = data_root / "exports" / "graph"
        export_dir.mkdir(parents=True, exist_ok=True)
        export_info = graph_proc.export_dataset(
            "demo_graphs",
            output_dir=export_dir,
            overwrite=True,
        )
        print(f"   Exported {len(export_info)} entities to: {export_dir}\n")
    except Exception as e:
        print(f"   Export: {e}\n")

    print("=== GraphProcessor Demo Complete ===")


if __name__ == "__main__":
    main()
