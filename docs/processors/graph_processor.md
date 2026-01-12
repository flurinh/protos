# GraphProcessor

Creates PyTorch Geometric compatible graphs from protein structures.

**Location:** `protos.processing.graph.graph_processor`

**Processor Type:** `graph`

**Dependencies:** torch, torch_geometric (optional)

```bash
pip install torch torch_geometric
```

## Overview

The GraphProcessor provides:
- Conversion of protein structures to graph representations
- Atom-level and residue-level graphs
- Distance-based edge construction
- Filtering by chain, GRN positions, or ligand proximity
- PyTorch Geometric Data format output

---

## API Reference

### Core Entity Methods

| Method | Description |
|--------|-------------|
| `load_entity(name)` | Load saved graph by name. Returns dict or PyG Data. |
| `save_entity(name, graph, metadata)` | Save graph and register. |
| `load_graph(entity_name)` | Load graph (alias for load_entity). |
| `list_graphs(dataset)` | List registered graph entities. |

### Graph Generation

| Method | Description |
|--------|-------------|
| `generate_graph(structure_id, cutoff, level, chain, grn_positions, protein_only)` | Convert single structure to graph. |
| `generate_graphs_from_dataset(dataset_name, cutoff, level, output_dataset)` | Batch convert structures to graphs. |
| `to_pyg(graph_dict)` | Convert dict representation to PyTorch Geometric Data. |

### Structure Filtering

| Method | Description |
|--------|-------------|
| `_filter_structure(df, chain, grn_positions, near_hetatm, protein_only)` | Filter structure before graph conversion. |
| `_filter_near_hetatm(df, hetatm_name, distance)` | Filter to residues near ligand. |
| `_normalize_columns(df)` | Normalize column names to canonical form. |

### Properties

| Property | Description |
|----------|-------------|
| `struct_proc` | Lazy-loaded StructureProcessor instance. |
| `dependencies_available` | Whether PyTorch Geometric is installed. |
| `graphs_dir` | Directory for saved graphs. |
| `datasets_dir` | Directory for dataset definitions. |
| `default_cutoff` | Default edge distance cutoff (Angstroms). |
| `default_level` | Default graph level ("atom" or "residue"). |

---

## Graph Levels

| Level | Nodes | Edges | Use Case |
|-------|-------|-------|----------|
| `atom` | One per atom | Atoms within cutoff distance | Fine-grained analysis |
| `residue` | One per residue (CA) | Residues within cutoff | Coarse-grained, faster |

---

## Graph Output Format

### Dictionary Format

```python
{
    'x': np.ndarray,           # Node features (N, F)
    'pos': np.ndarray,         # 3D coordinates (N, 3)
    'edge_index': np.ndarray,  # Edge connectivity (2, E)
    'edge_attr': np.ndarray,   # Edge features (E, Fe) [optional]
    'num_nodes': int,
    'structure_id': str,
    'metadata': dict
}
```

### PyTorch Geometric Format

```python
Data(
    x=[N, F],           # Node features
    pos=[N, 3],         # Coordinates
    edge_index=[2, E],  # Edge connectivity
    edge_attr=[E, Fe],  # Edge features
)
```

---

## Usage Examples

### Basic Graph Generation

```python
from protos.processing.graph import GraphProcessor

proc = GraphProcessor(
    default_cutoff=5.0,
    default_level="atom"
)

# Generate graph from structure
graph = proc.generate_graph("1ubq")

print(f"Nodes: {graph['num_nodes']}")
print(f"Edges: {graph['edge_index'].shape[1]}")
print(f"Node features: {graph['x'].shape}")
```

### Residue-Level Graphs

```python
# Coarser graph using CA atoms
graph = proc.generate_graph(
    "1ubq",
    cutoff=8.0,
    level="residue"
)

# Much fewer nodes than atom-level
print(f"Residue nodes: {graph['num_nodes']}")
```

### Filtering Options

```python
# Filter to specific chain
graph = proc.generate_graph(
    "3sn6",
    chain="R",
    protein_only=True
)

# Filter to specific GRN positions
graph = proc.generate_graph(
    "3sn6",
    grn_positions=["3.50x50", "6.50x50", "7.50x50"]
)

# Include HETATM records
graph = proc.generate_graph(
    "1ubq",
    protein_only=False
)
```

### Batch Processing

```python
# Convert entire dataset
proc.generate_graphs_from_dataset(
    "my_structures",
    cutoff=5.0,
    level="atom",
    output_dataset="my_graphs"
)

# List generated graphs
graphs = proc.list_graphs("my_graphs")
```

### PyTorch Geometric Integration

```python
# Generate and convert to PyG
graph_dict = proc.generate_graph("1ubq", cutoff=5.0)

if proc.dependencies_available:
    pyg_data = proc.to_pyg(graph_dict)

    print(f"PyG Data: {pyg_data}")
    print(f"Nodes: {pyg_data.num_nodes}")
    print(f"Edges: {pyg_data.num_edges}")

    # Use with PyG DataLoader
    from torch_geometric.loader import DataLoader

    graphs = [proc.generate_graph(s) for s in ["s1", "s2", "s3"]]
    pyg_graphs = [proc.to_pyg(g) for g in graphs]

    loader = DataLoader(pyg_graphs, batch_size=2)
    for batch in loader:
        print(f"Batch: {batch.num_graphs} graphs")
```

### Saving and Loading

```python
# Save graph
graph = proc.generate_graph("1ubq", cutoff=5.0)
proc.save_entity("1ubq_graph", graph, metadata={
    "cutoff": 5.0,
    "level": "atom"
})

# Load graph
loaded = proc.load_entity("1ubq_graph")

# Convert loaded dict to PyG
pyg = proc.to_pyg(loaded)
```

### Graph Neural Network Pipeline

```python
from protos import StructureLoader
from protos.processing.graph import GraphProcessor

# Load structures
loader = StructureLoader()
loader.download_batch(["1ubq", "3sn6", "4lde"], dataset_name="study")

# Create graphs
graph_proc = GraphProcessor(default_cutoff=8.0, default_level="residue")

graph_proc.generate_graphs_from_dataset(
    "study",
    output_dataset="study_graphs"
)

# Load for training
if graph_proc.dependencies_available:
    from torch_geometric.loader import DataLoader

    graphs = []
    for name in graph_proc.list_graphs("study_graphs"):
        g = graph_proc.load_entity(name)
        graphs.append(graph_proc.to_pyg(g))

    train_loader = DataLoader(graphs, batch_size=2, shuffle=True)
```

### Cutoff Distance Effects

```python
# Sparse graph (few edges)
sparse = proc.generate_graph("1ubq", cutoff=4.0)
print(f"Sparse edges: {sparse['edge_index'].shape[1]}")

# Dense graph (many edges)
dense = proc.generate_graph("1ubq", cutoff=10.0)
print(f"Dense edges: {dense['edge_index'].shape[1]}")
```

### Fallback Without PyG

```python
proc = GraphProcessor()

if not proc.dependencies_available:
    print("PyG not installed - using numpy dict format")

# Graphs saved as numpy dicts regardless
graph = proc.generate_graph("1ubq")
proc.save_entity("1ubq_graph", graph)

# Can still load and use numpy arrays
loaded = proc.load_graph("1ubq_graph")
print(f"Node features: {loaded['x'].shape}")
```
