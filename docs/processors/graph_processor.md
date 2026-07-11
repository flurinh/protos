# GraphProcessor

Import: `from protos.processing.graph import GraphProcessor`

`GraphProcessor` converts a structure DataFrame to a distance graph and stores
the graph as a pickled dictionary of NumPy arrays. PyTorch Geometric is not
required for graph generation or storage; it is required only for `to_pyg()`.

```python
import pandas as pd
from protos.processing.graph import GraphProcessor

atoms = pd.DataFrame(
    {
        "group": ["ATOM", "ATOM", "ATOM"],
        "atom_name": ["CA", "CA", "CA"],
        "element": ["C", "C", "C"],
        "auth_chain_id": ["A", "A", "A"],
        "auth_seq_id": [1, 2, 3],
        "res_name": ["ALA", "GLY", "SER"],
        "x": [0.0, 1.0, 4.0],
        "y": [0.0, 0.0, 0.0],
        "z": [0.0, 0.0, 0.0],
    }
)

graphs = GraphProcessor()
name = graphs.generate_graph(
    "toy_structure", structure=atoms, cutoff=1.5, entity_name="toy_graph"
)
graph = graphs.load_graph(name)

assert graph["graph_metadata"]["node_count"] == 3
assert graph["edge_index"].shape == (2, 2)
```

Edges are directed pairs for all distinct nodes whose Euclidean distance is at
most the cutoff; therefore one undirected contact produces two columns in
`edge_index`.

`generate_graph()` supports atom- or residue-level nodes, chain and GRN
filters, optional ligand proximity filtering, HETATM filtering, and an optional
graph dataset name. `generate_graphs_from_dataset()` applies those options to a
managed structure dataset. `load_entity()` catches load errors and returns
`None`, whereas `load_graph()` raises when an entity or file is missing.

Install PyTorch and PyTorch Geometric using versions appropriate for the local
platform before calling `to_pyg()`. Check `dependencies_available` at runtime;
ProtOS falls back to the NumPy dictionary representation when they are absent.
