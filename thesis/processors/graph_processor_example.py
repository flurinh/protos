#!/usr/bin/env python3
"""GraphProcessor Example: Binding pocket graph generation.

Demonstrates ProtOS capabilities:
- StructureLoader: Download structure from PDB
- StructureProcessor: Analyze ligand interactions
- GraphProcessor: Generate residue contact graphs

Question: "How do residues around retinal interact in bacteriorhodopsin?"
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.graph import GraphProcessor
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "graph"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

PDB_ID = "1C3W"
STRUCTURE_NAME = f"bacteriorhodopsin_{PDB_ID}"
CHAIN = "A"
LIGAND = "RET"
POCKET_DISTANCE = 7.0
EDGE_CUTOFF = 6.0


def main() -> int:
    """Run the GraphProcessor example."""
    print("=" * 60)
    print("GRAPH PROCESSOR EXAMPLE")
    print("=" * 60)
    print(f"Structure: {PDB_ID} (bacteriorhodopsin)")
    print(f"Ligand: {LIGAND}, Pocket: {POCKET_DISTANCE}A, Edges: {EDGE_CUTOFF}A")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Download structure using ProtOS StructureLoader
    # -------------------------------------------------------------------------
    print("\n[1] Downloading structure...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    # Download and register structure
    loader.download_and_register(f"pdb:{PDB_ID}", name=STRUCTURE_NAME)

    # Verify structure loaded
    entity = struct_proc.load_entity(STRUCTURE_NAME)
    n_atoms = len(entity) if entity is not None else 0
    print(f"  Loaded: {n_atoms} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Analyze binding pocket using ProtOS StructureProcessor
    # -------------------------------------------------------------------------
    print("\n[2] Analyzing binding pocket...")

    # ProtOS lists ligands in structure
    ligands = struct_proc.list_ligands(STRUCTURE_NAME)
    print(f"  Ligands found: {[l['res_name'] for l in ligands]}")

    # ProtOS provides built-in ligand interaction analysis
    interactions = struct_proc.get_ligand_interactions(
        STRUCTURE_NAME,
        ligand_id=LIGAND,
        chain_id=CHAIN,
        cutoff=POCKET_DISTANCE,
    )

    binding_residues = interactions.get("binding_residues", [])
    print(f"  Binding pocket residues: {len(binding_residues)}")
    if binding_residues:
        res_names = [f"{r['res_name']}{r['res_id']}" for r in binding_residues[:5]]
        print(f"  Closest: {', '.join(res_names)}...")

    # -------------------------------------------------------------------------
    # Step 3: Generate graph using ProtOS GraphProcessor
    # -------------------------------------------------------------------------
    print("\n[3] Generating residue contact graph...")
    graph_proc = GraphProcessor(default_cutoff=EDGE_CUTOFF, default_level="residue")

    # ProtOS generates graph with built-in filtering
    graph_name = graph_proc.generate_graph(
        STRUCTURE_NAME,
        chain=CHAIN,
        near_hetatm=(LIGAND, POCKET_DISTANCE),
        protein_only=True,
        level="residue",
        cutoff=EDGE_CUTOFF,
    )

    # Load graph through processor
    graph_data = graph_proc.load_entity(graph_name)
    if graph_data:
        pos = graph_data.get("node_positions", np.array([]))
        n_nodes = len(pos) if pos.size > 0 else 0
        edge_index = graph_data.get("edge_index", np.array([[], []]))
        n_edges = edge_index.shape[1] if edge_index.ndim == 2 else 0
        edge_weight = graph_data.get("edge_weight", np.array([]))
        mean_dist = np.mean(edge_weight) if edge_weight.size > 0 else 0

        print(f"  Nodes: {n_nodes}, Edges: {n_edges}")
        print(f"  Mean edge distance: {mean_dist:.2f}A")

    # -------------------------------------------------------------------------
    # Step 4: Generate graphs at multiple cutoffs
    # -------------------------------------------------------------------------
    print("\n[4] Comparing cutoff distances...")

    for cutoff in [4.0, 5.0, 6.0, 8.0]:
        name = graph_proc.generate_graph(
            STRUCTURE_NAME, chain=CHAIN, near_hetatm=(LIGAND, POCKET_DISTANCE),
            protein_only=True, level="residue", cutoff=cutoff,
        )
        g = graph_proc.load_entity(name)
        if g:
            n_e = len(g.get("edge_index", [[]])[0])
            print(f"  {cutoff}A: {n_e} edges")

    # -------------------------------------------------------------------------
    # Step 5: Visualize binding pocket graph (3D)
    # -------------------------------------------------------------------------
    print("\n[5] Creating 3D visualization...")
    import plotly.graph_objects as go

    fig = go.Figure()

    if graph_data:
        pos = graph_data.get("node_positions", np.array([]))
        edge_index = graph_data.get("edge_index", np.array([[], []]))
        node_metadata = graph_data.get("node_metadata", [])
        node_names = [f"{m['residue_name']}{m['residue_number']}" for m in node_metadata] if node_metadata else [f"R{i}" for i in range(len(pos))]

        # Draw edges and nodes only if we have valid position data
        if pos.ndim == 2 and pos.shape[0] > 0:
            # Find nodes with edges (filter isolated nodes)
            nodes_with_edges = set()
            if edge_index.ndim == 2 and edge_index.shape[1] > 0:
                nodes_with_edges = set(edge_index[0].tolist()) | set(edge_index[1].tolist())

            # Draw edges
            if edge_index.ndim == 2 and edge_index.shape[1] > 0:
                for i in range(edge_index.shape[1]):
                    src, tgt = int(edge_index[0, i]), int(edge_index[1, i])
                    if src < len(pos) and tgt < len(pos):
                        fig.add_trace(go.Scatter3d(
                            x=[pos[src, 0], pos[tgt, 0]], y=[pos[src, 1], pos[tgt, 1]], z=[pos[src, 2], pos[tgt, 2]],
                            mode="lines", line=dict(color="gray", width=2), showlegend=False,
                        ))

            # Draw only nodes that have edges (no text labels)
            connected_indices = sorted(nodes_with_edges)
            if connected_indices:
                connected_pos = pos[connected_indices]
                connected_names = [node_names[i] for i in connected_indices]
                fig.add_trace(go.Scatter3d(
                    x=connected_pos[:, 0], y=connected_pos[:, 1], z=connected_pos[:, 2],
                    mode="markers",
                    marker=dict(size=6, color="#1f77b4"),
                    text=connected_names,
                    hovertemplate="%{text}<extra></extra>",
                    name="Residues",
                ))

    # Draw ligand atoms and bonds from structure
    df = struct_proc.load_entity(STRUCTURE_NAME)
    if df is not None:
        ligand_df = df[(df["res_name"] == LIGAND) & (~df["element"].isin(["H"]))].copy()
        if len(ligand_df) > 0:
            ligand_coords = ligand_df[["x", "y", "z"]].values

            # Draw chemical bonds (atoms within ~1.8A are bonded)
            bond_cutoff = 1.8
            for i in range(len(ligand_coords)):
                for j in range(i + 1, len(ligand_coords)):
                    dist = np.linalg.norm(ligand_coords[i] - ligand_coords[j])
                    if dist < bond_cutoff:
                        fig.add_trace(go.Scatter3d(
                            x=[ligand_coords[i, 0], ligand_coords[j, 0]],
                            y=[ligand_coords[i, 1], ligand_coords[j, 1]],
                            z=[ligand_coords[i, 2], ligand_coords[j, 2]],
                            mode="lines", line=dict(color="#d62728", width=3),
                            showlegend=False,
                        ))

            # Draw ligand atoms
            fig.add_trace(go.Scatter3d(
                x=ligand_df["x"], y=ligand_df["y"], z=ligand_df["z"],
                mode="markers", marker=dict(size=5, color="#d62728"),
                name="Retinal",
            ))

    fig.update_layout(
        title=f"Binding Pocket Graph - 3D ({PDB_ID})",
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            aspectmode="data",
        ),
        width=800, height=700,
    )
    fig.write_html(str(FIGURES_DIR / "graph_binding_pocket_3d.html"))
    print(f"  Saved: {FIGURES_DIR / 'graph_binding_pocket_3d.html'}")

    # -------------------------------------------------------------------------
    # Step 6: Create 2D contact network visualization
    # -------------------------------------------------------------------------
    print("\n[6] Creating 2D contact network...")
    import networkx as nx

    fig2d = go.Figure()

    if graph_data:
        pos = graph_data.get("node_positions", np.array([]))
        edge_index = graph_data.get("edge_index", np.array([[], []]))
        edge_weight = graph_data.get("edge_weight", [])
        node_metadata = graph_data.get("node_metadata", [])
        node_names = [f"{m['residue_name']}{m['residue_number']}" for m in node_metadata] if node_metadata else [f"R{i}" for i in range(len(pos))]

        if pos.ndim == 2 and pos.shape[0] > 0:
            # Build networkx graph for 2D layout
            G = nx.Graph()
            for i, name in enumerate(node_names):
                G.add_node(i, label=name)

            if edge_index.ndim == 2 and edge_index.shape[1] > 0:
                for i in range(edge_index.shape[1]):
                    src, tgt = int(edge_index[0, i]), int(edge_index[1, i])
                    if src < len(pos) and tgt < len(pos):
                        weight = edge_weight[i] if i < len(edge_weight) else 1.0
                        G.add_edge(src, tgt, weight=weight)

            # Use spring layout for 2D positioning
            pos_2d = nx.spring_layout(G, seed=42, k=2.0, iterations=100)

            # Draw edges with distance-based color
            if edge_index.ndim == 2 and edge_index.shape[1] > 0:
                for i in range(edge_index.shape[1]):
                    src, tgt = int(edge_index[0, i]), int(edge_index[1, i])
                    if src in pos_2d and tgt in pos_2d:
                        dist = edge_weight[i] if i < len(edge_weight) else 5.0
                        # Color based on distance: closer = darker blue
                        intensity = max(0, min(1, 1 - (dist - 3) / 5))
                        color = f"rgba(31, 119, 180, {0.3 + 0.5 * intensity})"
                        fig2d.add_trace(go.Scatter(
                            x=[pos_2d[src][0], pos_2d[tgt][0]],
                            y=[pos_2d[src][1], pos_2d[tgt][1]],
                            mode="lines",
                            line=dict(color=color, width=1 + 2 * intensity),
                            showlegend=False,
                            hoverinfo="skip",
                        ))

            # Draw nodes
            node_x = [pos_2d[i][0] for i in range(len(node_names)) if i in pos_2d]
            node_y = [pos_2d[i][1] for i in range(len(node_names)) if i in pos_2d]
            node_text = [node_names[i] for i in range(len(node_names)) if i in pos_2d]

            # Color nodes by degree (number of contacts)
            degrees = [G.degree(i) for i in range(len(node_names)) if i in pos_2d]

            fig2d.add_trace(go.Scatter(
                x=node_x, y=node_y, mode="markers+text",
                marker=dict(
                    size=12,
                    color=degrees,
                    colorscale="Blues",
                    showscale=True,
                    colorbar=dict(title="Contacts"),
                    line=dict(width=1, color="black"),
                ),
                text=node_text,
                textposition="top center",
                textfont=dict(size=9),
                name="Residues",
                hovertemplate="%{text}<br>Contacts: %{marker.color}<extra></extra>",
            ))

    fig2d.update_layout(
        title=f"Binding Pocket Contact Network - 2D ({PDB_ID})",
        xaxis=dict(visible=False, showgrid=False, zeroline=False),
        yaxis=dict(visible=False, showgrid=False, zeroline=False),
        plot_bgcolor="white",
        paper_bgcolor="white",
        width=800,
        height=700,
        showlegend=False,
    )
    fig2d.write_html(str(FIGURES_DIR / "graph_contact_network_2d.html"))
    print(f"  Saved: {FIGURES_DIR / 'graph_contact_network_2d.html'}")

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Outputs: {OUTPUT_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
