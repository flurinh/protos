#!/usr/bin/env python3
"""GraphProcessor Example: Binding pocket graph from a crystal structure.

Demonstrates ProtOS capabilities:
- StructureProcessor: Identify ligand binding pocket residues
- GraphProcessor: Convert pocket residues into a residue contact graph

Question: "What does a binding pocket graph look like?"

This example generates a residue contact graph around the retinal binding pocket
of bovine rhodopsin (1U19). Nodes are pocket residues; edges encode spatial
contacts within a distance cutoff. The resulting graph is the data format that
downstream analyses (e.g., spectral prediction in Chapter 4) consume.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))
sys.path.insert(0, str(THESIS_DIR / "shared"))

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.graph import GraphProcessor
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Load Color Scheme
# =============================================================================
with open(THESIS_DIR / "colorscales.yaml") as f:
    COLORS = yaml.safe_load(f)


# =============================================================================
# Configuration
# =============================================================================
FIGURES_DIR = Path(__file__).resolve().parent

PDB_ID = "1U19"
STRUCTURE_INFO = {
    "name": "Bovine Rhodopsin",
    "type": "Type II (Animal/GPCR)",
    "chain": "A",
    "ligand": "RET",
    "color": COLORS["structures"]["1U19"]["hex"],
}

LIGAND_COLOR = COLORS["molecules"]["retinal"]["hex"]
EDGE_COLOR = COLORS["graphs"]["edges"]["hex"]

POCKET_DISTANCE = 7.0  # Angstroms — residues within this distance of ligand
EDGE_CUTOFF = 6.0      # Angstroms — edges between residue pairs closer than this


# =============================================================================
# 3D pocket visualization
# =============================================================================
def create_3d_pocket_visualization(
    graph_data: dict,
    struct_df,
    output_path: Path,
) -> None:
    """Create 3D visualization of the binding pocket graph."""
    import plotly.graph_objects as go

    fig = go.Figure()
    color = STRUCTURE_INFO["color"]

    if graph_data:
        pos = graph_data.get("node_positions", np.array([]))
        edge_index = graph_data.get("edge_index", np.array([[], []]))
        node_metadata = graph_data.get("node_metadata", [])
        node_names = (
            [f"{m['residue_name']}{m['residue_number']}" for m in node_metadata]
            if node_metadata
            else [f"R{i}" for i in range(len(pos))]
        )

        if pos.ndim == 2 and pos.shape[0] > 0:
            nodes_with_edges = set()
            if edge_index.ndim == 2 and edge_index.shape[1] > 0:
                nodes_with_edges = set(edge_index[0].tolist()) | set(edge_index[1].tolist())

                # Draw edges
                for i in range(edge_index.shape[1]):
                    src, tgt = int(edge_index[0, i]), int(edge_index[1, i])
                    if src < len(pos) and tgt < len(pos):
                        fig.add_trace(go.Scatter3d(
                            x=[pos[src, 0], pos[tgt, 0]],
                            y=[pos[src, 1], pos[tgt, 1]],
                            z=[pos[src, 2], pos[tgt, 2]],
                            mode="lines",
                            line=dict(color=EDGE_COLOR, width=2),
                            showlegend=False,
                        ))

            # Draw nodes
            connected_indices = sorted(nodes_with_edges)
            if connected_indices:
                connected_pos = pos[connected_indices]
                connected_names = [node_names[i] for i in connected_indices]
                fig.add_trace(go.Scatter3d(
                    x=connected_pos[:, 0],
                    y=connected_pos[:, 1],
                    z=connected_pos[:, 2],
                    mode="markers",
                    marker=dict(size=6, color=color),
                    text=connected_names,
                    hovertemplate="%{text}<extra></extra>",
                    name="Pocket residues",
                ))

    # Draw ligand
    chain = STRUCTURE_INFO["chain"]
    ligand = STRUCTURE_INFO["ligand"]
    if struct_df is not None:
        ligand_df = struct_df[
            (struct_df["res_name"] == ligand)
            & (struct_df["auth_chain_id"] == chain)
            & (~struct_df["element"].isin(["H"]))
        ].copy()
        if len(ligand_df) > 0:
            ligand_coords = ligand_df[["x", "y", "z"]].values
            bond_cutoff = 1.8
            for i in range(len(ligand_coords)):
                for j in range(i + 1, len(ligand_coords)):
                    dist = np.linalg.norm(ligand_coords[i] - ligand_coords[j])
                    if dist < bond_cutoff:
                        fig.add_trace(go.Scatter3d(
                            x=[ligand_coords[i, 0], ligand_coords[j, 0]],
                            y=[ligand_coords[i, 1], ligand_coords[j, 1]],
                            z=[ligand_coords[i, 2], ligand_coords[j, 2]],
                            mode="lines",
                            line=dict(color=LIGAND_COLOR, width=3),
                            showlegend=False,
                        ))
            fig.add_trace(go.Scatter3d(
                x=ligand_df["x"],
                y=ligand_df["y"],
                z=ligand_df["z"],
                mode="markers",
                marker=dict(size=5, color=LIGAND_COLOR),
                name="Retinal",
            ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
        ),
        width=800,
        height=700,
        margin=dict(t=30, b=30),
    )
    fig.write_html(str(output_path))
    fig.write_image(str(output_path.with_suffix(".png")), scale=2)


# =============================================================================
# Main
# =============================================================================
def main() -> int:
    """Run the GraphProcessor example."""
    print("=" * 70)
    print("GRAPH PROCESSOR EXAMPLE")
    print(f"Binding Pocket Graph: {PDB_ID} ({STRUCTURE_INFO['name']})")
    print("=" * 70)
    print(f"\nParameters: Pocket={POCKET_DISTANCE}A, Edge cutoff={EDGE_CUTOFF}A")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Load structure
    # -------------------------------------------------------------------------
    print("\n[1] Loading structure...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    df = struct_proc.load_entity(PDB_ID)
    if df is None:
        try:
            loader.download_and_register(f"pdb:{PDB_ID}", name=PDB_ID)
            df = struct_proc.load_entity(PDB_ID)
        except Exception as e:
            print(f"  Failed to load {PDB_ID}: {e}")
            return 1

    if df is None:
        print(f"  {PDB_ID} not available")
        return 1
    print(f"  {PDB_ID}: {len(df)} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Identify binding pocket
    # -------------------------------------------------------------------------
    print("\n[2] Identifying retinal binding pocket...")

    interactions = struct_proc.get_ligand_interactions(
        PDB_ID,
        ligand_id=STRUCTURE_INFO["ligand"],
        chain_id=STRUCTURE_INFO["chain"],
        cutoff=POCKET_DISTANCE,
    )

    binding_residues = interactions.get("binding_residues", [])
    print(f"  Binding residues within {POCKET_DISTANCE}A: {len(binding_residues)}")

    # -------------------------------------------------------------------------
    # Step 3: Generate residue contact graph
    # -------------------------------------------------------------------------
    print("\n[3] Generating binding pocket graph...")
    graph_proc = GraphProcessor(default_cutoff=EDGE_CUTOFF, default_level="residue")

    graph_name = graph_proc.generate_graph(
        PDB_ID,
        chain=STRUCTURE_INFO["chain"],
        near_hetatm=(STRUCTURE_INFO["ligand"], POCKET_DISTANCE),
        protein_only=True,
        level="residue",
        cutoff=EDGE_CUTOFF,
    )

    graph_data = graph_proc.load_entity(graph_name)
    if not graph_data:
        print("  Graph generation failed")
        return 1

    pos = graph_data.get("node_positions", np.array([]))
    n_nodes = len(pos) if pos.size > 0 else 0
    edge_index = graph_data.get("edge_index", np.array([[], []]))
    n_edges = edge_index.shape[1] if edge_index.ndim == 2 else 0
    edge_weight = graph_data.get("edge_weight", np.array([]))
    mean_dist = np.mean(edge_weight) if edge_weight.size > 0 else 0

    print(f"  {n_nodes} nodes, {n_edges} edges, mean distance={mean_dist:.2f}A")

    # -------------------------------------------------------------------------
    # Step 4: Create 3D visualization
    # -------------------------------------------------------------------------
    print("\n[4] Creating 3D visualization...")

    output_path = FIGURES_DIR / f"{PDB_ID.lower()}_graph_pocket.html"
    create_3d_pocket_visualization(graph_data, df, output_path)
    print(f"  Saved: {output_path.name}")
    print(f"  Saved: {output_path.with_suffix('.png').name}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"\n{PDB_ID} ({STRUCTURE_INFO['name']}):")
    print(f"  Pocket residues: {len(binding_residues)}")
    print(f"  Graph: {n_nodes} nodes, {n_edges} edges")
    print(f"  Mean contact distance: {mean_dist:.2f}A")
    print(f"\nThis graph format is the input for LAMBDA spectral prediction (Ch. 4)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
