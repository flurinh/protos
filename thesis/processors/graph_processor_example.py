#!/usr/bin/env python3
"""GraphProcessor Example: Binding pocket graph generation for Type I and Type II opsins.

Demonstrates ProtOS capabilities:
- StructureLoader: Download structures from PDB
- StructureProcessor: Analyze ligand interactions
- GraphProcessor: Generate residue contact graphs

Question: "How do binding pocket graphs compare between Type I and Type II opsins?"

KEY INSIGHT: Despite completely different protein folds, both Type I (microbial)
and Type II (animal) opsins create similar graph representations around retinal.
This demonstrates that LAMBDA's graph-based approach is fold-agnostic.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import yaml

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
# Load Color Scheme
# =============================================================================
with open(THESIS_DIR / "colorscales.yaml") as f:
    COLORS = yaml.safe_load(f)


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "graph"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Two structures: Type I (microbial) and Type II (animal)
STRUCTURES = {
    "1C3W": {
        "name": "Bacteriorhodopsin",
        "type": "Type I (Microbial)",
        "chain": "A",
        "ligand": "RET",
        "color": COLORS["structures"]["1C3W"]["hex"],
        "function": "Proton pump",
    },
    "1U19": {
        "name": "Bovine Rhodopsin",
        "type": "Type II (Animal/GPCR)",
        "chain": "A",
        "ligand": "RET",
        "color": COLORS["structures"]["1U19"]["hex"],
        "function": "Dim light vision",
    },
}

# Ligand color from colorscales
LIGAND_COLOR = COLORS["molecules"]["retinal"]["hex"]
EDGE_COLOR = COLORS["graphs"]["edges"]["hex"]

POCKET_DISTANCE = 7.0
EDGE_CUTOFF = 6.0


def create_3d_pocket_visualization(
    graph_data: dict,
    struct_df,
    pdb_id: str,
    info: dict,
    output_path: Path,
) -> None:
    """Create 3D visualization for a single binding pocket."""
    import plotly.graph_objects as go

    fig = go.Figure()
    color = info["color"]

    if graph_data:
        pos = graph_data.get("node_positions", np.array([]))
        edge_index = graph_data.get("edge_index", np.array([[], []]))
        node_metadata = graph_data.get("node_metadata", [])
        node_names = [f"{m['residue_name']}{m['residue_number']}" for m in node_metadata] if node_metadata else [f"R{i}" for i in range(len(pos))]

        if pos.ndim == 2 and pos.shape[0] > 0:
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
                            mode="lines", line=dict(color=EDGE_COLOR, width=2), showlegend=False,
                        ))

            # Draw nodes
            connected_indices = sorted(nodes_with_edges)
            if connected_indices:
                connected_pos = pos[connected_indices]
                connected_names = [node_names[i] for i in connected_indices]
                fig.add_trace(go.Scatter3d(
                    x=connected_pos[:, 0], y=connected_pos[:, 1], z=connected_pos[:, 2],
                    mode="markers",
                    marker=dict(size=6, color=color),
                    text=connected_names,
                    hovertemplate="%{text}<extra></extra>",
                    name="Residues",
                ))

    # Draw ligand (only from specified chain)
    ligand = info["ligand"]
    chain = info["chain"]
    if struct_df is not None:
        ligand_df = struct_df[
            (struct_df["res_name"] == ligand) &
            (struct_df["auth_chain_id"] == chain) &
            (~struct_df["element"].isin(["H"]))
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
                            mode="lines", line=dict(color=LIGAND_COLOR, width=3),
                            showlegend=False,
                        ))
            fig.add_trace(go.Scatter3d(
                x=ligand_df["x"], y=ligand_df["y"], z=ligand_df["z"],
                mode="markers", marker=dict(size=5, color=LIGAND_COLOR),
                name="Retinal",
            ))

    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            aspectmode="data",
        ),
        width=800, height=700,
        margin=dict(t=30, b=30),
    )
    fig.write_html(str(output_path))


def main() -> int:
    """Run the GraphProcessor example."""
    print("=" * 70)
    print("GRAPH PROCESSOR EXAMPLE")
    print("Type I vs Type II Opsin Binding Pocket Comparison")
    print("=" * 70)
    print("\nKEY QUESTION: Are binding pocket graphs comparable despite different folds?")
    print(f"\nStructures to compare:")
    for pdb_id, info in STRUCTURES.items():
        print(f"  {pdb_id}: {info['name']} - {info['type']}")
    print(f"\nParameters: Pocket={POCKET_DISTANCE}Å, Edge cutoff={EDGE_CUTOFF}Å")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Load both structures (download if needed)
    # -------------------------------------------------------------------------
    print("\n[1] Loading structures...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    structure_dfs = {}
    for pdb_id, info in STRUCTURES.items():
        # Try loading by PDB ID first (most common case)
        df = struct_proc.load_entity(pdb_id)

        if df is None:
            # Try downloading if not found
            structure_name = f"{info['name'].lower().replace(' ', '_')}_{pdb_id}"
            try:
                loader.download_and_register(f"pdb:{pdb_id}", name=pdb_id)
                df = struct_proc.load_entity(pdb_id)
            except Exception as e:
                print(f"  {pdb_id}: Failed to load - {e}")
                continue

        if df is not None:
            structure_dfs[pdb_id] = (pdb_id, df)
            print(f"  {pdb_id}: {len(df)} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Analyze binding pockets for both structures
    # -------------------------------------------------------------------------
    print("\n[2] Analyzing binding pockets...")

    pocket_stats = {}
    for pdb_id, info in STRUCTURES.items():
        if pdb_id not in structure_dfs:
            continue

        structure_name, df = structure_dfs[pdb_id]
        chain = info["chain"]
        ligand = info["ligand"]

        interactions = struct_proc.get_ligand_interactions(
            structure_name,
            ligand_id=ligand,
            chain_id=chain,
            cutoff=POCKET_DISTANCE,
        )

        binding_residues = interactions.get("binding_residues", [])
        pocket_stats[pdb_id] = {"n_residues": len(binding_residues)}
        print(f"  {pdb_id} ({info['type']}): {len(binding_residues)} binding residues")

    # -------------------------------------------------------------------------
    # Step 3: Generate graphs for both structures
    # -------------------------------------------------------------------------
    print("\n[3] Generating binding pocket graphs...")
    graph_proc = GraphProcessor(default_cutoff=EDGE_CUTOFF, default_level="residue")

    graph_data_dict = {}
    for pdb_id, info in STRUCTURES.items():
        if pdb_id not in structure_dfs:
            continue

        structure_name, df = structure_dfs[pdb_id]
        chain = info["chain"]
        ligand = info["ligand"]

        graph_name = graph_proc.generate_graph(
            structure_name,
            chain=chain,
            near_hetatm=(ligand, POCKET_DISTANCE),
            protein_only=True,
            level="residue",
            cutoff=EDGE_CUTOFF,
        )

        graph_data = graph_proc.load_entity(graph_name)
        if graph_data:
            graph_data_dict[pdb_id] = graph_data
            pos = graph_data.get("node_positions", np.array([]))
            n_nodes = len(pos) if pos.size > 0 else 0
            edge_index = graph_data.get("edge_index", np.array([[], []]))
            n_edges = edge_index.shape[1] if edge_index.ndim == 2 else 0
            edge_weight = graph_data.get("edge_weight", np.array([]))
            mean_dist = np.mean(edge_weight) if edge_weight.size > 0 else 0

            pocket_stats[pdb_id].update({
                "n_nodes": n_nodes,
                "n_edges": n_edges,
                "mean_edge_dist": mean_dist,
                "edge_density": n_edges / (n_nodes * (n_nodes - 1) / 2) if n_nodes > 1 else 0,
            })
            print(f"  {pdb_id}: {n_nodes} nodes, {n_edges} edges, mean dist={mean_dist:.2f}Å")

    # -------------------------------------------------------------------------
    # Step 4: Create individual 3D visualizations
    # -------------------------------------------------------------------------
    print("\n[4] Creating 3D visualizations...")

    for pdb_id, info in STRUCTURES.items():
        if pdb_id not in graph_data_dict:
            continue

        structure_name, df = structure_dfs[pdb_id]
        output_path = FIGURES_DIR / f"graph_binding_pocket_3d_{pdb_id}.html"

        create_3d_pocket_visualization(
            graph_data_dict[pdb_id],
            df,
            pdb_id,
            info,
            output_path,
        )
        print(f"  Saved: {output_path.name}")

    # -------------------------------------------------------------------------
    # Step 5: Create comparison statistics figure (PNG)
    # -------------------------------------------------------------------------
    print("\n[5] Creating comparison statistics figure...")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))

    pdb_ids = list(STRUCTURES.keys())
    colors = [STRUCTURES[p]["color"] for p in pdb_ids]
    labels = [f"{p}" for p in pdb_ids]

    # Plot 1: Node count (no subplot title)
    nodes = [pocket_stats[p]["n_nodes"] for p in pdb_ids]
    axes[0].bar(labels, nodes, color=colors, edgecolor="black", linewidth=1.2)
    axes[0].set_ylabel("Residues", fontsize=10)
    for i, v in enumerate(nodes):
        axes[0].text(i, v + 1, str(v), ha="center", fontsize=10, fontweight="bold")
    axes[0].tick_params(labelsize=9)

    # Plot 2: Edge count (no subplot title)
    edges = [pocket_stats[p]["n_edges"] for p in pdb_ids]
    axes[1].bar(labels, edges, color=colors, edgecolor="black", linewidth=1.2)
    axes[1].set_ylabel("Contacts", fontsize=10)
    for i, v in enumerate(edges):
        axes[1].text(i, v + 3, str(v), ha="center", fontsize=10, fontweight="bold")
    axes[1].tick_params(labelsize=9)

    # Plot 3: Edge density (no subplot title)
    density = [pocket_stats[p]["edge_density"] * 100 for p in pdb_ids]
    axes[2].bar(labels, density, color=colors, edgecolor="black", linewidth=1.2)
    axes[2].set_ylabel("Density (%)", fontsize=10)
    for i, v in enumerate(density):
        axes[2].text(i, v + 0.5, f"{v:.1f}%", ha="center", fontsize=10, fontweight="bold")
    axes[2].tick_params(labelsize=9)

    # No suptitle - clean for thesis
    plt.tight_layout()

    comparison_path = FIGURES_DIR / "graph_comparison.png"
    plt.savefig(comparison_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Saved: {comparison_path.name}")

    # -------------------------------------------------------------------------
    # Step 6: Print summary table
    # -------------------------------------------------------------------------
    print("\n[6] Summary Comparison:")
    print("-" * 70)
    print(f"{'Structure':<12} {'Type':<20} {'Nodes':<8} {'Edges':<8} {'Density':<10} {'Mean Dist':<10}")
    print("-" * 70)
    for pdb_id in pdb_ids:
        info = STRUCTURES[pdb_id]
        stats = pocket_stats[pdb_id]
        print(f"{pdb_id:<12} {info['type']:<20} {stats['n_nodes']:<8} {stats['n_edges']:<8} "
              f"{stats['edge_density']*100:.1f}%{'':<5} {stats['mean_edge_dist']:.2f}Å")
    print("-" * 70)

    print("\n" + "=" * 70)
    print("KEY INSIGHT")
    print("=" * 70)
    print("Despite completely different protein folds:")
    print(f"  • Type I (1C3W): {pocket_stats['1C3W']['n_nodes']} nodes, {pocket_stats['1C3W']['n_edges']} edges")
    print(f"  • Type II (1U19): {pocket_stats['1U19']['n_nodes']} nodes, {pocket_stats['1U19']['n_edges']} edges")
    print("\n→ Both create COMPARABLE graph representations around retinal")
    print("→ This enables LAMBDA's fold-agnostic binding pocket analysis")

    print("\n" + "=" * 70)
    print("COMPLETE")
    print("=" * 70)
    print(f"Outputs: {OUTPUT_DIR}")
    print(f"Figures:")
    print(f"  • graph_binding_pocket_3d_1C3W.html (Type I)")
    print(f"  • graph_binding_pocket_3d_1U19.html (Type II)")
    print(f"  • graph_comparison.png (statistics)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
