#!/usr/bin/env python3
"""pLM-Enriched Graph Workflow.

Demonstrates cross-processor composition:
- StructureProcessor: Load structure, extract sequence, annotate with GRN
- EmbeddingProcessor: Generate per-residue pLM embeddings
- GraphProcessor: Build binding pocket graph with embedding features

Question: "Can we enrich a binding pocket graph with pLM embeddings?"
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

OUTPUT_DIR = THESIS_DIR / "outputs" / "plm_graph"
FIGURES_DIR = THESIS_DIR / "workflows" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.graph import GraphProcessor
from protos.io.ingest.structure_loader import StructureLoader


# =============================================================================
# Configuration
# =============================================================================
PDB_ID = "1U19"  # Bovine rhodopsin - canonical GPCR structure
CHAIN = "A"
LIGAND = "RET"  # Retinal
POCKET_DISTANCE = 7.0
EDGE_CUTOFF = 6.0
EMBEDDING_MODEL = "ankh_large"  # Same model used by LAMBDA


def main() -> int:
    """Run the pLM-enriched graph workflow."""
    print("=" * 60)
    print("pLM-ENRICHED GRAPH WORKFLOW")
    print("=" * 60)
    print(f"Structure: {PDB_ID} (bovine rhodopsin)")
    print(f"Ligand: {LIGAND}, Pocket: {POCKET_DISTANCE}A")
    print(f"Embedding model: {EMBEDDING_MODEL}")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Load structure using ProtOS
    # -------------------------------------------------------------------------
    print("\n[1] Loading structure...")
    struct_proc = StructureProcessor()
    loader = StructureLoader(processor=struct_proc)

    # Download and register
    loader.download_and_register(PDB_ID, name=PDB_ID.lower())

    df = struct_proc.load_entity(PDB_ID.lower())
    if df is None:
        print("  ERROR: Failed to load structure")
        return 1
    print(f"  Loaded: {len(df)} atoms")

    # -------------------------------------------------------------------------
    # Step 2: Extract chain sequence using ProtOS StructureProcessor
    # -------------------------------------------------------------------------
    print("\n[2] Extracting chain sequence...")
    sequence = struct_proc.get_sequence(PDB_ID.lower(), chain_id=CHAIN)
    if not sequence:
        print("  ERROR: Failed to extract sequence")
        return 1
    print(f"  Chain {CHAIN}: {len(sequence)} residues")
    print(f"  Sequence: {sequence[:50]}...")

    # Register sequence with SequenceProcessor
    seq_proc = SequenceProcessor()
    seq_name = f"{PDB_ID}_{CHAIN}"
    seq_proc.save_entity(seq_name, sequence, metadata={"pdb_id": PDB_ID, "chain": CHAIN})

    # -------------------------------------------------------------------------
    # Step 3: Annotate with GRN using ProtOS StructureProcessor
    # -------------------------------------------------------------------------
    print("\n[3] Annotating with GRN (GPCRdb numbering)...")
    grn_df = struct_proc.annotate_with_grn(PDB_ID.lower(), chains=[CHAIN])

    grn_annotations = {}
    if grn_df is not None and hasattr(grn_df, 'empty') and not grn_df.empty:
        # Build GRN lookup from DataFrame
        if "grn" in grn_df.columns and "residue_number" in grn_df.columns:
            for _, row in grn_df.iterrows():
                res_num = row.get("residue_number")
                grn = row.get("grn", "")
                chain = row.get("chain_id", CHAIN)
                if res_num and grn:
                    grn_annotations[f"{chain}:{res_num}"] = grn

        print(f"  Annotated residues: {len(grn_annotations)}")
        examples = list(grn_annotations.items())[:5]
        for res_key, grn in examples:
            print(f"    {res_key}: {grn}")
    else:
        print("  GRN annotation not available (non-GPCR or service unavailable)")

    # -------------------------------------------------------------------------
    # Step 4: Generate pLM embeddings using ProtOS EmbeddingProcessor
    # -------------------------------------------------------------------------
    print("\n[4] Generating pLM embeddings...")
    # Ankh-large is ~1.5GB, use CPU to avoid OOM on smaller GPUs
    emb_proc = EmbeddingProcessor(model_name=EMBEDDING_MODEL, device="cpu")

    if not emb_proc.dependencies_available:
        print("  WARNING: Embedding dependencies not available")
        print("  Install with: pip install torch transformers")
        embeddings = None
    else:
        # Generate per-residue embeddings
        embeddings = emb_proc.embed_sequences(
            {seq_name: sequence},
            embedding_type="per_residue",
            save_dataset=f"{seq_name}__{EMBEDDING_MODEL}__per_residue",
            register_entities=True,
        )

        if seq_name in embeddings:
            emb = embeddings[seq_name]
            # Remove special tokens (CLS and EOS)
            residue_emb = emb_proc.get_residue_embeddings(emb)
            print(f"  Embedding shape: {residue_emb.shape}")
            print(f"  Embedding dim: {emb_proc.get_embedding_dim()}")

    # -------------------------------------------------------------------------
    # Step 5: Generate binding pocket graph using ProtOS GraphProcessor
    # -------------------------------------------------------------------------
    print("\n[5] Generating binding pocket graph...")
    graph_proc = GraphProcessor(default_cutoff=EDGE_CUTOFF, default_level="residue")

    graph_name = graph_proc.generate_graph(
        PDB_ID.lower(),
        chain=CHAIN,
        near_hetatm=(LIGAND, POCKET_DISTANCE),
        protein_only=True,
        level="residue",
        cutoff=EDGE_CUTOFF,
    )

    graph_data = graph_proc.load_entity(graph_name)
    if not graph_data:
        print("  ERROR: Failed to generate graph")
        return 1

    pos = graph_data.get("node_positions", np.array([]))
    edge_index = graph_data.get("edge_index", np.array([[], []]))
    node_metadata = graph_data.get("node_metadata", [])

    n_nodes = len(pos) if pos.size > 0 else 0
    n_edges = edge_index.shape[1] if edge_index.ndim == 2 else 0
    print(f"  Nodes: {n_nodes}, Edges: {n_edges}")

    # -------------------------------------------------------------------------
    # Step 6: Enrich graph nodes with embeddings
    # -------------------------------------------------------------------------
    print("\n[6] Enriching graph with pLM embeddings...")

    node_embeddings = []
    # grn_annotations already defined in Step 3

    if embeddings and seq_name in embeddings:
        emb = embeddings[seq_name]
        residue_emb = emb_proc.get_residue_embeddings(emb)
        residue_emb_np = residue_emb.numpy() if hasattr(residue_emb, 'numpy') else np.array(residue_emb)

        for meta in node_metadata:
            res_num = meta.get("residue_number", 0)
            res_name = meta.get("residue_name", "UNK")

            # Get embedding for this residue (0-indexed in embedding, 1-indexed in structure)
            emb_idx = res_num - 1
            if 0 <= emb_idx < len(residue_emb_np):
                node_emb = residue_emb_np[emb_idx]
            else:
                node_emb = np.zeros(emb_proc.get_embedding_dim())

            # Get GRN if available
            res_key = f"{CHAIN}:{res_num}"
            grn = grn_annotations.get(res_key, "")

            node_embeddings.append({
                "residue_number": res_num,
                "residue_name": res_name,
                "grn": grn,
                "embedding": node_emb,
                "embedding_norm": float(np.linalg.norm(node_emb)),
            })

        print(f"  Enriched {len(node_embeddings)} nodes with embeddings")

        # Compute embedding similarity between connected nodes
        edge_similarities = []
        for i in range(n_edges):
            src, tgt = int(edge_index[0, i]), int(edge_index[1, i])
            if src < len(node_embeddings) and tgt < len(node_embeddings):
                emb_src = node_embeddings[src]["embedding"]
                emb_tgt = node_embeddings[tgt]["embedding"]
                # Cosine similarity
                norm_src = np.linalg.norm(emb_src)
                norm_tgt = np.linalg.norm(emb_tgt)
                if norm_src > 0 and norm_tgt > 0:
                    sim = np.dot(emb_src, emb_tgt) / (norm_src * norm_tgt)
                else:
                    sim = 0.0
                edge_similarities.append(sim)

        if edge_similarities:
            print(f"  Edge embedding similarity: mean={np.mean(edge_similarities):.3f}, std={np.std(edge_similarities):.3f}")
    else:
        print("  Skipping embedding enrichment (embeddings not available)")

    # -------------------------------------------------------------------------
    # Step 7: Visualize enriched graph
    # -------------------------------------------------------------------------
    print("\n[7] Creating visualization...")
    import plotly.graph_objects as go

    fig = go.Figure()

    if pos.ndim == 2 and pos.shape[0] > 0:
        # Find nodes with edges
        nodes_with_edges = set()
        if edge_index.ndim == 2 and edge_index.shape[1] > 0:
            nodes_with_edges = set(edge_index[0].tolist()) | set(edge_index[1].tolist())

        # Draw edges colored by embedding similarity
        if edge_index.ndim == 2 and edge_index.shape[1] > 0 and edge_similarities:
            for i in range(edge_index.shape[1]):
                src, tgt = int(edge_index[0, i]), int(edge_index[1, i])
                if src < len(pos) and tgt < len(pos):
                    sim = edge_similarities[i] if i < len(edge_similarities) else 0.5
                    # Color: low similarity = red, high similarity = blue
                    r = int(255 * (1 - sim))
                    b = int(255 * sim)
                    color = f"rgb({r}, 100, {b})"
                    fig.add_trace(go.Scatter3d(
                        x=[pos[src, 0], pos[tgt, 0]],
                        y=[pos[src, 1], pos[tgt, 1]],
                        z=[pos[src, 2], pos[tgt, 2]],
                        mode="lines",
                        line=dict(color=color, width=2),
                        showlegend=False,
                        hoverinfo="skip",
                    ))

        # Draw nodes colored by embedding norm
        connected_indices = sorted(nodes_with_edges)
        if connected_indices and node_embeddings:
            connected_pos = pos[connected_indices]
            connected_norms = [node_embeddings[i]["embedding_norm"] for i in connected_indices]
            connected_names = []
            for i in connected_indices:
                meta = node_embeddings[i]
                grn_str = f" ({meta['grn']})" if meta['grn'] else ""
                connected_names.append(f"{meta['residue_name']}{meta['residue_number']}{grn_str}")

            fig.add_trace(go.Scatter3d(
                x=connected_pos[:, 0],
                y=connected_pos[:, 1],
                z=connected_pos[:, 2],
                mode="markers",
                marker=dict(
                    size=8,
                    color=connected_norms,
                    colorscale="Viridis",
                    showscale=True,
                    colorbar=dict(title="Emb Norm"),
                ),
                text=connected_names,
                hovertemplate="%{text}<br>Norm: %{marker.color:.2f}<extra></extra>",
                name="Residues",
            ))

    # Draw ligand with bonds
    if df is not None:
        ligand_df = df[(df["res_name"] == LIGAND) & (~df["element"].isin(["H"]))].copy()
        if len(ligand_df) > 0:
            ligand_coords = ligand_df[["x", "y", "z"]].values

            # Draw bonds
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
                            line=dict(color="#d62728", width=3),
                            showlegend=False,
                        ))

            fig.add_trace(go.Scatter3d(
                x=ligand_df["x"], y=ligand_df["y"], z=ligand_df["z"],
                mode="markers",
                marker=dict(size=5, color="#d62728"),
                name="Retinal",
            ))

    fig.update_layout(
        title=f"pLM-Enriched Binding Pocket Graph ({PDB_ID})",
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
        ),
        width=900, height=750,
    )
    fig.write_html(str(FIGURES_DIR / "plm_enriched_graph.html"))
    print(f"  Saved: {FIGURES_DIR / 'plm_enriched_graph.html'}")

    # -------------------------------------------------------------------------
    # Step 8: Save enriched graph data
    # -------------------------------------------------------------------------
    print("\n[8] Saving enriched graph data...")
    import json

    # Prepare serializable output
    output_data = {
        "pdb_id": PDB_ID,
        "chain": CHAIN,
        "ligand": LIGAND,
        "n_nodes": int(n_nodes),
        "n_edges": int(n_edges),
        "embedding_model": EMBEDDING_MODEL,
        "nodes": [
            {
                "residue_number": int(ne["residue_number"]),
                "residue_name": ne["residue_name"],
                "grn": ne["grn"],
                "embedding_norm": float(ne["embedding_norm"]),
            }
            for ne in node_embeddings
        ] if node_embeddings else [],
        "edge_similarities": {
            "mean": float(np.mean(edge_similarities)) if edge_similarities else None,
            "std": float(np.std(edge_similarities)) if edge_similarities else None,
        },
    }

    with open(OUTPUT_DIR / "plm_enriched_graph.json", "w") as f:
        json.dump(output_data, f, indent=2)
    print(f"  Saved: {OUTPUT_DIR / 'plm_enriched_graph.json'}")

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)
    print(f"Outputs: {OUTPUT_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
