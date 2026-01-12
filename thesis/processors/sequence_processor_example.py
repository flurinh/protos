#!/usr/bin/env python3
"""SequenceProcessor Example: Building an opsin diversity dataset.

Demonstrates ProtOS capabilities:
- SequenceLoader: Load sequences from UniProt
- NCBILoader: Run BLAST homology searches and batch downloads
- Built-in sequence registration and dataset management

Question: "What is the sequence diversity across cone opsin types?"

This example builds a comprehensive opsin dataset by searching for homologs
of short-wave (SW/blue) and long-wave (LW/red) cone opsins. We use these two
types because SW and LW are spectrally distinct (~420nm vs ~560nm), while
MW (~530nm) overlaps significantly with LW. The resulting 200-sequence dataset
feeds into the EmbeddingProcessor example for comparative analysis.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.io.ingest.sequence_loader import SequenceLoader
from protos.io.ingest.ncbi_loader import NCBILoader


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "sequence"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

BLAST_DATABASE = "swissprot"
HITS_PER_QUERY = 100  # Get 100 hits per opsin type -> 200 total (SW + LW)

# Query opsins: Human cone opsins representing SW and LW sensitivity
# Note: MW omitted because it overlaps heavily with LW (96% sequence identity)
QUERY_OPSINS = {
    "OPN1SW": {
        "uniprot": "P03999",
        "name": "OPSB_HUMAN",
        "type": "short_wave",
        "lambda_max": 420,
        "description": "Blue cone opsin (S-cone)",
    },
    "OPN1LW": {
        "uniprot": "P04000",
        "name": "OPSR_HUMAN",
        "type": "long_wave",
        "lambda_max": 560,
        "description": "Red cone opsin (L-cone)",
    },
}

DATASET_NAME = "cone_opsin_diversity"


def main() -> int:
    """Run the SequenceProcessor example."""
    print("=" * 70)
    print("SEQUENCE PROCESSOR EXAMPLE")
    print("Building Cone Opsin Diversity Dataset")
    print("=" * 70)
    print(f"\nQuery opsins:")
    for gene, info in QUERY_OPSINS.items():
        print(f"  {gene}: {info['description']} (λmax={info['lambda_max']} nm)")
    print(f"\nDatabase: {BLAST_DATABASE}")
    print(f"Hits per query: {HITS_PER_QUERY}")

    # Initialize ProtOS
    protos.set_data_path(str(REPO_ROOT / "data"))

    # -------------------------------------------------------------------------
    # Step 1: Download query sequences using ProtOS SequenceLoader
    # -------------------------------------------------------------------------
    print("\n[1] Downloading query sequences...")
    seq_proc = SequenceProcessor()
    loader = SequenceLoader(processor=seq_proc)

    query_sequences = {}
    for gene, info in QUERY_OPSINS.items():
        loader.download_and_register(
            f"uniprot:{info['uniprot']}",
            name=info["name"],
            materialize_entities=True,
        )
        seq = seq_proc.load_entity(info["name"])
        if seq:
            query_sequences[gene] = seq
            print(f"  {info['name']}: {len(seq)} aa ({info['type']})")

    # -------------------------------------------------------------------------
    # Step 2: Run BLAST for each opsin type using ProtOS NCBILoader
    # -------------------------------------------------------------------------
    print("\n[2] Running NCBI BLAST searches...")
    ncbi_loader = NCBILoader(processor=seq_proc)

    all_hits: Dict[str, pd.DataFrame] = {}

    for gene, info in QUERY_OPSINS.items():
        if gene not in query_sequences:
            continue

        print(f"\n  Searching for {gene} ({info['type']}) homologs...")

        result = ncbi_loader.blast_search(
            sequence=query_sequences[gene],
            query_id=info["name"],
            program="blastp",
            database=BLAST_DATABASE,
            hitlist_size=HITS_PER_QUERY,
            expect=0.001,
        )

        if result:
            hits_df = ncbi_loader.to_dataframe(result)
            hits_df["query_gene"] = gene
            hits_df["opsin_type"] = info["type"]
            hits_df["query_lambda_max"] = info["lambda_max"]
            all_hits[gene] = hits_df
            print(f"    Found {len(hits_df)} hits")
            print(f"    Identity range: {hits_df['identity_percent'].min():.1f}% - {hits_df['identity_percent'].max():.1f}%")
        else:
            print(f"    BLAST search failed for {gene}")

    # Combine all hits
    if not all_hits:
        print("\nNo BLAST results obtained!")
        return 1

    combined_df = pd.concat(all_hits.values(), ignore_index=True)
    print(f"\n  Total hits: {len(combined_df)}")

    # Remove duplicates (same sequence found by multiple queries)
    combined_df = combined_df.drop_duplicates(subset=["hit_accession"])
    print(f"  Unique sequences: {len(combined_df)}")

    # Save combined results
    combined_df.to_csv(OUTPUT_DIR / "blast_all_opsins.csv", index=False)

    # -------------------------------------------------------------------------
    # Step 3: Fetch homolog sequences using ProtOS NCBILoader
    # -------------------------------------------------------------------------
    print("\n[3] Fetching homolog sequences...")

    # Get all unique accessions
    accessions = combined_df["hit_accession"].tolist()

    # ProtOS batch download with dataset registration
    result = ncbi_loader.download_batch(
        accessions=accessions,
        dataset_name=DATASET_NAME,
        metadata={
            "source": "ncbi_blast",
            "queries": list(QUERY_OPSINS.keys()),
            "description": "Cone opsin diversity dataset for embedding analysis",
        },
    )

    success_count = len(result.get("success", []))
    fail_count = len(result.get("failed", []))
    print(f"  Downloaded: {success_count}/{len(accessions)} sequences")
    if fail_count:
        print(f"  Failed: {fail_count}")

    # -------------------------------------------------------------------------
    # Step 4: Annotate sequences with opsin type
    # -------------------------------------------------------------------------
    print("\n[4] Annotating sequences with opsin type...")

    # Create annotation mapping from BLAST results
    seq_annotations = {}
    for _, row in combined_df.iterrows():
        seq_annotations[row["hit_accession"]] = {
            "opsin_type": row["opsin_type"],
            "query_gene": row["query_gene"],
            "identity_to_query": row["identity_percent"],
            "query_lambda_max": row["query_lambda_max"],
        }

    # Save annotation mapping for EmbeddingProcessor
    import json
    with open(OUTPUT_DIR / "opsin_annotations.json", "w") as f:
        json.dump(seq_annotations, f, indent=2)
    print(f"  Saved annotations for {len(seq_annotations)} sequences")

    # -------------------------------------------------------------------------
    # Step 5: Visualize dataset composition
    # -------------------------------------------------------------------------
    print("\n[5] Creating visualizations...")
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    # Type colors matching spectral sensitivity
    type_colors = {
        "short_wave": "#4169E1",   # Blue
        "long_wave": "#DC143C",    # Red
    }

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Sequences by Opsin Type", "Identity Distribution by Type"),
        specs=[[{"type": "pie"}, {"type": "box"}]],
    )

    # Pie chart of opsin types
    type_counts = combined_df["opsin_type"].value_counts()
    fig.add_trace(go.Pie(
        labels=[t.replace("_", " ").title() for t in type_counts.index],
        values=type_counts.values,
        marker_colors=[type_colors.get(t, "#888") for t in type_counts.index],
        hole=0.4,
        textinfo="label+value",
    ), row=1, col=1)

    # Box plot of identity by type
    for opsin_type in ["short_wave", "long_wave"]:
        type_df = combined_df[combined_df["opsin_type"] == opsin_type]
        if not type_df.empty:
            fig.add_trace(go.Box(
                y=type_df["identity_percent"],
                name=opsin_type.replace("_", " ").title(),
                marker_color=type_colors.get(opsin_type, "#888"),
                boxpoints="outliers",
            ), row=1, col=2)

    fig.update_layout(
        title="Cone Opsin Diversity Dataset",
        height=450,
        width=950,
        showlegend=False,
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    fig.update_yaxes(title_text="Identity to Query (%)", row=1, col=2)
    fig.update_yaxes(showgrid=False)
    fig.write_image(str(FIGURES_DIR / "sequence_opsin_diversity.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'sequence_opsin_diversity.png'}")

    # Identity scatter plot
    fig = go.Figure()

    for opsin_type in ["short_wave", "long_wave"]:
        type_df = combined_df[combined_df["opsin_type"] == opsin_type]
        if not type_df.empty:
            fig.add_trace(go.Scatter(
                x=type_df["identity_percent"],
                y=type_df["evalue"],
                mode="markers",
                name=opsin_type.replace("_", " ").title(),
                marker=dict(
                    size=8,
                    color=type_colors.get(opsin_type, "#888"),
                    opacity=0.7,
                ),
                text=type_df["hit_accession"],
                hovertemplate="%{text}<br>Identity: %{x:.1f}%<br>E-value: %{y:.2e}<extra></extra>",
            ))

    fig.update_layout(
        title="BLAST Hits: Identity vs E-value by Opsin Type",
        xaxis_title="Sequence Identity (%)",
        yaxis_title="E-value",
        height=500,
        width=800,
        paper_bgcolor="white",
        plot_bgcolor="white",
        legend=dict(x=0.02, y=0.98),
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(type="log", showgrid=False)
    fig.write_image(str(FIGURES_DIR / "sequence_blast_scatter.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'sequence_blast_scatter.png'}")

    # -------------------------------------------------------------------------
    # Step 6: All-vs-all similarity and phylogenetic tree
    # -------------------------------------------------------------------------
    print("\n[6] Computing all-vs-all sequence similarity (MMseqs2)...")

    from protos.processing.sequence.seq_alignment import mmseqs2_align2

    # Load sequences from the dataset we just created
    downloaded_seqs = seq_proc.load_dataset(DATASET_NAME)

    if not downloaded_seqs:
        print("  Failed to load dataset for alignment")
    elif len(downloaded_seqs) < 3:
        print(f"  Loaded {len(downloaded_seqs)} sequences - not enough for phylogenetic analysis")
    else:
        print(f"  Loaded {len(downloaded_seqs)} sequences for alignment")
        # Run MMseqs2 all-vs-all alignment
        temp_dir = str(OUTPUT_DIR / "mmseqs_temp")
        alignment_df = mmseqs2_align2(downloaded_seqs, downloaded_seqs, temp_folder=temp_dir)

        if alignment_df is not None and not alignment_df.empty:
            print(f"  Computed {len(alignment_df)} pairwise alignments")

            # Build distance matrix from sequence identities
            # Distance = 1 - identity (identity is 0-1)
            seq_names = list(downloaded_seqs.keys())
            n_seqs = len(seq_names)
            name_to_idx = {name: i for i, name in enumerate(seq_names)}

            # Initialize with max distance (for missing pairs)
            dist_matrix = np.ones((n_seqs, n_seqs))
            np.fill_diagonal(dist_matrix, 0)

            # Fill in distances from alignment results
            for _, row in alignment_df.iterrows():
                q_id = row["query_id"]
                t_id = row["target_id"]
                identity = row["sequence_identity"] / 100.0  # Convert to 0-1

                if q_id in name_to_idx and t_id in name_to_idx:
                    i, j = name_to_idx[q_id], name_to_idx[t_id]
                    dist = 1.0 - identity
                    dist_matrix[i, j] = min(dist_matrix[i, j], dist)
                    dist_matrix[j, i] = min(dist_matrix[j, i], dist)

            # Save distance matrix
            dist_df = pd.DataFrame(dist_matrix, index=seq_names, columns=seq_names)
            dist_df.to_csv(OUTPUT_DIR / "sequence_distance_matrix.csv")
            print(f"  Saved distance matrix: {OUTPUT_DIR / 'sequence_distance_matrix.csv'}")

            # Hierarchical clustering for phylogenetic tree
            print("\n  Building phylogenetic tree...")
            from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
            from scipy.spatial.distance import squareform

            # Convert to condensed distance matrix
            condensed_dist = squareform(dist_matrix)

            # Hierarchical clustering (UPGMA-like)
            linkage_matrix = linkage(condensed_dist, method="average")

            # Helper to extract accession from sequence name
            def get_accession(name: str) -> str:
                if "|" in name:
                    parts = name.split("|")
                    if len(parts) >= 2:
                        return parts[1].split(".")[0]
                return name

            # Get opsin types for coloring (using accession lookup)
            def get_opsin_type(name):
                acc = get_accession(name)
                return seq_annotations.get(acc, {}).get("opsin_type", "unknown")

            # Create circular dendrogram using plotly
            # First, create a polar/radial tree visualization
            fig_tree = go.Figure()

            # Get leaf order from dendrogram (for coloring)
            from scipy.cluster.hierarchy import leaves_list
            leaf_order = leaves_list(linkage_matrix)

            # Count sequences by type
            type_counts_tree = {"short_wave": 0, "long_wave": 0}
            for idx in leaf_order:
                t = get_opsin_type(seq_names[idx])
                if t in type_counts_tree:
                    type_counts_tree[t] += 1

            # Create a radial layout for the tree
            n_leaves = len(leaf_order)
            angles = np.linspace(0, 2 * np.pi, n_leaves, endpoint=False)

            # Plot leaves as colored points in a circle
            for opsin_type, color in type_colors.items():
                indices = [i for i, idx in enumerate(leaf_order) if get_opsin_type(seq_names[idx]) == opsin_type]
                if indices:
                    theta = angles[indices]
                    r = np.ones(len(indices))
                    x = r * np.cos(theta)
                    y = r * np.sin(theta)
                    fig_tree.add_trace(go.Scatter(
                        x=x, y=y,
                        mode="markers",
                        name=f"{opsin_type.replace('_', ' ').title()} (n={type_counts_tree.get(opsin_type, 0)})",
                        marker=dict(size=8, color=color, opacity=0.8),
                        hoverinfo="skip",
                    ))

            # Draw arcs connecting clusters at different levels
            # Sample a few key linkage levels to show structure
            max_dist = linkage_matrix[:, 2].max()
            for level_frac in [0.3, 0.5, 0.7, 0.9]:
                level_dist = max_dist * level_frac
                clusters = fcluster(linkage_matrix, t=level_dist, criterion="distance")

                # Draw circle at this distance level
                circle_r = 1 - level_frac * 0.6
                theta_circle = np.linspace(0, 2 * np.pi, 100)
                fig_tree.add_trace(go.Scatter(
                    x=circle_r * np.cos(theta_circle),
                    y=circle_r * np.sin(theta_circle),
                    mode="lines",
                    line=dict(color="rgba(150,150,150,0.3)", width=1),
                    showlegend=False,
                    hoverinfo="skip",
                ))

            fig_tree.update_layout(
                title=f"Cone Opsin Clustering (n={n_seqs})",
                width=700,
                height=700,
                showlegend=True,
                legend=dict(x=0.02, y=0.98),
                xaxis=dict(showgrid=False, showticklabels=False, zeroline=False, scaleanchor="y"),
                yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
                plot_bgcolor="white",
                paper_bgcolor="white",
            )
            fig_tree.write_image(str(FIGURES_DIR / "sequence_phylogenetic_tree.png"), scale=2)
            print(f"  Saved: {FIGURES_DIR / 'sequence_phylogenetic_tree.png'}")

        else:
            print("  MMseqs2 alignment failed - skipping phylogenetic tree")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPLETE - Dataset Ready for Embedding Analysis")
    print("=" * 70)
    print(f"\nDataset: {DATASET_NAME}")
    print(f"Total sequences: {len(combined_df)}")
    print(f"\nBy opsin type:")
    for opsin_type, count in type_counts.items():
        print(f"  {opsin_type.replace('_', ' ').title()}: {count}")
    print(f"\nOutputs: {OUTPUT_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    print(f"\n→ This dataset feeds into the EmbeddingProcessor example")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
