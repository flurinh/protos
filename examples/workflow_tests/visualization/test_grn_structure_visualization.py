#!/usr/bin/env python3
"""Visualize GRN-annotated structures with interactive 3D alignment and GRN position slider."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.processing.structure import StructureProcessor


def ensure_data_root() -> Path:
    """Set up data root directory."""
    data_root = REPO_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def load_grn_annotated_structures() -> Tuple[StructureProcessor, Dict[str, pd.DataFrame]]:
    """Load GRN-annotated structures from saved pickles."""
    struct_proc = StructureProcessor()
    annotated_dir = Path("data/structure/annotated")
    
    structures = {}
    
    # Load all GRN-annotated structures
    for pkl_file in annotated_dir.glob("grn_annotated_*.pkl"):
        struct_id = pkl_file.stem.replace("grn_annotated_", "")
        try:
            df = pd.read_pickle(pkl_file)
            structures[struct_id] = df
            struct_proc.frames[struct_id] = df
            print(f"✓ Loaded GRN-annotated structure: {struct_id}")
        except Exception as e:
            print(f"✗ Failed to load {pkl_file}: {e}")
    
    return struct_proc, structures


def get_gpcr_chains(structures: Dict[str, pd.DataFrame]) -> Dict[str, List[str]]:
    """Identify GPCR chains (those with GRN annotations)."""
    gpcr_chains = {}
    
    for struct_id, df in structures.items():
        if 'grn' not in df.columns:
            continue
            
        # Find chains with GRN annotations
        chains_with_grn = df[df['grn'] != '-']['auth_chain_id'].unique()
        if len(chains_with_grn) > 0:
            gpcr_chains[struct_id] = list(chains_with_grn)
            print(f"  {struct_id}: chains {', '.join(chains_with_grn)} have GRN annotations")
    
    return gpcr_chains


def extract_ca_with_grn(df: pd.DataFrame, chain_id: str) -> pd.DataFrame:
    """Extract CA atoms with GRN annotations for a specific chain."""
    # Ensure grn column exists
    if 'grn' not in df.columns:
        return pd.DataFrame()
    
    ca_atoms = df[
        (df['auth_chain_id'] == chain_id) & 
        (df['atom_name'] == 'CA') &
        (df['grn'].notna()) &
        (df['grn'] != '-') &
        (df['grn'] != '')
    ].copy()
    
    # Sort by residue number
    ca_atoms = ca_atoms.sort_values(['auth_seq_id', 'insertion']).reset_index(drop=True)
    
    return ca_atoms


def align_structures_by_grn(
    struct_proc: StructureProcessor,
    structures: Dict[str, pd.DataFrame],
    gpcr_chains: Dict[str, List[str]],
    reference_struct: str = "5d5a"
) -> Dict[str, pd.DataFrame]:
    """Align all structures to reference, extracting aligned GPCR chains."""
    aligned_chains = {}
    
    # Get reference chain
    ref_chains = gpcr_chains.get(reference_struct, [])
    if not ref_chains:
        print(f"No GPCR chains found in reference structure {reference_struct}")
        return aligned_chains
    
    ref_chain = ref_chains[0]  # Use first GPCR chain as reference
    print(f"\nUsing {reference_struct} chain {ref_chain} as reference")
    
    # Add reference structure
    ref_ca = extract_ca_with_grn(structures[reference_struct], ref_chain)
    aligned_chains[f"{reference_struct}_{ref_chain}"] = ref_ca
    
    # Align other structures
    mobile_ids = [sid for sid in structures.keys() if sid != reference_struct]
    
    if mobile_ids:
        print("\nPerforming structural alignments...")
        rmsd_map, results = struct_proc.align_structures(
            structure_ids=mobile_ids,
            reference_id=reference_struct,
            apply_transform=True,
        )
        
        # Print RMSD summary
        print("\nAlignment RMSD (Å):")
        for target, refs in rmsd_map.items():
            for ref, rmsd in refs.items():
                if not np.isnan(rmsd):
                    print(f"  {target} -> {ref}: {rmsd:.3f}")
    
    # Extract aligned CA atoms for each structure
    for struct_id in mobile_ids:
        if struct_id not in gpcr_chains:
            continue
            
        # Get aligned structure (with suffix)
        aligned_id = f"{struct_id}_aligned"
        aligned_df = struct_proc.frames.get(aligned_id)
        
        if aligned_df is None:
            # Use original if alignment failed
            aligned_df = structures[struct_id]
            print(f"⚠ Using original coordinates for {struct_id}")
        
        # Extract CA atoms for each GPCR chain
        for chain_id in gpcr_chains[struct_id]:
            ca_atoms = extract_ca_with_grn(aligned_df, chain_id)
            if len(ca_atoms) > 0:
                aligned_chains[f"{struct_id}_{chain_id}"] = ca_atoms
    
    return aligned_chains


def get_unique_grn_positions(aligned_chains: Dict[str, pd.DataFrame]) -> List[str]:
    """Get all unique GRN positions across all structures."""
    all_grns = set()
    
    for chain_data in aligned_chains.values():
        grns = chain_data['grn'].unique()
        # Filter out invalid GRN values
        valid_grns = [g for g in grns if g and g != '-' and g != '' and isinstance(g, str)]
        all_grns.update(valid_grns)
    
    # Sort GRN positions properly
    from protos.processing.grn.grn_utils import parse_grn_str2float, sort_grns
    
    grn_floats = []
    for grn in all_grns:
        if not grn or grn == '-' or grn == '':
            continue
        try:
            grn_floats.append(parse_grn_str2float(grn))
        except Exception as e:
            print(f"Warning: Could not parse GRN '{grn}': {e}")
            continue
    
    sorted_floats = sort_grns(grn_floats)
    sorted_grns = []
    
    for gf in sorted_floats:
        # Find original string representation
        for grn in all_grns:
            try:
                if abs(parse_grn_str2float(grn) - gf) < 0.0001:
                    sorted_grns.append(grn)
                    break
            except:
                continue
    
    return sorted_grns


def create_grn_visualization(aligned_chains: Dict[str, pd.DataFrame]) -> go.Figure:
    """Create interactive 3D visualization with GRN position slider."""
    
    # Get all unique GRN positions
    grn_positions = get_unique_grn_positions(aligned_chains)
    if not grn_positions:
        print("No GRN positions found!")
        return None
    
    print(f"\nFound {len(grn_positions)} unique GRN positions")
    
    # Create figure with slider
    fig = go.Figure()
    
    # Color palette for different structures
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
    
    # Create traces for each structure/chain
    for idx, (chain_name, ca_atoms) in enumerate(aligned_chains.items()):
        color = colors[idx % len(colors)]
        
        # Create scatter plot for CA atoms
        fig.add_trace(
            go.Scatter3d(
                x=ca_atoms['x'],
                y=ca_atoms['y'],
                z=ca_atoms['z'],
                mode='markers+lines',
                marker=dict(
                    size=5,
                    color=color,
                    opacity=0.6,
                ),
                line=dict(
                    color=color,
                    width=2,
                ),
                name=chain_name,
                text=[f"{row['res_name']}{row['auth_seq_id']} (GRN: {row['grn']})" 
                      for _, row in ca_atoms.iterrows()],
                hoverinfo='text+name',
                visible=True,
                opacity=0.8,
            )
        )
    
    # Create traces for highlighted positions (one per structure, initially hidden)
    num_structures = len(aligned_chains)
    
    for grn_idx, grn_pos in enumerate(grn_positions):
        for struct_idx, (chain_name, ca_atoms) in enumerate(aligned_chains.items()):
            # Find atoms with this GRN position
            grn_atoms = ca_atoms[ca_atoms['grn'] == grn_pos]
            
            if len(grn_atoms) > 0:
                # Add highlighted position
                fig.add_trace(
                    go.Scatter3d(
                        x=grn_atoms['x'],
                        y=grn_atoms['y'],
                        z=grn_atoms['z'],
                        mode='markers',
                        marker=dict(
                            size=15,
                            color='red',
                            symbol='diamond',
                            line=dict(
                                color='darkred',
                                width=2
                            ),
                        ),
                        name=f"{chain_name} - GRN {grn_pos}",
                        text=[f"{row['res_name']}{row['auth_seq_id']} (GRN: {grn_pos})" 
                              for _, row in grn_atoms.iterrows()],
                        hoverinfo='text',
                        visible=False,  # Initially hidden
                        showlegend=False,
                    )
                )
            else:
                # Add empty trace to maintain indexing
                fig.add_trace(
                    go.Scatter3d(
                        x=[None],
                        y=[None],
                        z=[None],
                        mode='markers',
                        visible=False,
                        showlegend=False,
                    )
                )
    
    # Create slider steps
    steps = []
    
    # First step: show all structures without highlighting
    step = dict(
        method="update",
        args=[{"visible": [True] * num_structures + [False] * (len(fig.data) - num_structures)}],
        label="None"
    )
    steps.append(step)
    
    # Steps for each GRN position
    for grn_idx, grn_pos in enumerate(grn_positions):
        # Visibility array: show structures + highlights for this GRN
        visible = [True] * num_structures  # Always show structures
        
        # Show highlights for this GRN position
        for i in range(len(grn_positions) * num_structures):
            if i >= num_structures:  # This is a highlight trace
                highlight_grn_idx = (i - num_structures) // num_structures
                visible.append(highlight_grn_idx == grn_idx)
        
        step = dict(
            method="update",
            args=[{"visible": visible}],
            label=grn_pos
        )
        steps.append(step)
    
    # Create slider
    sliders = [dict(
        active=0,
        yanchor="top",
        xanchor="left",
        currentvalue=dict(
            prefix="GRN Position: ",
            visible=True,
            xanchor="right"
        ),
        transition=dict(duration=300),
        pad=dict(b=10, t=50),
        len=0.9,
        x=0.1,
        y=0,
        steps=steps
    )]
    
    # Update layout
    fig.update_layout(
        title="GRN-Annotated GPCR Structures (Aligned)",
        scene=dict(
            xaxis_title="X (Å)",
            yaxis_title="Y (Å)",
            zaxis_title="Z (Å)",
            aspectmode='data',
        ),
        showlegend=True,
        legend=dict(
            x=0.7,
            y=0.9,
            bgcolor='rgba(255, 255, 255, 0.8)',
        ),
        sliders=sliders,
        height=800,
    )
    
    # Add annotation explaining the slider
    fig.add_annotation(
        text="Use the slider below to highlight specific GRN positions across all structures",
        xref="paper", yref="paper",
        x=0.5, y=1.05,
        showarrow=False,
        font=dict(size=12),
    )
    
    return fig


def save_visualization(fig: go.Figure, filename: str = "grn_structure_alignment.html") -> Path:
    """Save the interactive visualization."""
    output_dir = Path("data/structure/visualizations")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    output_path = output_dir / filename
    fig.write_html(str(output_path))
    
    return output_path


def create_grn_heatmap(aligned_chains: Dict[str, pd.DataFrame]) -> go.Figure:
    """Create a heatmap showing GRN conservation across structures."""
    # Get all unique GRN positions
    grn_positions = get_unique_grn_positions(aligned_chains)
    
    # Create matrix: structures x GRN positions
    structure_names = list(aligned_chains.keys())
    matrix = []
    
    for struct_name in structure_names:
        ca_atoms = aligned_chains[struct_name]
        row = []
        for grn in grn_positions:
            # Check if this GRN position exists in this structure
            if grn in ca_atoms['grn'].values:
                residue = ca_atoms[ca_atoms['grn'] == grn].iloc[0]
                # Use 1 for present, could also encode amino acid type
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=grn_positions,
        y=structure_names,
        colorscale='Blues',
        showscale=True,
        hovertemplate='Structure: %{y}<br>GRN: %{x}<br>Present: %{z}<extra></extra>',
    ))
    
    fig.update_layout(
        title="GRN Position Coverage Across Structures",
        xaxis_title="GRN Position",
        yaxis_title="Structure/Chain",
        height=400 + len(structure_names) * 30,
        xaxis=dict(
            tickmode='linear',
            tick0=0,
            dtick=10,  # Show every 10th GRN position
            tickangle=-45,
        ),
    )
    
    return fig


def main() -> None:
    """Main workflow for GRN structure visualization."""
    # Initialize
    data_root = ensure_data_root()
    print(f"Data root: {data_root}\n")
    
    # Load GRN-annotated structures
    print("Loading GRN-annotated structures...")
    struct_proc, structures = load_grn_annotated_structures()
    
    if not structures:
        print("No GRN-annotated structures found!")
        print("Please run test_structure_grn_annotation.py first.")
        return
    
    # Identify GPCR chains
    print("\nIdentifying GPCR chains with GRN annotations...")
    gpcr_chains = get_gpcr_chains(structures)
    
    if not gpcr_chains:
        print("No chains with GRN annotations found!")
        return
    
    # Align structures
    aligned_chains = align_structures_by_grn(struct_proc, structures, gpcr_chains)
    
    if not aligned_chains:
        print("No aligned chains to visualize!")
        return
    
    # Create main visualization
    print("\nCreating interactive 3D visualization...")
    fig = create_grn_visualization(aligned_chains)
    
    if fig:
        output_path = save_visualization(fig)
        print(f"✓ Saved interactive visualization to: {output_path}")
        
        # Create conservation heatmap
        print("\nCreating GRN conservation heatmap...")
        heatmap_fig = create_grn_heatmap(aligned_chains)
        heatmap_path = save_visualization(heatmap_fig, "grn_conservation_heatmap.html")
        print(f"✓ Saved heatmap to: {heatmap_path}")
        
        # Print summary
        print("\n=== Visualization Summary ===")
        print(f"Structures visualized: {len(aligned_chains)}")
        grn_positions = get_unique_grn_positions(aligned_chains)
        print(f"Total GRN positions: {len(grn_positions)}")
        print("\nKey conserved positions included:")
        key_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
        for grn in key_grns:
            if grn in grn_positions:
                print(f"  ✓ GRN {grn}")
        
        print(f"\nOpen {output_path} in your browser to interact with the visualization!")


if __name__ == "__main__":
    main()