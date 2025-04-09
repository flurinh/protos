"""
Visualization tools for GRN assignments.

This module provides functions for visualizing GRN assignments,
including sequence alignments with GRN annotations.
"""

import os
import logging
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from typing import Dict, List, Tuple, Union, Optional, Any

from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.schema.grn_utils_updated import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string,
    sort_grns
)

# Configure logging
logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)

def visualize_alignment_with_grns(
    query_seq: str,
    ref_seq: str,
    alignment: Tuple[str, str, str],
    query_grns: Dict[str, str],
    ref_grns: Dict[str, str],
    output_path: Optional[str] = None,
    title: str = "Sequence Alignment with GRN Annotations",
    dpi: int = 300
) -> plt.Figure:
    """
    Visualize a sequence alignment with GRN annotations.
    
    Args:
        query_seq: Query amino acid sequence
        ref_seq: Reference amino acid sequence
        alignment: Alignment tuple from format_alignment
        query_grns: Dictionary mapping GRNs to query sequence positions
        ref_grns: Dictionary mapping GRNs to reference sequence positions
        output_path: Path to save the visualization
        title: Plot title
        dpi: DPI for saved figure
        
    Returns:
        Matplotlib figure object
    """
    query_aligned, match_line, ref_aligned = alignment
    
    # Create a figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 8), gridspec_kw={'height_ratios': [1, 0.5, 1]})
    
    # Invert query_grns and ref_grns for position lookup
    query_pos_to_grn = {v: k for k, v in query_grns.items()}
    ref_pos_to_grn = {v: k for k, v in ref_grns.items()}
    
    # Process the aligned sequences
    query_positions = []
    ref_positions = []
    q_idx = 0
    r_idx = 0
    
    for i in range(len(match_line)):
        if query_aligned[i] != '-':
            q_idx += 1
            query_positions.append(q_idx)
        else:
            query_positions.append(None)
            
        if ref_aligned[i] != '-':
            r_idx += 1
            ref_positions.append(r_idx)
        else:
            ref_positions.append(None)
    
    # Setup colors for different GRNs
    helix_colors = {
        1: 'tab:blue',
        2: 'tab:orange',
        3: 'tab:green',
        4: 'tab:red',
        5: 'tab:purple',
        6: 'tab:brown',
        7: 'tab:pink',
        8: 'tab:olive'
    }
    
    loop_color = 'tab:gray'
    n_term_color = 'tab:cyan'
    c_term_color = 'tab:cyan'
    
    # Helper function to get color for a GRN
    def get_grn_color(grn):
        if grn is None:
            return 'white'
        
        normalized = normalize_grn_format(grn)
        if 'x' in normalized:
            helix = int(normalized.split('x')[0])
            return helix_colors.get(helix, 'white')
        elif normalized.startswith('n.'):
            return n_term_color
        elif normalized.startswith('c.'):
            return c_term_color
        elif '.' in normalized and normalized[0].isdigit() and normalized[1].isdigit():
            return loop_color
        else:
            return 'white'
    
    # Plot the reference sequence
    for i, (aa, pos) in enumerate(zip(ref_aligned, ref_positions)):
        if pos is not None:
            pos_str = f"{aa}{pos}"
            grn = ref_pos_to_grn.get(pos_str)
            color = get_grn_color(grn)
            ax1.text(i, 0, aa, ha='center', va='center', fontsize=10, 
                    bbox=dict(facecolor=color, alpha=0.3, edgecolor='none'))
            if grn:
                ax1.text(i, 0.5, grn, ha='center', va='center', fontsize=8, rotation=90)
        else:
            ax1.text(i, 0, '-', ha='center', va='center', fontsize=10)
    
    # Plot the match line
    for i, c in enumerate(match_line):
        ax2.text(i, 0, c, ha='center', va='center', fontsize=10)
    
    # Plot the query sequence
    for i, (aa, pos) in enumerate(zip(query_aligned, query_positions)):
        if pos is not None:
            pos_str = f"{aa}{pos}"
            grn = query_pos_to_grn.get(pos_str)
            color = get_grn_color(grn)
            ax3.text(i, 0, aa, ha='center', va='center', fontsize=10, 
                    bbox=dict(facecolor=color, alpha=0.3, edgecolor='none'))
            if grn:
                ax3.text(i, -0.5, grn, ha='center', va='center', fontsize=8, rotation=90)
        else:
            ax3.text(i, 0, '-', ha='center', va='center', fontsize=10)
    
    # Remove axis ticks
    for ax in [ax1, ax2, ax3]:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-0.5, len(alignment[0]) - 0.5)
        ax.set_ylim(-1, 1)
    
    # Set labels
    ax1.set_title('Reference Sequence')
    ax3.set_title('Query Sequence')
    fig.suptitle(title, fontsize=16)
    
    # Add legend for helices
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=f'Helix {helix}')
        for helix, color in helix_colors.items()
    ]
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=loop_color, markersize=10, label='Loop Region'))
    legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=n_term_color, markersize=10, label='N/C Termini'))
    
    fig.legend(handles=legend_elements, loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0), fontsize=10)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15)
    
    # Save if output path is provided
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    
    return fig

def visualize_grn_table(
    grn_table_path: str,
    protein_ids: Optional[List[str]] = None,
    grn_subset: Optional[List[str]] = None,
    normalize_formats: bool = True,
    output_path: Optional[str] = None,
    title: str = "GRN Table Visualization",
    figsize: Tuple[int, int] = (15, 10),
    dpi: int = 300
) -> plt.Figure:
    """
    Visualize a GRN table as a heatmap.
    
    Args:
        grn_table_path: Path to GRN table CSV file
        protein_ids: Optional list of protein IDs to include (all if None)
        grn_subset: Optional list of GRNs to include (all if None)
        normalize_formats: Whether to normalize GRN formats
        output_path: Path to save the visualization
        title: Plot title
        figsize: Figure size
        dpi: DPI for saved figure
        
    Returns:
        Matplotlib figure object
    """
    # Load the GRN table
    df = pd.read_csv(grn_table_path, index_col=0)
    
    # Filter proteins if specified
    if protein_ids:
        df = df.loc[df.index.isin(protein_ids)]
    
    # Normalize column names if requested
    if normalize_formats:
        df.columns = [normalize_grn_format(col) for col in df.columns]
    
    # Sort columns by GRN order
    df = df.loc[:, sort_grns(df.columns.tolist())]
    
    # Filter GRNs if specified
    if grn_subset:
        normalized_subset = [normalize_grn_format(grn) for grn in grn_subset]
        df = df.loc[:, df.columns.isin(normalized_subset)]
    
    # Create amino acid mapping
    aa_list = list("ACDEFGHIKLMNPQRSTVWY-")
    aa_to_num = {aa: i for i, aa in enumerate(aa_list)}
    
    # Replace amino acids with numbers
    numeric_data = df.applymap(lambda x: aa_to_num.get(x, -1) if isinstance(x, str) else -1)
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create custom colormap for amino acids
    cmap = plt.cm.tab20
    bounds = np.arange(-0.5, len(aa_list) + 0.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    # Plot heatmap
    im = ax.imshow(numeric_data, cmap=cmap, norm=norm)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, ticks=range(len(aa_list)))
    cbar.ax.set_yticklabels(aa_list)
    cbar.set_label('Amino Acid')
    
    # Add labels
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, rotation=90)
    
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index)
    
    # Add grid
    ax.set_xticks(np.arange(-0.5, len(df.columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(df.index), 1), minor=True)
    ax.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
    
    # Color-code the helices in the x-axis labels
    for i, grn in enumerate(df.columns):
        helix_color = 'black'
        if 'x' in grn:
            helix = int(grn.split('x')[0])
            helix_color = helix_colors.get(helix, 'black')
        ax.get_xticklabels()[i].set_color(helix_color)
    
    # Add title
    plt.title(title)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save if output path is provided
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    
    return fig

def visualize_grn_distribution(
    grn_table_path: str,
    grn_subset: Optional[List[str]] = None,
    normalize_formats: bool = True,
    output_path: Optional[str] = None,
    title: str = "GRN Distribution",
    figsize: Tuple[int, int] = (15, 10),
    dpi: int = 300
) -> plt.Figure:
    """
    Visualize the distribution of amino acids at each GRN position.
    
    Args:
        grn_table_path: Path to GRN table CSV file
        grn_subset: Optional list of GRNs to include (all if None)
        normalize_formats: Whether to normalize GRN formats
        output_path: Path to save the visualization
        title: Plot title
        figsize: Figure size
        dpi: DPI for saved figure
        
    Returns:
        Matplotlib figure object
    """
    # Load the GRN table
    df = pd.read_csv(grn_table_path, index_col=0)
    
    # Normalize column names if requested
    if normalize_formats:
        df.columns = [normalize_grn_format(col) for col in df.columns]
    
    # Sort columns by GRN order
    df = df.loc[:, sort_grns(df.columns.tolist())]
    
    # Filter GRNs if specified
    if grn_subset:
        normalized_subset = [normalize_grn_format(grn) for grn in grn_subset]
        df = df.loc[:, df.columns.isin(normalized_subset)]
    
    # Create figure
    fig, axes = plt.subplots(len(df.columns), 1, figsize=figsize)
    
    # Create a dictionary to store amino acid frequencies
    aa_list = list("ACDEFGHIKLMNPQRSTVWY-")
    
    # Plot distribution for each GRN
    for i, grn in enumerate(df.columns):
        ax = axes[i] if len(df.columns) > 1 else axes
        
        # Calculate amino acid frequencies
        aa_counts = df[grn].value_counts().reindex(aa_list, fill_value=0)
        
        # Plot barplot
        bars = ax.bar(aa_counts.index, aa_counts.values)
        
        # Add GRN label
        ax.set_title(f"GRN: {grn}")
        ax.set_xlabel("Amino Acid")
        ax.set_ylabel("Frequency")
        
        # Set y-axis limit
        ax.set_ylim(0, max(aa_counts.values) * 1.1)
        
        # Add value labels
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                        f"{height:.0f}", ha='center', va='bottom')
    
    # Add overall title
    fig.suptitle(title, fontsize=16)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    
    # Save if output path is provided
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    
    return fig

def main():
    """Main function when run as script."""
    parser = argparse.ArgumentParser(description='Visualize GRN assignments')
    parser.add_argument('-t', '--table', required=True, 
                        help='Path to GRN table CSV file')
    parser.add_argument('-o', '--output', required=True,
                        help='Path to save visualization')
    parser.add_argument('-v', '--vis_type', default='heatmap',
                        choices=['heatmap', 'distribution'],
                        help='Type of visualization')
    parser.add_argument('-g', '--grns', 
                        help='Comma-separated list of GRNs to include')
    parser.add_argument('-p', '--proteins', 
                        help='Comma-separated list of protein IDs to include')
    parser.add_argument('--title', default=None,
                        help='Plot title')
    parser.add_argument('--dpi', type=int, default=300,
                        help='DPI for saved figure')
    
    args = parser.parse_args()
    
    # Parse GRNs and proteins if provided
    grn_subset = args.grns.split(',') if args.grns else None
    protein_ids = args.proteins.split(',') if args.proteins else None
    
    # Generate title if not provided
    title = args.title or f"GRN {args.vis_type.capitalize()} Visualization"
    
    # Create visualization
    if args.vis_type == 'heatmap':
        visualize_grn_table(
            grn_table_path=args.table,
            protein_ids=protein_ids,
            grn_subset=grn_subset,
            output_path=args.output,
            title=title,
            dpi=args.dpi
        )
    elif args.vis_type == 'distribution':
        visualize_grn_distribution(
            grn_table_path=args.table,
            grn_subset=grn_subset,
            output_path=args.output,
            title=title,
            dpi=args.dpi
        )
    
    print(f"Visualization saved to {args.output}")
    
if __name__ == '__main__':
    main()