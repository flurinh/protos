#!/usr/bin/env python
"""
Visualize microbial opsin structures with GRN annotations.
This script properly uses the built-in assign_grns method which includes expand_annotation.
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# Setup paths following Protos principles
project_root = Path(__file__).parent.absolute()
data_dir = project_root / "data"
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# Protos imports
from protos.io.paths import ProtosPaths
from protos.processing.structure.structure_processor import StructureProcessor
from protos.processing.grn.grn_utils import sort_grns_str

def main():
    """Main function to create GRN visualizations for microbial opsins."""
    
    # Initialize processor - let it handle everything internally
    struct_proc = StructureProcessor(name="mo_viz")
    
    # Load structure
    print("Microbial Opsin GRN Visualization (Correct Implementation)")
    print("="*60)
    print("\nLoading microbial opsin structure...")
    
    # Try to load common microbial opsin structures
    pdb_ids = ["1UAZ", "6RMK", "1M0L", "3UG9", "4HYJ", "4PXK", "1VGO", "1XIO", "1JGJ"]
    
    loaded_pdb_id = None
    for pdb_id in pdb_ids:
        try:
            print(f"Trying to load {pdb_id}...")
            struct_proc.load_structure(pdb_id)
            struct_proc.update_pdb_ids()
            loaded_pdb_id = pdb_id
            print(f"Successfully loaded {pdb_id}")
            break
        except Exception as e:
            print(f"Could not load {pdb_id}: {e}")
            continue
    
    if not loaded_pdb_id:
        print("ERROR: Could not load any microbial opsin structure!")
        return
    
    pdb_id = loaded_pdb_id
    print(f"\nUsing {pdb_id} with {len(struct_proc.data)} atoms")
    
    # Use the built-in assign_grns method which properly:
    # 1. Loads the reference table using processor methods
    # 2. Uses expand_annotation functionality
    # 3. Handles all the complex mapping correctly
    print("\nAssigning GRN positions for microbial opsins...")
    print("(Using built-in assign_grns with expand_annotation)")
    
    grn_assignments = struct_proc.assign_grns(
        protein_family='microbial_opsins',
        reference_table='mo_ref',  # This will be loaded by the processor
        similarity_threshold=0.20,  # 20% as used in annotate_microbial_opsins_proper.py
        verbose=True
    )
    
    print(f"\nAssigned GRNs to {len(grn_assignments)} chains")
    
    # Find the chain with most GRN assignments
    grn_data = struct_proc.data[struct_proc.data['grn'].notna()]
    if grn_data.empty:
        print("ERROR: No GRN assignments found!")
        return
    
    chain_counts = grn_data.groupby('auth_chain_id').size()
    opsin_chain = chain_counts.idxmax()
    print(f"\nUsing chain {opsin_chain} (has {chain_counts[opsin_chain]} residues with GRN)")
    
    # Show GRN coverage statistics
    print("\nGRN Assignment Statistics:")
    unique_grns = sort_grns_str(grn_data['grn'].unique().tolist())
    print(f"  Total unique GRN positions: {len(unique_grns)}")
    print(f"  GRN range: {unique_grns[0]} to {unique_grns[-1]}")
    
    # Count by helix
    helix_counts = {}
    for grn in unique_grns:
        helix = grn.split('.')[0]
        helix_counts[helix] = helix_counts.get(helix, 0) + 1
    
    print("\n  Positions per helix:")
    for helix in sorted(helix_counts.keys()):
        print(f"    Helix {helix}: {helix_counts[helix]} positions")
    
    # Verify key positions
    print("\nVerifying key microbial opsin positions:")
    verify_key_positions(struct_proc, pdb_id, opsin_chain)
    
    # Create visualizations
    print("\nCreating visualizations...")
    create_helix_wheel_view(struct_proc, pdb_id, opsin_chain)
    create_snake_plot(struct_proc, pdb_id, opsin_chain)
    create_3d_structure_view(struct_proc, pdb_id, opsin_chain)
    create_interactive_3d_with_slider(struct_proc, pdb_id, opsin_chain)
    create_grn_heatmap(struct_proc, pdb_id, opsin_chain)
    
    print("\n" + "="*60)
    print("Visualization files created:")
    print("- mo_helix_wheel_correct.html")
    print("- mo_snake_plot_correct.html")
    print("- mo_3d_structure_correct.html")
    print("- mo_3d_interactive_slider_correct.html")
    print("- mo_grn_heatmap_correct.html")


def verify_key_positions(struct_proc, pdb_id, chain_id):
    """Verify key microbial opsin positions are correctly mapped."""
    chain_mask = (struct_proc.data['pdb_id'] == pdb_id) & (struct_proc.data['auth_chain_id'] == chain_id)
    grn_data = struct_proc.data[chain_mask & struct_proc.data['grn'].notna()]
    
    # Key x.50 positions in microbial opsins
    key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    
    print("\nKey conserved positions (x.50) found:")
    for grn in key_positions:
        residues = grn_data[(grn_data['grn'] == grn) & (grn_data['res_atom_name'] == 'CA')]
        if not residues.empty:
            res = residues.iloc[0]
            print(f"  {grn}: {res['res_name3l']}{res['auth_seq_id']}")
        else:
            print(f"  {grn}: Not found")


def create_helix_wheel_view(struct_proc, pdb_id, chain_id):
    """Create a 2D helix wheel projection view for microbial opsins."""
    print("  Creating helix wheel view...")
    
    chain_mask = (struct_proc.data['pdb_id'] == pdb_id) & (struct_proc.data['auth_chain_id'] == chain_id)
    grn_data = struct_proc.data[chain_mask & struct_proc.data['grn'].notna() & 
                                (struct_proc.data['res_atom_name'] == 'CA')]
    
    fig = make_subplots(
        rows=2, cols=4,
        subplot_titles=[f'Helix {i}' for i in range(1, 8)] + [''],
        specs=[[{'type': 'polar'}]*4, [{'type': 'polar'}]*4]
    )
    
    helix_colors = ['#FF6B6B', '#4ECDC4', '#45B7D1', '#96CEB4', '#FECA57', '#FF9FF3', '#54A0FF']
    
    for helix_num in range(1, 8):
        helix_data = grn_data[grn_data['grn'].str.startswith(f'{helix_num}.')]
        
        if not helix_data.empty:
            # Calculate angles for helix wheel (100 degrees per residue for alpha helix)
            angles = []
            radii = []
            texts = []
            colors = []
            
            sorted_data = helix_data.sort_values('auth_seq_id')
            for i, (_, res) in enumerate(sorted_data.iterrows()):
                angle = (i * 100) % 360
                angles.append(angle)
                radii.append(1)
                texts.append(f"{res['res_name1l']}{res['auth_seq_id']}<br>{res['grn']}")
                
                # Highlight conserved x.50 positions
                if res['grn'].endswith('.50'):
                    colors.append('orange')
                else:
                    colors.append(helix_colors[helix_num-1])
            
            row = (helix_num - 1) // 4 + 1
            col = (helix_num - 1) % 4 + 1
            
            fig.add_trace(
                go.Scatterpolar(
                    r=radii,
                    theta=angles,
                    mode='markers+text',
                    marker=dict(size=20, color=colors),
                    text=texts,
                    textposition='top center',
                    name=f'Helix {helix_num}',
                    showlegend=False
                ),
                row=row, col=col
            )
    
    fig.update_layout(
        title="Microbial Opsin Helix Wheel Projections (with expand_annotation)",
        height=800,
        showlegend=False
    )
    
    fig.write_html("mo_helix_wheel_correct.html")


def create_snake_plot(struct_proc, pdb_id, chain_id):
    """Create a 2D snake-like diagram of the microbial opsin."""
    print("  Creating snake plot...")
    
    chain_mask = (struct_proc.data['pdb_id'] == pdb_id) & (struct_proc.data['auth_chain_id'] == chain_id)
    grn_data = struct_proc.data[chain_mask & struct_proc.data['grn'].notna() & 
                                (struct_proc.data['res_atom_name'] == 'CA')]
    
    # Define helix positions for snake plot
    helix_positions = {
        1: {'x': 0, 'y': 0, 'direction': 'up'},
        2: {'x': 1, 'y': 0, 'direction': 'down'},
        3: {'x': 2, 'y': 0, 'direction': 'up'},
        4: {'x': 3, 'y': 0, 'direction': 'down'},
        5: {'x': 4, 'y': 0, 'direction': 'up'},
        6: {'x': 5, 'y': 0, 'direction': 'down'},
        7: {'x': 6, 'y': 0, 'direction': 'up'},
    }
    
    fig = go.Figure()
    
    # Plot each helix
    for helix_num in range(1, 8):
        helix_data = grn_data[grn_data['grn'].str.startswith(f'{helix_num}.')].sort_values('auth_seq_id')
        
        if not helix_data.empty:
            pos_info = helix_positions[helix_num]
            x_base = pos_info['x'] * 100
            
            for i, (_, res) in enumerate(helix_data.iterrows()):
                if pos_info['direction'] == 'up':
                    y = i * 10
                else:
                    y = -i * 10
                
                # Color based on properties
                if res['grn'].endswith('.50'):  # Conserved x.50 positions
                    color = 'orange'
                    size = 15
                elif res['res_name1l'] in 'RKH':  # Basic
                    color = 'blue'
                    size = 10
                elif res['res_name1l'] in 'DE':  # Acidic
                    color = 'red'
                    size = 10
                elif res['res_name1l'] in 'AVLIMFYW':  # Hydrophobic
                    color = 'gray'
                    size = 10
                else:
                    color = 'lightgray'
                    size = 10
                
                fig.add_trace(go.Scatter(
                    x=[x_base],
                    y=[y],
                    mode='markers+text',
                    marker=dict(size=size, color=color),
                    text=f"{res['res_name1l']}<br>{res['grn']}",
                    textposition='middle right',
                    showlegend=False,
                    hovertemplate=f"{res['res_name3l']}{res['auth_seq_id']}<br>GRN: {res['grn']}<extra></extra>"
                ))
    
    # Add helix labels
    for helix_num, pos_info in helix_positions.items():
        fig.add_annotation(
            x=pos_info['x'] * 100,
            y=-50 if pos_info['direction'] == 'up' else 50,
            text=f"TM{helix_num}",
            showarrow=False,
            font=dict(size=14, color='black')
        )
    
    # Add legend for key residues
    fig.add_annotation(
        x=700, y=100,
        text="Key: Orange=x.50 conserved positions, Blue=Basic, Red=Acidic, Gray=Hydrophobic",
        showarrow=False,
        font=dict(size=12)
    )
    
    fig.update_layout(
        title="Microbial Opsin Snake Plot - 2D Topology (with expand_annotation)",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        height=800,
        width=1200,
        plot_bgcolor='white'
    )
    
    fig.write_html("mo_snake_plot_correct.html")


def create_3d_structure_view(struct_proc, pdb_id, chain_id):
    """Create a 3D structure view highlighting functional regions."""
    print("  Creating 3D structure view...")
    
    chain_mask = (struct_proc.data['pdb_id'] == pdb_id) & (struct_proc.data['auth_chain_id'] == chain_id)
    chain_data = struct_proc.data[chain_mask].copy()
    
    # Get CA atoms for backbone
    ca_atoms = chain_data[chain_data['res_atom_name'] == 'CA'].sort_values('auth_seq_id')
    
    # Get GRN data
    grn_data = chain_data[chain_data['grn'].notna()]
    unique_grns = sort_grns_str(grn_data['grn'].unique().tolist())
    
    # Create figure
    fig = go.Figure()
    
    # Add backbone trace
    fig.add_trace(
        go.Scatter3d(
            x=ca_atoms['x'].values,
            y=ca_atoms['y'].values,
            z=ca_atoms['z'].values,
            mode='lines',
            line=dict(color='lightgray', width=5),
            name='Backbone',
            hoverinfo='skip'
        )
    )
    
    # Add all GRN positions with different colors per helix
    helix_colors = {
        '1': '#FF6B6B', '2': '#4ECDC4', '3': '#45B7D1', '4': '#96CEB4',
        '5': '#FECA57', '6': '#FF9FF3', '7': '#54A0FF'
    }
    
    # Group by helix
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        helix_grns = [grn for grn in unique_grns if grn.startswith(f'{helix_num}.')]
        if helix_grns:
            helix_mask = grn_data['grn'].isin(helix_grns) & (grn_data['res_atom_name'] == 'CA')
            helix_atoms = grn_data[helix_mask]
            
            if not helix_atoms.empty:
                hover_text = []
                for _, res in helix_atoms.iterrows():
                    hover_text.append(
                        f"{res['res_name3l']}{res['auth_seq_id']}<br>"
                        f"GRN: {res['grn']}<br>"
                        f"Helix {helix_num}"
                    )
                
                fig.add_trace(
                    go.Scatter3d(
                        x=helix_atoms['x'].values,
                        y=helix_atoms['y'].values,
                        z=helix_atoms['z'].values,
                        mode='markers',
                        marker=dict(
                            size=6,
                            color=helix_colors.get(helix_num, 'gray'),
                            line=dict(width=1, color='black')
                        ),
                        text=hover_text,
                        hovertemplate='%{text}<extra></extra>',
                        name=f'Helix {helix_num}',
                        legendgroup=f'helix{helix_num}'
                    )
                )
    
    # Highlight conserved x.50 positions
    key_positions = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    key_mask = grn_data['grn'].isin(key_positions) & (grn_data['res_atom_name'] == 'CA')
    key_atoms = grn_data[key_mask]
    
    if not key_atoms.empty:
        hover_text = []
        for _, res in key_atoms.iterrows():
            hover_text.append(
                f"CONSERVED POSITION<br>"
                f"{res['res_name3l']}{res['auth_seq_id']}<br>"
                f"GRN: {res['grn']}"
            )
        
        fig.add_trace(
            go.Scatter3d(
                x=key_atoms['x'].values,
                y=key_atoms['y'].values,
                z=key_atoms['z'].values,
                mode='markers',
                marker=dict(
                    size=12,
                    color='orange',
                    symbol='square',
                    line=dict(width=2, color='darkorange')
                ),
                text=hover_text,
                hovertemplate='%{text}<extra></extra>',
                name='Conserved Positions (x.50)'
            )
        )
    
    # Update layout
    fig.update_layout(
        title=f"Microbial Opsin Structure {pdb_id} Chain {chain_id} - GRN Positions (with expand_annotation)",
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),
                center=dict(x=0, y=0, z=0)
            ),
            aspectmode='data'
        ),
        height=800,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )
    
    fig.write_html("mo_3d_structure_correct.html")


def create_interactive_3d_with_slider(struct_proc, pdb_id, chain_id):
    """Create an interactive 3D structure view with GRN position slider."""
    print("  Creating interactive 3D structure with GRN slider...")
    
    chain_mask = (struct_proc.data['pdb_id'] == pdb_id) & (struct_proc.data['auth_chain_id'] == chain_id)
    chain_data = struct_proc.data[chain_mask].copy()
    
    # Get CA atoms for backbone
    ca_atoms = chain_data[chain_data['res_atom_name'] == 'CA'].sort_values('auth_seq_id')
    
    # Get GRN data
    grn_data = chain_data[chain_data['grn'].notna()]
    unique_grns = sort_grns_str(grn_data['grn'].unique().tolist())
    
    print(f"    Found {len(unique_grns)} unique GRN positions for slider")
    
    # Create figure
    fig = go.Figure()
    
    # Add backbone trace (always visible)
    fig.add_trace(
        go.Scatter3d(
            x=ca_atoms['x'].values,
            y=ca_atoms['y'].values,
            z=ca_atoms['z'].values,
            mode='lines+markers',
            line=dict(color='lightgray', width=4),
            marker=dict(size=3, color='lightgray'),
            name='Backbone',
            hoverinfo='skip'
        )
    )
    
    # Add initial empty trace for highlighted GRN position
    fig.add_trace(
        go.Scatter3d(
            x=[],
            y=[],
            z=[],
            mode='markers',
            marker=dict(size=15, color='red', line=dict(color='darkred', width=2)),
            name='GRN Position',
            hovertemplate='%{text}<extra></extra>'
        )
    )
    
    # Create frames for each GRN position
    frames = []
    
    for grn_pos in unique_grns:
        # Get all atoms for this GRN position
        grn_residues = grn_data[grn_data['grn'] == grn_pos]
        grn_ca = grn_residues[grn_residues['res_atom_name'] == 'CA']
        
        if not grn_ca.empty:
            # Create hover text
            hover_text_ca = []
            for _, res in grn_ca.iterrows():
                special_label = ""
                if grn_pos.endswith('.50'):
                    special_label = " (CONSERVED)"
                
                hover_text_ca.append(
                    f"GRN: {grn_pos}{special_label}<br>"
                    f"{res['res_name3l']}{res['auth_seq_id']}<br>"
                    f"Helix {grn_pos.split('.')[0]}"
                )
            
            # Determine color based on conservation
            if grn_pos.endswith('.50'):
                color = 'orange'
                size = 15
            else:
                color = 'lightblue'
                size = 12
            
            # Create frame with both backbone and highlighted residue
            frame_data = [
                # Backbone (same for all frames)
                go.Scatter3d(
                    x=ca_atoms['x'].values,
                    y=ca_atoms['y'].values,
                    z=ca_atoms['z'].values,
                    mode='lines+markers',
                    line=dict(color='lightgray', width=4),
                    marker=dict(size=3, color='lightgray'),
                    name='Backbone',
                    hoverinfo='skip'
                ),
                # Highlighted GRN position
                go.Scatter3d(
                    x=grn_ca['x'].values,
                    y=grn_ca['y'].values,
                    z=grn_ca['z'].values,
                    mode='markers',
                    marker=dict(
                        size=size,
                        color=color,
                        line=dict(color='black', width=2)
                    ),
                    text=hover_text_ca,
                    hovertemplate='%{text}<extra></extra>',
                    name=f'GRN {grn_pos}'
                )
            ]
            
            
            frame = go.Frame(
                data=frame_data,
                name=str(grn_pos)
            )
            frames.append(frame)
    
    # Add frames to figure
    fig.frames = frames
    
    # Create slider
    sliders = [{
        'active': 0,
        'yanchor': 'top',
        'xanchor': 'left',
        'currentvalue': {
            'font': {'size': 16},
            'prefix': 'GRN Position: ',
            'visible': True,
            'xanchor': 'right'
        },
        'transition': {'duration': 300, 'easing': 'cubic-in-out'},
        'pad': {'b': 10, 't': 50},
        'len': 0.9,
        'x': 0.1,
        'y': 0,
        'steps': []
    }]
    
    # Add steps for each GRN position
    for grn_pos in unique_grns:
        step = {
            'args': [[str(grn_pos)], {
                'frame': {'duration': 300, 'redraw': True},
                'mode': 'immediate',
                'transition': {'duration': 300}
            }],
            'label': str(grn_pos),
            'method': 'animate'
        }
        sliders[0]['steps'].append(step)
    
    # Set initial data to first GRN position
    if unique_grns:
        first_grn = unique_grns[0]
        grn_residues = grn_data[grn_data['grn'] == first_grn]
        grn_ca = grn_residues[grn_residues['res_atom_name'] == 'CA']
        
        if not grn_ca.empty:
            hover_text = []
            for _, res in grn_ca.iterrows():
                special_label = ""
                if first_grn.endswith('.50'):
                    special_label = " (CONSERVED)"
                
                hover_text.append(
                    f"GRN: {first_grn}{special_label}<br>"
                    f"{res['res_name3l']}{res['auth_seq_id']}<br>"
                    f"Helix {first_grn.split('.')[0]}"
                )
            
            # Update the second trace with the first GRN position data
            fig.data[1].x = grn_ca['x'].values
            fig.data[1].y = grn_ca['y'].values
            fig.data[1].z = grn_ca['z'].values
            fig.data[1].text = hover_text
            
            if first_grn.endswith('.50'):
                fig.data[1].marker.color = 'orange'
                fig.data[1].marker.size = 15
            
            fig.data[1].name = f'GRN {first_grn}'
    
    # Update layout
    fig.update_layout(
        title=f"Microbial Opsin {pdb_id} Chain {chain_id} - Interactive GRN Navigator (with expand_annotation)",
        scene=dict(
            xaxis_title='X (Å)',
            yaxis_title='Y (Å)',
            zaxis_title='Z (Å)',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.5),
                center=dict(x=0, y=0, z=0)
            ),
            aspectmode='data',
            bgcolor='white'
        ),
        sliders=sliders,
        updatemenus=[{
            'buttons': [
                {
                    'args': [None, {'frame': {'duration': 500, 'redraw': True},
                                   'fromcurrent': True,
                                   'transition': {'duration': 300, 'easing': 'quadratic-in-out'}}],
                    'label': 'Play',
                    'method': 'animate'
                },
                {
                    'args': [[None], {'frame': {'duration': 0, 'redraw': True},
                                     'mode': 'immediate',
                                     'transition': {'duration': 0}}],
                    'label': 'Pause',
                    'method': 'animate'
                }
            ],
            'direction': 'left',
            'pad': {'r': 10, 't': 87},
            'showactive': False,
            'type': 'buttons',
            'x': 0.1,
            'xanchor': 'right',
            'y': 0,
            'yanchor': 'top'
        }],
        height=800,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=0.99
        )
    )
    
    # Add annotations for key information
    fig.add_annotation(
        text="Use slider to navigate GRN positions<br>Orange = x.50 conserved positions",
        xref="paper", yref="paper",
        x=0.5, y=0.95,
        showarrow=False,
        font=dict(size=12)
    )
    
    fig.write_html("mo_3d_interactive_slider_correct.html")


def create_grn_heatmap(struct_proc, pdb_id, chain_id):
    """Create a heatmap of GRN positions and properties for microbial opsins."""
    print("  Creating GRN heatmap...")
    
    chain_mask = (struct_proc.data['pdb_id'] == pdb_id) & (struct_proc.data['auth_chain_id'] == chain_id)
    grn_data = struct_proc.data[chain_mask & struct_proc.data['grn'].notna() & 
                                (struct_proc.data['res_atom_name'] == 'CA')]
    
    # Create a matrix of properties for each helix
    helices = ['1', '2', '3', '4', '5', '6', '7']
    properties = {
        'Hydrophobic': lambda x: x in 'AVLIMFYW',
        'Aromatic': lambda x: x in 'FYW',
        'Positive': lambda x: x in 'RKH',
        'Negative': lambda x: x in 'DE',
        'Polar': lambda x: x in 'STNQ',
        'Small': lambda x: x in 'AGST',
        'Proline': lambda x: x == 'P',
        'Glycine': lambda x: x == 'G',
        'Cysteine': lambda x: x == 'C'
    }
    
    # Build data for heatmap
    heatmap_data = []
    grn_labels = []
    
    for helix in helices:
        helix_data = grn_data[grn_data['grn'].str.startswith(f'{helix}.')].sort_values('grn')
        
        for _, res in helix_data.iterrows():
            row_data = []
            label = f"{res['grn']} ({res['res_name1l']}{res['auth_seq_id']})"
            
            # Add special marker for conserved positions
            if res['grn'].endswith('.50'):
                label += " *"
            
            grn_labels.append(label)
            
            for prop_name, prop_func in properties.items():
                if prop_func(res['res_name1l']):
                    row_data.append(1)
                else:
                    row_data.append(0)
            
            heatmap_data.append(row_data)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=list(properties.keys()),
        y=grn_labels,
        colorscale='Viridis',
        showscale=True,
        hovertemplate='GRN: %{y}<br>Property: %{x}<br>Value: %{z}<extra></extra>'
    ))
    
    # Add helix separators
    helix_boundaries = []
    current_helix = None
    for i, label in enumerate(grn_labels):
        helix_num = label.split('.')[0]
        if helix_num != current_helix:
            if current_helix is not None:
                helix_boundaries.append(i - 0.5)
            current_helix = helix_num
    
    for boundary in helix_boundaries:
        fig.add_shape(
            type="line",
            x0=-0.5, x1=len(properties)-0.5,
            y0=boundary, y1=boundary,
            line=dict(color="white", width=2)
        )
    
    fig.update_layout(
        title="Microbial Opsin GRN Position Properties (with expand_annotation)<br>* = x.50 conserved positions",
        xaxis_title="Residue Properties",
        yaxis_title="GRN Positions",
        height=1200,
        width=800
    )
    
    fig.write_html("mo_grn_heatmap_correct.html")


if __name__ == "__main__":
    main()