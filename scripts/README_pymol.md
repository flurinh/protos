# PyMOL GPCR Visualization

## Files

- `gpcr_visualization.pml` - Main script: loads, aligns, and visualizes 8 adrenergic receptor structures
- `gpcr_visualization_pocket1.pml` - ADRB1 binding pocket view (2VT4, 2Y02, 2Y04, 2Y00)
- `gpcr_visualization_pocket2.pml` - ADRB2 binding pocket view (2RH1, 3NY9, 3SN6, 4LDO)
- `gpcr_visualization_pocket3.pml` - Active state binding pocket (3SN6, 4LDO, 2Y02)
- `gpcr_visualization_pocket4.pml` - Inactive state binding pocket (2RH1, 3NY9, 2VT4)
- `gpcr_visualization_pocket5.pml` - Agonist-bound binding pocket (all agonists)
- `gpcr_visualization_pocket6.pml` - Inverse agonist binding pocket (2RH1, 3NY9)
- `gpcr_render_all.pml` - Batch render all views to PNG
- `prepare_pymol_grn_labels.py` - Python script to generate GRN annotation files

## Generated Files (in `data/visualizations/`)

### Rendered Images
- `overall_structures.png` - All 8 structures aligned
- `pocket1_adrb1.png` - ADRB1 binding pocket with waters
- `pocket2_adrb2.png` - ADRB2 binding pocket with waters
- `pocket3_active.png` - Active state pocket with waters
- `pocket4_inactive.png` - Inactive state pocket with waters
- `pocket5_agonist.png` - Agonist-bound pocket with waters
- `pocket6_inverse_agonist.png` - Inverse agonist pocket with waters

### GRN Annotation Scripts
- `grn_selections.pml` - GRN position selections (3.28, 3.32, 5.43, etc.)
- `grn_colors.pml` - Color scheme by TM helix
- `water_networks.pml` - Water-mediated H-bond distances
- `state_comparison.pml` - Active vs inactive coloring
- `gpcr_master.pml` - Master script that loads all annotations

## Running PyMOL

### Option 1: Conda/Mamba (Recommended)

```bash
# Create a separate environment for PyMOL
mamba create -n pymol-env python=3.11 pymol-open-source -c conda-forge
mamba activate pymol-env

# Run from project directory
cd /data/fast/projects/protos
pymol scripts/gpcr_visualization.pml
```

### Option 2: Docker

```bash
# Run PyMOL in Docker
docker run -it --rm \
    -v /data/fast/projects/protos:/workspace \
    -e DISPLAY=$DISPLAY \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    schrodinger/pymol \
    pymol /workspace/scripts/gpcr_visualization.pml
```

### Option 3: PyMOL AppImage

Download from: https://pymol.org/2/
```bash
chmod +x PyMOL-*.AppImage
./PyMOL-*.AppImage scripts/gpcr_visualization.pml
```

## Usage

### Interactive Mode

```bash
# Load main visualization
pymol scripts/gpcr_visualization.pml

# Then in PyMOL, load a specific pocket view:
@scripts/gpcr_visualization_pocket1.pml  # ADRB1
@scripts/gpcr_visualization_pocket2.pml  # ADRB2
@scripts/gpcr_visualization_pocket3.pml  # Active states
@scripts/gpcr_visualization_pocket4.pml  # Inactive states
@scripts/gpcr_visualization_pocket5.pml  # Agonist-bound
@scripts/gpcr_visualization_pocket6.pml  # Inverse agonist
```

### Batch Rendering

```bash
# Render all 7 views to PNG
pymol -cq scripts/gpcr_visualization.pml scripts/gpcr_render_all.pml
# Output: /tmp/gpcr_output/*.png (also copied to data/visualizations/)
```

### PyMOL Commands

```
# Switch views
view pocket_view, recall
view overall_view, recall

# Save/render
save_session
render_current my_figure
```

## Structures

| PDB | Receptor | Ligand | Type | State | Color |
|-----|----------|--------|------|-------|-------|
| 2RH1 | β2-AR | carazolol | inverse_agonist | inactive | blue |
| 3NY9 | β2-AR | ICI-118551 | inverse_agonist | inactive | slate |
| 3SN6 | β2-AR | BI-167107 | full_agonist | active | green |
| 4LDO | β2-AR | adrenaline | full_agonist | active | lime |
| 2VT4 | β1-AR | cyanopindolol | antagonist | inactive | teal |
| 2Y02 | β1-AR | isoprenaline | full_agonist | active_like | chartreuse |
| 2Y04 | β1-AR | salbutamol | partial_agonist | intermediate | yellow |
| 2Y00 | β1-AR | dobutamine | partial_agonist | intermediate | orange |

## Key Observations

1. **TM6 Movement**: Active states show ~14Å outward movement of TM6
2. **Water Networks**: Inactive states have many more water networks (2RH1: 30 vs 3SN6: 0)
3. **Binding Pocket**: Agonists make closer contacts at positions 3.32, 5.43, 7.38
