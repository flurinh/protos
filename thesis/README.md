# ProtOS Thesis Figures

This folder contains scripts and assets for generating thesis figures using ProtOS.

## Folder Structure

```
thesis/
├── shared/                          # Shared settings (use in all scripts)
│   ├── colorscales.yaml             # Color definitions
│   └── thesis_pymol_settings.pml    # PyMOL settings
│
├── ch1_introduction/                # Chapter 1: Introduction
│   ├── fig_1.1_type_comparison/     # Figure 1.1: Type I/II structural comparison
│   └── fig_1.2_spectral_tuning/     # Figure 1.2: Spectral tuning problem
│
├── ch2_protos/                      # Chapter 2: ProtOS
│   ├── fig_2.2_opsin_diversity/     # Figure 2.2: Cone opsin diversity
│   ├── fig_2.3_br_binding_pocket/   # Figure 2.3: BR binding pocket
│   ├── fig_2.4_pocket_graphs/       # Figure 2.4: Binding pocket graphs
│   ├── fig_2.5_grn_comparison/      # Figure 2.5: GRN comparison
│   └── fig_2.6_embeddings/          # Figure 2.6: t-SNE + heatmap
│
├── ch5_applications/                # Chapter 5: Applications
│   └── fig_5.2_rhodozyme/           # Figure 5.2: Rhodozyme design
│
├── outputs/                         # Data outputs (JSON, CSV, etc.)
│   ├── binding_pocket/
│   ├── embedding/
│   ├── graph/
│   ├── grn/
│   ├── molecule/
│   ├── property/
│   ├── sequence/
│   ├── structure/
│   └── structure_grn/
│
├── _archive/                        # Deprecated/archived scripts
│   ├── gpcr_analysis/               # Old GPCR analysis (not used)
│   ├── old_scripts/                 # Old processor examples
│   ├── html_visualizations/         # Interactive HTML figures
│   └── tiger/                       # Tiger image for vision demo
│
└── logs/                            # Execution logs
```

## Figure Status

| Figure | Description | Status | Folder |
|--------|-------------|--------|--------|
| 1.1 | Type I/II structural comparison | ☐ not done | `ch1_introduction/fig_1.1_type_comparison/` |
| 1.2 | Spectral tuning: (A) Type I/II reps, (B) Retinal + key determinants | ☐ not done | `ch1_introduction/fig_1.2_spectral_tuning/` |
| 2.2 | Cone opsin diversity | ☑ done | `ch2_protos/fig_2.2_opsin_diversity/` |
| 2.3 | BR binding pocket | ☑ done | `ch2_protos/fig_2.3_br_binding_pocket/` |
| 2.4 | Binding pocket graphs | ☑ done | `ch2_protos/fig_2.4_pocket_graphs/` |
| 2.5 | GRN comparison | ☑ done | `ch2_protos/fig_2.5_grn_comparison/` |
| 2.6 | t-SNE + heatmap | ☑ done | `ch2_protos/fig_2.6_embeddings/` |
| 5.2 | Rhodozyme design | ☐ not done | `ch5_applications/fig_5.2_rhodozyme/` |

## PyMOL Usage

Always load shared settings at the start of PyMOL scripts:

```
@../../shared/thesis_pymol_settings.pml
```

Or from the protos root:
```
@thesis_pymol_settings.pml
```

### Key PyMOL Colors

| Name | Hex | Use |
|------|-----|-----|
| `type1_slate` | #3d5a80 | Type I microbial opsins |
| `type2_terracotta` | #c1666b | Type II animal opsins |
| `retinal_rust` | #bc6c25 | Retinal chromophore |
| `conserved_olive` | #6a994e | Functional positions |

### Export Commands

```
ray_pub        # 2400x1800 publication quality
ray_preview    # 800x600 quick preview
png figure.png, dpi=300
```

## Related Projects

- **lambda**: Figures 4.1-4.12 (spectral data, scatter plots)
- **MOGRN**: Figure 3.1 (MO numbering)
- **Protos_MCP**: Figures 6.1, 6.3 (vision examples)
- **Ghostwriter**: Diagrams 2.1, 4.2, 5.1, 6.2 (PowerPoint)
