# Figure Coloring Reference

All colors reference `colorscales.yaml`. All scripts load colors through shared infrastructure — never hardcoded.

**Shared files:**
- `shared/colorscales_pymol.pml` — PyMOL color + render settings (all chapters)
- `shared/thesis_style.py` — matplotlib rcParams + plotly layout defaults + color constants
- `ch5_applications/fig_5.2_rhodozyme/scripts/rhodozyme_colors.pml` — Ch5 superset (includes TM spectral colors)

---

## Fig 2.2 — Opsin Diversity (Sequence Processor)

**Palette: Sequences (Sage/Ochre)**

| Component | Color | Hex | Source |
|-----------|-------|-----|--------|
| SW opsins (donut, boxplot, scatter, tree leaves) | Sage | `#457b6b` | `TYPE_COLORS["short_wave"]` |
| LW opsins (donut, boxplot, scatter, tree leaves) | Ochre | `#d4a03c` | `TYPE_COLORS["long_wave"]` |
| Dendrogram links | Gray | `#6c757d` | `GRAY` |
| Bar edges | Dark Gray | `#343a40` | `DARK_GRAY` |

**Files**: `sequence_opsin_diversity.png`, `sequence_blast_scatter.png`, `sequence_phylogenetic_tree.png`

---

## Fig 2.3 — BR Binding Pocket (Structure Processor)

**Palette: Structures (Slate/Terracotta) + Molecules (Rust)**

| Component | Color | Hex | Source |
|-----------|-------|-----|--------|
| 1C3W backbone (Bacteriorhodopsin) | Slate | `#3d5a80` | `COLORS["structures"]["1C3W"]` |
| 3UG9 backbone (Channelrhodopsin) | Slate Light | `#98c1d9` | `COLORS["structures"]["3UG9"]` |
| 1U19 backbone (Bovine Rhodopsin) | Terracotta | `#c1666b` | `COLORS["structures"]["1U19"]` |
| 2Z73 backbone (Squid Rhodopsin) | Terracotta Light | `#e4c1c1` | `COLORS["structures"]["2Z73"]` |
| Retinal (all structures) | Rust | `#bc6c25` | `COLORS["molecules"]["retinal"]` |

**Files**: `opsin_alignment.html`, 4x `*_aligned.cif`

---

## Fig 2.4 — Pocket Graphs (Graph Processor)

**Palette: Structures (Slate/Terracotta) + Molecules (Rust) + Graphs (Gray)**

| Component | Color | Hex | Source |
|-----------|-------|-----|--------|
| 1C3W residue nodes + bar | Slate | `#3d5a80` | `COLORS["structures"]["1C3W"]` |
| 1U19 residue nodes + bar | Terracotta | `#c1666b` | `COLORS["structures"]["1U19"]` |
| Graph edges | Gray | `#6c757d` | `COLORS["graphs"]["edges"]` |
| Retinal atoms + bonds | Rust | `#bc6c25` | `COLORS["molecules"]["retinal"]` |
| Bar edges + labels | Dark Gray | `#343a40` | `DARK_GRAY` |

**Files**: `graph_binding_pocket_3d_1C3W.html`, `graph_binding_pocket_3d_1U19.html`, `graph_comparison.png`

---

## Fig 2.5 — GRN Comparison (GRN Processor)

**Palette: Structures (Terracotta variants) + Molecules (Rust)**

Both structures are Type II, so both use terracotta variants. PyMOL scripts source `shared/colorscales_pymol.pml`.

| Component | PyMOL color name | Hex | Source |
|-----------|-----------------|-----|--------|
| Bovine (1U19) cartoon + sticks | `color_1U19` | `#c1666b` | `colorscales_pymol.pml` |
| Squid (2Z73) cartoon + sticks | `color_2Z73` | `#e4c1c1` | `colorscales_pymol.pml` |
| Retinal sticks | `retinal_rust` | `#bc6c25` | `colorscales_pymol.pml` |
| Labels | `dark_gray` | `#343a40` | `colorscales_pymol.pml` |

**Files**: `grn_structural_alignment.pml`, `grn_spectral_tuning.pml`, `1u19_grn_aligned.cif`, `2z73_grn_aligned.cif`

---

## Fig 2.6 — Embeddings (Embedding Processor)

**Palette: Sequences/Embeddings (Sage/Ochre) + Neutral**

| Component | Color | Hex | Source |
|-----------|-------|-----|--------|
| SW points (t-SNE scatter) | Sage | `#457b6b` | `TYPE_COLORS["short_wave"]` |
| LW points (t-SNE scatter) | Ochre | `#d4a03c` | `TYPE_COLORS["long_wave"]` |
| "Same" box (similarity analysis) | Sage | `#457b6b` | `SAGE` |
| "Diff" box (similarity analysis) | Ochre | `#d4a03c` | `OCHRE` |
| Histogram bars | Gray | `#6c757d` | `GRAY` |
| Heatmap colorscale | RdBu | — | plotly built-in |
| Heatmap SW/LW tick marks | Sage/Ochre | — | `TYPE_COLORS` |

**Files**: `embedding_tsne.png`, `embedding_similarity_analysis.png`, `embedding_heatmap.png`

---

## Color Logic Summary

| Data domain | Palette | Primary | Secondary |
|-------------|---------|---------|-----------|
| Structures / Graphs | Slate / Terracotta | `#3d5a80` | `#c1666b` |
| Sequences / Embeddings | Sage / Ochre | `#457b6b` | `#d4a03c` |
| Molecules (retinal) | Rust | `#bc6c25` | — |
| Edges / neutral | Gray | `#6c757d` | `#adb5bd` |
| Per-structure (light variants) | 3UG9: `#98c1d9`, 2Z73: `#e4c1c1` | | |

## Matplotlib Style (all figures)

Applied via `thesis_style.apply_style()`:
- White background, no top/right spines
- Axis/tick/label color: Dark Gray (`#343a40`)
- Grid: `#e9ecef` when enabled
- Font: sans-serif, 10pt base
- DPI: 300 (save), 150 (display)
