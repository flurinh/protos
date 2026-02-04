# ProtOS Figures

This project generates PyMOL structure visualizations, binding pocket graphs, and related figures.

**Color Reference:** `colorscales.yaml`

---

## Checklist

| Fig | Description | Status |
|-----|-------------|--------|
| 1.1 | Type I/II structural comparison | ☐ not done |
| 2.2 | Cone opsin diversity (boxplot + dendrogram) | ☑ done |
| 2.3 | BR binding pocket | ☑ done |
| 2.4 | Binding pocket graphs (3 images) | ☑ done |
| 2.5 | Bovine/squid GRN comparison | ☑ done |
| 2.6 | t-SNE + heatmap (2 images) | ☑ done |
| 5.2 | Rhodozyme design candidates | ☐ not done |

**Summary:** 5 done / 2 not done

---

## Figures to Generate

### Figure 1.1 - Type I/II Structural Comparison (Chapter 1)
**Priority:** HIGH | **Status:** TODO

**Components:**
- (A) Bacteriorhodopsin 1C3W ribbon
- (B) Bovine rhodopsin 1U19 ribbon
- (C) Binding pocket overlay with retinal + counterion

| Element | Color Key | Hex | PyMOL |
|---------|-----------|-----|-------|
| Type I ribbon (1C3W) | `structures.type_i` | #3d5a80 | `set_color type1_slate, [0.239, 0.353, 0.502]` |
| Type II ribbon (1U19) | `structures.type_ii` | #c1666b | `set_color type2_terracotta, [0.757, 0.400, 0.420]` |
| Unfocused backbone | `figure_elements.backbone_unfocused` | #adb5bd | `set_color backbone_gray, [0.678, 0.710, 0.741]` |
| Retinal | `molecules.retinal` | #bc6c25 | `set_color retinal_rust, [0.737, 0.424, 0.145]` |
| Counterion | `functional.conserved` | #6a994e | `set_color conserved_olive, [0.416, 0.600, 0.306]` |

---

### Figure 2.2 - Cone Opsin Diversity (Chapter 2)
**Priority:** PRESENT | **Status:** Update colors

**Components:**
- Boxplot showing sequence diversity
- Dendrogram showing clustering

| Element | Color Key | Hex |
|---------|-----------|-----|
| Short-wave opsins | `sequences.short_wave` | #457b6b |
| Long-wave opsins | `sequences.long_wave` | #d4a03c |
| Dendrogram branches | `neutral.gray` | #6c757d |

---

### Figure 2.3 - BR Binding Pocket (Chapter 2)
**Priority:** PRESENT | **Status:** Update colors

**Components:**
- Bacteriorhodopsin structure
- Binding pocket close-up with retinal

| Element | Color Key | Hex | PyMOL |
|---------|-----------|-----|-------|
| BR protein | `structures.type_i` | #3d5a80 | `set_color type1_slate, [0.239, 0.353, 0.502]` |
| Retinal | `molecules.retinal` | #bc6c25 | `set_color retinal_rust, [0.737, 0.424, 0.145]` |
| Binding pocket residues | `neutral.gray` | #6c757d | `set_color pocket_gray, [0.424, 0.459, 0.490]` |

---

### Figure 2.4 - Binding Pocket Graphs (Chapter 2)
**Priority:** PRESENT | **Status:** Update colors

**Components:**
- 1C3W binding pocket graph
- 1U19 binding pocket graph
- Comparison overlay

| Element | Color Key | Hex |
|---------|-----------|-----|
| Type I graph nodes | `graphs.type_i` | #3d5a80 |
| Type II graph nodes | `graphs.type_ii` | #c1666b |
| Graph edges | `graphs.edges` | #6c757d |

---

### Figure 2.5 - Bovine/Squid GRN Comparison (Chapter 2)
**Priority:** PRESENT | **Status:** Update colors

**Components:**
- Bovine rhodopsin structure
- Squid rhodopsin structure
- GRN position annotations

| Element | Color Key | Hex | PyMOL |
|---------|-----------|-----|-------|
| Bovine rhodopsin | `structures.type_ii` | #c1666b | `set_color type2_terracotta, [0.757, 0.400, 0.420]` |
| Squid rhodopsin | `structures.type_ii_light` | #e4c1c1 | `set_color type2_light, [0.894, 0.757, 0.757]` |
| Conserved positions | `functional.conserved` | #6a994e | `set_color conserved_olive, [0.416, 0.600, 0.306]` |

---

### Figure 2.6 - t-SNE + Heatmap (Chapter 2)
**Priority:** PRESENT | **Status:** Update colors

**Components:**
- t-SNE projection of cone opsin embeddings
- Similarity heatmap

| Element | Color Key | Hex |
|---------|-----------|-----|
| t-SNE points (group A) | `embeddings.group_a` | #457b6b |
| t-SNE points (group B) | `embeddings.group_b` | #d4a03c |
| Heatmap | `continuous.diverging_type` | #3d5a80 ↔ #f7f5f3 ↔ #c1666b |

---

### Figure 5.2 - Rhodozyme Design Candidates (Chapter 5)
**Priority:** LOW (work in progress) | **Status:** TODO

**Components:**
- RFdiffusion2 outputs
- Predicted structures
- LAMBDA predictions

| Element | Color Key | Hex |
|---------|-----------|-----|
| TBD | — | — |

---

## PyMOL Standard Settings

**Master settings file:** `thesis_pymol_settings.pml`

Load at the start of every PyMOL session:
```
@thesis_pymol_settings.pml
```

### Defined Colors

| Name | Hex | RGB norm | Use |
|------|-----|----------|-----|
| `type1_slate` | #3d5a80 | [0.239, 0.353, 0.502] | Type I structures |
| `type1_light` | #98c1d9 | [0.596, 0.757, 0.851] | Type I secondary |
| `type2_terracotta` | #c1666b | [0.757, 0.400, 0.420] | Type II structures |
| `type2_light` | #e4c1c1 | [0.894, 0.757, 0.757] | Type II secondary |
| `retinal_rust` | #bc6c25 | [0.737, 0.424, 0.145] | Retinal (11-cis) |
| `retinal_light` | #ddb892 | [0.867, 0.722, 0.573] | Retinal (all-trans) |
| `retinal_deprot` | #7a5980 | [0.478, 0.349, 0.502] | Deprotonated retinal |
| `conserved_olive` | #6a994e | [0.416, 0.600, 0.306] | Conserved/functional |
| `neutral_gray` | #6c757d | [0.424, 0.459, 0.490] | Binding pocket |
| `neutral_light` | #adb5bd | [0.678, 0.710, 0.741] | Unfocused backbone |

### Rendering Settings (from thesis_pymol_settings.pml)

```
# Background
bg_color white

# Quality
set antialias, 2
set ray_shadows, 0
set orthoscopic, 1
set depth_cue, 0

# Cartoon
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_transparency, 0.0

# Sticks/Spheres
set stick_radius, 0.15
set sphere_scale, 0.25

# Labels
set label_size, 14
set label_color, neutral_dark
```

### Quick Commands (aliases)

| Command | Action |
|---------|--------|
| `ray_preview` | Quick 800x600 render |
| `ray_pub` | Publication 2400x1800 render |
| `ray_poster` | High-res 4800x3600 render |
| `show_retinal` | Highlight retinal as sticks |
| `show_pocket` | Show binding pocket residues |
| `cartoon_trans` | Set cartoon transparency to 0.7 |

### Standard Export Workflow

```
# After setting up view
ray_pub
png figure_name.png, dpi=300
```

---

## Quick Reference

| Data Type | Primary | Light | PyMOL color name |
|-----------|---------|-------|------------------|
| Type I | #3d5a80 | #98c1d9 | `type1_slate` / `type1_light` |
| Type II | #c1666b | #e4c1c1 | `type2_terracotta` / `type2_light` |
| Retinal | #bc6c25 | #ddb892 | `retinal_rust` / `retinal_light` |
| Conserved | #6a994e | #a7c494 | `conserved_olive` / `conserved_light` |
| Neutral | #6c757d | #adb5bd | `neutral_gray` / `neutral_light` |

---

## Existing PyMOL Scripts

| Script | Purpose | Status |
|--------|---------|--------|
| `thesis/processors/figures/colorscales_pymol.pml` | Old color definitions | **Deprecated** - use `thesis_pymol_settings.pml` |
| `thesis/processors/figures/structure_type1_vs_type2.pml` | Type I/II comparison | Needs color update |
| `thesis/processors/figures/grn_structural_alignment.pml` | GRN alignment | Needs color update |
| `thesis/processors/figures/grn_spectral_tuning.pml` | Spectral tuning sites | Needs color update |
