# Chapter 5 — Rhodozyme Figures

All colors reference `colorscales.yaml`. Convention across all panels: **gray = from crystal structure (locked/known), color = from model (designed/generated)**.

---

## Global Color Assignments

| Role | Color | Hex | RGB norm | Source |
|------|-------|-----|----------|--------|
| Locked backbone (crystal structure) | Light Gray | `#adb5bd` | `[0.678, 0.710, 0.741]` | `figure_elements.backbone_unfocused` |
| Designed backbone (model output) | Terracotta | `#c1666b` | `[0.757, 0.400, 0.420]` | `structures.type_ii` |
| Retinal | Rust | `#bc6c25` | `[0.737, 0.424, 0.145]` | `molecules.retinal` |
| Substrate (succinyl-AAPR) | Ochre | `#d4a03c` | `[0.831, 0.627, 0.235]` | `sequences.long_wave` |
| Theozyme residues (Ser, His, Asp) | Conserved Green | `#6a994e` | `[0.416, 0.600, 0.306]` | `functional.conserved` |
| Theozyme geometry (triangle, vectors) | Dark Gray | `#343a40` | `[0.204, 0.227, 0.251]` | `neutral.dark_gray` |
| GRN candidate region | Sage | `#457b6b` | `[0.271, 0.482, 0.420]` | `sequences.short_wave` |
| Edges / dashed lines / H-bonds | Gray | `#6c757d` | `[0.424, 0.459, 0.490]` | `neutral.gray` |
| Dark-state rhodopsin (superposition) | Gray | `#6c757d` | `[0.424, 0.459, 0.490]` | `neutral.gray` |
| Labels / annotations | Dark Gray | `#343a40` | `[0.204, 0.227, 0.251]` | `figure_elements.annotation` |
| Sequence identity (alignment) | Light Gray | `#adb5bd` | `[0.678, 0.710, 0.741]` | `neutral.light_gray` |
| Sequence mutations (alignment) | Terracotta | `#c1666b` | `[0.757, 0.400, 0.420]` | `structures.type_ii` |
| Background | White | `#ffffff` | — | PyMOL `bg_color white` |

PyMOL color definitions (copy-paste):
```
set_color locked_gray, [0.678, 0.710, 0.741]
set_color designed_terracotta, [0.757, 0.400, 0.420]
set_color retinal_rust, [0.737, 0.424, 0.145]
set_color substrate_ochre, [0.831, 0.627, 0.235]
set_color theozyme_green, [0.416, 0.600, 0.306]
set_color geo_dark, [0.204, 0.227, 0.251]
set_color candidate_sage, [0.271, 0.482, 0.420]
set_color hbond_gray, [0.424, 0.459, 0.490]
set_color darkstate_gray, [0.424, 0.459, 0.490]
bg_color white
```

---

## Fig 5.1 — Input Structures and the Design Premise

**Two panels. The point: show the physical space that opens upon activation, and the chemistry to be transplanted.**

### Panel A — Rhodopsin dark vs. active superposition

| Component | Representation | Color |
|-----------|---------------|-------|
| Dark-state rhodopsin (inactive, e.g. 1U19 or modeled) | Cartoon | Gray `#6c757d` |
| Active-state rhodopsin (3PQR) | Cartoon | Terracotta `#c1666b` |
| Retinal (3PQR) | Sticks | Rust `#bc6c25` |
| Intracellular cavity (active state only) | Transparent surface (50% opacity) | Terracotta Light `#e4c1c1` |
| TM6 displacement arrow | Arrow annotation | Dark Gray `#343a40` |
| TM6 label, distance label (10–14 Å) | Text | Dark Gray `#343a40` |

**View**: Side view (membrane plane horizontal). Dark state behind, active state in front. The cavity surface at the intracellular face should be the visual focus — it is the "real estate" for the active site.

### Panel B — Trypsin catalytic triad

| Component | Representation | Color |
|-----------|---------------|-------|
| Trypsin backbone (2AGE) | Cartoon | Light Gray `#adb5bd` |
| Catalytic triad (Ser195, His57, Asp102) | Sticks | Theozyme Green `#6a994e` |
| Substrate (succinyl-AAPR) | Sticks | Ochre `#d4a03c` |
| H-bond distances (Ser–His, His–Asp) | Dashed lines + distance labels | Gray `#6c757d` |

**View**: Close-up of the active site. Backbone recedes (light gray), triad and substrate are the focus. Label each residue. Label H-bond distances in Å.

---

## Fig 5.2 — Theozyme Extraction

**Single panel. The point: show what is extracted from the enzyme and what is discarded.**

| Component | Representation | Color |
|-----------|---------------|-------|
| Trypsin backbone (2AGE) | Cartoon (faded) | Light Gray `#adb5bd` |
| Catalytic triad sidechains | Sticks | Theozyme Green `#6a994e` |
| Cα positions (Ser195, His57, Asp102) | Spheres | Dark Gray `#343a40` |
| Cα→Cβ direction vectors | Arrows (CGO) | Dark Gray `#343a40` |
| Pairwise Cα–Cα distances | Dashed lines + distance labels | Gray `#6c757d` |
| Residue labels | Text | Dark Gray `#343a40` |

**View**: Same orientation as Fig 5.1B for continuity. The trypsin backbone is present but faded — the geometric abstraction (spheres, arrows, triangle) is drawn on top. Caption states: "This triangle and these vectors are the input to the placement search — everything else about trypsin is discarded."

---

## Fig 5.3 — Theozyme Placement

**Two panels. The point: where on the rhodopsin does the theozyme land, and is there room for a substrate?**

### Panel A — Candidate region

| Component | Representation | Color |
|-----------|---------------|-------|
| Rhodopsin TM helices (non-candidate) | Cartoon | Light Gray `#adb5bd` |
| Candidate region (TM helix ends, ICL1–3, H8) | Cartoon | Sage `#457b6b` |
| Retinal | Sticks | Rust `#bc6c25` |
| Intracellular cavity | Transparent surface (50% opacity) | Sage Light `#a3c4b8` |

**View**: From cytoplasm, looking up into the intracellular face. TM helices form a ring; the candidate region is in the center and on the helix ends. The transparent surface shows the available volume.

### Panel B — Winning triplet placed

| Component | Representation | Color |
|-----------|---------------|-------|
| Rhodopsin backbone | Cartoon | Light Gray `#adb5bd` |
| Candidate region | Cartoon | Sage `#457b6b` |
| Matched Cα positions | Spheres | Theozyme Green `#6a994e` |
| Cβ direction vectors | Arrows (CGO) | Dark Gray `#343a40` |
| Theozyme triangle (Cα–Cα distances) | Dashed lines | Gray `#6c757d` |
| Cavity surface | Transparent surface | Sage Light `#a3c4b8` |
| Residue on TM5/6 label | Text + annotation | Dark Gray `#343a40` |

**View**: Same cytoplasmic view. The theozyme triangle sits inside the cavity. Annotate which position is on TM5 or TM6 — "in the dark state, this distance breaks."

---

## Fig 5.4 — RFdiffusion Mask and Output

**Two panels. The point: what is from the crystal structure (trusted) and what is from the model (evaluate critically). This is the most important figure for a skeptical reader.**

### Panel A — The mask

| Component | Representation | Color |
|-----------|---------------|-------|
| Locked regions: TM helices | Cartoon | Light Gray `#adb5bd` |
| Locked regions: theozyme Cα positions | Spheres | Theozyme Green `#6a994e` |
| Designable regions: ICL loops | Cartoon (dashed or distinct style) | Terracotta `#c1666b` |
| Retinal | Sticks | Rust `#bc6c25` |

**View**: Side view. The locked/designable boundary should be visually unambiguous. Caption opens with: "Gray regions are locked — they come directly from the crystal structure (3PQR) and are not modified."

### Panel B — Backbone designs

| Component | Representation | Color |
|-----------|---------------|-------|
| Locked scaffold (shared across designs) | Cartoon | Light Gray `#adb5bd` |
| Design 1 — loops | Cartoon | Terracotta `#c1666b` |
| Design 2 — loops | Cartoon | Terracotta Light `#e4c1c1` |
| Design 3 — loops | Cartoon | Ochre `#d4a03c` |
| Theozyme positions (all designs) | Spheres | Theozyme Green `#6a994e` |
| Retinal | Sticks | Rust `#bc6c25` |

**View**: Same orientation as panel A. Three designs superimposed on the same locked scaffold. Loop diversity is visible; the gray scaffold and green theozyme positions are identical in all three. The retinal pocket at the top is untouched.

---

## Fig 5.5 — Sequence Design

**Two panels. The point: how different is the designed sequence from wild-type rhodopsin? Mutations should cluster in the loops, not the TM helices.**

### Panel A — Sequence alignment (matplotlib or text figure)

| Component | Color | Hex |
|-----------|-------|-----|
| Identical positions (WT = designed) | Light Gray | `#adb5bd` |
| Mutated positions (designed ≠ WT) | Terracotta | `#c1666b` |
| Theozyme positions (Ser, His, Asp) | Conserved Green | `#6a994e` |
| TM helix region markers | Slate | `#3d5a80` |
| ICL region markers | Sage | `#457b6b` |
| Background | White | — |

**Format**: Horizontal alignment strip. One row = WT (3PQR), one row = designed sequence. Color each column by identity/mutation. Region annotations (TM1, ICL1, TM2, ...) above. The visual pattern should show gray blocks (TM helices, conserved) interrupted by terracotta blocks (loops, redesigned).

### Panel B — Sequence logo (matplotlib)

| Component | Color | Hex |
|-----------|-------|-----|
| Variable positions (low info content) | Terracotta | `#c1666b` |
| Conserved positions (high info content) | Sage | `#457b6b` |
| Fixed theozyme positions | Conserved Green | `#6a994e` |
| Axis, ticks, labels | Dark Gray | `#343a40` |

**Format**: Standard sequence logo for the designed region across 8 sampled sequences for one backbone. Theozyme positions annotated.

---

## Fig 5.6 — Boltz2 Evaluation

**Three panels. The point: does the predicted structure match the parent rhodopsin (fold) and the trypsin triad (catalysis)?**

### Panel A — Structural overlay

| Component | Representation | Color |
|-----------|---------------|-------|
| Parent rhodopsin (3PQR) | Cartoon | Light Gray `#adb5bd` |
| Predicted structure — locked regions | Cartoon | Light Gray `#adb5bd` |
| Predicted structure — designed loops | Cartoon | Terracotta `#c1666b` |
| Predicted theozyme residues | Sticks | Theozyme Green `#6a994e` |
| Retinal (predicted) | Sticks | Rust `#bc6c25` |
| Substrate (predicted) | Sticks | Ochre `#d4a03c` |
| RMSD label (overall backbone) | Text | Dark Gray `#343a40` |

**View**: Side view. The 3PQR reference and the prediction overlap. Designed loops in terracotta diverge slightly from the gray reference — the deviation is the model's contribution. Label overall backbone RMSD.

### Panel B — Catalytic geometry comparison

| Component | Representation | Color |
|-----------|---------------|-------|
| Reference trypsin triad (2AGE) | Sticks + Cα spheres | Light Gray `#adb5bd` |
| Predicted theozyme | Sticks + Cα spheres | Theozyme Green `#6a994e` |
| Cα–Cα distance lines (reference) | Dashed lines | Gray `#6c757d` |
| Cα–Cα distance lines (predicted) | Solid lines | Theozyme Green `#6a994e` |
| RMSD label (triad Cα) | Text | Dark Gray `#343a40` |

**View**: Close-up. Reference triad (gray) and predicted triad (green) superimposed after Kabsch alignment. Both triangles visible. RMSD value prominent.

### Panel C — pLDDT per residue (matplotlib)

| Component | Color | Hex |
|-----------|-------|-----|
| pLDDT — locked regions (from crystal structure) | Light Gray | `#adb5bd` |
| pLDDT — designed regions (from model) | Terracotta | `#c1666b` |
| pLDDT — theozyme positions | Conserved Green | `#6a994e` |
| Confidence threshold line (e.g. 70) | Dashed | Dark Gray `#343a40` |
| Region annotations (TM1, ICL1, ...) | Text | Dark Gray `#343a40` |
| Axis, ticks | Dark Gray | `#343a40` |

**Format**: Line or bar plot, residue number on x-axis, pLDDT on y-axis. Color each point/bar by its category (locked vs. designed vs. theozyme). Region boundaries annotated. The reader should see that locked regions have high confidence (expected — they come from a crystal structure) and evaluate whether designed regions also score well.

---

## Color Logic Summary — Chapter 5

The palette uses two visual channels:

1. **Gray ↔ Color** encodes **provenance**: gray means the coordinates come from an experimentally determined crystal structure; color means they were generated by a model (RFdiffusion, LigandMPNN, or Boltz2).

2. **Which color** encodes **role**:
   - Terracotta `#c1666b` — designed backbone (the model's contribution)
   - Theozyme Green `#6a994e` — catalytic residues (the constraint carried across all steps)
   - Rust `#bc6c25` — retinal (the light-sensitivity element, always locked)
   - Ochre `#d4a03c` — substrate (the target of catalysis)
   - Sage `#457b6b` — GRN candidate region (the search space, Step 3 only)

This convention is consistent across all six figures. A reader who understands it for one figure can read the rest without re-learning the color code.
