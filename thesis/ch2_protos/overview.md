# Chapter 2 — ProtOS: Figures Overview

This document describes all figures in Chapter 2, which demonstrates the ProtOS
framework through six figures. Each figure pair follows the same pattern:
a didactic code example showing the ProtOS API, followed by a scaled-up
application to real research data from the LAMBDA opsin atlas (27,640 sequences,
13 Type II opsin subfamilies).

---

## Figure 2.1 — ProtOS Architecture

**Location**: `fig_2.1_architecture/` (diagram, not script-generated)

**Content**: System architecture diagram showing the ProtOS data management
framework. Central database surrounded by six processor modules
(SequenceProcessor, StructureProcessor, GraphProcessor, GRNProcessor,
EmbeddingProcessor, PropertyProcessor), external data sources (PDB, UniProt,
NCBI), and downstream analysis pipelines.

**Key message**: ProtOS provides a unified interface for multi-modal protein data
— sequences, structures, graphs, annotations, and embeddings — under a single
entity-centric data model.

**Figure caption**: *ProtOS system architecture. Protein entities are registered
in a central database and accessed through specialized processors that handle
sequences, 3D structures, residue contact graphs, generic residue number (GRN)
annotations, and embeddings. External databases (PDB, UniProt, NCBI) are
accessed through loader modules. All processors share a common entity namespace,
enabling seamless cross-modal analysis.*

---

## Figure 2.2 — Opsin Sequence Diversity

**Location**: `fig_2.2_opsin_diversity/`

### 2.2 Code example: `sequence_processor_example.py`

Demonstrates the SequenceProcessor API by building a small opsin dataset:
loading a reference sequence from UniProt, running a BLAST homology search via
NCBILoader, and registering results as a ProtOS dataset. This is the code
snippet shown in the thesis text.

### 2.2 Atlas figure: `atlas_sequence_diversity.py` → `atlas_identity_ridge.png`

**Content**: Per-subfamily sequence identity distributions for the full Type II
opsin atlas. A ridge plot (joy plot) showing one KDE curve per subfamily,
vertically stacked in descending order of atlas representation. Each curve is
colored by subfamily and labeled with sequence count.

**Data source**: LAMBDA atlas property table (27,640 sequences, 13 subfamilies).

**Key message**: ProtOS scales seamlessly from a single BLAST search to managing
a 27,640-sequence atlas spanning 13 subfamilies. The ridge plot reveals that
sequence identity distributions vary substantially across subfamilies — rod
opsins form a tight cluster while non-visual opsins are more divergent.

**Figure caption**: *Type II opsin sequence diversity across 13 subfamilies.
Kernel density estimates of pairwise sequence identity (%) to each subfamily's
query sequence. Rod opsins (n = 12,543) show a narrow identity distribution
centered at ~75%, reflecting high conservation among vertebrate scotopic
receptors. Cone opsins (LWS, SWS1, MWS) show broader distributions. Non-visual
opsins (melanopsin, neuropsin, encephalopsin) are more divergent with lower
median identity. Data: LAMBDA Type II opsin atlas, constructed using ProtOS
SequenceProcessor.*

---

## Figure 2.3 — Retinal Binding Pocket Structure

**Location**: `fig_2.3_br_binding_pocket/`

### 2.3: `structure_processor_example.py` → PyMOL render

**Content**: Structural comparison of Type I (microbial) and Type II (animal)
opsins aligned on the retinal binding site. Shows bacteriorhodopsin (1C3W) and
bovine rhodopsin (1U19) superimposed around the retinal chromophore, revealing
their different 7TM topologies despite converging on the same ligand-binding
function.

**Key message**: ProtOS StructureProcessor handles loading, filtering, ligand
extraction, and structural alignment — the building blocks for any structure-
based analysis. Type I and Type II opsins bind the same retinal chromophore
in topologically distinct folds, motivating fold-agnostic representations
(Chapter 4).

**Figure caption**: *Structural comparison of microbial and animal opsins around
the retinal binding site. Bacteriorhodopsin (1C3W, blue) and bovine rhodopsin
(1U19, terracotta) aligned on retinal atom positions (rust). Despite binding
the same chromophore, the two protein families adopt different seven-
transmembrane topologies — Type I threads N-terminus extracellularly while
Type II threads it intracellularly. Structures loaded and aligned using ProtOS
StructureProcessor.*

---

## Figure 2.4 — Binding Pocket Graph

**Location**: `fig_2.4_pocket_graphs/`

### 2.4: `graph_processor_example.py` → `1u19_graph_pocket.png`

**Content**: 3D interactive visualization of the residue contact graph around
the retinal binding pocket of bovine rhodopsin (1U19). Nodes represent pocket
residues (within 7.0 Å of retinal); edges connect residue pairs whose
Cα atoms are within 6.0 Å. The retinal chromophore is shown as a stick model
at the center of the graph.

**Key message**: ProtOS GraphProcessor converts a 3D binding pocket into a graph
representation — the data format consumed by downstream machine learning models.
This graph captures the spatial arrangement of pocket residues without requiring
explicit sequence alignment, enabling cross-family comparisons (Chapter 4).

**Figure caption**: *Residue contact graph of the retinal binding pocket in
bovine rhodopsin (1U19). Nodes (terracotta spheres) are protein residues within
7.0 Å of the retinal chromophore (rust sticks). Edges (gray lines) connect
residue pairs with Cα–Cα distance < 6.0 Å. This graph representation
encodes the three-dimensional architecture of the binding pocket in a format
amenable to graph neural networks (Chapter 4). Generated by ProtOS
GraphProcessor.*

---

## Figure 2.5 — Generic Residue Numbers (GRN)

**Location**: `fig_2.5_grn_comparison/`

### 2.5 Code example: `grn_processor_example.py`

Demonstrates the GRN annotation workflow: extracting sequences from crystal
structures, annotating with Ballesteros–Weinstein generic residue numbers via
`annotate_with_grn()`, comparing binding pocket positions across species, and
generating PyMOL alignment scripts. Uses bovine (1U19) and squid (2Z73)
rhodopsin as examples.

### 2.5a: `grn_alignment_block.py` → `grn_alignment_block.png`

**Content**: Table of amino acid identities at key GRN positions across the 10
Type II opsin reference sequences used to build the LAMBDA atlas. Rows are
reference opsins (one per subfamily query); columns are 15 functionally important
GRN positions grouped into spectral tuning sites (3.28, 3.32, 3.36) and
microswitch positions (E/DRY, PIF, Pro kink, CWxP, F-switch, Schiff base,
NPxxY). Cells show the one-letter amino acid code, background-colored by
physicochemical property: positive (slate), negative (terracotta), aromatic
(ochre), polar (sage), nonpolar (light gray). Row labels include a subfamily-
colored dot. TM helix brackets are shown below the table.

**Data source**: GRN annotations from LAMBDA atlas construction
(`opsin_phylogeny/type_ii_grns.csv`, 27,752 sequences). The 10 reference
sequences are the query accessions used in `step1_collect_opsins.py`.

**Key message**: GRN provides a structure-independent coordinate system that
transfers across all Type II opsins by sequence alignment alone. Key functional
positions (Schiff base K at 7.43, counterion E at 3.28, NPxxY motif) are fully
conserved across all 10 references, while spectral tuning positions (3.32, 3.36)
vary — reflecting the molecular basis of wavelength sensitivity. The amino acid
property coloring reveals that microswitches maintain strict physicochemical
constraints (e.g., the E/DRY arginine is always positively charged).

**Figure caption**: *GRN microswitch and spectral tuning table across 10 atlas
reference sequences. Columns correspond to 15 Ballesteros–Weinstein positions
spanning TM helices 3, 5, 6, and 7, divided into spectral tuning sites (amber
labels) and microswitch positions. Amino acids are colored by physicochemical
property: positive (slate), negative (terracotta), aromatic (ochre), polar
(sage), nonpolar (gray). The Schiff base lysine (K, 7.43), E/DRY arginine (R,
3.50), and NPxxY motif (N/P/Y, 7.49–7.53) are fully conserved. Spectral tuning
positions (3.32, 3.36) show subfamily-specific variation. GRN annotations
generated during LAMBDA atlas construction using ProtOS GRNProcessor.*

### 2.5b: `grn_microswitches.pml` → PyMOL render

**Content**: PyMOL visualization of the GPCR signal transduction microswitch
network on bovine rhodopsin (1U19). Seven microswitch positions are highlighted
as colored sticks on a transparent cartoon backbone, with GRN labels and dashed
lines indicating functional interactions. The pathway runs from the retinal
chromophore (extracellular side) through the transmembrane bundle to the
cytoplasmic G-protein interface.

**Microswitch positions**:
| GRN | Residue | Function | Color |
|-----|---------|----------|-------|
| 7.42 | K296 | Schiff base (retinal anchor) | Rust |
| 3.28 | E113 | Counterion (Schiff base stabilization) | Terracotta |
| 6.48 | W265 | CWxP toggle switch | Green |
| 6.44 | F261 | Transmission switch | Ochre |
| 5.50/3.40 | P215/I123 | PIF connector | Sage |
| 3.49–3.51 | E134–R135–Y136 | E/DRY ionic lock | Slate |
| 7.49–7.53 | N302–P303–Y306 | NPxxY activation motif | Mauve |

**Key message**: GRN positions are not merely an indexing convenience — they
label functionally critical residues involved in GPCR signal transduction. The
microswitch network shows how conformational changes propagate from retinal
isomerization through the toggle switch and transmission switch to the
cytoplasmic ionic lock. These positions are the basis for structure–function
analysis in Chapter 4.

**Figure caption**: *Microswitch signal transduction network in bovine rhodopsin
(1U19), annotated with Ballesteros–Weinstein generic residue numbers. Upon
photoisomerization of retinal (rust), the Schiff base K296 (7.42) and
counterion E113 (3.28) trigger conformational changes that propagate through the
CWxP toggle switch W265 (6.48, green) and transmission switch F261 (6.44,
ochre). The PIF connector (sage) links helices 3, 5, and 6. On the cytoplasmic
side, the E/DRY ionic lock (E134–R135–Y136, slate) and NPxxY motif (N302–P303–
Y306, mauve) gate G-protein coupling. Dashed lines indicate key inter-residue
contacts. These conserved microswitch positions, identified by GRN, form the
mechanistic basis for spectral tuning analysis in Chapter 4.*

---

## Figure 2.6 — Opsin Embedding Space

**Location**: `fig_2.6_embeddings/`

### 2.6 Atlas figure: `atlas_embedding_space.py` → `atlas_umap_embedding.png`

**Content**: UMAP 2D projection of 27,639 mean-pooled Ankh Large protein
language model embeddings, colored by subfamily. Each point represents one opsin
sequence from the atlas, positioned by its embedding similarity to all others.
The 9 query/reference sequences used to build the atlas are highlighted as star
markers with dark outlines, each colored by its subfamily.

**Data source**: Ankh Large per-residue embeddings (1,536 dimensions), mean-
pooled to sequence-level vectors. UMAP with n_neighbors=30, min_dist=0.3,
cosine metric. Query status from atlas property table (`is_query` column).

**Key message**: Protein language model embeddings capture both evolutionary and
functional relationships without requiring alignment. The UMAP projection
reveals clear subfamily clustering — rod opsins form the largest dense cluster,
cone opsins separate by spectral class, and non-visual opsins occupy distinct
regions. The query reference sequences (stars) sit at or near the center of
their respective subfamily clusters, confirming they are representative
exemplars. This learned representation is the foundation for the embedding-based
analyses in Chapter 4.

**Figure caption**: *UMAP projection of Ankh Large protein language model
embeddings for 27,639 Type II opsins. Each point represents one sequence,
colored by opsin subfamily; star markers denote the 9 query sequences used to
construct the atlas. Rod opsins (blue, n = 12,543) form a dense central cluster.
Cone opsins separate into spectral classes: LWS (red), SWS1 (purple), and MWS
(green). Non-visual opsins (melanopsin, neuropsin, encephalopsin) occupy distinct
peripheral regions. Query sequences cluster centrally within their subfamilies,
validating them as representative exemplars. The embedding space captures
functional relationships without sequence alignment. Embeddings computed using
ProtOS EmbeddingProcessor.*

---

## Color Palette

All figures use the thesis color scheme defined in `thesis/colorscales.yaml`:

- **Subfamily colors**: 13 spectrally motivated colors for the opsin atlas
  (defined in `colorscales.yaml` → `subfamilies` section)
- **Structure colors**: Slate (Type I) / Terracotta (Type II)
- **Molecule colors**: Rust (retinal) / Mauve (carazolol)
- **Functional**: Green (conserved positions)
- **Neutral**: Dark gray (text), Light gray (backbone), Warm white (background)

Shared styling: `thesis/shared/thesis_style.py` (matplotlib rcParams) and
`thesis/shared/colorscales_pymol.pml` (PyMOL colors and render settings).

---

## Reproduction

All programmatic figures can be regenerated from the `thesis/` directory:

```bash
# Fig 2.2 — Atlas sequence diversity
python ch2_protos/fig_2.2_opsin_diversity/atlas_sequence_diversity.py

# Fig 2.4 — Binding pocket graph
python ch2_protos/fig_2.4_pocket_graphs/graph_processor_example.py

# Fig 2.5a — GRN alignment block
python ch2_protos/fig_2.5_grn_comparison/grn_alignment_block.py

# Fig 2.5b — Microswitches (open in PyMOL)
pymol ch2_protos/fig_2.5_grn_comparison/grn_microswitches.pml

# Fig 2.6 — Atlas embedding space
python ch2_protos/fig_2.6_embeddings/atlas_embedding_space.py
```

Figures 2.1, 2.3, and 2.5b require manual rendering (diagram tool / PyMOL).
