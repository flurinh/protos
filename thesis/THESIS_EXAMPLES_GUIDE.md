# ProtOS Thesis Examples - Comprehensive Guide

**Last Updated:** 2026-01-12
**Status:** 8 Processor Examples + 4 Advanced Workflows = COMPLETE

This document serves as the complete reference for all thesis examples, their outputs,
figures, and the story they tell about ProtOS capabilities.

---

## The Big Picture: Information Flow

The examples demonstrate how data flows through ProtOS processors to answer increasingly
complex biological questions. Each processor transforms or enriches data, and outputs
from one processor become inputs to another.

```
                                    ┌─────────────────────────────────┐
                                    │     ADVANCED WORKFLOWS          │
                                    │  (Multi-processor pipelines)    │
                                    └─────────────────────────────────┘
                                               ▲
         ┌─────────────────────────────────────┼─────────────────────────────────────┐
         │                                     │                                     │
         ▼                                     ▼                                     ▼
┌─────────────────┐              ┌─────────────────┐              ┌─────────────────┐
│  PropertyProc   │◄────────────►│  EmbeddingProc  │◄────────────►│   GraphProc     │
│  (annotations)  │              │  (pLM vectors)  │              │  (contacts)     │
└────────┬────────┘              └────────┬────────┘              └────────┬────────┘
         │                                │                                │
         │         ┌──────────────────────┼──────────────────────┐        │
         │         │                      │                      │        │
         ▼         ▼                      ▼                      ▼        ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              GRNProcessor                                        │
│                     (Generic Residue Numbers - cross-family mapping)             │
└─────────────────────────────────────────────────────────────────────────────────┘
                                          ▲
                    ┌─────────────────────┼─────────────────────┐
                    │                     │                     │
                    ▼                     ▼                     ▼
         ┌─────────────────┐   ┌─────────────────┐   ┌─────────────────┐
         │  SequenceProc   │   │  StructureProc  │   │  MoleculeProc   │
         │  (FASTA, BLAST) │   │  (PDB, mmCIF)   │   │  (SMILES, FP)   │
         └─────────────────┘   └─────────────────┘   └─────────────────┘
                    ▲                     ▲                     ▲
                    │                     │                     │
                    └─────────────────────┼─────────────────────┘
                                          │
                              ┌───────────────────────┐
                              │     ModelManager      │
                              │ (External compute:    │
                              │  Boltz2, ESM2, etc.)  │
                              └───────────────────────┘
```

**Key Insight:** Understanding at one level enables applications at another:
- Sequence similarity → Embedding clustering → Functional prediction
- Structure analysis → GRN annotation → Cross-family comparison
- Mechanism understanding → Ligand design → Drug discovery

---

## Story Thread

The examples follow a logical progression:

1. **Basic Data Processors** - Load, store, and manage fundamental biological entities
   - Sequence, Structure, Molecule

2. **Derived Data Processors** - Build on basic data with specialized annotations
   - GRN (requires Sequence), Property (links to entities), Embedding (requires Sequence), Graph (requires Structure)

3. **ModelManager** - Orchestrate computational jobs for external compute

4. **Advanced Workflows** - Combine processors to answer real biological questions:
   - **GPCR Mechanism** → Understanding enables drug design
   - **pLM-Enriched Graphs** → Data structures for AI models
   - **Prey Vision Enhancement** → Evolution of color vision
   - **Light-Controlled Chemistry** → Spatiotemporal reaction control

---

## Folder Structure

```
thesis/
├── THESIS_EXAMPLES_GUIDE.md    # This file
├── processors/                  # Single-processor examples (8 total)
│   ├── figures/                 # Generated figures (PNG, HTML)
│   ├── sequence_processor_example.py
│   ├── structure_processor_example.py
│   ├── molecule_processor_example.py
│   ├── grn_processor_example.py
│   ├── property_processor_example.py
│   ├── embedding_processor_example.py
│   ├── graph_processor_example.py
│   └── model_manager_example.py
│
├── workflows/                   # Multi-processor workflows (4 advanced)
│   ├── figures/                 # Workflow figures + PyMOL scripts
│   ├── gpcr_binding_pocket_workflow.py
│   ├── plm_graph_workflow.py
│   ├── redshift_mutation_workflow.py
│   ├── rhodozyme_design_workflow.py
│   └── pymol_rhodozyme_visualization.py
│
└── outputs/                     # Data outputs by processor/workflow
```

---

# PART I: BASIC DATA PROCESSORS

These processors handle the fundamental biological data types.

---

## 1. SequenceProcessor
**File:** `processors/sequence_processor_example.py`

**Question:** "What is the sequence diversity across cone opsin spectral types?"

**Connection to Other Processors:** This example creates the `cone_opsin_diversity` dataset
(200 sequences) that feeds directly into the EmbeddingProcessor example, demonstrating
cross-processor data flow. The same dataset is used by PropertyProcessor to associate
experimental spectral measurements.

**What It Does:**
- Downloads human cone opsin query sequences (SW and LW)
- Runs NCBI BLAST searches for each spectral type (100 hits each)
- Fetches 200 homolog sequences and annotates by opsin type
- Computes all-vs-all sequence similarity using MMseqs2
- Creates phylogenetic tree showing spectral type clustering

**Model System:** 2 query opsins → 200 total sequences
| Gene | Type | λmax | Hits |
|------|------|------|------|
| OPN1SW | short_wave | 420 nm | 100 |
| OPN1LW | long_wave | 560 nm | 100 |

**Note:** Medium wave (OPN1MW) was removed because LW and MW are difficult to distinguish
spectrally and phylogenetically - this is central to the Prey Vision workflow!

**ProtOS Capabilities:**
- `SequenceLoader.download_and_register()` - UniProt download
- `NCBILoader.blast_search()` - Remote BLAST
- `NCBILoader.download_batch()` - Batch sequence fetch
- `mmseqs2_align2()` - All-vs-all sequence alignment
- Automatic entity registration and dataset management

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `sequence_opsin_diversity.png` | Pie chart by type + identity box plots |
| `sequence_blast_scatter.png` | Identity vs E-value by opsin type |
| `sequence_phylogenetic_tree.png` | Dendrogram colored by opsin type |

**Key Results:**
- 100 sequences per opsin type from SwissProt
- Clear separation of SW (blue) and LW (red) clades
- Dataset annotated with opsin_type for EmbeddingProcessor

**Data Flow:**
```
SequenceProcessor → cone_opsin_diversity dataset
                         ↓
              EmbeddingProcessor (clustering)
                         ↓
              PropertyProcessor (spectral values)
```

---

## 2. StructureProcessor
**File:** `processors/structure_processor_example.py`

**Question:** "How do microbial and animal opsins compare structurally around retinal?"

**What It Does:**
- Downloads opsin structures from PDB (microbial + animal)
- Analyzes retinal binding pockets
- Aligns structures using Kabsch algorithm
- Exports aligned structures to CIF

**Model System:** 4 structures
| PDB | Name | Type |
|-----|------|------|
| 1C3W | Bacteriorhodopsin | microbial |
| 3UG9 | Channelrhodopsin-2 | microbial |
| 1U19 | Bovine Rhodopsin | animal |
| 4ZWJ | Squid Rhodopsin | animal |

**ProtOS Capabilities:**
- `StructureLoader.download_batch()` - PDB download
- `StructureProcessor.get_ligand_interactions()` - Binding pocket analysis
- `write_cif_file()` - CIF export

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `opsin_alignment.html` | Interactive 3D alignment visualization |

**Key Results:**
- RMSD: 1C3W=0.0 (ref), 3UG9=2.8A, 1U19=3.1A, 4ZWJ=3.5A
- Retinal binding pocket conserved despite sequence divergence

---

## 3. MoleculeProcessor
**File:** `processors/molecule_processor_example.py`

**Question:** "Which compounds are similar to carazolol (the inverse agonist in 2RH1)?"

**Connection to GPCR Workflow:** Carazolol is the crystallographic inverse agonist in 2RH1,
analyzed in the GPCR binding workflow. Once we understand *why* carazolol acts as an
inverse agonist (binding closer to W6.48), we can use MoleculeProcessor to find
structurally similar compounds that might share this mechanism - or deliberately
search for compounds with different binding modes.

**What It Does:**
- Registers query ligand (carazolol from 2RH1 structure)
- Registers ~70 compounds from adrenergic drug database
- Calculates Tanimoto similarity using Morgan fingerprints
- Ranks and classifies hits by drug class

**Model System:** ~70 compounds across adrenergic drug classes
- Non-selective beta-blockers (same class as carazolol)
- Beta-1 selective blockers (cardioselective)
- Beta-2 agonists (opposing pharmacology)
- Natural catecholamines (endogenous ligands)
- Alpha-1/Alpha-2 adrenergic ligands
- Carbazole derivatives (similar scaffold)

**ProtOS Capabilities:**
- `MoleculeProcessor.save_entity()` - Register compounds with SMILES
- `MoleculeProcessor.calculate_properties()` - RDKit property calculation

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `molecule_similarity_distribution.png` | Histogram with similarity thresholds |
| `molecule_top_hits.png` | Top 30 hits bar chart, colored by class |
| `molecule_similarity_mw.png` | Scatter: similarity vs MW by drug class |
| `molecule_2d_structures.png` | RDKit 2D structures: query + top 3 hits |

**Key Results:**
- Top hits: carvedilol, propranolol (non-selective beta-blockers)
- Carbazole derivatives show high similarity (shared scaffold)
- Beta-blockers cluster together, agonists separate
- Clear pharmacological class separation by structure

**Insight:** Structure similarity predicts pharmacological class, but mechanism
understanding (from GPCR workflow) reveals *why* - the binding pose matters more
than the scaffold.

---

# PART II: DERIVED DATA PROCESSORS

These processors build on basic data with specialized annotations.

---

## 4. GRNProcessor
**File:** `processors/grn_processor_example.py`

**Question:** "What are the conserved positions in beta-adrenergic receptors?"

**Why GRN Matters:** Generic Residue Numbers enable comparison across GPCRs with
<30% sequence identity. Position "6.48" in any Class A GPCR is the same functional
position - the "toggle switch" tryptophan. This is essential for the GPCR workflow
where we compare 8 structures across 2 receptor subtypes.

**What It Does:**
- Downloads 6 beta-AR sequences from UniProt
- Annotates with GPCRdb generic residue numbers (GRN)
- Queries key conserved positions (x.50)

**Model System:** 6 beta-adrenergic receptors (ADRB1, ADRB2, ADRB3)

**ProtOS Capabilities:**
- `SequenceProcessor.annotate_with_grn()` - GPCRdb numbering
- `GRNProcessor.load_table()` - Query GRN annotations

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `grn_coverage.png` | Coverage bar chart + key positions heatmap |

**Key Results (x.50 positions, 100% conserved):**
- 1.50: N (Asn), 2.50: D (Asp), 3.50: R (Arg)
- 4.50: W (Trp), 5.50: P (Pro), 6.50: P (Pro), 7.50: P (Pro)

---

## 5. PropertyProcessor
**File:** `processors/property_processor_example.py`

**Question:** "What is the absorption maximum (lambda_max) for cone opsins in our dataset?"

**Design Principle:** Store ANY annotation/property for ANY registered entity.
The PropertyProcessor creates a unified interface for experimental measurements,
computed predictions, and custom metadata - all linked to entities by name.

**Connection to Workflows:** In the Prey Vision workflow, predicted lambda_max
values from the LAMBDA model are stored using PropertyProcessor. This means
predictions become *associated data*, not just outputs - enabling queries like
"which mutations shift lambda_max by more than 5 nm?"

**What It Does:**
- Loads annotations from the cone opsin dataset
- Matches known spectral data to sequences in the dataset
- Records properties linked to sequence entities
- Creates spectral sensitivity visualizations

**Model System:** Subset of 200 cone opsins with known lambda_max values
- Short wave (blue): ~360-440 nm
- Long wave (red): ~555-565 nm

**ProtOS Capabilities:**
- `PropertyProcessor.record_properties()` - Store entity-linked properties
- `PropertyProcessor.filter_by_property()` - Query by value
- Zero-configuration: no schema changes for new property types

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `property_sensitivity.png` | Gaussian curves on visible spectrum, colored by type |

**Key Results:**
- Properties successfully linked to registered sequence entities
- Spectral types clearly separated by lambda_max
- Demonstrates partial coverage (not all sequences have measurements)

**Data Flow:**
```
PropertyProcessor stores:
  - Experimental lambda_max (from literature)
  - Predicted lambda_max (from LAMBDA model in Workflow 3)
  - Any future annotations (binding affinity, stability, etc.)
```

---

## 6. EmbeddingProcessor
**File:** `processors/embedding_processor_example.py`

**Question:** "Do cone opsins cluster by spectral type in embedding space?"

**Connection to SequenceProcessor:** Loads the `cone_opsin_diversity` dataset
(200 sequences: 100 SW + 100 LW) created by the SequenceProcessor example.

**Connection to Workflow 2:** The pLM-enriched graph workflow uses per-residue
embeddings to create node features for structural graphs. This example shows
sequence-level (mean-pooled) embeddings; the workflow extends to residue-level.

**Insight:** Sequences with similar function cluster in embedding space,
demonstrating that pLMs capture functional relationships beyond sequence identity.

**What It Does:**
- Loads 200 cone opsin sequences from SequenceProcessor dataset
- Generates ESM2 embeddings (mean pooled per sequence)
- Creates t-SNE visualization colored by opsin type
- Creates similarity heatmap with 2x2 block structure showing SW vs LW clustering

**Model System:** 200 cone opsins (100 short_wave + 100 long_wave)

**ProtOS Capabilities:**
- `SequenceProcessor.load_dataset()` - Load dataset from another processor
- `EmbeddingProcessor.embed_sequences()` - ESM2 embedding generation
- `EmbeddingProcessor.load_embeddings()` - Retrieve from registry

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `embedding_tsne.png` | t-SNE colored by opsin type (blue/red) |
| `embedding_heatmap.png` | Pairwise similarity with 2x2 block structure, colored axis ticks |

**Key Results:**
- Clear 2x2 block structure in similarity heatmap (SW-SW, SW-LW, LW-SW, LW-LW)
- t-SNE shows clear clustering by spectral sensitivity
- ESM2 embeddings capture functional relationships

---

## 7. GraphProcessor
**File:** `processors/graph_processor_example.py`

**Question:** "How do residues around retinal interact in bacteriorhodopsin?"

**Connection to Workflow 2:** The pLM-enriched graph workflow combines this
structural contact information with per-residue embeddings, creating multi-modal
representations suitable for graph neural networks.

**What It Does:**
- Downloads bacteriorhodopsin structure (1C3W)
- Identifies retinal binding pocket residues
- Generates residue contact graph at multiple cutoffs
- Creates 3D and 2D network visualizations

**ProtOS Capabilities:**
- `StructureProcessor.list_ligands()` - Find ligands
- `GraphProcessor.generate_graph()` - Build contact graph
- `GraphProcessor.load_entity()` - Retrieve graph data

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `graph_binding_pocket_3d.html` | Interactive 3D with retinal bonds |
| `graph_contact_network_2d.html` | 2D spring layout, nodes by degree |

**Key Results:**
- 46 binding pocket residues within 7A of retinal
- Edge counts: 4A=12, 5A=72, 6A=162, 8A=324

---

# PART III: MODEL ORCHESTRATION

---

## 8. ModelManager
**File:** `processors/model_manager_example.py`

**Question:** "How do I prepare Boltz2 structure predictions for cluster submission?"

**Focus:** Wrapping processor data into compute-ready jobs for external resources.

**What It Does:**
- Loads sequences via SequenceProcessor
- Prepares Boltz2 jobs (wild-type + mutants)
- Generates job artifacts (config.yaml, sequences.fasta, metadata.json)
- Creates batch submission manifest

**Model System:** 3 rhodopsins + 2 mutant variants = 5 jobs

**ProtOS Capabilities:**
- `ModelManager.prepare_input()` - Generate job configuration
- `ModelManager.cards` - List available models (12 registered)
- Zero-configuration job packaging

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `model_manager_jobs.png` | Pie chart + Sankey workflow diagram |

**Available Models:**
- Structure: boltz2, boltzgen, ligand_mpnn
- Embeddings: esm2_t12_35m, esm2_t33_650m, ankh_large
- Docking: unidock, equibind
- Affinity: pocketdta, graphscoredta
- Generation: pocket2mol

---

# PART IV: ADVANCED WORKFLOWS

Each workflow combines multiple processors to answer a biological question.
The workflows are interconnected - understanding from one enables applications in another.

---

## Workflow 1: GPCR Binding Pocket Mechanism
**File:** `workflows/gpcr_binding_pocket_workflow.py`

### The Question
"What is the molecular basis for agonist vs inverse agonist vs antagonist action?"

The same receptor protein (ADRB2) can be activated (agonist), inhibited below basal
activity (inverse agonist), or blocked without changing basal activity (antagonist).
What structural features determine these different outcomes?

### Real-World Impact
Understanding GPCR mechanisms enables:
- **Drug design:** Design ligands with specific signaling profiles
- **Biased agonism:** Create drugs that activate only beneficial pathways
- **Side effect prediction:** Understand why similar compounds have different effects

### Connection to MoleculeProcessor
Once we understand the mechanism (inverse agonists bind closer to W6.48), we can
return to MoleculeProcessor to search for compounds that match or avoid this pattern.
The carazolol similarity search gains new meaning - compounds with similar scaffolds
but different binding poses may have opposite effects.

### Processors Used
```
StructureProcessor → Load 8 GPCR structures
        ↓
annotate_with_grn() → Map residues to GRN positions
        ↓
LigandInteractionAnalyzer → Find binding residues
        ↓
Hypothesis Testing → Compare distances by ligand type
        ↓
PropertyProcessor → Store results for future queries
```

### Model System: 8 adrenergic receptor structures
| PDB | Receptor | Ligand | Type | State |
|-----|----------|--------|------|-------|
| 3SN6 | ADRB2 | BI-167107 | full_agonist | active |
| 4LDO | ADRB2 | Adrenaline | full_agonist | active |
| 2Y02 | ADRB1 | Isoprenaline | full_agonist | active_like |
| 2Y04 | ADRB1 | Salbutamol | partial_agonist | intermediate |
| 2Y00 | ADRB1 | Dobutamine | partial_agonist | intermediate |
| 2RH1 | ADRB2 | Carazolol | inverse_agonist | inactive |
| 3NY9 | ADRB2 | ICI 118,551 | inverse_agonist | inactive |
| 2VT4 | ADRB1 | Cyanopindolol | antagonist | inactive |

### Hypotheses Tested
- **H1:** Agonists bind CLOSER to S5.43 than inverse agonists
- **H2:** Inverse agonists bind CLOSER to W6.48 (toggle switch) than agonists
- **H3:** Water at N6.55 is EXCLUSIVE to agonist-bound active structures

### Key Results
| Hypothesis | Result |
|------------|--------|
| H1 | SUPPORTED: Agonists (3.05A) closer to S5.43 than inverse agonists (3.45A) |
| H2 | SUPPORTED: Inverse agonists (3.43A) closer to W6.48 than agonists (4.34A) |
| H3 | SUPPORTED: Water at N6.55 only in active states (4LDO, 2Y02) |

### Cycling Back: From Mechanism to Molecules
```
GPCR Mechanism Understanding
        ↓
"Inverse agonists bind close to W6.48"
        ↓
MoleculeProcessor: Search for compounds that:
  - Similar scaffold but different W6.48 distance → may be agonists
  - Maximize W6.48 contact → potent inverse agonists
  - No W6.48 contact → neutral antagonists
```

### Key Figures
| Figure | Description |
|--------|-------------|
| `hypothesis_analysis.html` | Interactive bar charts for H1, H2, H3 |
| `hypothesis_analysis.png` | Static version for publication |
| `gpcr_mechanism_analysis.pml` | PyMOL script with GRN selections |

---

## Workflow 2: pLM-Enriched Structural Graphs
**File:** `workflows/plm_graph_workflow.py`

### The Question
"How do we create data structures suitable for AI models?"

Modern machine learning on proteins requires multi-modal representations: sequence
information (what the protein is), structural information (how residues are arranged),
and functional annotations (what positions matter).

### Real-World Impact
Enriched graph representations enable:
- **Property prediction:** Train GNNs to predict binding affinity, stability, activity
- **Transfer learning:** Pre-computed embeddings accelerate model training
- **Multi-task learning:** Same representation for multiple prediction tasks

### Connection to Other Workflows
This workflow creates the data structures used by:
- **Workflow 3 (Prey Vision):** LAMBDA model uses similar sequence+structure features
- **PropertyProcessor:** Any property can be linked to graph nodes for supervised learning

### The Key Insight
```
Traditional: sequence → model → prediction
Enriched:    sequence + structure + embeddings + properties → model → better prediction
```

### Processors Used
```
StructureProcessor → Load rhodopsin structure
        ↓
SequenceProcessor → Extract chain sequence
        ↓
GRNProcessor → Annotate functional positions
        ↓
EmbeddingProcessor → Generate per-residue pLM embeddings
        ↓
GraphProcessor → Build residue contact graph
        ↓
Combine → Enriched graph with node features
```

### What Gets Created
Each node (residue) in the graph has:
- **Position features:** 3D coordinates, secondary structure
- **Sequence features:** Amino acid identity, conservation
- **Embedding features:** 1536-dim pLM vector (ANKH-large)
- **Functional annotations:** GRN position, binding pocket membership
- **Properties:** Lambda_max (if known), any other annotations

Each edge (contact) has:
- **Distance:** CA-CA distance
- **Embedding similarity:** Cosine similarity of node embeddings

### Extending to New Properties
```python
# Any property can become a training target:
PropertyProcessor.record_properties("binding_affinity", [
    {"entity_name": "1U19_A", "Kd": 1.2, "source": "experiment"},
    {"entity_name": "mutant_K296A", "Kd": 15.3, "source": "predicted"},
])

# Graph nodes automatically gain this annotation
# Train a GNN to predict Kd from enriched graph features
```

### Key Figures
| Figure | Description |
|--------|-------------|
| `plm_enriched_graph.html` | 3D graph with embedding-colored nodes |

---

## Workflow 3: Prey Vision Enhancement (Red-Shift Mutation Screen)
**File:** `workflows/redshift_mutation_workflow.py`

### The Biological Story

**Why are tigers orange?**

Tigers are orange because their main prey (deer, wild boar) are dichromats - they have
only two types of cone opsins (SW and LW), not three like humans. To a dichromatic
animal, orange and green are indistinguishable. The tiger's orange coat is
*camouflaged* against green foliage from the prey's perspective.

### The Visualization (3-panel figure)

| Panel | Description |
|-------|-------------|
| **Human view** | Tiger in grass - orange clearly visible against green |
| **Prey view (current)** | Same image through dichromatic filter - tiger blends with grass |
| **Prey view (enhanced)** | With red-shifted opsin - tiger becomes visible |

### The Question
"Can we engineer a mutation that would give prey animals the ability to see tigers?"

### The Approach
1. Start with prey's LW opsin (absorbs ~560 nm)
2. Screen mutations that red-shift absorption toward green sensitivity
3. A duplicated, shifted opsin would create MW-like sensitivity
4. Dichromat → Trichromat-like = ability to distinguish orange from green

### Processors Used
```
StructureProcessor → Load rhodopsin structure
        ↓
GRNProcessor → Identify binding pocket positions
        ↓
Generate mutations → In silico mutagenesis
        ↓
EmbeddingProcessor → Per-residue features
        ↓
LAMBDA model → Predict lambda_max shift
        ↓
PropertyProcessor → Store predictions
        ↓
Rank candidates → Best mutations for enhanced vision
```

### Key Results
| Mutation | λmax (nm) | Shift | Progress to Green |
|----------|-----------|-------|-------------------|
| L114S | 451.1 | +8.7 nm | 11% |
| L114G | 447.6 | +5.3 nm | 7% |
| F93I | 445.4 | +3.1 nm | 4% |
| Wild-type | 442.3 | - | 0% |

**Biological Interpretation:** Single mutations achieve modest shifts. Evolution
of trichromacy likely required multiple mutations accumulated over time.

### Real-World Connection
This workflow demonstrates:
- **Gene duplication + divergence:** How new sensory capabilities evolve
- **Structure-function relationships:** Why specific positions affect absorption
- **Computational screening:** Test hypotheses before expensive experiments

### Key Figures
| Figure | Description |
|--------|-------------|
| `prey_vision_enhancement.html` | Mutation screen results |
| `tiger_trichromat.png` | Human view (to be added) |
| `tiger_dichromat.png` | Prey view - camouflaged (to be added) |
| `tiger_enhanced.png` | Enhanced prey view (to be added) |

---

## Workflow 4: Rhodozyme - Light-Activated Enzyme Design
**File:** `workflows/rhodozyme_design_workflow.py`
**Visualization:** `workflows/pymol_rhodozyme_visualization.py`

### The Vision: Light-Controlled Chemistry

Imagine a single reaction vessel where you control multiple chemical reactions
with light. Different wavelengths activate different enzymes. Reactions happen
only where and when you shine the light.

**Spatiotemporal control:**
- Step 1: Blue light → Protease cleaves substrate
- Step 2: Green light → Ligase joins fragments
- Step 3: Red light → Kinase adds phosphate

All in one container, controlled by light patterns.

### The Question
"Can we design an enzyme that is active only when illuminated?"

### The Mechanism
GPCRs undergo large conformational changes upon activation:
- TM5 moves ~3Å
- TM6 moves ~14Å outward
- These movements are triggered by light (in rhodopsin) or ligand binding

**The Rhodozyme Concept:**
Place a catalytic triad on GPCR helices such that:
- In dark state: Triad geometry is disrupted → inactive
- In light state: Helices move, triad aligns → active

### Processors Used
```
StructureProcessor → Load active rhodopsin (3PQR)
        ↓
GRNProcessor → Map helix positions
        ↓
For each enzyme (trypsin, chymotrypsin, papain, subtilisin):
        ↓
  Extract catalytic triad geometry
        ↓
  Screen GPCR surface for matching geometry
        ↓
  Require ≥1 residue on TM5 or TM6 (dynamic helices)
        ↓
  Rank by geometric fit + placement
        ↓
ModelManager → Prepare Boltz2 jobs for structure prediction
```

### Design Constraint
At least one catalytic residue must be on TM5 or TM6 (the helices that move
during activation). This ensures the catalytic geometry changes with light.

### Key Results
| Enzyme | Best Design | Score | Dynamic Helix Residues |
|--------|-------------|-------|------------------------|
| Trypsin | Rhodozyme-TRP_01 | 4.61 | F230(TM5), Q246(TM6) |
| Chymotrypsin | Rhodozyme-CHY_01 | 4.69 | F230(TM5), Q246(TM6) |
| Papain | Rhodozyme-PAP_01 | 4.37 | V139(TM3), F230(TM5) |
| Subtilisin | Rhodozyme-SUB_01 | 4.84 | V139(TM3), F230(TM5) |

**31 total designs** across 4 enzyme types, all with residues on dynamic helices.

### Real-World Applications
- **Prodrug activation:** Light-triggered drug release in tumors
- **Biosensors:** Enzyme activity reports on illumination
- **Synthetic biology:** Light-controlled metabolic pathways
- **Biomanufacturing:** Precise reaction timing without chemical additives

### Next Steps
1. Boltz2 structure prediction of chimeras
2. Verify catalytic geometry in both states
3. Molecular dynamics of light-induced conformational change
4. Experimental validation

### Key Figures
| Figure | Description |
|--------|-------------|
| `rhodozyme_designs_combined.html` | All designs with scores |
| PyMOL script | Active vs inactive triad geometry |

---

## Running All Examples

To run all processor examples:
```bash
cd thesis/processors
for f in *_example.py; do python "$f"; done
```

To run all workflows:
```bash
cd thesis/workflows
python gpcr_binding_pocket_workflow.py
python plm_graph_workflow.py
python redshift_mutation_workflow.py
python rhodozyme_design_workflow.py
```

---

## Information Flow Summary

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                              THE CYCLE OF UNDERSTANDING                          │
└─────────────────────────────────────────────────────────────────────────────────┘

     STRUCTURE                    MECHANISM                     APPLICATION
         │                            │                              │
         ▼                            ▼                              ▼
  Load structures ──────────► Analyze binding ──────────► Design new molecules
  (StructureProc)              (GPCR workflow)            (MoleculeProc)
         │                            │                              │
         │                            ▼                              │
         │                    Test hypotheses                        │
         │                    (H1, H2, H3)                          │
         │                            │                              │
         └────────────────────────────┼──────────────────────────────┘
                                      │
                                      ▼
                          PREDICT NEW PROPERTIES
                             (PropertyProc)
                                      │
                    ┌─────────────────┼─────────────────┐
                    │                 │                 │
                    ▼                 ▼                 ▼
             Spectral shift    Binding affinity    Enzyme activity
             (Workflow 3)      (future work)       (Workflow 4)
```

### Key Connections

1. **SequenceProcessor → EmbeddingProcessor → PropertyProcessor**
   - Sequences get embeddings, embeddings predict properties

2. **StructureProcessor → GRNProcessor → Cross-family comparison**
   - GRN enables comparison of structures with <30% sequence identity

3. **GPCR Mechanism → MoleculeProcessor → Drug Design**
   - Understanding binding modes guides molecular similarity searches

4. **GraphProcessor + EmbeddingProcessor → AI Training Data**
   - Enriched representations enable machine learning

5. **All Workflows → PropertyProcessor → Queryable Results**
   - Every prediction becomes associated data for future queries

---

## Design Principles

1. **Zero-Configuration Data Management**
   - All processors auto-register entities
   - Datasets managed through unified API
   - Cross-processor references via entity names

2. **Composable Workflows**
   - Each processor does one thing well
   - Workflows chain processors naturally
   - ModelManager bridges to external compute

3. **Biological Relevance**
   - Examples answer real questions
   - Model systems chosen for insight
   - Results connect to experimental data

4. **ProtOS Utilities Throughout**
   - Use ProtOS functions for each analytical step
   - Avoid custom implementations when ProtOS provides the functionality
   - GRN annotation, ligand analysis, etc. all via ProtOS

5. **Information Flows, Understanding Accumulates**
   - Each workflow builds on previous understanding
   - Results become inputs to future analyses
   - The whole is greater than the sum of parts
