# ProtOS Thesis Examples - Comprehensive Guide

**Last Updated:** 2026-01-11
**Status:** 8 Processor Examples + 4 Advanced Workflows = COMPLETE

This document serves as the complete reference for all thesis examples, their outputs,
figures, and the story they tell about ProtOS capabilities.

---

## Story Thread

The examples follow a logical progression:

1. **Basic Data Processors** - Load, store, and manage fundamental biological entities
   - Sequence, Structure, Molecule

2. **Derived Data Processors** - Build on basic data with specialized annotations
   - GRN (requires Sequence), Property (links to entities), Embedding (requires Sequence), Graph (requires Structure)

3. **ModelManager** - Orchestrate computational jobs for external compute

4. **Advanced Workflows** - Combine processors to answer real biological questions:
   - GPCR binding mechanisms
   - pLM-enriched structural graphs
   - Spectral tuning mutations
   - Light-activated enzyme design

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
│   └── rhodozyme_design_workflow.py
│
└── outputs/                     # Data outputs by processor/workflow
```

---

# PART I: BASIC DATA PROCESSORS

These processors handle the fundamental biological data types.

---

## 1. SequenceProcessor
**File:** `processors/sequence_processor_example.py`

**Question:** "What is the sequence diversity across cone opsin types?"

**Connection to Other Processors:** This example creates the `cone_opsin_diversity` dataset
that feeds directly into the EmbeddingProcessor example, demonstrating cross-processor data flow.

**What It Does:**
- Downloads human cone opsin query sequences (SW/MW/LW)
- Runs separate NCBI BLAST searches for each spectral type
- Fetches ~150 homolog sequences and annotates by opsin type
- Creates dataset with spectral type annotations for downstream analysis

**Model System:** 3 query opsins → ~150 total sequences
| Gene | Type | λmax | Description |
|------|------|------|-------------|
| OPN1SW | short_wave | 420 nm | Blue cone opsin |
| OPN1MW | medium_wave | 530 nm | Green cone opsin |
| OPN1LW | long_wave | 560 nm | Red cone opsin |

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
| `sequence_phylogenetic_tree.png` | Dendrogram colored by opsin type (no leaf labels) |

**Key Results:**
- ~50 sequences per opsin type from SwissProt
- MMseqs2 all-vs-all similarity → distance matrix → phylogenetic tree
- Tree shows clear separation of SW/MW/LW clades
- Dataset annotated with opsin_type for EmbeddingProcessor

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

**Connection to Workflows:** Carazolol is the crystallographic inverse agonist in 2RH1,
the key structure from the GPCR binding workflow. This search finds structurally related
beta-blockers and adrenergic ligands.

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

---

# PART II: DERIVED DATA PROCESSORS

These processors build on basic data with specialized annotations.

---

## 4. GRNProcessor
**File:** `processors/grn_processor_example.py`

**Question:** "What are the conserved positions in beta-adrenergic receptors?"

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

**Question:** "What is the absorption maximum (lambda_max) for these opsins?"

**Design Principle:** Store ANY annotation/property for ANY registered entity.
The PropertyProcessor creates a unified interface for experimental measurements,
computed predictions, and custom metadata - all linked to entities by name.

**Connection to Workflows:** In the Red-Shift Mutation workflow, predicted lambda_max
values are stored using PropertyProcessor - demonstrating that predictions become
*associated data*, not just outputs. This enables querying predicted vs experimental values.

**What It Does:**
- Records experimental spectral properties for 12 opsins
- Links properties to sequence entities by name
- Filters by property values (UV-sensitive, long-wavelength)
- Creates spectral sensitivity visualizations and property table

**Model System:** 12 opsins with lambda_max values (360-568 nm)

**ProtOS Capabilities:**
- `PropertyProcessor.record_properties()` - Store entity-linked properties
- `PropertyProcessor.filter_by_property()` - Query by value
- Zero-configuration: no schema changes for new property types

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `property_spectral.png` | Lambda_max bar + box plot by type |
| `property_sensitivity.png` | Gaussian curves on visible spectrum |
| `property_table.png` | Property table with wavelength-colored cells |

**Key Results:**
- UV-sensitive: OPSB_MOUSE (360 nm), OPN5_HUMAN (380 nm)
- Long-wavelength: OPSR_HUMAN (560 nm), BRHOD_HALSA (568 nm)

---

## 6. EmbeddingProcessor
**File:** `processors/embedding_processor_example.py`

**Question:** "Do cone opsins cluster by spectral type in embedding space?"

**Connection to SequenceProcessor:** This example loads the `cone_opsin_diversity` dataset
created by the SequenceProcessor example, demonstrating cross-processor data flow.

**Insight:** Sequences with similar function may cluster in embedding space even when
sequence identity is not the highest predictor - pLMs capture functional relationships.

**What It Does:**
- Loads ~150 cone opsin sequences from SequenceProcessor dataset
- Generates ESM2 embeddings for all sequences
- Compares sequence identity vs embedding similarity by opsin type
- Creates t-SNE visualization showing spectral type clustering

**Model System:** ~150 cone opsins from SequenceProcessor (short/medium/long wave)

**ProtOS Capabilities:**
- `SequenceProcessor.load_dataset()` - Load dataset from another processor
- `EmbeddingProcessor.embed_sequences()` - ESM2 embedding generation
- `EmbeddingProcessor.load_embeddings()` - Retrieve from registry

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `embedding_tsne.png` | t-SNE colored by opsin type (blue/green/red) |
| `embedding_similarity_analysis.png` | Same-type vs different-type similarity |
| `embedding_heatmap.png` | Pairwise similarity heatmap |

**Key Results:**
- Same-type pairs have higher embedding similarity than different-type pairs
- t-SNE shows clear clustering by spectral sensitivity (SW/MW/LW)
- ESM2 embeddings capture functional relationships beyond sequence identity

---

## 7. GraphProcessor
**File:** `processors/graph_processor_example.py`

**Question:** "How do residues around retinal interact in bacteriorhodopsin?"

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
Overview figures (showing processor connections) will be provided separately.

---

## Workflow 1: GPCR Binding Pocket Mechanism
**File:** `workflows/gpcr_binding_pocket_workflow.py`

**Question:** "What differentiates agonist vs inverse agonist binding in histamine receptors?"

**Processors:** Structure -> GRN -> Analysis

**Story:** Uses GRN numbering to compare binding pockets across receptor states,
revealing how ligand type correlates with specific residue movements.

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `binding_pocket_comparison.html` | 3D pocket overlay |
| `grn_heatmap.png` | Distance matrix at key positions |
| PyMOL sessions for publication figures |

---

## Workflow 2: pLM-Enriched Structural Graphs
**File:** `workflows/plm_graph_workflow.py`

**Question:** "How do sequence embeddings correlate with structural contacts?"

**Processors:** Structure -> Sequence -> GRN -> Embedding -> Graph

**Story:** Demonstrates multi-modal data integration - combining structural
contact graphs with pLM embeddings to create enriched representations.

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `plm_graph_3d.html` | 3D graph with embedding-colored nodes |
| `embedding_contact_correlation.png` | Similarity vs structural distance |

---

## Workflow 3: Red-Shift Mutation Screen
**File:** `workflows/redshift_mutation_workflow.py`

**Question:** "Which single mutations could red-shift rhodopsin absorption?"

**Processors:** Structure -> GRN -> Sequence -> Property -> (Boltz2)

**Story:** Systematic mutation screen in the retinal binding pocket,
scoring candidates by structural and evolutionary constraints.

**Key Figures:**
| Figure | Description |
|--------|-------------|
| `redshift_candidates.png` | Ranked mutation scores |
| `binding_pocket_mutations.html` | 3D view of mutation sites |

---

## Workflow 4: Rhodozyme Light-Activated Enzyme Design
**File:** `workflows/rhodozyme_design_workflow.py`

**Question:** "Can we design a light-activated enzyme using rhodopsin's conformational change?"

**Processors:** Structure -> GRN -> Geometry matching -> Design ranking

**Story:** Places catalytic triads (trypsin, chymotrypsin, papain, subtilisin)
on GPCR helices that move during activation, creating "rhodozymes" that
should be active only in one conformational state.

**Key Innovation:** Requires at least one catalytic residue on TM5 or TM6
(the helices that move 3-14A during GPCR activation).

**Key Figures:**
| Figure | Description |
|--------|-------------|
| PyMOL visualization script | Active vs inactive triad geometry |
| Design summary tables | Per-enzyme ranked designs |

**Key Results:**
- 31 designs across 4 enzyme types
- Best: Rhodozyme-SUB_01 (subtilisin, score=4.84)
- All designs have residues on dynamic helices (TM5/TM6)

---

## Figure Quality Notes

Current figures are functional drafts. For publication:
- Increase font sizes
- Standardize color schemes
- Add proper axis labels and legends
- Consider PyMOL for structural figures

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
