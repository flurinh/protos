# Refactoring Processor Examples

This document specifies exactly how to update each processor example to align with the thesis narrative. All examples should use consistent structures throughout.

---

## Consistent Structure Set

Use these structures across ALL processor examples and workflows:

### Type I (Microbial Rhodopsins)
| PDB | Name | Function | Notes |
|-----|------|----------|-------|
| **1C3W** | Bacteriorhodopsin | Proton pump | Primary Type I reference |
| 3UG9 | Channelrhodopsin-2 | Cation channel | Secondary Type I (optogenetics) |

### Type II (Animal Opsins / GPCRs)
| PDB | Name | Function | Notes |
|-----|------|----------|-------|
| **1U19** | Bovine rhodopsin | Dim light vision | Primary Type II reference |
| 2Z73 | Squid rhodopsin | Vision | Secondary Type II (invertebrate) |

### GPCR (Non-opsin)
| PDB | Name | Function | Notes |
|-----|------|----------|-------|
| 2RH1 | β2-adrenergic receptor | Adrenergic signaling | For Molecule/GPCR workflow only |

---

## Example Refactoring Details

### 1. SequenceProcessor (`sequence_processor_example.py`)

**Status:** ✅ No changes needed

**Current:** Builds cone opsin diversity dataset (200 sequences, SW vs LW)
- Query: Human SW opsin (P03999), Human LW opsin (P04000)
- BLAST search against SwissProt
- Creates `cone_opsin_diversity` dataset

**Figures:**
- `sequence_opsin_diversity.png` ✅ Exists
- `sequence_phylogenetic_tree.png` ✅ Exists

---

### 2. StructureProcessor (`structure_processor_example.py`)

**Status:** ✅ Updated - Type I vs Type II emphasis with colorscales

**TODO:** Aligns 4 structures (1C3W, 3UG9, 1U19, 2Z73). First the MO and AOs on each other to show the conserved 'fold'.
Here retinal is in the same binding pocket.
Then on retinal between 1c3w and 1u19 to show that their 'binding modes' are completely different.

**Required changes:**
1. Update visualization title to emphasize "Type I vs Type II - Different Folds"
2. Add text output comparing fold topology
3. Emphasize that binding is NOT homologous despite same ligand

**Structures to use:**
- Reference: **1C3W** (bacteriorhodopsin, Type I)
- Mobile: **1U19** (bovine rhodopsin, Type II)
- Optional additional: 3UG9, 2Z73

**Key insight to demonstrate:**
- Same ligand (retinal)
- Different 7TM topology
- Binding poses are NOT homologous
- Yet both tune λmax through binding pocket

**Figures:**
- `opsin_alignment.html` ✅ Exists - update title/labels

---

### 3. GraphProcessor (`graph_processor_example.py`)

**Status:** ✅ Completed - Both Type I and Type II with comparison figure

**Current:** Compares binding pocket graphs for 1C3W (Type I) and 1U19 (Type II)

**Required changes:**
1. Generate binding pocket graph for BOTH structures:
   - **1C3W** (bacteriorhodopsin, Type I)
   - **1U19** (bovine rhodopsin, Type II)
2. Create side-by-side comparison visualization
3. Show graph statistics for both (nodes, edges)
4. Demonstrate that topology is comparable despite fold difference

**Structures to use:**
- **1C3W** (Type I) - extract 7Å pocket around RET
- **1U19** (Type II) - extract 7Å pocket around RET

**New figures needed:**
- `graph_binding_pocket_3d_1C3W.html` - Type I pocket
- `graph_binding_pocket_3d_1U19.html` - Type II pocket
- `graph_comparison.png` - Side-by-side statistics

**Key insight to demonstrate:**
- Both pockets have similar node count (~40-50 residues)
- Both have similar edge density
- Graphs ARE comparable despite fold difference
- This is LAMBDA's input representation

---

### 4. GRNProcessor (`grn_processor_example.py`)

**Status:** ✅ Completed - Compares bovine vs squid rhodopsin with offset correction

**Current:** Compares bovine rhodopsin (1U19) vs squid rhodopsin (2Z73) binding pocket positions

**Implemented changes:**
1. Changed comparison to TWO ANIMAL OPSINS:
   - **1U19** (bovine rhodopsin)
   - **2Z73** (squid rhodopsin)
2. Focus on BINDING POCKET positions with x.50 alignment:
   - Counterion position (3.28)
   - Spectral tuning sites
   - Schiff base lysine (7.42, NOT 7.43)
3. Includes MOGRN gap discussion in output
4. Fixed residue numbering offset for squid (PDB starts at residue 9)

**Structures used:**
- **1U19** (bovine rhodopsin, Type II) - reference
- **2Z73** (squid rhodopsin, Type II) - mobile

**GRN positions (verified from GRN annotation):**
| GRN | Function | Bovine (1U19) | Squid (2Z73) |
|-----|----------|---------------|--------------|
| 3.28 | Counterion | E113 | Y111 |
| 3.32 | Tuning site | A117 | G115 |
| 6.48 | Rotamer switch | W265 | W274 |
| 7.42 | Schiff base K | K296 | K305 |

**Note:** Squid has Y (Tyr) at 3.28 counterion position instead of E (Glu) - interesting evolutionary difference!

**Figures:**
- `grn_spectral_tuning.pml` - Update to show two opsins
- `grn_structural_alignment.pml` - Update to bovine vs squid

**Key insight to demonstrate:**
- GRN enables systematic comparison of binding pocket positions
- Same GRN position = same functional role
- Works for animal opsins (GPCRs with Ballesteros-Weinstein)
- BUT microbial rhodopsins don't have this → MOGRN gap

---

### 5. EmbeddingProcessor (`embedding_processor_example.py`)

**Status:** ✅ No changes needed

**Current:**
- Uses `cone_opsin_diversity` dataset from SequenceProcessor
- Generates ESM2 embeddings
- Shows SW vs LW clustering

**Figures:**
- `embedding_tsne.png` ✅ Exists
- `embedding_heatmap.png` ✅ Exists
- `embedding_similarity_analysis.png` ✅ Exists

---

### 6. PropertyProcessor (`property_processor_example.py`)

**Status:** ✅ Minor update only

**Current:** Associates λmax values with cone opsin sequences

**Optional improvements:**
1. Add clearer explanation of property association concept
2. Show both literature values AND predicted values (if available)

**Figures:**
- `property_sensitivity.png` ✅ Exists

---

### 7. MoleculeProcessor (`molecule_processor_example.py`)

**Status:** ✅ No changes needed

**Current:** Carazolol similarity search with beta-blockers/agonists

**Note:** This example is intentionally separate from opsin narrative. It connects to GPCR Binding Pocket workflow.

**Figures:**
- `molecule_similarity_distribution.png` ✅ Exists
- `molecule_top_hits.png` ✅ Exists
- `molecule_2d_structures.png` ✅ Exists

---

## Workflow Consistency Check

### Binding Pocket Graphs Workflow
**Must use:**
- 1C3W (Type I) + 1U19 (Type II) aligned on retinal
- Extract binding pocket graphs from both
- Should reuse outputs from StructureProcessor and GraphProcessor examples

### Prey Vision Enhancement Workflow
**Must use:**
- Cone opsin sequences from SequenceProcessor dataset
- GRN annotation via GRNProcessor
- Embeddings from EmbeddingProcessor

### GPCR Binding Pocket Workflow
**Must use:**
- 2RH1 (β2AR) from MoleculeProcessor example
- Carazolol and other ligands
- GRN for position annotation

### Rhodozyme Design Workflow
**Uses different structures** (enzyme-specific, not part of opsin set)

---

## Priority Order for Updates

1. **GraphProcessor** - ✅ Completed (both Type I and Type II structures)
2. **GRNProcessor** - ✅ Completed (two opsins, binding pocket focus, offset fix)
3. **StructureProcessor** - ✅ Completed (Type I vs Type II emphasis)
4. Others are already aligned

---

## Structure Consistency Verification

After refactoring, verify that:

| Example | 1C3W | 1U19 | 2Z73 | 3UG9 | 2RH1 |
|---------|------|------|------|------|------|
| SequenceProcessor | - | - | - | - | - |
| StructureProcessor | ✓ | ✓ | ✓ | ✓ | - |
| GraphProcessor | ✓ | ✓ | - | - | - |
| GRNProcessor | - | ✓ | ✓ | - | - |
| EmbeddingProcessor | - | - | - | - | - |
| PropertyProcessor | - | - | - | - | - |
| MoleculeProcessor | - | - | - | - | ✓ |
| Binding Pocket WF | ✓ | ✓ | - | - | - |
| GPCR Binding WF | - | - | - | - | ✓ |

---

## Figure Numbering (Final)

| Section | Figure | File | Status |
|---------|--------|------|--------|
| 2.2a | SequenceProcessor - diversity | `sequence_opsin_diversity.png` | ✅ |
| 2.2b | SequenceProcessor - tree | `sequence_phylogenetic_tree.png` | ✅ |
| 2.3 | StructureProcessor - alignment | `opsin_alignment.html` | ✅ Updated |
| 2.4a | GraphProcessor - Type I pocket | `graph_binding_pocket_3d_1C3W.html` | ✅ Created |
| 2.4b | GraphProcessor - Type II pocket | `graph_binding_pocket_3d_1U19.html` | ✅ Created |
| 2.4c | GraphProcessor - comparison | `graph_comparison.png` | ✅ Created |
| 2.5a | GRNProcessor - tuning sites | `grn_spectral_tuning.pml` | ✅ Updated |
| 2.5b | GRNProcessor - alignment | `grn_structural_alignment.pml` | ✅ Updated |
| 2.6a | EmbeddingProcessor - t-SNE | `embedding_tsne.png` | ✅ |
| 2.6b | EmbeddingProcessor - heatmap | `embedding_heatmap.png` | ✅ |
| 2.6c | EmbeddingProcessor - analysis | `embedding_similarity_analysis.png` | ✅ |
| 2.7 | PropertyProcessor - sensitivity | `property_sensitivity.png` | ✅ |
| 2.8a | MoleculeProcessor - distribution | `molecule_similarity_distribution.png` | ✅ |
| 2.8b | MoleculeProcessor - top hits | `molecule_top_hits.png` | ✅ |
| 2.8c | MoleculeProcessor - structures | `molecule_2d_structures.png` | ✅ |


  UPDATE: Complete Figure Inventory
  ┌──────────────┬──────────────────────────────────────┬──────┬────────────────────────┐
  │  Processor   │             Figure File              │ Type │      Description       │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Sequence     │ sequence_opsin_diversity.png         │ PNG  │ Taxonomy diversity     │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Sequence     │ sequence_phylogenetic_tree.png       │ PNG  │ Phylogenetic tree      │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Sequence     │ sequence_blast_scatter.png           │ PNG  │ BLAST hit scatter      │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Structure    │ opsin_alignment.html                 │ HTML │ 3D alignment viewer    │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Structure    │ structure_type1_vs_type2.pml         │ PML  │ PyMOL script           │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Graph        │ graph_binding_pocket_3d_1C3W.html    │ HTML │ Type I pocket 3D       │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Graph        │ graph_binding_pocket_3d_1U19.html    │ HTML │ Type II pocket 3D      │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Graph        │ graph_comparison.png                 │ PNG  │ Statistics comparison  │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ GRN          │ grn_structural_alignment.pml         │ PML  │ Alignment + GRN labels │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ GRN          │ grn_spectral_tuning.pml              │ PML  │ Tuning sites view      │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Embedding    │ embedding_tsne.png                   │ PNG  │ t-SNE clustering       │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Embedding    │ embedding_heatmap.png                │ PNG  │ Similarity heatmap     │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Embedding    │ embedding_similarity_analysis.png    │ PNG  │ Class similarity       │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Property     │ property_sensitivity.png             │ PNG  │ λmax distribution      │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Molecule     │ molecule_similarity_distribution.png │ PNG  │ Tanimoto histogram     │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Molecule     │ molecule_top_hits.png                │ PNG  │ Top hits bar chart     │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Molecule     │ molecule_similarity_mw.png           │ PNG  │ Similarity vs MW       │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ Molecule     │ molecule_2d_structures.png           │ PNG  │ 2D structures grid     │
  ├──────────────┼──────────────────────────────────────┼──────┼────────────────────────┤
  │ ModelManager │ model_manager_jobs.png               │ PNG  │ Job status chart       │