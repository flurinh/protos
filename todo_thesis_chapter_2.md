# TODO: Thesis Chapter 2 - ProtOS

## Chapter Location
**Thesis:** `/data/fast/ghostwriter/PKB7mgfarp/chapters/02-protos/`
**Project:** `/data/fast/projects/protos/`

---

## Current Chapter Draft (index.md)

```markdown
# ProtOS

The previous chapter identified three gaps in opsin research: microbial rhodopsins lack the standardized positional annotation that enabled systematic GPCR comparison, no method predicts spectral properties across the Type I / Type II divide, and existing tools do not compose into workflows that link sequences, structures, and functional properties while maintaining consistent identity. This chapter addresses the third gap through ProtOS, a framework for protein data access and management that also provides the infrastructure underlying the research contributions that follow.

Protein research requires data from multiple sources, and each source uses its own conventions. Bacteriorhodopsin appears as P02945 in UniProt, 1C3W in the Protein Data Bank (one of over fifty experimental structures), and AF-P02945-F1 in AlphaFold DB. Assembling data for even a single protein means navigating separate databases, downloading files in different formats, and tracking which data belongs together. The burden grows with scope. Comparing binding pockets across 129 microbial rhodopsin structures, the validation set for MOGRN, requires hundreds of downloads. Collecting the 1,799 sequences with measured absorption maxima for LAMBDA training means reconciling identifiers across sources and linking measurements scattered through decades of literature. Each step is simple; collectively, they dominate the work.

ProtOS addresses this through two connected capabilities. Unified data access provides a single interface to major protein databases, routing requests by identifier format and caching results locally. The entity system maintains identity across data types, so that a protein loaded from the PDB, annotated with standardized positions, enriched with embeddings, and linked to measured properties remains a single coherent object rather than a collection of disconnected files.

Seven processors handle the data types protein research requires: sequences, structures, residue contact graphs, standardized position annotations, embeddings, properties, and molecules. Each processor manages one data modality. Outputs from one become inputs to another through persistent identity, enabling workflows that span modalities without the format conversion and identifier reconciliation that currently fragment computational biology.

The following sections present the framework through examples drawn from opsin research. These examples are not demonstrations of ProtOS capabilities in the abstract; they reveal pieces of the spectral tuning puzzle, building toward the research contributions in subsequent chapters.

![Figure 2.1: ProtOS architecture overview...]

The framework has boundaries. ProtOS requires Python for programmatic use, though ProtOS-MCP provides natural language access for researchers without programming experience. Database coverage includes UniProt, PDB, AlphaFold DB, and NCBI, but not every biological database; extending to additional sources requires implementing new loaders. GRN annotation depends on having a numbering system defined for the protein family of interest—GPCRs have Ballesteros-Weinstein, microbial rhodopsins now have MOGRN, but families without standardized numbering cannot use the GRN Processor until such systems are developed. These limitations define where ProtOS applies and where future development is needed.
```

## Subchapters

| # | Title | File |
|---|-------|------|
| 2.1 | Unified Data Access | `01-unified-data-access/index.md` |
| 2.2 | Entity Registry and Datasets | `02-entity-registry-datasets/index.md` |
| 2.3 | Zero Configuration | `03-zero-configuration/index.md` |
| 2.4 | Processors | `04-processors/index.md` |
| 2.5 | Model Manager | `05-model-manager/index.md` |

### Processor Subchapters (2.4.x)
- 2.4.1 Sequence Processor
- 2.4.2 Structure Processor
- 2.4.3 Graph Processor
- 2.4.4 GRN Processor
- 2.4.5 Embedding Processor
- 2.4.6 Property Processor
- 2.4.7 Molecule Processor

---

## Project Context (from /data/fast/projects/protos/)

### Available Resources
- **README.md** (13.5 KB) - Complete API reference
- **docs/** - Comprehensive processor documentation
- **thesis/processors/** - 8 processor example scripts with outputs
- **thesis/workflows/** - 4 advanced workflow demonstrations
- **tests/** - 92+ test files validating features
- **workflow_report.md** - 20 workflow tests, all passed

### Key Demonstrations Available
1. **Cone opsin diversity** (sequence_processor_example.py) - 200 sequences, BLAST, clustering
2. **Type I vs II comparison** (structure_processor_example.py) - Retinal binding alignment
3. **Ligand similarity** (molecule_processor_example.py) - Carazolol analogs, drug classes
4. **GRN positions** (grn_processor_example.py) - β-adrenergic conserved sites
5. **Spectral properties** (property_processor_example.py) - Experimental data linking
6. **Embedding clustering** (embedding_processor_example.py) - ESM2 functional separation
7. **Binding pocket graphs** (graph_processor_example.py) - 1C3W vs 1U19 topology
8. **Structure prediction** (model_manager_example.py) - Boltz2 integration

### Rhodozyme Application
Located in `thesis/workflows/rhodozyme_design_workflow.py`:
- Light-activated enzyme engineering
- Catalytic triad grafting onto rhodopsin
- 31 candidate designs identified
- PyMOL visualization session (1.9 MB)

---

## Issues to Address

### HIGH PRIORITY

#### H1. Chapter is too long / diluted focus
- **Problem**: 7 processors + Model Manager creates very long chapter
- **Impact**: Reader fatigue, molecules processor admits opsins don't use it
- **Action**: Consider:
  - [ ] Move Molecule Processor to supplementary (opsins don't need it)
  - [ ] Condense processor sections to essential capabilities
  - [ ] Focus examples strictly on opsin-relevant use cases

#### H2. Missing explicit connection to MOGRN/LAMBDA
- **Problem**: Chapter ends without previewing how ProtOS enables subsequent chapters
- **Impact**: Reader doesn't understand ProtOS as foundation for MOGRN/LAMBDA
- **Action**:
  - [ ] Add transition paragraph connecting to Chapter 3
  - [ ] Show how GRNProcessor specifically enables MOGRN integration
  - [ ] Preview LAMBDA's use of processors (embedding + graph + property)

#### H3. Weak transitions between processors
- **Problem**: Each processor section feels standalone
- **Impact**: Reader loses sense of how they compose
- **Action**:
  - [ ] Add bridging sentences between processors
  - [ ] Show a mini-workflow threading multiple processors
  - [ ] Graph Processor should preview its role in LAMBDA

### MEDIUM PRIORITY

#### M1. No error handling discussion
- **Problem**: None of the processors discuss malformed data, conflicts
- **Impact**: Real-world users will hit these issues
- **Action**:
  - [ ] Add brief error handling notes per processor
  - [ ] Document what happens with identifier collisions

#### M2. No performance/scalability discussion
- **Problem**: No mention of registry lookup time, disk space, batch size limits
- **Impact**: Users don't know system limits
- **Action**:
  - [ ] Add performance notes (what's been tested at scale)
  - [ ] Mention 129 structures for MOGRN, 40K sequences for atlas

#### M3. Database coverage not quantified
- **Problem**: Says "UniProt, PDB, AlphaFold DB, NCBI" without versions/dates
- **Impact**: Reproducibility concern
- **Action**:
  - [ ] Document database versions used
  - [ ] Note update frequency policy

#### M4. Code examples are pseudocode
- **Problem**: All examples illustrative, not directly executable
- **Impact**: Users can't copy-paste to learn
- **Action**:
  - [ ] Reference the actual example scripts in thesis/processors/
  - [ ] Consider code blocks from working examples

### LOW PRIORITY

#### L1. Entity Registry collision handling vague
- **Problem**: Says "raises an error" for ambiguous aliases
- **Impact**: Users don't know how to resolve
- **Action**:
  - [ ] Document resolution strategies
  - [ ] Show example of disambiguation

#### L2. Zero-configuration claim overstates
- **Problem**: Environment variable still needed for non-default paths
- **Impact**: Not truly zero-config for all users
- **Action**:
  - [ ] Clarify defaults vs configuration
  - [ ] Show when configuration IS needed

---

## Content Available in Project

### For enriching thesis content:
1. **Validated workflows** - 20 tests passed, can cite specific results
2. **Performance data** - workflow_report.md has timing/size info
3. **Rhodozyme designs** - 31 candidates, concrete application
4. **Test coverage** - 92+ tests demonstrate robustness

### Figures available:
- Architecture diagram exists (thesis_overview_protos_v1.png)
- Cone opsin clustering outputs
- Structure alignment visualizations
- Graph comparison figures

### Missing from project (needs creation):
- [ ] Performance benchmark table
- [ ] Cross-processor workflow diagram
- [ ] Error handling documentation

---

## Recommended Chapter Restructure

```
2. ProtOS
   2.1 The Data Fragmentation Problem (expand introduction)
   2.2 Unified Data Access
   2.3 Entity System (combine registry + datasets)
   2.4 Zero Configuration
   2.5 Processors (condensed to essentials)
       - Sequence
       - Structure
       - GRN (emphasize MOGRN connection)
       - Embedding
       - Graph (emphasize LAMBDA connection)
       - Property
   2.6 Model Manager
   2.7 Toward MOGRN and LAMBDA (new transition section)

   [Move Molecule Processor to Supplementary]
```

---

## Action Items

- [ ] Read and assess each processor subchapter
- [ ] Identify redundant content to cut
- [ ] Add explicit MOGRN/LAMBDA connections
- [ ] Include performance notes from project
- [ ] Reference working examples in thesis/processors/
- [ ] Consider moving Molecule Processor
- [ ] Add cross-processor workflow example
