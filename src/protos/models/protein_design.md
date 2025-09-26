# Protein Design: State-of-the-Art Workflows

## Core Design Loop (Practical & Proven)

### 1. Backbone/Scaffold Generation (or Inpainting)

Use **RFdiffusion** to generate whole backbones de novo, to inpaint around a fixed motif/active site, or to build symmetric oligomers and interfaces. It's a generative structure model (diffusion fine-tuned from RoseTTAFold) that natively supports constraints:
- Hotspots/epitopes
- Symmetry constraints
- Distance/angle restraints

*Watson et al., Nature, 2023*

### 2. Sequence Design on Fixed Backbone

Use **ProteinMPNN** (inverse folding) to assign sequences that are designable for your backbones. For cofactors/ligands/nucleic acids, use **LigandMPNN** so the sequence is conditioned on non-protein atoms. These are fast, robust, and experimentally validated upgrades over classic Rosetta-only flows.

*Dauparas et al., Science, 2022*

### 3. In-Silico Validation & Triage

Re-predict structures (AF2/ESMFold/RoseTTAFold) and filter on:
- **Structure quality**: pLDDT/PAE/interface metrics
- **Physical checks**: clash/packing analysis
- **Functional validation**: 
  - For binders: epitope/contact goals
  - For enzymes: catalytic geometry (distances/angles/planarity)

This validation step is the workhorse in all recent experimental papers. *Watson et al., Nature, 2023*

### 4. Targeted Refinement

- **Partial diffusion** ("diffusion local edits") to tweak loops/interfaces while keeping the rest fixed
- Run short MPNN→AF2 cycles to climb local optima
- For binders: re-guide with explicit hotspot residue constraints to hit desired epitopes with shape complementarity

*Watson et al., Nature, 2023*

### 5. Wet-Lab Screening

- **Binders**: Yeast display/SPR/BLI
- **Enzymes**: Activity assays
- Iterate the loop on the most promising designs

## Binders: What Works Best (2023-2025)

### Diffusion-First Interfaces

RFdiffusion can generate high-quality binder backbones conditioned on:
- Target surface
- Guide residues ("hotspots")

This is now the standard entry point for de novo binders. Follow with ProteinMPNN for sequences and AF2 triage; success rates have jumped vs. purely energy-based protocols.

*Watson et al., Nature, 2023*

### Epitope Targeting & Partial Diffusion

For flexible/peptidic epitopes (e.g., helical peptides, pMHC), studies demonstrate de novo high-specificity binders by:
- Steering the contact map during diffusion
- Polish with targeted refinement

*Watson et al., Nature, 2023*

## Enzymes: Current Best Practices

### Functional-Site Scaffolding

Two families of methods:

1. **Constrained hallucination / site completion** (pre-diffusion era)
   - Optimize a sequence so the predicted fold contains the catalytic motif
   - Geometry is part of the loss function
   - *Various early enzyme design papers, PMC*

2. **Diffusion-based scaffolding** (RFdiffusion) - **Modern Default**
   - Condition generation on an active-site motif ("theozyme") at atomic resolution
   - Then MPNN→AF2 loop to keep the motif intact and encodable
   - Generalized across multiple enzyme classes
   - *Watson et al., Nature, 2023*

### Emerging Trends

Pushing to atomic-level active-site fidelity (emerging RFdiffusion2 work) to:
- Reduce gap between designed and catalytic geometries
- Improve hit rates with far fewer screened sequences

*RFdiffusion2 preprint, 2024*

## Newer Directions to Watch

### Joint Sequence-Structure Generation

No hard "backbone→sequence" split. Example: **ProteinGenerator** (RF-based diffusion in sequence space with coupled structure prediction) for multistate/functional objectives. Useful when sequence-specific properties need optimization during generation.

*Nature Biotechnology, 2024*

### Data-Efficient Fitness Models

Local sequence optimization around a lead using:
- Supervised regressors over DMS or task-specific data
- Combine with generative backbone/sequence pipeline for exploitation around hits

*Recent reviews in Cell and other journals*

## Concrete Production-Ready Pipeline

### 1. Specification

**Binders**:
- Target structure
- Epitope/hotspots
- Symmetry (if any)

**Enzymes**:
- Catalytic motif (theozyme)
- Residue identities
- Geometry constraints

### 2. Backbone Generation / Inpainting

```python
# RFdiffusion with constraints
backbone_configs = {
    "epitope_mask": epitope_atoms,      # For binders
    "active_site": theozyme_atoms,      # For enzymes
    "n_samples": 1000,                  # Sample broadly
    "symmetry": "C3",                   # If applicable
}
```

*Watson et al., Nature, 2023*

### 3. Sequence Design

```python
# ProteinMPNN for standard cases
# LigandMPNN for cofactors/substrates
sequence_configs = {
    "temperature": [0.1, 0.2, 0.3],     # Multiple seeds
    "n_seqs_per_backbone": 8,
    "ligand_atoms": cofactor_coords,    # If using LigandMPNN
}
```

Filter by:
- MPNN log-likelihood/entropy
- Basic packing metrics

*Dauparas et al., Science, 2022*

### 4. Structure Re-prediction & Filtering

**AF2/ESMFold validation criteria**:
- ✓ Folds to intended structure/interface (high pLDDT, low PAE)
- ✓ Preserves geometric constraints (dist/angle RMSD vs. target motif)
- ✓ Clean stereochemistry

### 5. Interface/Active-Site Tightening

- Local (partial) diffusion on problematic regions
- Re-run MPNN → Re-predict with AF2
- **Binders**: Reweight hotspot contacts
- **Enzymes**: Penalize deviations in catalytic geometry

*Watson et al., Nature, 2023*

### 6. Selection & Experimental Plan

- Diversity pick (cluster by interface geometry/sequence)
- Add stabilizing mutations if needed
- Plan expression & assays

## Integration with Protos: Minimal Friction Implementation

### 1. **RFdiffusion Integration**
```python
# Default backbone/inpainting generator with constraint hooks
rf_adapter = RFdiffusionAdapter()
rf_adapter.generate(
    epitope_atoms=epitope_coords,      # For binders
    theozyme_atoms=active_site_coords, # For enzymes
    constraints=user_constraints
)
```

### 2. **Sequence Design Backend**
- Keep ProteinMPNN/LigandMPNN as sequence backends
- Standardized interface for both

### 3. **Validation Metrics as First-Class Outputs**
```python
validation_metrics = {
    "structure_quality": {
        "pLDDT_threshold": 80,
        "PAE_threshold": 5,
        "interface_ΔPAE": 3
    },
    "functional_metrics": {
        "hotspot_contact_recovery": 0.8,
        "catalytic_distance_error": 0.5,  # Angstroms
        "catalytic_angle_error": 5         # Degrees
    }
}
```

### 4. **Partial-Diffusion Tools**
- Provide tools for local redesign loops
- Enable iterative refinement workflows

### 5. **Workflow Presets**
```python
# Encapsulate steps into callable workflows
binder_preset = ProteinDesignWorkflow(
    mode="binder",
    target_structure=target_pdb,
    epitope_residues=[45, 46, 47, 89, 90],
    validation_strictness="high"
)

enzyme_preset = ProteinDesignWorkflow(
    mode="enzyme",
    theozyme=catalytic_motif,
    substrate=substrate_mol,
    geometry_tolerance="tight"
)
```

This architecture provides a modern, production-ready protein design pipeline that integrates seamlessly with Protos while maintaining flexibility for both standard and custom workflows.

## Key References

### Core Methods
- **RFdiffusion**: Watson JL, et al. De novo design of protein structure and function with RFdiffusion. Nature. 2023 Jul;620(7976):1089-1100.
- **ProteinMPNN**: Dauparas J, et al. Robust deep learning based protein sequence design using ProteinMPNN. Science. 2022 Sep 15;378(6615):49-56.
- **LigandMPNN**: Dauparas J, et al. Atomic context-conditioned protein sequence design using LigandMPNN. Nature Methods. 2025.

### Related Tools
- **ProteinGenerator**: Multistate and functional protein design using RoseTTAFold sequence space diffusion. Nature Biotechnology. 2024.
- **RFdiffusion2**: Computational enzyme design by catalytic motif scaffolding. bioRxiv preprint. 2024.
- **AlphaFold2**: Jumper J, et al. Highly accurate protein structure prediction with AlphaFold. Nature. 2021.

### Reviews and Perspectives
- Building Enzymes through Design and Evolution. ACS Catalysis. 2023.
- Generative models for protein sequence modeling: recent advances and future directions. Briefings in Bioinformatics. 2023.
- De novo protein design—From new structures to programmable functions. Cell. 2024.