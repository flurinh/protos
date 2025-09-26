# Migration Status

## Structure Processor
- [x] Alignment engine wraps CEalign/Kabsch/simple, in-place transforms.
- [x] Loader/exporter share processor, no ad-hoc instances.
- [x] CIF write/read use ProtosPaths only; outputs parse correctly.
- [x] Workflow scripts (`test_structure_*`) align/export using the dataset model.
- [x] DataFrame utilities: column filters, residue/ligand addition, atom reindex helpers while preserving canonical schema.
- [x] Chain sequence extraction opt-in flow: collect chain metadata, only register sequences when requested.
- [x] Registry bridge: on registration emit `derived_from` relationships with chain metadata (chain ID, residue span, semantic role).
- [x] Relationship lookup helpers exposed via BaseProcessor to surface related sequence entities for loaded structures/datasets.

## Sequence Processor
- [x] Paths split between `fasta/entities` and `fasta/datasets` with ProtosPaths.
- [x] Loader can ingest local or UniProt sequences, with lazy entity registration.
- [x] Alignment logic moved under `analysis.sequence`; engine exposes Biopython & MMseqs.
- [x] Exporter filters by sequence IDs; scripts demonstrate dataset + alignment flow.
- [x] Relationship lookup helpers mirroring structure side to resolve source structures for registered sequences.
- [x] Mutant library generation plus conservation/linkage analytics for engineering workflows.
- [ ] Example workflow: after chain registration, align sequences and classify chains, annotating source structures with results.

## GRN Processor
- [ ] Integrate sequence-driven GRN assignment so SequenceProcessor can align datasets to reference tables and hand recordings to GRNProcessor without leaving the registry workflow.
- [ ] Expose GRN utilities (alignment, annotation, cleaning) through SequenceProcessor APIs to replace legacy scripts in `analysis/grn/`.
- [ ] Ensure GRN tables are stored as datasets with CSV artifacts and relationships tying tables back to the originating sequence entities for discovery via the registry.

### GRN Annotation Usage (current state)

**Prerequisites**
- Install and expose `mmseqs` on PATH. When it is missing (e.g. Windows native shells), the workflow falls back to BioPython pairwise scoring and logs warnings. Until we ship an auto-detection shim, users must either install `mmseqs` or expect the slower fallback path.
- Keep the GRN reference tables under `data/grn/reference/` managed by ProtosPaths (e.g. `gpcrdb_ref.csv`). The reference name and the protein family (`gpcr_a`, `mo`, …) must be supplied together.
- Sequence datasets must already be registered via `SequenceProcessor` (e.g. `gpcr_sequences`). Chain IDs should remain stable so relationships back to structures resolve correctly.

**Workflow**
1. Create a `SequenceProcessor` and (optionally) load a reference sequence to classify query chains before annotation. In our GPCR demo we compare all chains against `5d5a_chain_A` and mark anything with a normalised score ≥0.35 as `gpcr_like`.
2. Call `SequenceProcessor.annotate_with_grn(dataset_name="gpcr_sequences", reference_table="gpcrdb_ref", protein_family="gpcr_a", output_table="gpcr_grn")`. The helper handles dataset loading + GRNProcessor annotation so callers stay in the registry workflow.
3. Internally we first attempt an MMseqs search to pick the closest reference row, then fall back to BioPython if MMseqs is unavailable. The annotation step still uses the legacy `expand_annotation` helpers; coverage metadata (`assigned_positions`, `coverage_fraction`) is recorded per sequence for visibility.
4. Results are written through `GRNProcessor.record_table`, producing a CSV artefact under `data/grn/datasets/<table>.csv` plus dataset metadata that links back to the source sequence dataset.

**Known gaps**
- `expand_annotation` throws `list index out of range` for several GPCR references; the processor falls back to seed mappings. We need to upstream fixes (loop handling, missing GRN columns) before calling the workflow “complete”.
- Windows shells without WSL `mmseqs` binaries emit noisy “[WinError 2]” logs. Add an early capability probe that downgrades to Biopython without printing the stack.
- Provide utilities that summarise coverage/quality so users can quickly identify sequences annotated only by the seed fallback.

## Graph Processor
- [ ] Introduce a `GraphProcessor` (BaseProcessor-derived) that persists PyTorch Geometric graph artefacts under `graph/graphs` and registers them as entities with provenance metadata (source structure ID, featurisation settings, PyG version).
- [ ] Implement structure→graph translation utilities (atom- and residue-level) with pluggable neighbour rules (distance cutoffs, k-NN) while recording those parameters alongside each entity.
- [ ] Mirror the registry bridge: register a `derived_from` relationship from each graph to its originating structure so datasets can be traced back.
- [ ] Provide dataset helpers to bulk-generate graphs from a StructureProcessor dataset, saving node/edge tensors plus summary statistics (node count, edge density) for health checks.
- [ ] Add a workflow script (guarded when PyG is missing) that demonstrates loading a structure dataset, building graphs, and enumerating basic graph metrics.

## Ligand Processor
- [ ] Finalise RDKit integration for descriptor generation, fingerprint caching, and SDF round-tripping (currently optional and warning-heavy).
- [ ] Expand structure-ligand extraction to classify ligand types (cofactors, ions) and support residue-level filtering.
- [ ] Extend `compute_interactions` beyond simple distance cut-offs to classify interaction types (H-bonds, π–π, hydrophobic) and record them via `PropertyProcessor`.
- [ ] Wire structure loading pipeline so StructureProcessor can request ligand extraction eg. during complex annotation workflows.
- [ ] Add loaders for external datasets (QM9, ChEMBL, Enamine) that register ligands through the new API while capturing provenance metadata.
- [ ] Plan downstream QSAR integration: property tables from LigandProcessor feeding into StructureProcessor analyses or dedicated modelling pipelines.

## Open Tasks
1. Implement registry-mediated exchange (lazy chain registration + relationship helpers) so processors can discover related entities without tight coupling.
2. Define standard alignment artifact formats (`.msa`, `.alm`) and integrate them into exporter outputs.
3. Expand tests for cross-processor interactions once registry flow is established (structure ↔ sequence roundtrip, classification annotations).
