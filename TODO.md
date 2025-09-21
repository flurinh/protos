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

## Open Tasks
1. Implement registry-mediated exchange (lazy chain registration + relationship helpers) so processors can discover related entities without tight coupling.
2. Define standard alignment artifact formats (`.msa`, `.alm`) and integrate them into exporter outputs.
3. Expand tests for cross-processor interactions once registry flow is established (structure ↔ sequence roundtrip, classification annotations).
