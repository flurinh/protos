# ProtOS documentation

These pages describe the API implemented on `master`. They were reduced to
behavior that can be traced to the current source and exercised without hidden
state. Model inference and network-backed acquisition are marked explicitly.

## Data management

- [Acquisition and processing](01_unified_data_access.md)
- [Entities, registry, and datasets](02_entity_registry_datasets.md)
- [Data-root configuration](03_zero_configuration.md)

## Processors

- [Processor overview](processors/index.md)
- [StructureProcessor](processors/structure_processor.md)
- [SequenceProcessor](processors/sequence_processor.md)
- [GRNProcessor](processors/grn_processor.md)
- [PropertyProcessor](processors/property_processor.md)
- [MoleculeProcessor](processors/molecule_processor.md)
- [GraphProcessor](processors/graph_processor.md)
- [EmbeddingProcessor](processors/embedding_processor.md)

## Models

- [ModelManager](05_model_manager.md)
- [Current protein-design integrations](../src/protos/models/protein_design.md)
- [Job server and client](../src/protos/models/remote.md)

The model implementations under `src/protos/models/boltz`, `boltzgen`, and
`lambda` are pinned external repositories. Their own documentation is not part
of the ProtOS API documentation.
