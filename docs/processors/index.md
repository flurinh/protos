# Processors

A processor owns one managed data type. Every processor receives the shared
`ProtosPaths`, registry, and a type-specific `DatasetManager` from
`BaseProcessor`.

| Processor | Managed representation | Optional runtime dependency |
| --- | --- | --- |
| [StructureProcessor](structure_processor.md) | canonical pandas DataFrame in PKL | none |
| [SequenceProcessor](sequence_processor.md) | FASTA | MMseqs2 for MMseqs operations |
| [GRNProcessor](grn_processor.md) | CSV GRN tables | none |
| [PropertyProcessor](property_processor.md) | CSV property tables and JSON index | none |
| [MoleculeProcessor](molecule_processor.md) | JSON descriptor records | RDKit for chemistry calculations |
| [GraphProcessor](graph_processor.md) | pickled NumPy graph dictionaries | PyTorch + PyG only for `to_pyg()` |
| [EmbeddingProcessor](embedding_processor.md) | pickled arrays and dataset metadata | PyTorch + Transformers for inference |

The common implemented operations are `list_entities`, `entity_exists`,
`delete_entity`, `list_datasets`, `create_dataset`, `add_to_dataset`,
`remove_from_dataset`, and `get_dataset_info`. `load_entity`, `save_entity`,
and sometimes `load_dataset` are specialized by each processor, so their exact
return values are documented on the processor pages.

Processors do not fetch remote data implicitly. Use a loader under
`protos.io.ingest` when acquisition is required.
