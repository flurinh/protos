# Entities, registry, and datasets

## Entity registry

`get_registry()` returns the process-wide `EntityRegistry`. Public lookups use
human-readable names. The persisted registry assigns an internal UUID so a
dataset can continue to resolve an entity after a rename.

An entity may have several formats. Each format record contains its own path
and metadata. Relationships connect entity UUIDs internally but are created and
queried with names.

```python
from protos.io.core import get_registry
from protos.processing.sequence import SequenceProcessor

sequences = SequenceProcessor()
sequences.save_entity("demo_protein", "ACDEFGHIK")

registry = get_registry()
assert registry.entity_exists("demo_protein", "sequence")
assert registry.get_entity_formats("demo_protein") == ["sequence"]
```

The main registry operations are:

| Method | Result |
| --- | --- |
| `register_entity(name, format_type, file_path, metadata=None)` | register or update one format |
| `find_entity(name, format_type=None)` | return `EntityInfo` or `None` |
| `list_entities(format_type=None)` | list names |
| `get_entity_formats(name)` | list registered formats |
| `add_alias(name, alias)` | add another lookup name |
| `rename_entity(old_name, new_name)` | rename while retaining the UUID |
| `add_relationship(source_name, target_name, rel_type, metadata=None)` | add a typed edge |
| `get_relationships(name, rel_type=None, direction="both")` | query relationship payloads |

The defined relationship types are `derived_from`, `subset_of`, `merged_from`,
`version_of`, `aligned_to`, and `annotated_by`. Both endpoints must already be
registered.

## Datasets

A dataset is a JSON definition under
`<data_root>/<processor>/datasets/<name>.json`. It stores entity names and UUIDs,
processor type, metadata, and timestamps. Creating a dataset warns about and
omits names that are not registered.

Processors expose the common dataset operations from `BaseProcessor`:

```python
from protos.processing.sequence import SequenceProcessor

sequences = SequenceProcessor()
sequences.save_entity("seq_a", "ACDE")
sequences.save_entity("seq_b", "FGHIK")
sequences.create_dataset("demo_sequences", ["seq_a", "seq_b"])

loaded = sequences.load_dataset("demo_sequences")
assert loaded == {"seq_a": "ACDE", "seq_b": "FGHIK"}

sequences.remove_from_dataset("demo_sequences", ["seq_b"])
info = sequences.get_dataset_info("demo_sequences")
assert [entity["name"] for entity in info["entities"]] == ["seq_a"]
```

`DatasetManager.load_dataset()` returns the raw JSON definition.
Processor-level `load_dataset()` methods return loaded data and may specialize
the return type. In particular, `StructureProcessor.load_dataset()` returns a
stacked DataFrame by default and accepts `return_format="dict"` for a mapping.
