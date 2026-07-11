# PropertyProcessor

Import: `from protos.processing.property import PropertyProcessor`

A property table is a CSV whose rows contain a required `scope`. A scope is a
non-empty list of `{format, name}` mappings identifying the entities annotated
by the row.

```python
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor

sequences = SequenceProcessor()
sequences.save_entity("opsin", "ACDEFGHIK")

properties = PropertyProcessor()
properties.record_properties(
    "measurements",
    [
        {
            "scope": [{"format": "sequence", "name": "opsin"}],
            "value": 1.25,
            "unit": "a.u.",
        }
    ],
)

rows = properties.load_dataset_rows(
    "measurements", "opsin", format_type="sequence"
)
assert rows.loc[0, "value"] == 1.25
```

By default, `record_properties()` registers one dataset-level property entity
and maintains a JSON row index. `materialize_entries=True` instead registers a
separate property entity for every row. Referenced entities must exist unless
`allow_create=True` is selected; that option creates placeholders.

| Method | Current behavior |
| --- | --- |
| `record_properties(table_name, rows, ...)` | append validated rows and update registry/index/dataset |
| `load_table(table_name)` | return a copy of the full DataFrame |
| `load_dataset_rows(table_name, entity_name=None, format_type=None)` | use the row index, with a scope-scan fallback |
| `get_properties(entity_name, table_name=None)` | find rows by table or incoming `annotated_by` relationships |
| `list_tables()` | list property dataset names |

`load_entity()` is meaningful for materialized row entities. With the default
dataset-level representation, query the table or dataset rows instead.
