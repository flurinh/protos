# GRNProcessor

Import: `from protos.processing.grn import GRNProcessor`

`GRNProcessor` reads bundled reference tables, stores generated GRN tables, and
projects reference positions onto sequences with pairwise sequence alignment.

## Reference table format

- rows are sequence/entity identifiers;
- columns are ProtOS dot-notation GRN labels such as `1.50`; and
- a populated cell is an amino-acid letter plus a 1-based sequence position,
  such as `A123`; `-` denotes no assignment.

```python
from protos.processing.grn import GRNProcessor

grn = GRNProcessor()
reference = grn.load_reference_table("gpcrdb_ref")

assert reference.index.name == "entity_name"
assert len(reference) > 0
assert GRNProcessor.parse_grn_value("A123") == ("A", 123)
assert GRNProcessor.parse_grn_value("-") is None
```

`load_reference_table(name)` reads
`<data_root>/grn/reference/<name>.csv`. It raises `FileNotFoundError` for an
unknown reference; it does not download one.

## Generated tables

| Method | Current behavior |
| --- | --- |
| `annotate_sequences(sequences, reference_table=..., protein_family=..., search="pairwise")` | return `(annotation_dataframe, summary)` |
| `record_table(name, table, ..., allow_create=False, link_entities=True)` | write CSV, register a table entity/dataset, and optionally link rows to existing entities |
| `load_table(name)` | load a generated table; return an empty DataFrame if absent |
| `list_tables()` | list generated table datasets, not bundled references |
| `build_grn_to_seq_index(table, sequence_id=...)` | map GRN labels to 1-based sequence indices |

With `link_entities=True`, row entities must already exist unless
`allow_create=True` is explicitly selected. The `protein_family` and `search`
values are recorded in the annotation summary; current matching is pairwise
alignment against all sequences in the selected reference table.

Bundle provenance and checksums are documented in
`../../src/protos/reference_data/grn/README.md` and `manifest.json`.
