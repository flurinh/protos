# SequenceProcessor

Import: `from protos.processing.sequence import SequenceProcessor`

Single sequence entities are stored as FASTA files. `save_entity()` accepts a
sequence string or a one-entry mapping; a mapping with several entries raises
`ValueError` and should instead be passed to `save_sequences()`.

```python
from protos.processing.sequence import SequenceProcessor

sp = SequenceProcessor()
sp.save_entity("wild_type", "ACDEFGHIK")

assert sp.load_entity("wild_type") == "ACDEFGHIK"

sp.create_dataset("demo_sequences", ["wild_type"])
assert sp.load_dataset("demo_sequences") == {"wild_type": "ACDEFGHIK"}

mutant = sp.mutate_sequence("ACDEFGHIK", ["A1V"], "variant")
assert mutant == "VCDEFGHIK"
```

## Implemented operations

| Group | Methods |
| --- | --- |
| entity I/O | `load_entity`, `save_entity`, `export_entity`, `get_sequence` |
| multi-FASTA and datasets | `load_sequences`, `save_sequences`, `load_dataset`, `export_dataset` |
| pairwise alignment | `align_sequences`, `find_best_match`, `align_and_record` |
| mutation | `mutate_sequence`, `generate_variants`, `create_mutant_library` |
| analysis | `compute_conservation`, `compute_linkage` |
| GRN | `annotate_with_grn`, `generate_mutants_for_sequence` |
| MMseqs2 | `align_sequences_mmseqs`, `create_mmseqs_database`, `list_mmseqs_databases` |

`align_sequences()` uses Biopython and returns `(score, formatted_alignment)`.
MMseqs methods require an `mmseqs` executable and may raise
`MMseqsUnavailableError`; ProtOS does not install MMseqs2.

Use `SequenceLoader` for local FASTA or UniProt acquisition. A UniProt request
is explicit and network-backed; `load_entity()` never queries UniProt.
