# GRN sequence fixtures

`uniprot_grn_sequences.fasta` contains full-length wild-type sequences fetched
through ProtOS `SequenceLoader` from UniProt on 2026-07-12. FASTA identifiers
are `<accession>|<entry_name>` so every fixture remains accession-traceable.

The deterministic GRN tests use these frozen records without network access.
The separately marked network suite fetches the same accessions through ProtOS,
registers them, and annotates them against the bundled tables.
