# GPCRdb, CGN, and CAN reference-data pipeline

ProtOS does not use the GPCRdb API at runtime. The bundled tables are generated
offline from a pinned checkout of
[`protwis/gpcrdb_data`](https://github.com/protwis/gpcrdb_data), then installed
from the package by `setup_project_data.py`.

## Exact inputs

The three numbering systems have different upstream representations:

| Numbering | Raw Git inputs | ProtOS output |
| --- | --- | --- |
| GPCR | `protein_data/proteins_and_families.txt`, `residue_data/generic_numbers/`, and `residue_data/reference_positions/` plus an explicit generated receptor export | One table per top-level class and the retained all-receptor `gpcrdb_ref.csv` |
| CGN | `g_protein_data/PDB_UNIPROT_ENSEMBLE_ALL.txt` and `CGN_lookup.csv` | All 16 human G-alpha proteins and Gs, Gi/o, Gq/11, and G12/13 subsets |
| CAN | `arrestin_data/CAN_aln.csv` | The four human arrestins |

GPCRdb does **not** commit an all-receptor residue CSV to `gpcrdb_data`.
Protwis builds that database table from sequences, curated reference-position
YAML files, numbering maps, and alignments. Consequently the importer requires
`--receptor-export`; it never mislabels the generated export as a raw Git file.
Its checksum and the raw-data commit are recorded separately in
`gpcrdb_provenance.json`.

The current receptor export is human. Class D1 is fungal, so its bundled human
table has a header and zero rows. That is explicit in the provenance statistics.

## Rebuild

Use a clean checkout at the commit recorded in `gpcrdb_provenance.json`:

```bash
python scripts/update_gpcrdb_references.py \
  --gpcrdb-data /path/to/gpcrdb_data \
  --receptor-export src/protos/reference_data/grn/reference/gpcrdb_ref.csv \
  --output-dir src/protos/reference_data/grn/reference \
  --provenance src/protos/reference_data/grn/gpcrdb_provenance.json \
  --manifest src/protos/reference_data/grn/manifest.json
```

The importer contains no HTTP client. It classifies every receptor row against
the raw hierarchy, refreshes CGN labels from `CGN_lookup.csv` by `sortColumn`
(as Protwis does), converts CGN/CAN labels to two-digit ProtOS dot notation,
checks residue/sequence-position integrity, and writes SHA-256 checksums.

## Annotation-algorithm comparison

GPCRdb/Protwis and ProtOS solve related but different problems:

| Behaviour | GPCRdb/Protwis | Current ProtOS |
| --- | --- | --- |
| Reference selection | Uses curated per-protein anchors; missing anchors are projected from a related-family template | Aligns the query against every row in the caller-selected reference table and chooses the best normalized score |
| Alignment | Clustal Omega projects template anchor positions | `SequenceAlignmentEngine` performs pairwise alignment with caller-configurable gap penalties and projects every reference-table label |
| Family knowledge | Template search walks the GPCR family hierarchy | The table name constrains candidates; `protein_family` is currently reporting metadata only |
| Structural anomalies | Scheme maps and explicit bulge/constriction handling adjust structure-based labels | Anomalies are preserved only when already encoded in the reference table |
| Confidence | Curated anchors or family-template provenance | Reports score and coverage and supports caller thresholds; benchmark-derived defaults and a best-vs-second-best margin remain future work |
| CGN/CAN | Imports curated residue mappings/alignments | Uses the same curated mappings as reference tables, then applies the generic projection algorithm to new sequences |

The safe improvements implemented now are class-specific GPCR tables, explicit
source provenance, current CGN lookup correction, arbitrary segment identifiers,
reference-derived segment ordering, native CGN/CAN validation, and explicit
insertion/deletion diagnostics with conservative insertion assignment.
Before changing scientific annotation behaviour, the next comparison should use
a held-out benchmark to choose family filtering, score/coverage thresholds,
best-hit margins, and explicit bulge/constriction tests. This avoids silently
changing thesis-era results without measured evidence.

The complete current ProtOS projection and numeric-expansion pseudocode,
including insertion/deletion handling and exceptions, is in
[the GRN annotation algorithm](grn-annotation-algorithm.md).
