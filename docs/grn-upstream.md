# GPCRdb, opsin, CGN, and CAN reference-data pipeline

ProtOS does not use the GPCRdb API at runtime. The bundled tables are generated
offline from a pinned checkout of
[`protwis/gpcrdb_data`](https://github.com/protwis/gpcrdb_data), then installed
from the package by `setup_project_data.py`.

## Exact inputs

The three numbering systems have different upstream representations:

| Numbering | Raw Git inputs | ProtOS output |
| --- | --- | --- |
| GPCR | `protein_data/proteins_and_families.txt`, `residue_data/generic_numbers/`, and `residue_data/reference_positions/` plus an explicit generated receptor export | One table per top-level class and normalized all-receptor `gpcrdb_ref.csv`; upstream loop motifs and 0/9 terminal pseudo-segments are removed |
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
  --receptor-export /path/to/generated/unmodified_gpcrdb_export.csv \
  --output-dir src/protos/reference_data/grn/reference \
  --provenance src/protos/reference_data/grn/gpcrdb_provenance.json \
  --manifest src/protos/reference_data/grn/manifest.json
```

The importer contains no HTTP client. It classifies every receptor row against
the raw hierarchy, refreshes CGN labels from `CGN_lookup.csv` by `sortColumn`
(as Protwis does), converts CGN/CAN labels to two-digit ProtOS dot notation,
checks residue/sequence-position integrity, and writes SHA-256 checksums.

The importer rejects GPCRdb loop-motif labels such as `45.50` and terminal
pseudo-segments `0.*`/`9.*`. Those labels do not have ProtOS structural
semantics. ProtOS instead generates directional loop coordinates from the two
flanking structural segments and assigns `n.<distance>`/`c.<distance>` from the
nearest terminal structural segment.

## Opsin tables

`scripts/update_opsin_references.py` normalizes the type-I table and selects the
type-II reference from an already annotated full VPOD table. The type-II rows
are selected by exact ungapped sequence equality against VPOD's WT FASTA. The
matching VPOD metadata must contain unique accessions and sequences and no
non-empty `Mutations` values. The accessions replace the source table's opaque
row hashes. The current bundle uses the VPOD 1.3 WT FASTA and metadata; these
contain the same 364 WT protein sequences as VPOD 1.2.

The command requires an unmodified, full type-II GRN source table because the
distributed table is intentionally WT-only:

```bash
python scripts/update_opsin_references.py \
  --type-i src/protos/reference_data/grn/reference/type_I_opsins.csv \
  --type-ii-all /path/to/full_type_II_opsins.csv \
  --type-ii-wt-fasta /path/to/VPOD_1.3/wt_aligned_VPOD_1.3_het.fasta \
  --type-ii-wt-metadata /path/to/VPOD_1.3/wt_meta.tsv \
  --source-version VPOD_1.3 \
  --source-revision <vpod-git-commit> \
  --output-dir src/protos/reference_data/grn/reference \
  --provenance src/protos/reference_data/grn/opsin_provenance.json
```

The provenance records checksums for the full GRN source table, WT FASTA,
metadata, and both output tables. GPCR-style directional loop columns are
removed from type-II for the same reason they are removed from GPCRdb imports.

## Annotation-algorithm comparison

GPCRdb/Protwis and ProtOS solve related but different problems:

| Behaviour | GPCRdb/Protwis | Current ProtOS |
| --- | --- | --- |
| Reference selection | Uses curated per-protein anchors; missing anchors are projected from a related-family template | Aligns the query against every row in the caller-selected reference table and chooses the best normalized score |
| Alignment | Clustal Omega projects template anchor positions | `SequenceAlignmentEngine` performs pairwise alignment with caller-configurable gap penalties and projects every reference-table label |
| Family knowledge | Template search walks the GPCR family hierarchy | The table name constrains candidates; `protein_family` is currently reporting metadata only |
| Structural anomalies | Scheme maps and explicit bulge/constriction handling adjust structure-based labels | Encoded segment insertions are preserved; new insertions and flexible regions are detected from the alignment |
| Confidence | Curated anchors or family-template provenance | Reports score and coverage and supports caller thresholds; benchmark-derived defaults and a best-vs-second-best margin remain future work |
| CGN/CAN | Imports curated residue mappings/alignments | Uses the same curated mappings as reference tables, then applies the generic projection algorithm to new sequences |

The safe improvements implemented now are class-specific GPCR tables, explicit
source provenance, current CGN lookup correction, arbitrary segment identifiers,
reference-derived segment ordering, native CGN/CAN validation, ProtOS-native
loop/terminal coordinates, and explicit insertion/deletion diagnostics.
Before changing scientific annotation behaviour, the next comparison should use
a held-out benchmark to choose family filtering, score/coverage thresholds,
best-hit margins, and explicit bulge/constriction tests. This avoids silently
changing thesis-era results without measured evidence.

The complete current ProtOS projection and numeric-expansion pseudocode,
including insertion/deletion handling and exceptions, is in
[the GRN annotation algorithm](grn-annotation-algorithm.md).
