# Bundled GRN reference tables

This directory is the authoritative, package-distributed GRN data bundle for
ProtOS. It is copied into `<data_root>/grn/` during first initialization. On
later initialization, the installed bundle version and CSV checksums are
verified against `manifest.json`; missing, damaged, or outdated files are
repaired from the package without a network request.

Run `python setup_project_data.py --refresh` from a source checkout to replace
that checkout's `data/grn/` copy with this bundle. The script deliberately
targets `<checkout>/data`; it does not refresh an arbitrary configured data
root. Explicit refresh is authoritative and removes unbundled CSVs and stale
configuration files. Automatic repair preserves unrelated custom CSVs.

No GRN table is downloaded from Zenodo or GPCRdb at runtime. Upstream data must
first be normalized to the ProtOS dot-notation schema and reviewed, then added
to `reference/` and recorded in `manifest.json`.

The July 2026 bundle contains only the supported current tables: type-I opsins,
WT-only type-II opsins, the aggregate and per-class GPCRdb receptor tables,
human G-alpha CGN tables, and the human arrestin CAN table. Deprecated thesis
aliases and superseded tables are intentionally absent.

GPCRdb receptor imports keep only structural segment coordinates in ProtOS dot
notation. GPCRdb loop-motif coordinates such as `45.50` and terminal
pseudo-segments `0.*` and `9.*` are excluded. During annotation, ProtOS derives
flexible-region coordinates from proximity to the flanking structural segments
and terminal coordinates as `n.<distance>` and `c.<distance>`.

`manifest.json` records the bundle version, every distributed CSV filename and
SHA-256 checksum. `gpcrdb_provenance.json` records
the pinned upstream commit, raw inputs, transformations, coverage, and output
checksums. The reproducible command and the annotation-algorithm comparison are
documented in [the upstream GRN pipeline](../../../../docs/grn-upstream.md).
`opsin_provenance.json` records the exact-sequence selection of wild-type VPOD
type-II opsins.

Large future bundles may move to a versioned external release. Until a manifest
and downloader are implemented, initialization remains entirely package-local.
