# Bundled GRN reference tables

This directory is the authoritative, package-distributed GRN data bundle for
ProtOS. It is copied into `<data_root>/grn/` during first initialization.

Run `python setup_project_data.py --refresh` from a source checkout to replace
that checkout's `data/grn/` copy with this bundle. The script deliberately
targets `<checkout>/data`; it does not refresh an arbitrary configured data
root.

No GRN table is downloaded from Zenodo or GPCRdb at runtime. Upstream data must
first be normalized to the ProtOS dot-notation schema and reviewed, then added
to `reference/` and recorded in `manifest.json`.

The July 2026 bundle contains the latest curated tables available at the thesis
handoff. `type_II_opsins.csv` and `vpod1_2.csv` are identical; the latter is
retained as a compatibility alias for existing analyses.

`manifest.json` records the bundle version, every distributed CSV filename and
SHA-256 checksum, and compatibility aliases. It does not currently contain
upstream release identifiers or URLs; those must be added when the planned
GPCRdb multi-family import pipeline is implemented.

Large future bundles may move to a versioned external release. Until a manifest
and downloader are implemented, initialization remains entirely package-local.
