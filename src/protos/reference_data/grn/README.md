# Bundled GRN reference tables

This directory is the authoritative, package-distributed GRN data bundle for
ProtOS. It is copied into `<data_root>/grn/` during first initialization. Run
`python setup_project_data.py --refresh` from a source checkout to replace the
installed copy with the bundle from that checkout.

No GRN table is downloaded from Zenodo or GPCRdb at runtime. Upstream data must
first be normalized to the ProtOS dot-notation schema and reviewed, then added
to `reference/` and recorded in `manifest.json`.

The July 2026 bundle contains the latest curated tables available at the thesis
handoff. `type_II_opsins.csv` and `vpod1_2.csv` are identical; the latter is
retained as a compatibility alias for existing analyses.

Large future bundles may move to a versioned external release. If that happens,
the manifest should retain filenames, checksums, source versions, and download
URLs so installation remains deterministic.
