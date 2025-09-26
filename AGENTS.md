# Repository Guidelines

## Project Structure & Module Organization
Source code sits in `src/protos`: `core/` for processors scaffolding, `io/` for paths/registry/input, `processing/` for domain logic, `analysis/structure/` for shared analytics, and Typer CLIs under `cli/`. Data is in `resources/` and `demo_structures/`. Tests are in `tests/` with migration coverage in `migration_tests/`. Mirror this layout so ProtosPaths auto-discovers assets.

## Build, Test, and Development Commands
- `pip install -e .` — editable install.
- `pytest` — main test suite.
- `pytest migration_tests` — checks registry/dataset/ProtosPaths migrations.
- `pytest -m "not integration and not slow"` — quick smoke run.
- `protos init --path ./data` — scaffold zero-config data roots; pair with `protos clear` to reset.
- `python setup_project_data.py --refresh` — reload packaged reference data when fixtures drift.

## Entity & Data Workflow
ProtosPaths is the lone authority for layout. `EntityRegistry` issues UUIDs while exposing only human names, and `DatasetManager` keeps identifiers. Structure data is PKL-canonical—CIF exists only for ingest/export—and InputManager handles registration via `data/input/{pending,processed,rejected}`. GRN and embedding processors still rely on deprecated `GlobalRegistry` helpers; queue refactors per `PROCESSING_DEPRECATED_CODE_SUMMARY.md`.

## Coding Style & Naming Conventions
Use 4-space indentation and PEP 8 naming (`snake_case` modules/functions, `PascalCase` classes). Format with `black` (88 columns) and `isort` (`profile=black`). Mypy runs with `disallow_untyped_defs`, so type every new public function. Never hardcode paths—call ProtosPaths helpers—and store metadata via registries instead of stray JSON files. Document new CLI/processor APIs briefly and update `__all__` when exposing symbols.

## Testing Guidelines
Pytest drives validation. Place unit tests in `tests/<area>/test_<feature>.py`, reuse fixtures from `tests/conftest.py`, and flag slow or external cases with the existing markers. Touch `migration_tests/` whenever registry, dataset, or ProtosPaths logic changes to guard UUID/name invariants. Add workflow checks under `tests/test_workflows/`, and keep fixtures in `tests/test-data/` with names.

## Commit & Pull Request Guidelines
Commit history is terse; prefer clearer imperative subjects (`processor: add mmCIF cache purge`) and limit bodies to focused rationale. Group changes logically and reference issue IDs when available. Pull requests must describe the motivation, summarize behaviour changes, highlight data migrations, and link screenshots or CLI transcripts when relevant. Confirm tests in the PR description (`pytest` plus targeted markers) and note any follow-up tasks.

## Data & Configuration Tips
Let ProtosPaths create and lock the data root; override via `PROTOS_DATA_ROOT` instead of path literals. Use `protos init`/`protos clear` for fast environment resets. Register assets through loaders so InputManager validates and archives them. Follow the `data/` tree (`structure/cache/`, `sequence/fasta/`, `grn/tables/`, etc.) and purge leftovers in `data/temp/` so health checks stay clean. Large property tables should rely on dataset-level registration (`record_properties(..., materialize_entries=False)`) and hydrate rows via `load_dataset_rows()`. For GRNs, call `SequenceProcessor.annotate_with_grn(...)` so annotations are persisted under `data/grn/tables/` without inflating the registry.
