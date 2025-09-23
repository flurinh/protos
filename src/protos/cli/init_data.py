"""Backward-compatible data initialization helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import typer

from . import _data


def init_data_directory(force: bool = False, path: Optional[str] = None) -> Dict[str, object]:
    """Programmatic API used during package installation."""

    target = _data.resolve_target_path(Path(path) if path else None)
    if target.exists() and force:
        _data.remove_data_root(target)
    summary = _data.initialize_data_root(target)

    # Provide legacy keys expected by setup scripts
    legacy_stats = {
        "directories_created": summary.get("directory_count", 0),
        "registries_created": 1 if summary.get("reference_installed") else 0,
        "reference_files_copied": 1 if summary.get("reference_installed") else 0,
    }
    legacy_stats.update(summary)
    return legacy_stats


def _cli(force: bool = False, path: Optional[str] = None) -> None:
    summary = init_data_directory(force=force, path=path)
    typer.echo("Protos data directory initialized:")
    typer.echo(f"  data_root: {summary['data_root']}")
    typer.echo(f"  directories: {summary['directory_count']}")
    typer.echo(f"  files: {summary['file_count']}")
    typer.echo(
        "  reference data: "
        + ("installed" if summary.get("reference_installed") else "not installed")
    )


def init_data() -> None:
    """Console entry point registered as ``protos-init``."""
    typer.run(_cli)
