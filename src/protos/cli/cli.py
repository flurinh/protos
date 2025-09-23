"""Typer-based CLI for managing Protos data directories."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import typer

from . import _data

app = typer.Typer(help="Command-line utilities for Protos data management.")


def _reset_paths_singleton() -> None:
    """Reset the cached ProtosPaths singleton so a fresh instance is created."""
    _data.reset_paths_singleton()


def _resolve_target_path(path_option: Optional[Path]) -> Path:
    return _data.resolve_target_path(path_option)


def _collect_summary(data_root: Path) -> Dict[str, object]:
    return _data.collect_summary(data_root)


def _initialize_data_root(target_root: Path) -> Dict[str, object]:
    """Create (or refresh) the Protos data directory at ``target_root``."""
    return _data.initialize_data_root(target_root)


@app.command()
def init(
    path: Optional[Path] = typer.Option(
        None,
        "--path",
        "-p",
        help="Target directory for Protos data. Defaults to ~/protos_data or PROTOS_DATA_ROOT.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Remove existing data directory before initializing.",
    ),
) -> None:
    """Initialize a Protos data directory."""

    target_root = _resolve_target_path(path)

    if target_root.exists() and force:
        typer.echo(f"Removing existing data directory: {target_root}")
        _data.remove_data_root(target_root)

    if target_root.exists() and not force:
        typer.echo(f"Data directory already exists at {target_root}. Ensuring structure is complete.")
    else:
        typer.echo(f"Initializing data directory at {target_root}")

    summary = _initialize_data_root(target_root)

    typer.echo("\nInitialization summary:")
    typer.echo(f"  data_root: {summary['data_root']}")
    typer.echo(f"  processors: {', '.join(summary['processors']) if summary['processors'] else 'none'}")
    typer.echo(f"  directories: {summary['directory_count']}")
    typer.echo(f"  files: {summary['file_count']}")
    typer.echo(
        "  reference data: "
        + ("installed" if summary['reference_installed'] else "not installed (marker missing)")
    )


@app.command()
def clear(
    path: Optional[Path] = typer.Option(
        None,
        "--path",
        "-p",
        help="Data directory to clear. Defaults to the active Protos data root.",
    ),
    force: bool = typer.Option(
        False,
        "--force",
        help="Skip confirmation before deleting the data directory.",
    ),
    reinit: bool = typer.Option(
        True,
        "--reinit/--no-reinit",
        help="Reinitialize the directory after deletion (default: reinitialize).",
    ),
) -> None:
    """Delete the Protos data directory (and optionally reinitialize it)."""

    target_root = _resolve_target_path(path)

    if not target_root.exists():
        typer.echo(f"No data directory found at {target_root}.")
        if reinit:
            typer.echo("Creating a fresh data directory.")
            summary = _initialize_data_root(target_root)
            typer.echo(f"Reinitialized data at {summary['data_root']}")
        return

    if not force:
        confirm = typer.confirm(
            f"This will delete all data under {target_root}. Continue?", default=False
        )
        if not confirm:
            typer.echo("Aborted.")
            raise typer.Exit(code=1)

    summary = _data.clear_data_root(target_root, reinitialize=reinit)
    typer.echo(f"Removed data directory: {target_root}")

    if reinit:
        typer.echo("\nReinitialized data directory.")
        typer.echo(f"  data_root: {summary['data_root']}")
        typer.echo(f"  processors: {', '.join(summary['processors']) if summary['processors'] else 'none'}")
        typer.echo(f"  directories: {summary['directory_count']}")
        typer.echo(f"  files: {summary['file_count']}")
    else:
        typer.echo("Singleton reset. Future Protos sessions will recreate data on demand.")
