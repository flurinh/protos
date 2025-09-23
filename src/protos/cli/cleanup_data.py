"""Backward-compatible data cleanup helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from . import _data


def cleanup_data_directory(path: Optional[str] = None, reinitialize: bool = False) -> None:
    target = _data.resolve_target_path(Path(path) if path else None)
    summary = _data.clear_data_root(target, reinitialize=reinitialize)
    if reinitialize:
        typer.echo("Data directory cleared and reinitialized:")
        typer.echo(f"  data_root: {summary['data_root']}")
        typer.echo(f"  directories: {summary['directory_count']}")
        typer.echo(f"  files: {summary['file_count']}")
    else:
        typer.echo(f"Data directory removed: {summary['data_root']}")


def _cli(
    path: Optional[str] = typer.Option(None, "--path", "-p", help="Data directory to clear"),
    reinit: bool = typer.Option(False, "--reinit", help="Recreate the directory after deletion."),
    force: bool = typer.Option(False, "--force", help="Skip the confirmation prompt."),
) -> None:
    target = _data.resolve_target_path(Path(path) if path else None)
    if not force:
        confirm = typer.confirm(
            f"This will delete all data under {target}. Continue?", default=False
        )
        if not confirm:
            typer.echo("Aborted.")
            raise typer.Exit(code=1)
    cleanup_data_directory(path=path, reinitialize=reinit)


def cleanup_data() -> None:
    """Console entry point registered as ``protos-cleanup``."""
    typer.run(_cli)
