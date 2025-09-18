"""
CLI commands for safe file registration from input folder.
"""

import typer
from pathlib import Path
from typing import Optional
from enum import Enum

from protos.io.paths import ProtosPaths
from protos.io.input_manager import InputManager, ConflictResolutionStrategy


# Create a sub-app for input registration commands
register_input_app = typer.Typer(help="Register files from input folder")


class StrategyChoice(str, Enum):
    """Conflict resolution strategy choices for CLI."""
    skip = "skip"
    version = "version"
    replace = "replace"


@register_input_app.command("scan")
def scan_input():
    """Scan input folder and show files ready for registration."""
    manager = InputManager()
    input_files = manager.scan_input_folder()
    
    if not input_files:
        typer.echo("No files found in input folder")
        typer.echo(f"Input folder: {manager.input_dir}")
        return
    
    typer.echo(f"Found {len(input_files)} files in input folder:")
    typer.echo()
    
    # Group by file type
    by_type = {}
    for f in input_files:
        by_type.setdefault(f.file_type, []).append(f)
    
    for file_type, files in by_type.items():
        typer.echo(f"{file_type.upper()} files ({len(files)}):")
        for f in files:
            size_mb = f.size / (1024 * 1024)
            typer.echo(f"  - {f.path.name} ({size_mb:.1f} MB) -> {f.entity_name}")
        typer.echo()


@register_input_app.command("process")
def process_input(
    dry_run: bool = typer.Option(False, "--dry-run", help="Preview without making changes"),
    strategy: StrategyChoice = typer.Option(
        StrategyChoice.skip, 
        "--strategy", 
        help="Conflict resolution strategy"
    ),
    yes: bool = typer.Option(False, "-y", "--yes", help="Skip confirmation prompt")
):
    """Process and register files from input folder."""
    manager = InputManager()
    
    # Show what's in input folder
    input_files = manager.scan_input_folder()
    if not input_files:
        typer.echo("No files found in input folder")
        typer.echo(f"Input folder: {manager.input_dir}")
        return
    
    typer.echo(f"Found {len(input_files)} files to process:")
    for f in input_files:
        typer.echo(f"  - {f.path.name} ({f.file_type})")
    typer.echo()
    
    # Get confirmation unless --yes flag or dry run
    if not dry_run and not yes:
        if not typer.confirm("Proceed with registration?"):
            typer.echo("Aborted.")
            return
    
    if dry_run:
        typer.echo("DRY RUN - No changes will be made")
        typer.echo()
    
    # Process files
    conflict_strat = ConflictResolutionStrategy[strategy.value.upper()]
    report = manager.process_input_files(
        conflict_strategy=conflict_strat,
        dry_run=dry_run
    )
    
    # Show results
    report.display()
    
    if dry_run:
        typer.echo("\nThis was a DRY RUN - no files were moved or registered")
        typer.echo("Run without --dry-run to actually process files")


@register_input_app.command("show-folder")
def show_input_folder():
    """Show the path to the input folder."""
    paths = ProtosPaths()
    input_dir = Path(paths.data_root) / "input"
    
    typer.echo(f"Input folder: {input_dir}")
    
    if input_dir.exists():
        # Count files
        files = list(f for f in input_dir.iterdir() if f.is_file() and not f.name.startswith('.'))
        typer.echo(f"Contains {len(files)} files")
    else:
        typer.echo("Folder does not exist yet - it will be created when needed")


@register_input_app.command("clean")
def clean_processed(
    older_than_days: int = typer.Option(
        30, 
        "--older-than", 
        help="Remove processed files older than this many days"
    ),
    yes: bool = typer.Option(False, "-y", "--yes", help="Skip confirmation prompt")
):
    """Clean old files from processed and rejected folders."""
    import datetime
    
    paths = ProtosPaths()
    processed_dir = Path(paths.data_root) / "input" / "processed"
    rejected_dir = Path(paths.data_root) / "input" / "rejected"
    
    now = datetime.datetime.now()
    cutoff = now - datetime.timedelta(days=older_than_days)
    
    # Find old files
    old_files = []
    for dir_path in [processed_dir, rejected_dir]:
        if dir_path.exists():
            for f in dir_path.iterdir():
                if f.is_file():
                    mtime = datetime.datetime.fromtimestamp(f.stat().st_mtime)
                    if mtime < cutoff:
                        old_files.append((f, dir_path.name))
    
    if not old_files:
        typer.echo(f"No files older than {older_than_days} days found")
        return
    
    typer.echo(f"Found {len(old_files)} files older than {older_than_days} days:")
    for f, folder in old_files:
        typer.echo(f"  - {folder}/{f.name}")
    
    if not yes:
        if not typer.confirm(f"Delete {len(old_files)} old files?"):
            typer.echo("Aborted.")
            return
    
    # Delete files
    for f, _ in old_files:
        f.unlink()
    
    typer.echo(f"Deleted {len(old_files)} files")


# For backward compatibility, also provide a single command at top level
@typer.command()
def register_input(
    dry_run: bool = typer.Option(False, "--dry-run", help="Preview without making changes"),
    strategy: StrategyChoice = typer.Option(
        StrategyChoice.skip, 
        "--strategy", 
        help="Conflict resolution strategy"
    )
):
    """Process and register files from input folder."""
    # Delegate to the sub-command
    process_input(dry_run=dry_run, strategy=strategy, yes=False)