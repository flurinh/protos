#!/usr/bin/env python3
"""
Initialize Protos data directory structure and reference data.

This command sets up:
- Directory structure for all processors
- Empty registries
- Reference data from src/protos/reference_data
"""

import os
import shutil
import json
import click
from pathlib import Path
from typing import Optional, Dict, Any
import logging

from protos.io.paths.path_config import ProtosPaths


logger = logging.getLogger(__name__)

# protos/cli/init.py
import typer
from pathlib import Path

init_app = typer.Typer(help="Initialization commands for protos.")

@init_app.command()
def folders(path: Path = typer.Option("./data", "--path", "-p", help="Path where data folders will be created")):
    """Set up data folders."""
    typer.echo(f"Setting up data folders at: {path}")
    path.mkdir(parents=True, exist_ok=True)
    # Add your folder creation logic here
    typer.echo("Data folders initialized successfully!")

# Example: Another subcommand under init, with its own params
@init_app.command()
def config(force: bool = typer.Option(False, "--force", "-f", help="Force overwrite existing config")):
    """Initialize configuration."""
    typer.echo(f"Setting up config (force: {force})...")
    # Your logic here


def create_directory_structure(data_root: Path) -> int:
    """
    Create the complete Protos directory structure.
    
    Args:
        data_root: Root data directory
        
    Returns:
        Number of directories created
    """
    directories = [
        # Structure directories
        "structure",
        "structure/mmcif",
        "structure/pdb",
        "structure/alignments",
        "structure/structure_dataset",
        "structure/structure_dataset/standard",
        
        # Sequence directories
        "sequence",
        "sequence/fasta",
        "sequence/alignments",
        "sequence/metadata",
        
        # GRN directories
        "grn",
        "grn/ref",
        "grn/tables",
        "grn/configs",
        "grn/assignments",
        "grn/grn",  # For dataset management
        
        # Embedding directories
        "embedding",
        "embedding/datasets",
        "embedding/models",
        
        # Property directories
        "property",
        "property/annotations",
        "property/calculated",
        
        # Ligand directories
        "ligand",
        "ligand/structures",
        "ligand/binding_sites",
        
        # Graph directories
        "graph",
        
        # Other processor directories
        "literature",
        "seq",
        "seq/datasets",
    ]
    
    created = 0
    for dir_path in directories:
        full_path = data_root / dir_path
        if not full_path.exists():
            full_path.mkdir(parents=True, exist_ok=True)
            created += 1
    
    return created


def create_registries(data_root: Path) -> int:
    """
    Create empty registry files.
    
    Args:
        data_root: Root data directory
        
    Returns:
        Number of registries created
    """
    registries = [
        # Global registries
        ("entity_registry.json", {"entities": {}, "version": "1.0"}),
        ("global_registry.json", {"entities": {}, "datasets": {}, "version": "1.0"}),
        
        # Processor registries
        ("structure/registry.json", {"datasets": {}, "version": "1.0"}),
        ("sequence/registry.json", {"datasets": {}, "version": "1.0"}),
        ("grn/registry.json", {"datasets": {}, "version": "1.0"}),
        ("embedding/registry.json", {"datasets": {}, "version": "1.0"}),
        ("property/registry.json", {"datasets": {}, "version": "1.0"}),
        ("ligand/registry.json", {"datasets": {}, "version": "1.0"}),
    ]
    
    created = 0
    for registry_path, content in registries:
        full_path = data_root / registry_path
        full_path.parent.mkdir(parents=True, exist_ok=True)
        
        if not full_path.exists() or full_path.stat().st_size == 0:
            with open(full_path, 'w') as f:
                json.dump(content, f, indent=2)
            created += 1
    
    return created


def copy_reference_data(ref_data_dir: Path, target_data_dir: Path) -> int:
    """
    Copy reference data files to target directory.
    
    Args:
        ref_data_dir: Source reference data directory
        target_data_dir: Target data directory
        
    Returns:
        Number of files copied
    """
    if not ref_data_dir.exists():
        logger.warning(f"Reference data directory not found: {ref_data_dir}")
        return 0
    
    # Define all reference files to copy
    files_to_copy = []
    
    # Walk through reference data and collect all files
    for root, dirs, files in os.walk(ref_data_dir):
        for file in files:
            src_path = Path(root) / file
            # Calculate relative path from ref_data_dir
            rel_path = src_path.relative_to(ref_data_dir)
            files_to_copy.append((str(rel_path), str(rel_path)))
    
    copied_count = 0
    for src_rel, dst_rel in files_to_copy:
        src_path = ref_data_dir / src_rel
        dst_path = target_data_dir / dst_rel
        
        if src_path.exists():
            # Create parent directories
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            
            try:
                shutil.copy2(src_path, dst_path)
                copied_count += 1
            except Exception as e:
                logger.warning(f"Failed to copy {src_rel}: {e}")
    
    return copied_count


def init_data_directory(
    data_root: Optional[str] = None,
    skip_reference: bool = False,
    force: bool = False
) -> Dict[str, Any]:
    """
    Initialize Protos data directory with structure and reference data.
    
    Args:
        data_root: Root data directory (uses ProtosPaths default if None)
        skip_reference: If True, don't copy reference data
        force: If True, overwrite existing files
        
    Returns:
        Dictionary with initialization statistics
    """
    # Get data root from ProtosPaths
    if data_root is None:
        paths = ProtosPaths()
        data_root = paths.data_root
        logger.info(f"Using ProtosPaths data root: {data_root}")
    
    data_dir = Path(data_root)
    
    stats = {
        "directories_created": 0,
        "registries_created": 0,
        "reference_files_copied": 0,
        "errors": []
    }
    
    click.echo(f"Initializing Protos data directory at: {data_dir}")
    
    # Check if directory exists and has content
    if data_dir.exists() and any(data_dir.iterdir()) and not force:
        click.echo(f"\n⚠️  Directory already exists and contains files: {data_dir}")
        if not click.confirm("Do you want to continue?"):
            click.echo("Aborted.")
            return stats
    
    # 1. Create directory structure
    click.echo("\n1. Creating directory structure...")
    data_dir.mkdir(parents=True, exist_ok=True)
    created = create_directory_structure(data_dir)
    stats["directories_created"] = created
    click.echo(f"   ✓ Created {created} directories")
    
    # 2. Create registries
    click.echo("\n2. Creating registry files...")
    created = create_registries(data_dir)
    stats["registries_created"] = created
    click.echo(f"   ✓ Created {created} registries")
    
    # 3. Copy reference data
    if not skip_reference:
        click.echo("\n3. Copying reference data...")
        
        # Find reference data in the source tree
        # Reference data is part of the source code, not in the data directory
        import protos
        protos_module_path = Path(protos.__file__).parent
        
        ref_paths = [
            protos_module_path / "reference_data",  # Standard location in package
            Path(__file__).parent.parent.parent / "reference_data",  # From CLI location
            Path.cwd() / "src" / "protos" / "reference_data",  # From project root
        ]
        
        ref_data_dir = None
        for path in ref_paths:
            if path.exists():
                ref_data_dir = path
                break
        
        if ref_data_dir:
            copied = copy_reference_data(ref_data_dir, data_dir)
            stats["reference_files_copied"] = copied
            click.echo(f"   ✓ Copied {copied} reference files")
        else:
            msg = "Reference data directory not found"
            stats["errors"].append(msg)
            click.echo(f"   ⚠️  {msg}")
    
    # 4. Summary
    click.echo("\n=== Initialization Summary ===")
    click.echo(f"Directories created: {stats['directories_created']}")
    click.echo(f"Registries created: {stats['registries_created']}")
    if not skip_reference:
        click.echo(f"Reference files copied: {stats['reference_files_copied']}")
    
    if stats["errors"]:
        click.echo("\nWarnings:")
        for error in stats["errors"]:
            click.echo(f"  - {error}")
    
    return stats


@click.command()
@click.option(
    '--data-root',
    type=click.Path(file_okay=False, dir_okay=True),
    help='Root data directory (uses default if not specified)'
)
@click.option(
    '--skip-reference',
    is_flag=True,
    help='Skip copying reference data'
)
@click.option(
    '--force',
    '-f',
    is_flag=True,
    help='Overwrite existing files'
)
def init_data(data_root, skip_reference, force):
    """
    Initialize Protos data directory structure and reference data.
    
    This command will:
    - Create all necessary directories
    - Initialize empty registry files
    - Copy reference data from src/protos/reference_data
    
    Use this command when:
    - Setting up Protos for the first time
    - Creating a new data directory
    - Restoring directory structure
    """
    try:
        stats = init_data_directory(
            data_root=data_root,
            skip_reference=skip_reference,
            force=force
        )
        
        click.echo("\n✅ Initialization completed successfully!")
        click.echo("\nYour Protos data directory is ready to use.")
        
    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        raise click.Abort()


if __name__ == "__main__":
    init_data()