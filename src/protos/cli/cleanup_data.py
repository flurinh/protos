#!/usr/bin/env python3
"""
Cleanup command for Protos data.

This command:
- Deletes all data in the Protos data directory
- Reinitializes the directory structure
- Copies reference data (same as init)

Essentially: cleanup = delete everything + init
"""

import os
import shutil
import click
from pathlib import Path
from typing import Optional, Dict, Any
import logging

from protos.io.paths.path_config import ProtosPaths
from protos.cli.init_data import init_data_directory


logger = logging.getLogger(__name__)


def clean_all_data(data_root: Path) -> int:
    """
    Remove all data files and directories.
    
    Args:
        data_root: Root data directory
        
    Returns:
        Number of items removed
    """
    if not data_root.exists():
        return 0
    
    removed_count = 0
    
    # Remove everything except the root directory itself
    for item in data_root.iterdir():
        try:
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
            removed_count += 1
        except Exception as e:
            logger.warning(f"Failed to remove {item}: {e}")
    
    return removed_count


def cleanup_and_reinit(
    data_root: Optional[str] = None,
    skip_reference: bool = False,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Clean up all Protos data and reinitialize.
    
    This is essentially: delete everything + init
    
    Args:
        data_root: Root data directory (uses ProtosPaths default if None)
        skip_reference: If True, don't copy reference data
        dry_run: If True, only show what would be done
        
    Returns:
        Dictionary with cleanup and init statistics
    """
    # Get data root
    if data_root is None:
        paths = ProtosPaths()
        data_root = paths.data_root
    
    data_dir = Path(data_root)
    
    stats = {
        "items_removed": 0,
        "init_stats": {},
        "errors": []
    }
    
    if dry_run:
        click.echo("DRY RUN - No changes will be made\n")
    
    # 1. Clean everything
    click.echo(f"Cleaning Protos data at: {data_dir}")
    click.echo("\n1. Removing all data...")
    
    if data_dir.exists():
        if dry_run:
            # Count what would be removed
            item_count = sum(1 for _ in data_dir.iterdir())
            file_count = sum(1 for _ in data_dir.rglob('*') if _.is_file())
            click.echo(f"   Would remove {item_count} top-level items ({file_count} total files)")
        else:
            removed = clean_all_data(data_dir)
            stats["items_removed"] = removed
            click.echo(f"   ✓ Removed {removed} items")
    else:
        click.echo("   Directory doesn't exist, nothing to clean")
    
    # 2. Reinitialize
    click.echo("\n2. Reinitializing data directory...")
    
    if not dry_run:
        # Actually run init
        init_stats = init_data_directory(
            data_root=str(data_dir),
            skip_reference=skip_reference,
            force=True  # Force since we just cleaned
        )
        stats["init_stats"] = init_stats
    else:
        click.echo("   Would run initialization:")
        click.echo("   - Create directory structure")
        click.echo("   - Initialize registries")
        if not skip_reference:
            click.echo("   - Copy reference data")
    
    # 3. Summary
    if not dry_run:
        click.echo("\n=== Cleanup Summary ===")
        click.echo(f"Items removed: {stats['items_removed']}")
        if stats["init_stats"]:
            click.echo(f"Directories created: {stats['init_stats'].get('directories_created', 0)}")
            click.echo(f"Registries created: {stats['init_stats'].get('registries_created', 0)}")
            if not skip_reference:
                click.echo(f"Reference files copied: {stats['init_stats'].get('reference_files_copied', 0)}")
    
    return stats


@click.command()
@click.option(
    '--data-root',
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help='Root data directory (uses default if not specified)'
)
@click.option(
    '--skip-reference',
    is_flag=True,
    help='Skip copying reference data after cleanup'
)
@click.option(
    '--dry-run',
    is_flag=True,
    help='Show what would be done without making changes'
)
@click.option(
    '--force',
    '-f',
    is_flag=True,
    help='Skip all confirmation prompts'
)
def cleanup_data(data_root, skip_reference, dry_run, force):
    """
    Clean up all Protos data and reinitialize.
    
    This command will:
    1. DELETE all data in the Protos data directory
    2. Recreate the directory structure
    3. Initialize empty registries
    4. Copy reference data (unless --skip-reference)
    
    This is equivalent to: delete everything + init
    
    Use --dry-run to see what would be done without making changes.
    """
    try:
        # Get data root
        if data_root is None:
            paths = ProtosPaths()
            data_root = paths.data_root
        
        data_dir = Path(data_root)
        
        # Show warning and get confirmation
        if not dry_run and not force:
            click.echo("\n" + "="*60)
            click.echo("⚠️  WARNING: DESTRUCTIVE OPERATION")
            click.echo("="*60)
            click.echo(f"\nThis will DELETE ALL DATA in: {data_dir}")
            click.echo("\nThis includes:")
            click.echo("  • All structure files (PDB/CIF)")
            click.echo("  • All sequence files (FASTA)")
            click.echo("  • All GRN tables")
            click.echo("  • All embeddings")
            click.echo("  • All registries and metadata")
            click.echo("\nThe directory will then be reinitialized with:")
            click.echo("  • Empty directory structure")
            click.echo("  • Empty registries")
            if not skip_reference:
                click.echo("  • Reference data from src/protos/reference_data")
            click.echo("\nThis action CANNOT be undone!")
            click.echo("="*60)
            
            # First confirmation
            if not click.confirm("\nAre you sure you want to continue?", default=False):
                click.echo("Aborted.")
                return
            
            # Second confirmation for safety
            click.echo("\n⚠️  FINAL WARNING: All data will be permanently deleted!")
            confirm_text = click.prompt(
                '\nType "DELETE ALL DATA" to confirm', 
                type=str
            )
            
            if confirm_text != "DELETE ALL DATA":
                click.echo("\nConfirmation text did not match. Aborted.")
                return
        
        # Run cleanup and reinit
        stats = cleanup_and_reinit(
            data_root=data_root,
            skip_reference=skip_reference,
            dry_run=dry_run
        )
        
        if not dry_run:
            click.echo("\n✅ Cleanup and reinitialization completed successfully!")
            click.echo("\nYour Protos data directory has been reset to a clean state.")
        
    except KeyboardInterrupt:
        click.echo("\n\nOperation cancelled by user.")
        raise click.Abort()
    except Exception as e:
        click.echo(f"\n❌ Error: {e}", err=True)
        raise click.Abort()


# Make it directly callable as a function as well
def reset_protos_data(
    data_root: Optional[str] = None,
    force: bool = False,
    skip_reference: bool = False
) -> Dict[str, Any]:
    """
    Python function interface to reset Protos data.
    
    This runs: cleanup (delete all) + init
    
    Args:
        data_root: Root data directory (uses default if None)
        force: If True, skip confirmation prompts
        skip_reference: If True, don't copy reference data
        
    Returns:
        Statistics dictionary
    """
    if not force:
        response = input("\n⚠️  This will DELETE ALL PROTOS DATA and reinitialize. Are you sure? (yes/no): ")
        if response.lower() != 'yes':
            print("Aborted.")
            return {}
    
    return cleanup_and_reinit(
        data_root=data_root,
        skip_reference=skip_reference,
        dry_run=False
    )


if __name__ == "__main__":
    cleanup_data()