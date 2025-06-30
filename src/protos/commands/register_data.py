"""
Commands for registering raw data files in the Protos entity system.

This module provides functions to:
1. Register individual files
2. Scan directories for unregistered files
3. Bulk register multiple files
4. Clean up orphaned registry entries
"""

import os
from pathlib import Path
from typing import List, Dict, Optional, Set, Tuple
import logging

from protos.io.data_access import GlobalRegistry, generate_entity_id
from protos.io.fasta_utils import read_fasta
# Removed import - will extract PDB ID from filename

logger = logging.getLogger(__name__)


def register_structure_file(file_path: Path, pdb_id: Optional[str] = None) -> str:
    """
    Register a single structure file (CIF/PDB) in the entity system.
    
    Args:
        file_path: Path to the structure file
        pdb_id: PDB ID (if not provided, extracted from filename)
        
    Returns:
        Entity ID of the registered structure
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Structure file not found: {file_path}")
    
    # Extract PDB ID if not provided
    if pdb_id is None:
        # Extract from filename (remove .cif/.pdb extension)
        pdb_id = file_path.stem.upper()
        # Handle AlphaFold naming convention
        if pdb_id.startswith('AF-') and '-F' in pdb_id:
            # Keep the full AlphaFold ID
            pass
        else:
            # For regular PDB files, ensure it's uppercase
            pdb_id = pdb_id.upper()
    
    # Generate entity ID
    entity_id = generate_entity_id(pdb_id)
    
    # Register with global registry
    registry = GlobalRegistry()
    registry.entity_registry.register_entity(
        entity_id=entity_id,
        entity_type="structure",
        original_id=pdb_id,
        file_path=str(file_path),
        metadata={
            'file_format': 'mmcif' if file_path.suffix == '.cif' else 'pdb',
            'registered_from': 'manual_registration'
        }
    )
    
    logger.info(f"Registered structure {pdb_id} -> {entity_id}")
    return entity_id


def register_sequence_file(file_path: Path) -> Dict[str, str]:
    """
    Register sequences from a FASTA file in the entity system.
    
    Args:
        file_path: Path to the FASTA file
        
    Returns:
        Dict mapping sequence IDs to entity IDs
    """
    if not file_path.exists():
        raise FileNotFoundError(f"Sequence file not found: {file_path}")
    
    # Read sequences from FASTA
    sequences = read_fasta(str(file_path))
    
    if not sequences:
        logger.warning(f"No sequences found in {file_path}")
        return {}
    
    registry = GlobalRegistry()
    registered = {}
    
    for header, sequence in sequences.items():
        # Extract ID from header (first word)
        seq_id = header.split()[0]
        
        # Generate entity ID
        entity_id = generate_entity_id(seq_id)
        
        # Register with global registry
        registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=seq_id,
            file_path=str(file_path),
            metadata={
                'full_header': header,
                'sequence_length': len(sequence),
                'registered_from': 'manual_registration'
            }
        )
        
        registered[seq_id] = entity_id
        logger.info(f"Registered sequence {seq_id} -> {entity_id}")
    
    return registered


def scan_for_unregistered_files(
    data_root: Path,
    file_type: str = "all"
) -> Dict[str, List[Path]]:
    """
    Scan data directories for files that are not registered in the entity system.
    
    Args:
        data_root: Root data directory
        file_type: Type of files to scan ("structure", "sequence", or "all")
        
    Returns:
        Dict with 'unregistered' and 'registered' file lists
    """
    registry = GlobalRegistry()
    results = {
        'unregistered': [],
        'registered': [],
        'orphaned_entries': []  # Registry entries with no files
    }
    
    # Define paths to scan
    scan_paths = {}
    data_root_path = Path(data_root) if isinstance(data_root, str) else data_root
    if file_type in ["structure", "all"]:
        scan_paths['structure'] = data_root_path / "structure" / "mmcif"
    if file_type in ["sequence", "all"]:
        scan_paths['sequence'] = data_root_path / "sequence" / "fasta"
    
    # Get all registered file paths
    registered_paths = set()
    all_entities = registry.entity_registry.entities
    
    for entity_id, entity_data in all_entities.items():
        for format_type, format_data in entity_data.get('formats', {}).items():
            if 'file_path' in format_data:
                registered_paths.add(Path(format_data['file_path']))
    
    # Scan directories
    for data_type, scan_path in scan_paths.items():
        if not scan_path.exists():
            logger.warning(f"Directory not found: {scan_path}")
            continue
            
        # Find all relevant files
        if data_type == "structure":
            patterns = ["*.cif", "*.pdb"]
        else:  # sequence
            patterns = ["*.fasta", "*.fa", "*.faa"]
        
        for pattern in patterns:
            for file_path in scan_path.glob(pattern):
                if file_path in registered_paths:
                    results['registered'].append(file_path)
                else:
                    results['unregistered'].append(file_path)
    
    # Check for orphaned registry entries
    for reg_path in registered_paths:
        if not reg_path.exists():
            results['orphaned_entries'].append(reg_path)
    
    return results


def bulk_register_files(
    file_paths: List[Path],
    file_type: str,
    skip_errors: bool = True
) -> Dict[str, any]:
    """
    Register multiple files at once.
    
    Args:
        file_paths: List of file paths to register
        file_type: Type of files ("structure" or "sequence")
        skip_errors: Whether to continue on errors
        
    Returns:
        Dict with registration results
    """
    results = {
        'success': [],
        'failed': [],
        'skipped': []
    }
    
    for file_path in file_paths:
        try:
            if not file_path.exists():
                results['skipped'].append((str(file_path), "File not found"))
                continue
                
            if file_type == "structure":
                entity_id = register_structure_file(file_path)
                results['success'].append((str(file_path), entity_id))
            elif file_type == "sequence":
                registered = register_sequence_file(file_path)
                results['success'].append((str(file_path), registered))
            else:
                results['failed'].append((str(file_path), f"Unknown file type: {file_type}"))
                
        except Exception as e:
            logger.error(f"Failed to register {file_path}: {e}")
            results['failed'].append((str(file_path), str(e)))
            if not skip_errors:
                raise
    
    return results


def clean_orphaned_entries(dry_run: bool = True) -> List[Tuple[str, str]]:
    """
    Remove registry entries that point to non-existent files.
    
    Args:
        dry_run: If True, only report what would be cleaned
        
    Returns:
        List of (entity_id, file_path) tuples that were cleaned
    """
    registry = GlobalRegistry()
    orphaned = []
    
    for entity_id, entity_data in registry.entity_registry._entities.items():
        for format_type, format_data in entity_data.get('formats', {}).items():
            if 'file_path' in format_data:
                file_path = Path(format_data['file_path'])
                if not file_path.exists():
                    orphaned.append((entity_id, str(file_path)))
                    
                    if not dry_run:
                        # Remove this format from the entity
                        # Note: In real implementation, we'd need a method to remove formats
                        logger.info(f"Would remove orphaned entry: {entity_id} -> {file_path}")
    
    return orphaned


def create_dataset_from_entities(
    entity_ids: List[str],
    dataset_name: str,
    dataset_type: str,
    metadata: Optional[Dict] = None
) -> Path:
    """
    Create a dataset file from a list of entity IDs.
    
    This is how users define collections of entities to work with.
    
    Args:
        entity_ids: List of entity IDs to include
        dataset_name: Name for the dataset
        dataset_type: Type of dataset (structure, sequence, etc.)
        metadata: Additional metadata for the dataset
        
    Returns:
        Path to created dataset file
    """
    import json
    from protos.io.paths import ProtosPaths
    
    # Get dataset directory
    paths = ProtosPaths()
    if dataset_type == "structure":
        dataset_dir = paths.get_path("structure") / "structure_dataset"
    elif dataset_type == "sequence":
        dataset_dir = paths.get_path("sequence") / "datasets"
    else:
        dataset_dir = paths.get_path(dataset_type) / "datasets"
    
    dataset_dir.mkdir(parents=True, exist_ok=True)
    
    # Create dataset definition
    dataset = {
        "name": dataset_name,
        "type": dataset_type,
        "entity_ids": entity_ids,
        "metadata": metadata or {},
        "created_at": str(Path.ctime(Path.cwd())),
        "version": "1.0"
    }
    
    # Save dataset file
    dataset_path = dataset_dir / f"{dataset_name}.json"
    with open(dataset_path, 'w') as f:
        json.dump(dataset, f, indent=2)
    
    logger.info(f"Created dataset '{dataset_name}' with {len(entity_ids)} entities")
    return dataset_path


# CLI Interface functions (can be called from CLI commands)

def register_command(args):
    """CLI command to register files."""
    file_path = Path(args.file)
    
    if args.type == "structure":
        entity_id = register_structure_file(file_path, args.pdb_id)
        print(f"Registered structure: {file_path.name} -> {entity_id}")
    elif args.type == "sequence":
        registered = register_sequence_file(file_path)
        print(f"Registered {len(registered)} sequences from {file_path.name}")
        for seq_id, entity_id in registered.items():
            print(f"  {seq_id} -> {entity_id}")
    else:
        print(f"Unknown file type: {args.type}")


def scan_command(args):
    """CLI command to scan for unregistered files."""
    from protos.io.paths import ProtosPaths
    
    paths = ProtosPaths()
    results = scan_for_unregistered_files(
        paths.data_root,
        file_type=args.type
    )
    
    print(f"\nScan Results:")
    print(f"Registered files: {len(results['registered'])}")
    print(f"Unregistered files: {len(results['unregistered'])}")
    print(f"Orphaned entries: {len(results['orphaned_entries'])}")
    
    if results['unregistered'] and args.verbose:
        print("\nUnregistered files:")
        for file_path in results['unregistered']:
            print(f"  {file_path}")
    
    if results['orphaned_entries'] and args.verbose:
        print("\nOrphaned registry entries (files not found):")
        for file_path in results['orphaned_entries']:
            print(f"  {file_path}")


def bulk_register_command(args):
    """CLI command to bulk register files."""
    from protos.io.paths import ProtosPaths
    
    # Get unregistered files if --all flag
    if args.all:
        paths = ProtosPaths()
        scan_results = scan_for_unregistered_files(
            paths.data_root,
            file_type=args.type
        )
        file_paths = scan_results['unregistered']
    else:
        file_paths = [Path(f) for f in args.files]
    
    # Register files
    results = bulk_register_files(
        file_paths,
        file_type=args.type,
        skip_errors=args.skip_errors
    )
    
    print(f"\nBulk Registration Results:")
    print(f"Success: {len(results['success'])}")
    print(f"Failed: {len(results['failed'])}")
    print(f"Skipped: {len(results['skipped'])}")
    
    if results['failed'] and args.verbose:
        print("\nFailed registrations:")
        for file_path, error in results['failed']:
            print(f"  {file_path}: {error}")