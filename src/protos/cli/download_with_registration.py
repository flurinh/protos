"""
Enhanced download functions that automatically register downloaded files.

These functions wrap the existing download functionality to ensure
all downloaded files are immediately registered in the entity system.
"""

import os
from pathlib import Path
from typing import List, Optional, Dict
import logging

from protos.processing.structure import StructureProcessor
from protos.commands.register_data import register_structure_file, register_sequence_file
from protos.io.data_access import GlobalRegistry

logger = logging.getLogger(__name__)


def download_and_register_structure(
    pdb_id: str, 
    save_dir: Optional[Path] = None,
    skip_if_registered: bool = True
) -> Optional[str]:
    """
    Download a structure file and automatically register it.
    
    Args:
        pdb_id: PDB ID to download
        save_dir: Directory to save file (uses default if None)
        skip_if_registered: Skip download if already registered
        
    Returns:
        Entity ID if successful, None if failed
    """
    # Check if already registered
    if skip_if_registered:
        registry = GlobalRegistry()
        existing_entity = registry.entity_registry.resolve_identifier(pdb_id)
        if existing_entity:
            logger.info(f"Structure {pdb_id} already registered as {existing_entity}")
            return existing_entity
    
    # Download the file
    if save_dir:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)
        success = StructureProcessor.download_cif(pdb_id, str(save_dir))
    else:
        success = StructureProcessor.download_cif(pdb_id, None)
    
    if not success:
        logger.error(f"Failed to download structure {pdb_id}")
        return None
    
    # Determine file path using the same method as download_cif
    if save_dir is not None:
        file_path = Path(save_dir) / f"{pdb_id}.cif"
    else:
        from protos.io.paths import get_structure_path
        file_path = Path(get_structure_path(pdb_id, create_if_missing=False))
    
    # Register the downloaded file
    try:
        entity_id = register_structure_file(file_path, pdb_id)
        logger.info(f"Downloaded and registered {pdb_id} -> {entity_id}")
        return entity_id
    except Exception as e:
        logger.error(f"Failed to register downloaded file {pdb_id}: {e}")
        return None


def bulk_download_structures(
    pdb_ids: List[str],
    save_dir: Optional[Path] = None,
    skip_existing: bool = True,
    max_workers: int = 4
) -> Dict[str, str]:
    """
    Download multiple structures and register them.
    
    Args:
        pdb_ids: List of PDB IDs to download
        save_dir: Directory to save files
        skip_existing: Skip already registered structures
        max_workers: Number of parallel downloads
        
    Returns:
        Dict mapping PDB IDs to entity IDs
    """
    import concurrent.futures
    
    results = {}
    
    # Use thread pool for parallel downloads
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit download tasks
        future_to_pdb = {
            executor.submit(
                download_and_register_structure, 
                pdb_id, 
                save_dir,
                skip_existing
            ): pdb_id 
            for pdb_id in pdb_ids
        }
        
        # Collect results
        for future in concurrent.futures.as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                entity_id = future.result()
                if entity_id:
                    results[pdb_id] = entity_id
                else:
                    results[pdb_id] = None
            except Exception as e:
                logger.error(f"Error downloading {pdb_id}: {e}")
                results[pdb_id] = None
    
    # Summary
    successful = sum(1 for v in results.values() if v is not None)
    logger.info(f"Downloaded and registered {successful}/{len(pdb_ids)} structures")
    
    return results


def download_sequence_database(
    database: str,
    sequence_ids: List[str],
    save_dir: Optional[Path] = None
) -> Dict[str, str]:
    """
    Download sequences from a database and register them.
    
    Args:
        database: Database name (e.g., "uniprot", "ncbi")
        sequence_ids: List of sequence IDs
        save_dir: Directory to save files
        
    Returns:
        Dict mapping sequence IDs to entity IDs
    """
    # This is a placeholder - actual implementation would depend on the database
    # For now, let's show the pattern
    
    if database.lower() == "uniprot":
        from protos.loaders.uniprot_utils import download_uniprot_fasta
        
        results = {}
        for seq_id in sequence_ids:
            try:
                # Download sequence
                file_path = download_uniprot_fasta(seq_id, save_dir)
                
                # Register it
                registered = register_sequence_file(Path(file_path))
                results.update(registered)
                
            except Exception as e:
                logger.error(f"Failed to download {seq_id}: {e}")
                results[seq_id] = None
        
        return results
    else:
        raise NotImplementedError(f"Database {database} not yet supported")


# CLI command implementations

def download_structure_command(args):
    """CLI command to download and register a structure."""
    entity_id = download_and_register_structure(
        args.pdb_id,
        Path(args.dir) if args.dir else None,
        skip_if_registered=not args.force
    )
    
    if entity_id:
        print(f"Successfully downloaded and registered {args.pdb_id} -> {entity_id}")
    else:
        print(f"Failed to download {args.pdb_id}")
        return 1
    return 0


def bulk_download_command(args):
    """CLI command to bulk download structures."""
    # Read PDB IDs from file or args
    if args.file:
        with open(args.file, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
    else:
        pdb_ids = args.pdb_ids
    
    results = bulk_download_structures(
        pdb_ids,
        Path(args.dir) if args.dir else None,
        skip_existing=not args.force,
        max_workers=args.workers
    )
    
    # Report results
    successful = sum(1 for v in results.values() if v is not None)
    print(f"\nDownloaded {successful}/{len(pdb_ids)} structures")
    
    if args.verbose:
        print("\nDetails:")
        for pdb_id, entity_id in results.items():
            if entity_id:
                print(f"  {pdb_id} -> {entity_id}")
            else:
                print(f"  {pdb_id} -> FAILED")