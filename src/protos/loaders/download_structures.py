from Bio.PDB import PDBList
import os
import shutil
import logging
from pathlib import Path

# Import ProtosPaths for centralized path management
try:
    from protos.io.paths.path_config import ProtosPaths
    _HAS_PROTOS_PATHS = True
except ImportError:
    _HAS_PROTOS_PATHS = False

logger = logging.getLogger(__name__)

def download_protein_structures(pdb_ids, target_folder=None, overwrite=False, processor=None):
    """
    Download protein structures using ProtosPaths for centralized path management.
    
    This function follows the core philosophy: "Users work with names, never paths. 
    Protos handles ALL file system complexity."
    
    Args:
        pdb_ids: List of PDB IDs to download
        target_folder: Legacy parameter for backward compatibility (optional)
        overwrite: Whether to overwrite existing files
        processor: StructureProcessor instance to use for path resolution (preferred)
        
    Returns:
        Tuple of (successful_downloads, failed_downloads)
    """
    pdbl = PDBList()
    
    # Determine target path using ProtosPaths (preferred) or legacy fallback
    if processor is not None:
        # Use processor's structure directory (follows core philosophy)
        target_path = processor.path_structure_dir
    elif _HAS_PROTOS_PATHS and target_folder is None:
        # Use ProtosPaths directly if no processor provided
        paths = ProtosPaths()
        target_path = Path(paths.get_subdir_path('structure', 'structure_dir'))
    else:
        # Legacy fallback for backward compatibility
        target_path = Path(target_folder or 'data/mmcif/')
        if not target_path.exists():
            target_path.mkdir(parents=True, exist_ok=True)

    successful = []
    failed = []
    
    for pdb_id in pdb_ids:
        try:
            # Ensure lowercase for consistency
            pdb_id = pdb_id.lower()
            
            # Check if file already exists
            expected_path = target_path / f"{pdb_id}.cif"
            if expected_path.exists() and not overwrite:
                logger.info(f"File already exists: {expected_path}")
                successful.append(pdb_id)
                continue
                
            # PDBList downloads to a subdirectory structure
            # The file will be saved as something like: target_folder/xy/1xyz.cif
            # where xy are the middle two characters of the PDB ID
            file_path = pdbl.retrieve_pdb_file(
                pdb_id, 
                pdir=str(target_folder), 
                file_format="mmCif",
                overwrite=overwrite
            )
            
            if file_path and os.path.exists(file_path):
                # Move the file to the expected location if needed
                if file_path != str(expected_path):
                    # Need to move the file
                    shutil.move(file_path, expected_path)
                    
                    # Clean up the subdirectory if empty
                    parent_dir = Path(file_path).parent
                    if parent_dir != target_path:
                        try:
                            parent_dir.rmdir()
                        except OSError:
                            # Directory not empty, that's fine
                            pass
                
                successful.append(pdb_id)
            else:
                logger.warning(f"Download seemed to succeed but file not found for {pdb_id}")
                failed.append(pdb_id)
                
        except Exception as e:
            logger.warning(f"Failed to download {pdb_id}: {e}")
            failed.append(pdb_id)
    
    if failed:
        logger.info(f"Successfully downloaded {len(successful)} structures, failed {len(failed)}: {failed}")
    else:
        logger.info(f"Successfully downloaded all {len(successful)} structures")
    
    return successful, failed


def download_structures_with_processor(pdb_ids, processor, overwrite=False):
    """
    Download protein structures using a processor for ProtosPaths integration.
    
    This is the preferred method that fully follows the core philosophy.
    
    Args:
        pdb_ids: List of PDB IDs to download
        processor: StructureProcessor instance (handles all path management)
        overwrite: Whether to overwrite existing files
        
    Returns:
        Tuple of (successful_downloads, failed_downloads)
    """
    return download_protein_structures(
        pdb_ids=pdb_ids,
        processor=processor,
        overwrite=overwrite
    )
