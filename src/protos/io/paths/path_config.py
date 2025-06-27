"""
Simplified path configuration for the Protos framework.

This module provides a simplified path management system using a single
data directory for all Protos data, eliminating the complexity of
separate reference and user data paths.
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple
from enum import Enum

from .path_constants import (
    ENV_DATA_ROOT,
    ENV_REF_DATA_ROOT,
    DEFAULT_PROCESSOR_DIRS,
    DEFAULT_STRUCTURE_SUBDIRS,
    DEFAULT_GRN_SUBDIRS,
    DEFAULT_SEQUENCE_SUBDIRS,
    DEFAULT_TEST_SUBDIRS,
    DEFAULT_REGISTRY_FILENAME,
    DEFAULT_GLOBAL_REGISTRY_FILENAME,
    join_path
)

# Configure logger
logger = logging.getLogger(__name__)


class DataSource(Enum):
    """Legacy enum for backward compatibility."""
    REFERENCE = "reference"
    USER = "user"
    AUTO = "auto"


def get_default_data_root() -> str:
    """
    Get the default data root directory.
    
    Order of precedence:
    1. PROTOS_DATA_ROOT environment variable
    2. ~/protos_data/ (home directory)
    
    Returns:
        Default data root directory path
    """
    # Check environment variable first
    env_root = os.environ.get(ENV_DATA_ROOT)
    if env_root:
        return os.path.expanduser(env_root)
    
    # Default to home directory
    return os.path.expanduser("~/protos_data")


class ProtosPaths:
    """
    Simplified path management for Protos using a single data directory.
    
    This class provides a standardized way to manage paths for different
    data types and processors in the Protos framework using a unified
    data directory approach.
    """
    
    # Class-level configuration
    _global_data_root: Optional[str] = None
    
    @classmethod
    def set_data_root(cls, data_root: Optional[str]) -> None:
        """
        Set the global data root for all ProtosPaths instances.
        
        This allows configuring the data directory once for the entire session,
        particularly useful for testing where all paths should point to test-data.
        
        Args:
            data_root: Path to set as global data root, or None to clear
        """
        cls._global_data_root = data_root
        if data_root:
            logger.info(f"Global data root set to: {data_root}")
        else:
            logger.info("Global data root cleared")
        
        # Reset the default resolver so it picks up the new global setting
        global _DEFAULT_PATH_RESOLVER
        _DEFAULT_PATH_RESOLVER = None
    
    @classmethod
    def get_global_data_root(cls) -> Optional[str]:
        """
        Get the globally configured data root.
        
        Returns:
            The global data root if set, None otherwise
        """
        return cls._global_data_root
    
    def __init__(self, 
                 data_root: Optional[str] = None,
                 user_data_root: Optional[str] = None,  # For backward compatibility
                 ref_data_root: Optional[str] = None,   # For backward compatibility
                 create_dirs: bool = True,
                 validate: bool = True):
        """
        Initialize the path manager with a single data directory.
        
        Args:
            data_root: Root directory for all data (default: ~/protos_data/)
            user_data_root: Deprecated - if provided, used as data_root
            ref_data_root: Deprecated - ignored with warning
            create_dirs: Whether to create directories that don't exist
            validate: Whether to validate path structure
        """
        # Handle backward compatibility
        if user_data_root or ref_data_root:
            logger.warning(
                "user_data_root and ref_data_root are deprecated. "
                "Using a single data_root instead."
            )
            # Use user_data_root if provided, otherwise data_root
            data_root = data_root or user_data_root
        
        # Set data root with priority: parameter > global setting > env var > default
        if data_root:
            self.data_root = data_root
        elif self._global_data_root:
            self.data_root = self._global_data_root
        else:
            self.data_root = get_default_data_root()
        
        self.data_root = os.path.expanduser(self.data_root)
        
        # Make path absolute if it's not already
        if not os.path.isabs(self.data_root):
            self.data_root = os.path.abspath(self.data_root)
        
        # For backward compatibility, set both old attributes to the same value
        self.user_data_root = self.data_root
        self.ref_data_root = self.data_root
        
        # Create data root directory if requested
        if create_dirs:
            os.makedirs(self.data_root, exist_ok=True)
            logger.info(f"Initialized Protos data directory at: {self.data_root}")
        
        # Initialize directory maps
        self.processor_dirs = DEFAULT_PROCESSOR_DIRS.copy()
        self.structure_dirs = DEFAULT_STRUCTURE_SUBDIRS.copy()
        self.grn_dirs = DEFAULT_GRN_SUBDIRS.copy()
        self.sequence_dirs = DEFAULT_SEQUENCE_SUBDIRS.copy()
        self.test_dirs = DEFAULT_TEST_SUBDIRS.copy()
        
        # Create standard directories if requested
        if create_dirs:
            self._create_standard_dirs(self.data_root)
            
        # Validate directory structure if requested
        if validate:
            self._validate_directory_structure()
    
    def _create_standard_dirs(self, root: str):
        """
        Create the standard directory structure in the specified root.
        
        Args:
            root: Root directory where to create the structure
        """
        # Create processor directories
        for processor_type, dir_name in self.processor_dirs.items():
            processor_path = join_path(root, dir_name)
            os.makedirs(processor_path, exist_ok=True)
            
            # Create subdirectories for each processor type
            if processor_type == 'structure':
                for subdir in self.structure_dirs.values():
                    os.makedirs(join_path(processor_path, subdir), exist_ok=True)
            elif processor_type == 'grn':
                for subdir in self.grn_dirs.values():
                    os.makedirs(join_path(processor_path, subdir), exist_ok=True)
                # Create ref subdirectory for reference tables
                os.makedirs(join_path(processor_path, "ref"), exist_ok=True)
            elif processor_type == 'sequence':
                for subdir in self.sequence_dirs.values():
                    os.makedirs(join_path(processor_path, subdir), exist_ok=True)
            elif processor_type in ['test', 'test_processor', 'simple']:
                for subdir in self.test_dirs.values():
                    os.makedirs(join_path(processor_path, subdir), exist_ok=True)
    
    def _validate_directory_structure(self):
        """
        Validate that the directory structure is as expected.
        
        Logs warnings if directories are missing but doesn't raise exceptions.
        """
        # Check data root
        if not os.path.exists(self.data_root):
            logger.warning(f"Data root directory does not exist: {self.data_root}")
            return
        
        # Check processor directories
        for processor_type, dir_name in self.processor_dirs.items():
            processor_path = join_path(self.data_root, dir_name)
            if not os.path.exists(processor_path):
                logger.debug(f"Processor directory does not exist: {processor_path}")
    
    def get_processor_path(self, 
                           processor_type: str, 
                           source: DataSource = DataSource.AUTO) -> str:
        """
        Get the path for a specific processor type.
        
        Args:
            processor_type: Type of processor ('structure', 'grn', etc.)
            source: Ignored - kept for backward compatibility
            
        Returns:
            Full path to the processor directory
        """
        if source != DataSource.AUTO:
            logger.debug(f"DataSource parameter '{source}' is deprecated and ignored")
            
        # If processor type is not in the mapping, use it directly as directory name
        if processor_type not in self.processor_dirs:
            logger.debug(f"Unknown processor type '{processor_type}', using as directory name")
            return join_path(self.data_root, processor_type)
            
        return join_path(self.data_root, self.processor_dirs[processor_type])
    
    def _resolve_data_root(self, source: DataSource) -> str:
        """
        Resolve data root - always returns the single data root.
        
        Args:
            source: Ignored - kept for backward compatibility
            
        Returns:
            Path to the data root
        """
        if source != DataSource.AUTO:
            logger.debug(f"DataSource parameter '{source}' is deprecated")
        return self.data_root
            
    def get_structure_subdir_path(self, 
                                  subdir_type: str, 
                                  source: DataSource = DataSource.AUTO) -> str:
        """
        Get the path for a structure subdirectory.
        
        Args:
            subdir_type: Type of subdirectory ('structure_dir', 'dataset_dir', etc.)
            source: Ignored - kept for backward compatibility
            
        Returns:
            Full path to the structure subdirectory
            
        Raises:
            ValueError: If subdirectory type is not recognized
        """
        if subdir_type not in self.structure_dirs:
            raise ValueError(f"Unknown structure subdirectory type: {subdir_type}")
            
        structure_path = self.get_processor_path('structure', source)
        return join_path(structure_path, self.structure_dirs[subdir_type])
    
    def get_grn_subdir_path(self, 
                           subdir_type: str, 
                           source: DataSource = DataSource.AUTO) -> str:
        """
        Get the path for a GRN subdirectory.
        
        Args:
            subdir_type: Type of subdirectory ('table_dir', 'configs_dir', etc.)
            source: Ignored - kept for backward compatibility
            
        Returns:
            Full path to the GRN subdirectory
            
        Raises:
            ValueError: If subdirectory type is not recognized
        """
        # Handle special case for reference tables
        if subdir_type == "ref" or subdir_type == "ref_dir":
            grn_path = self.get_processor_path('grn', source)
            return join_path(grn_path, "ref")
            
        if subdir_type not in self.grn_dirs:
            raise ValueError(f"Unknown GRN subdirectory type: {subdir_type}")
            
        grn_path = self.get_processor_path('grn', source)
        return join_path(grn_path, self.grn_dirs[subdir_type])
    
    def get_sequence_subdir_path(self, 
                                subdir_type: str, 
                                source: DataSource = DataSource.AUTO) -> str:
        """
        Get the path for a sequence subdirectory.
        
        Args:
            subdir_type: Type of subdirectory ('fasta_dir', 'alignment_dir', etc.)
            source: Ignored - kept for backward compatibility
            
        Returns:
            Full path to the sequence subdirectory
            
        Raises:
            ValueError: If subdirectory type is not recognized
        """
        if subdir_type not in self.sequence_dirs:
            raise ValueError(f"Unknown sequence subdirectory type: {subdir_type}")
            
        sequence_path = self.get_processor_path('sequence', source)
        return join_path(sequence_path, self.sequence_dirs[subdir_type])
    
    def get_registry_path(self, 
                         processor_type: str, 
                         source: DataSource = DataSource.USER) -> str:
        """
        Get the path for a registry file.
        
        Args:
            processor_type: Type of processor ('structure', 'grn', etc.)
            source: Ignored - kept for backward compatibility
            
        Returns:
            Full path to the registry file
        """
        processor_path = self.get_processor_path(processor_type, source)
        return join_path(processor_path, DEFAULT_REGISTRY_FILENAME)
    
    def get_global_registry_path(self) -> str:
        """
        Get the path for the global registry file.
        
        Returns:
            Full path to the global registry file
        """
        return join_path(self.data_root, DEFAULT_GLOBAL_REGISTRY_FILENAME)
    
    def get_dataset_path(self, 
                        processor_type: str, 
                        dataset_name: str,
                        source: DataSource = DataSource.AUTO,
                        file_extension: Optional[str] = None) -> str:
        """
        Get the path for a dataset file.
        
        Args:
            processor_type: Type of processor ('structure', 'grn', etc.)
            dataset_name: Name of the dataset
            source: Ignored - kept for backward compatibility
            file_extension: Optional file extension (with dot)
            
        Returns:
            Full path to the dataset file
        """
        processor_path = self.get_processor_path(processor_type, source)
        
        # Get appropriate dataset directory
        if processor_type == 'structure':
            dataset_dir = self.structure_dirs['dataset_dir']
        elif processor_type == 'grn':
            dataset_dir = self.grn_dirs['table_dir']
        elif processor_type == 'sequence':
            dataset_dir = self.sequence_dirs['metadata_dir']
        elif processor_type in ['test', 'test_processor', 'simple']:
            dataset_dir = self.test_dirs['dataset_dir']
        else:
            dataset_dir = 'datasets'
            
        # Add extension if provided
        filename = f"{dataset_name}{file_extension or ''}"
        
        return join_path(processor_path, dataset_dir, filename)
    
    def resolve_path(self, 
                    path: Optional[str], 
                    source: DataSource = DataSource.AUTO,
                    relative_to: Optional[str] = None) -> str:
        """
        Resolve a path, handling relative paths intelligently.
        
        Args:
            path: Path to resolve (absolute or relative)
            source: Ignored - kept for backward compatibility
            relative_to: Base directory for relative paths
            
        Returns:
            Resolved absolute path
        """
        if path is None:
            return self.data_root if relative_to is None else relative_to
        
        path = os.path.expanduser(path)
        
        if os.path.isabs(path):
            return path
        
        if relative_to is not None:
            return join_path(relative_to, path)
        
        # Use data root for relative paths
        return join_path(self.data_root, path)
    
    def exists(self, 
              path: str, 
              check_both_sources: bool = True) -> Tuple[bool, Optional[DataSource]]:
        """
        Check if a path exists.
        
        Args:
            path: Path to check
            check_both_sources: Ignored - kept for backward compatibility
            
        Returns:
            Tuple of (exists, source) where source is always USER if exists
        """
        # First check if the path is absolute and exists
        if os.path.isabs(path) and os.path.exists(path):
            return True, DataSource.USER
        
        # Check in data directory
        full_path = self.resolve_path(path)
        if os.path.exists(full_path):
            return True, DataSource.USER
        
        return False, None
    
    def update_paths(self, 
                    user_data_root: Optional[str] = None, 
                    ref_data_root: Optional[str] = None,
                    data_root: Optional[str] = None,
                    processor_dirs: Optional[Dict[str, str]] = None,
                    structure_dirs: Optional[Dict[str, str]] = None,
                    grn_dirs: Optional[Dict[str, str]] = None,
                    sequence_dirs: Optional[Dict[str, str]] = None):
        """
        Update path configurations.
        
        Args:
            user_data_root: Deprecated - use data_root
            ref_data_root: Deprecated - ignored
            data_root: New data root directory
            processor_dirs: New processor directory mapping
            structure_dirs: New structure subdirectory mapping
            grn_dirs: New GRN subdirectory mapping
            sequence_dirs: New sequence subdirectory mapping
        """
        if user_data_root or ref_data_root:
            logger.warning(
                "user_data_root and ref_data_root are deprecated. "
                "Use data_root parameter instead."
            )
            data_root = data_root or user_data_root
            
        if data_root is not None:
            self.data_root = os.path.expanduser(data_root)
            if not os.path.isabs(self.data_root):
                self.data_root = os.path.abspath(self.data_root)
            # Update backward compatibility attributes
            self.user_data_root = self.data_root
            self.ref_data_root = self.data_root
        
        if processor_dirs is not None:
            self.processor_dirs.update(processor_dirs)
            
        if structure_dirs is not None:
            self.structure_dirs.update(structure_dirs)
            
        if grn_dirs is not None:
            self.grn_dirs.update(grn_dirs)
            
        if sequence_dirs is not None:
            self.sequence_dirs.update(sequence_dirs)


# Global helper functions that delegate to the default instance

_DEFAULT_PATH_RESOLVER = None

def get_default_resolver():
    """Get or create the default path resolver."""
    global _DEFAULT_PATH_RESOLVER
    if _DEFAULT_PATH_RESOLVER is None:
        _DEFAULT_PATH_RESOLVER = ProtosPaths(create_dirs=True, validate=True)
    return _DEFAULT_PATH_RESOLVER

def resolve_path(path: Optional[str], 
                relative_to: Optional[str] = None,
                source: DataSource = DataSource.AUTO) -> str:
    """
    Resolve a path, handling relative paths intelligently.
    
    Args:
        path: Path to resolve (absolute or relative)
        relative_to: Base directory for relative paths
        source: Ignored - kept for backward compatibility
        
    Returns:
        Resolved absolute path
    """
    return get_default_resolver().resolve_path(path, source, relative_to)

def get_structure_path(pdb_id: str, 
                      structure_dir: Optional[str] = None,
                      source: DataSource = DataSource.AUTO,
                      create_if_missing: bool = False) -> str:
    """
    Get the path for a structure file.
    
    Args:
        pdb_id: PDB identifier
        structure_dir: Optional custom directory for structure files
        source: Ignored - kept for backward compatibility
        create_if_missing: Whether to create parent directories if they don't exist
        
    Returns:
        Path to the structure file
    """
    # Preserve original PDB ID to avoid scientific notation issues
    original_pdb_id = str(pdb_id)
    
    if structure_dir is not None:
        path = join_path(structure_dir, f"{original_pdb_id}.cif")
    else:
        structure_dir = get_default_resolver().get_structure_subdir_path('structure_dir', source)
        path = join_path(structure_dir, f"{original_pdb_id}.cif")
    
    # Create parent directory if requested
    if create_if_missing:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        
    return path

def get_grn_path(table_name: str, 
                table_dir: Optional[str] = None,
                source: DataSource = DataSource.AUTO) -> str:
    """
    Get the path for a GRN table file.
    
    Args:
        table_name: Name of the GRN table
        table_dir: Optional custom directory for GRN tables
        source: Ignored - kept for backward compatibility
        
    Returns:
        Path to the GRN table file
    """
    if table_dir is not None:
        return join_path(table_dir, f"{table_name}.csv")
    
    table_dir = get_default_resolver().get_grn_subdir_path('table_dir', source)
    return join_path(table_dir, f"{table_name}.csv")

def get_sequence_path(sequence_id: str, 
                     fasta_dir: Optional[str] = None,
                     source: DataSource = DataSource.AUTO) -> str:
    """
    Get the path for a sequence file.
    
    Args:
        sequence_id: Sequence identifier
        fasta_dir: Optional custom directory for FASTA files
        source: Ignored - kept for backward compatibility
        
    Returns:
        Path to the sequence file
    """
    if fasta_dir is not None:
        return join_path(fasta_dir, f"{sequence_id}.fasta")
    
    fasta_dir = get_default_resolver().get_sequence_subdir_path('fasta_dir', source)
    return join_path(fasta_dir, f"{sequence_id}.fasta")

def get_dataset_path(dataset_name: str, 
                    processor_type: str = 'structure',
                    file_extension: str = '.json',
                    source: DataSource = DataSource.AUTO,
                    create_if_missing: bool = False) -> str:
    """
    Get the path for a dataset file.
    
    Args:
        dataset_name: Name of the dataset
        processor_type: Type of processor ('structure', 'grn', etc.)
        file_extension: File extension with dot
        source: Ignored - kept for backward compatibility
        create_if_missing: Whether to create parent directories if they don't exist
        
    Returns:
        Path to the dataset file
    """
    path = get_default_resolver().get_dataset_path(
        processor_type, dataset_name, source, file_extension)
    
    # Create parent directory if requested
    if create_if_missing:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        
    return path


def get_data_root(source: DataSource = DataSource.USER) -> str:
    """
    Get the data root directory.
    
    Args:
        source: Ignored - kept for backward compatibility
        
    Returns:
        Path to the data root
    """
    return get_default_resolver().data_root


# Functions for backward compatibility
def get_user_data_root() -> str:
    """
    Get the user data root directory.
    Deprecated - use get_data_root() instead.
    
    Returns:
        Data root directory path
    """
    logger.warning("get_user_data_root() is deprecated. Use get_data_root() instead.")
    return get_data_root()


def get_reference_data_root() -> str:
    """
    Get the reference data root directory.
    Deprecated - use get_data_root() instead.
    
    Returns:
        Data root directory path
    """
    logger.warning("get_reference_data_root() is deprecated. Use get_data_root() instead.")
    return get_data_root()


def ensure_directory(directory: Union[str, Path]) -> str:
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        directory: Directory path to ensure
        
    Returns:
        Normalized absolute path to the directory
    """
    # Convert to Path object for reliable handling
    dir_path = Path(directory).expanduser().resolve()
    
    # Create if it doesn't exist
    os.makedirs(dir_path, exist_ok=True)
    
    # Return normalized string path
    return str(dir_path)


def is_package_resource(path: str) -> bool:
    """
    Check if a path is within the package resource directory.
    Deprecated - always returns False in simplified version.
    
    Args:
        path: Path to check
        
    Returns:
        False (no separate reference data in simplified version)
    """
    return False