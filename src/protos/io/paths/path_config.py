"""
Simplified path configuration for the Protos framework.
"""

import os
import json
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Union

from .path_constants import (
    ENV_DATA_ROOT,
    DEFAULT_PROCESSOR_DIRS,
    DEFAULT_STRUCTURE_SUBDIRS,
    DEFAULT_GRN_SUBDIRS,
    DEFAULT_SEQUENCE_SUBDIRS,
    DEFAULT_GLOBAL_REGISTRY_FILENAME,
    join_path,
    DEFAULT_PROPERTY_SUBDIRS,
    DEFAULT_EMBEDDING_SUBDIRS,
    DEFAULT_LIGAND_SUBDIRS,
    DEFAULT_GRAPH_SUBDIRS,
    DEFAULT_INPUT_SUBDIRS,
    DEFAULT_TEMP_SUBDIRS
)


def get_default_data_root() -> str:
    env_root = os.environ.get(ENV_DATA_ROOT)
    if env_root:
        return os.path.expanduser(env_root)
    return os.path.expanduser("~/protos_data")


class ProtosPaths:
    """
    Path management for Protos using a single data directory.
    Implements singleton pattern with lazy initialization.
    """
    _instance = None
    _initialized = False

    def __new__(cls, data_root: Optional[str] = None):
        if cls._instance is None or data_root is not None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self, data_root: Optional[str] = None):
        # Avoid re-initialization if not needed
        if self._initialized and data_root is None:
            return
            
        self.data_root = data_root or get_default_data_root()
        self.data_root = os.path.abspath(os.path.expanduser(self.data_root))

        self.processor_dirs = DEFAULT_PROCESSOR_DIRS.copy()

        self.subdirs = {
            'structure': DEFAULT_STRUCTURE_SUBDIRS.copy(),
            'grn': {**DEFAULT_GRN_SUBDIRS.copy(), 'reference': 'reference'},  # Maintain backward compatibility
            'sequence': DEFAULT_SEQUENCE_SUBDIRS.copy(),
            'property': DEFAULT_PROPERTY_SUBDIRS.copy(),
            'embedding': DEFAULT_EMBEDDING_SUBDIRS.copy(),
            'ligand': DEFAULT_LIGAND_SUBDIRS.copy(),
            'graph': DEFAULT_GRAPH_SUBDIRS.copy(),
            'input': DEFAULT_INPUT_SUBDIRS.copy(),
            'temp': DEFAULT_TEMP_SUBDIRS.copy(),
            'test': DEFAULT_STRUCTURE_SUBDIRS.copy(),  # Use structure subdirs for test compatibility
            'test_processor': DEFAULT_STRUCTURE_SUBDIRS.copy(),
            'simple': DEFAULT_STRUCTURE_SUBDIRS.copy()
        }
        
        # Mark as not yet initialized (directories not created)
        self._initialized = False

    def _ensure_initialized(self):
        """Lazy initialization - called on first actual use."""
        if not self._initialized:
            # Check if data directory exists
            data_path = Path(self.data_root)
            if not data_path.exists():
                self._initialize_directory_structure()
                self._initialize_registry()
                self._install_reference_data()
            else:
                # Directory exists, ensure it has all required subdirs
                self._ensure_complete_structure()
            
            self._initialized = True
    
    def _initialize_directory_structure(self):
        """Create all required directories."""
        # Create base directory
        Path(self.data_root).mkdir(parents=True, exist_ok=True)
        
        # Create all processor directories and their subdirectories
        for processor, subdirs in self.subdirs.items():
            if processor in ['test', 'test_processor', 'simple']:
                continue  # Skip test processors
            processor_path = Path(self.data_root) / self.processor_dirs.get(processor, processor)
            processor_path.mkdir(exist_ok=True)
            
            for subdir_name in subdirs.values():
                if subdir_name:  # Skip empty subdirs
                    (processor_path / subdir_name).mkdir(parents=True, exist_ok=True)
    
    def _initialize_registry(self):
        """Create global registry with proper structure."""
        registry_path = Path(self.get_global_registry_path())
        
        if not registry_path.exists():
            registry_path.parent.mkdir(parents=True, exist_ok=True)
            registry_data = {
                "entities": {},
                "name_index": {},
                "version": "2.0"  # UUID-based version
            }
            with open(registry_path, 'w') as f:
                json.dump(registry_data, f, indent=2)
    
    def _install_reference_data(self):
        """Copy reference data from package to data directory if not present."""
        # Only install if this is a fresh data directory
        marker_file = Path(self.data_root) / '.protos_initialized'
        if marker_file.exists():
            return  # Already initialized
        
        # Find reference data in package
        try:
            import protos
            package_dir = Path(protos.__file__).parent
            ref_data_src = package_dir / 'reference_data'
            
            if ref_data_src.exists():
                # Copy GRN reference data
                grn_src = ref_data_src / 'grn'
                if grn_src.exists():
                    grn_dest = Path(self.data_root) / 'grn'

                    configs_src = grn_src / 'configs'
                    if configs_src.exists():
                        shutil.copytree(configs_src, grn_dest / 'configs', dirs_exist_ok=True)

                    reference_dest = grn_dest / 'reference'
                    reference_candidates = [grn_src / 'reference', grn_src / 'ref']
                    copied_reference = False
                    for candidate in reference_candidates:
                        if candidate.exists():
                            shutil.copytree(candidate, reference_dest, dirs_exist_ok=True)
                            copied_reference = True

                    if not copied_reference:
                        reference_dest.mkdir(parents=True, exist_ok=True)
        except Exception:
            # Reference data not available, skip
            pass
        
        # Mark as initialized
        marker_file.write_text(f"Initialized on {datetime.now().isoformat()}\n")
    
    def _ensure_complete_structure(self):
        """Ensure existing directory has all required subdirectories."""
        # Check and create any missing directories
        for processor, subdirs in self.subdirs.items():
            if processor in ['test', 'test_processor', 'simple']:
                continue  # Skip test processors
            processor_path = Path(self.data_root) / self.processor_dirs.get(processor, processor)
            
            for subdir_name in subdirs.values():
                if subdir_name:  # Skip empty subdirs
                    subdir_path = processor_path / subdir_name
                    if not subdir_path.exists():
                        subdir_path.mkdir(parents=True, exist_ok=True)
        
        # Ensure registry exists
        registry_path = Path(self.get_global_registry_path())
        if not registry_path.exists():
            self._initialize_registry()
        
        # Install reference data if missing
        grn_config = Path(self.data_root) / 'grn' / 'configs' / 'config.json'
        if not grn_config.exists():
            self._install_reference_data()

    def get_processor_path(self, processor_type: str) -> str:
        """
        Get processor directory path - triggers initialization.
        """
        self._ensure_initialized()  # Lazy init on first use
        dir_name = self.processor_dirs.get(processor_type, processor_type)
        path = join_path(self.data_root, dir_name)
        os.makedirs(path, exist_ok=True)
        return path

    def get_subdir_path(self, processor_type: str, subdir_type: str) -> str:
        """
        Get subdirectory path for a processor - triggers initialization.

        Raises:
            ValueError: If processor or subdir type unknown.
        """
        self._ensure_initialized()  # Lazy init on first use
        
        if processor_type not in self.subdirs:
            raise ValueError(f"Unknown processor type: {processor_type}")

        subdirs = self.subdirs[processor_type]
        if subdir_type not in subdirs:
            raise ValueError(f"Unknown subdir type for {processor_type}: {subdir_type}")

        processor_path = self.get_processor_path(processor_type)
        path = join_path(processor_path, subdirs[subdir_type])
        os.makedirs(path, exist_ok=True)
        return path

    # DEPRECATED: _ensure_registry and get_registry_path removed
    # Use unified EntityRegistry from io.core instead

    def get_global_registry_path(self) -> str:
        """
        Get global registry file path.
        """
        path = join_path(self.data_root, DEFAULT_GLOBAL_REGISTRY_FILENAME)
        # Registry initialization handled by EntityRegistry
        return path

    def get_dataset_path(self,
                         processor_type: str,
                         dataset_name: str,
                         file_extension: str = '.json') -> str:
        """
        Get dataset definition file path.
        """
        processor_path = self.get_processor_path(processor_type)
        dataset_dir = join_path(processor_path, 'datasets')
        os.makedirs(dataset_dir, exist_ok=True)
        return join_path(dataset_dir, f"{dataset_name}{file_extension}")

    def resolve_path(self,
                     path: Optional[str],
                     relative_to: Optional[str] = None) -> str:
        """
        Resolve path to absolute.
        """
        if path is None:
            return self.data_root if relative_to is None else relative_to

        path = os.path.expanduser(path)

        if os.path.isabs(path):
            return path

        base = relative_to or self.data_root
        return join_path(base, path)

    def exists(self, path: str) -> bool:
        """
        Check if path exists.
        """
        full_path = self.resolve_path(path)
        return os.path.exists(full_path)

    def update(self,
               data_root: Optional[str] = None,
               processor_dirs: Optional[Dict[str, str]] = None,
               subdirs: Optional[Dict[str, Dict[str, str]]] = None):
        """
        Update configurations.
        """
        if data_root is not None:
            self.data_root = os.path.abspath(os.path.expanduser(data_root))

        if processor_dirs is not None:
            self.processor_dirs.update(processor_dirs)

        if subdirs is not None:
            for pt, sd in subdirs.items():
                if pt in self.subdirs:
                    self.subdirs[pt].update(sd)
                else:
                    self.subdirs[pt] = sd.copy()


def resolve_path(paths: ProtosPaths, path: Optional[str], relative_to: Optional[str] = None) -> str:
    return paths.resolve_path(path, relative_to)


def get_structure_path(paths: ProtosPaths, pdb_id: str,
                       structure_dir: Optional[str] = None,
                       create_if_missing: bool = False) -> str:
    original_pdb_id = str(pdb_id)

    if structure_dir is not None:
        path = join_path(structure_dir, f"{original_pdb_id}.cif")
    else:
        structure_dir = paths.get_subdir_path('structure', 'structure_dir')
        path = join_path(structure_dir, f"{original_pdb_id}.cif")

    if create_if_missing:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    return path


def get_grn_path(paths: ProtosPaths, table_name: str,
                 table_dir: Optional[str] = None) -> str:
    if table_dir is not None:
        return join_path(table_dir, f"{table_name}.csv")

    table_dir = paths.get_subdir_path('grn', 'table_dir')
    return join_path(table_dir, f"{table_name}.csv")


def get_sequence_path(paths: ProtosPaths, sequence_id: str,
                      fasta_dir: Optional[str] = None) -> str:
    if fasta_dir is not None:
        return join_path(fasta_dir, f"{sequence_id}.fasta")

    fasta_dir = paths.get_subdir_path('sequence', 'fasta_dir')
    return join_path(fasta_dir, f"{sequence_id}.fasta")


def get_dataset_path(paths: ProtosPaths, dataset_name: str,
                     processor_type: str = 'structure',
                     file_extension: str = '.json',
                     create_if_missing: bool = False) -> str:
    path = paths.get_dataset_path(processor_type, dataset_name, file_extension)

    if create_if_missing:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    return path


def get_data_root(paths: ProtosPaths) -> str:
    return paths.data_root


def ensure_directory(directory: Union[str, Path]) -> str:
    dir_path = Path(directory).expanduser().resolve()
    os.makedirs(dir_path, exist_ok=True)
    return str(dir_path)


# Singleton access functions
_paths_instance = None


def get_protos_paths(data_root: Optional[str] = None) -> ProtosPaths:
    """
    Get or create ProtosPaths singleton.
    
    IMPORTANT: If data_root is provided after initialization has occurred,
    this will raise an error! Path must be set BEFORE any processor usage.
    """
    global _paths_instance
    
    if _paths_instance is None:
        # First time - create instance
        _paths_instance = ProtosPaths(data_root)
    elif data_root is not None:
        # Trying to change path after initialization
        if _paths_instance._initialized:
            raise RuntimeError(
                "Cannot change ProtosPaths after initialization! "
                "Call set_data_path() BEFORE creating any processors."
            )
        # Not initialized yet, safe to change
        _paths_instance = ProtosPaths(data_root)
    
    return _paths_instance


def set_protos_paths(data_root: str) -> ProtosPaths:
    """
    Update the global ProtosPaths instance with a new location.
    
    Raises:
        RuntimeError: If paths have already been initialized.
    """
    return get_protos_paths(data_root)
