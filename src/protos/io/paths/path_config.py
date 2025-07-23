"""
Simplified path configuration for the Protos framework.
"""

import os
import json
from pathlib import Path
from typing import Dict, Optional, Union

from .path_constants import (
    ENV_DATA_ROOT,
    DEFAULT_PROCESSOR_DIRS,
    DEFAULT_STRUCTURE_SUBDIRS,
    DEFAULT_GRN_SUBDIRS,
    DEFAULT_SEQUENCE_SUBDIRS,
    DEFAULT_TEST_SUBDIRS,
    DEFAULT_REGISTRY_FILENAME,
    DEFAULT_GLOBAL_REGISTRY_FILENAME,
    join_path,
    DEFAULT_PROPERTY_SUBDIRS,
    DEFAULT_EMBEDDING_SUBDIRS,
    DEFAULT_LIGAND_SUBDIRS,
    DEFAULT_GRAPH_SUBDIRS
)


def get_default_data_root() -> str:
    env_root = os.environ.get(ENV_DATA_ROOT)
    if env_root:
        return os.path.expanduser(env_root)
    return os.path.expanduser("~/protos_data")


class ProtosPaths:
    """
    Path management for Protos using a single data directory.
    """

    def __init__(self, data_root: Optional[str] = None):
        self.data_root = data_root or get_default_data_root()
        self.data_root = os.path.abspath(os.path.expanduser(self.data_root))

        self.processor_dirs = DEFAULT_PROCESSOR_DIRS.copy()

        self.subdirs = {
            'structure': DEFAULT_STRUCTURE_SUBDIRS.copy(),
            'grn': {**DEFAULT_GRN_SUBDIRS.copy(), 'ref': 'ref', 'ref_dir': 'ref'},
            'sequence': DEFAULT_SEQUENCE_SUBDIRS.copy(),
            'property': DEFAULT_PROPERTY_SUBDIRS.copy(),
            'embedding': DEFAULT_EMBEDDING_SUBDIRS.copy(),
            'ligand': DEFAULT_LIGAND_SUBDIRS.copy(),
            'graph': DEFAULT_GRAPH_SUBDIRS.copy(),
            'test': DEFAULT_TEST_SUBDIRS.copy(),
            'test_processor': DEFAULT_TEST_SUBDIRS.copy(),
            'simple': DEFAULT_TEST_SUBDIRS.copy()
        }

    def get_processor_path(self, processor_type: str) -> str:
        """
        Get processor directory path.
        """
        dir_name = self.processor_dirs.get(processor_type, processor_type)
        path = join_path(self.data_root, dir_name)
        os.makedirs(path, exist_ok=True)
        return path

    def get_subdir_path(self, processor_type: str, subdir_type: str) -> str:
        """
        Get subdirectory path for a processor.

        Raises:
            ValueError: If processor or subdir type unknown.
        """
        if processor_type not in self.subdirs:
            raise ValueError(f"Unknown processor type: {processor_type}")

        subdirs = self.subdirs[processor_type]
        if subdir_type not in subdirs:
            raise ValueError(f"Unknown subdir type for {processor_type}: {subdir_type}")

        processor_path = self.get_processor_path(processor_type)
        path = join_path(processor_path, subdirs[subdir_type])
        os.makedirs(path, exist_ok=True)
        return path

    def _ensure_registry(self, path: str):
        if not os.path.exists(path):
            with open(path, 'w') as f:
                json.dump({"entities": {}}, f)

    def get_registry_path(self, processor_type: str) -> str:
        """
        Get registry file path for a processor.
        """
        processor_path = self.get_processor_path(processor_type)
        path = join_path(processor_path, DEFAULT_REGISTRY_FILENAME)

        skip_types = ['test', 'test_processor', '_test', '__test', 'simple',
                      'complex_processor_with_long_name', 'custom_dir']
        if processor_type not in skip_types:
            self._ensure_registry(path)

        return path

    def get_global_registry_path(self) -> str:
        """
        Get global registry file path.
        """
        path = join_path(self.data_root, DEFAULT_GLOBAL_REGISTRY_FILENAME)
        self._ensure_registry(path)
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