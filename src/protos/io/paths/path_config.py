"""
Simplified path configuration for the Protos framework.
"""

import os
import hashlib
import json
import re
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Union

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
    DEFAULT_MOLECULE_SUBDIRS,
    DEFAULT_GRAPH_SUBDIRS,
    DEFAULT_INPUT_SUBDIRS,
    DEFAULT_TEMP_SUBDIRS,
)


def get_default_data_root() -> str:
    env_root = os.environ.get(ENV_DATA_ROOT)
    if env_root:
        return os.path.expanduser(env_root)
    return os.path.expanduser("~/protos_data")


# ----------------------- Structure path helpers -----------------------

def get_structure_entity_path(paths: 'ProtosPaths', name: str, *, extension: str = 'cif') -> Path:
    """Compute default export path for a structure entity based on extension.

    - 'cif' -> structure/mmcif/{name}.cif
    - 'sdf' -> structure/sdf/{name}.sdf
    """
    ext = extension.lower().lstrip('.')
    if ext == 'cif':
        base = Path(paths.get_subdir_path('structure', 'structure_dir'))
    elif ext == 'sdf':
        base = Path(paths.get_subdir_path('structure', 'sdf_dir'))
    elif ext == 'pdb':
        base = Path(paths.get_subdir_path('structure', 'pdb_dir'))
    else:
        base = Path(paths.get_processor_path('structure'))
    return base / f"{name}.{ext}"


def get_structure_dataset_path(paths: 'ProtosPaths', dataset_name: str, *, extension: str = 'cif') -> Path:
    """Compute default export path for a structure dataset based on extension.

    - 'cif' -> structure/mmcif/{dataset_name}.cif (not typical; usually per-entity)
    - 'sdf' -> structure/sdf/{dataset_name}.sdf
    """
    ext = extension.lower().lstrip('.')
    if ext == 'sdf':
        base = Path(paths.get_subdir_path('structure', 'sdf_dir'))
    elif ext == 'pdb':
        base = Path(paths.get_subdir_path('structure', 'pdb_dir'))
    else:
        base = Path(paths.get_subdir_path('structure', 'structure_dir'))
    return base / f"{dataset_name}.{ext}"


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
            'molecule': DEFAULT_MOLECULE_SUBDIRS.copy(),
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

    def reinitialize(
        self,
        *,
        wipe: bool = False,
        reinstall_reference: bool = True,
    ) -> None:
        """Reset directory layout and registry state for the configured data root.

        Args:
            wipe: Remove the existing data root before recreating it. Use with care.
            reinstall_reference: When False, skip re-installing bundled reference
                data after clearing the root. The directory structure is still
                recreated.
        """

        root_path = Path(self.data_root)

        if wipe and root_path.exists():
            shutil.rmtree(root_path)

        # Mark the instance as uninitialized so subsequent calls rebuild layout
        self._initialized = False

        # Directory creation is idempotent. Do not call _ensure_complete_structure
        # here because that method repairs references, which would violate
        # reinstall_reference=False.
        self._initialize_directory_structure()

        self._initialize_registry()

        if reinstall_reference:
            self._install_reference_data(force=True)

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
    
    @staticmethod
    def _sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def _bundled_grn_path(self) -> Path:
        import protos

        path = Path(protos.__file__).parent / "reference_data" / "grn"
        if not path.is_dir():
            raise FileNotFoundError(f"Bundled GRN reference data not found: {path}")
        return path

    @staticmethod
    def _read_manifest(path: Path) -> Dict[str, object]:
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise ValueError(f"Invalid GRN reference manifest: {path}") from exc
        if not isinstance(manifest.get("files"), dict):
            raise ValueError(f"GRN reference manifest has no file map: {path}")
        return manifest

    def _reference_data_complete(self) -> bool:
        """Return whether data/grn exactly contains the current bundled release."""

        try:
            bundled = self._read_manifest(self._bundled_grn_path() / "manifest.json")
            installed_root = Path(self.data_root) / "grn"
            installed = self._read_manifest(installed_root / "manifest.json")
            if installed.get("bundle_version") != bundled.get("bundle_version"):
                return False
            if installed.get("files") != bundled.get("files"):
                return False
            for name, expected_hash in bundled["files"].items():
                destination = installed_root / "reference" / str(name)
                if not destination.is_file() or self._sha256(destination) != expected_hash:
                    return False
            return (installed_root / "configs" / "config.json").is_file()
        except (FileNotFoundError, TypeError, ValueError):
            return False

    def _install_reference_data(self, *, force: bool = False) -> None:
        """Install or repair the package-local, manifest-authenticated GRN bundle."""

        marker_file = Path(self.data_root) / ".protos_initialized"
        if not force and marker_file.exists() and self._reference_data_complete():
            return

        grn_src = self._bundled_grn_path()
        manifest_src = grn_src / "manifest.json"
        manifest = self._read_manifest(manifest_src)
        reference_src = grn_src / "reference"
        if not reference_src.is_dir():
            reference_src = grn_src / "ref"
        if not reference_src.is_dir():
            raise FileNotFoundError("Bundled GRN reference CSV directory is missing")

        expected_files = {str(name): value for name, value in manifest["files"].items()}
        bundled_names = {path.name for path in reference_src.glob("*.csv")}
        if bundled_names != set(expected_files):
            raise ValueError("Bundled GRN manifest and reference CSV directory disagree")
        for name, expected_hash in expected_files.items():
            if self._sha256(reference_src / name) != expected_hash:
                raise ValueError(f"Bundled GRN checksum mismatch: {name}")

        grn_dest = Path(self.data_root) / "grn"
        reference_dest = grn_dest / "reference"
        grn_dest.mkdir(parents=True, exist_ok=True)
        reference_dest.mkdir(parents=True, exist_ok=True)

        previous_files: set[str] = set()
        previous_manifest = grn_dest / "manifest.json"
        if previous_manifest.is_file():
            try:
                previous_files = set(self._read_manifest(previous_manifest)["files"])
            except (TypeError, ValueError):
                previous_files = set()
        stale_names = (
            {path.name for path in reference_dest.glob("*.csv")}
            if force
            else previous_files - bundled_names
        )
        for stale_name in stale_names:
            (reference_dest / stale_name).unlink(missing_ok=True)

        for metadata_name in (
            "README.md",
            "manifest.json",
            "gpcrdb_provenance.json",
            "opsin_provenance.json",
        ):
            metadata_src = grn_src / metadata_name
            if metadata_src.exists():
                shutil.copy2(metadata_src, grn_dest / metadata_name)

        configs_src = grn_src / "configs"
        if not configs_src.is_dir():
            raise FileNotFoundError("Bundled GRN configuration directory is missing")
        configs_dest = grn_dest / "configs"
        if force and configs_dest.exists():
            shutil.rmtree(configs_dest)
        shutil.copytree(configs_src, configs_dest, dirs_exist_ok=True)
        shutil.copytree(reference_src, reference_dest, dirs_exist_ok=True)

        if not self._reference_data_complete():
            raise RuntimeError("Installed GRN reference bundle failed verification")
        marker_file.write_text(
            f"bundle={manifest.get('bundle_version')}\n"
            f"installed={datetime.now().isoformat()}\n",
            encoding="utf-8",
        )
    
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
        
        # Install or repair reference data if missing, damaged, or outdated.
        if not self._reference_data_complete():
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


def sanitize_storage_name(name: str, *, default: str = "item") -> str:
    """Produce a filesystem-safe name for stored artifacts."""

    if not name:
        sanitized = default
    else:
        sanitized = re.sub(r"\s+", "_", name.strip())
        sanitized = re.sub(r"[^A-Za-z0-9._-]", "_", sanitized)

    # Collapse duplicate underscores and trim
    sanitized = re.sub(r"_+", "_", sanitized).strip("._")

    if not sanitized:
        sanitized = default

    return sanitized[:128]


def _sequence_dir(paths: "ProtosPaths", key: str) -> Path:
    return Path(paths.get_subdir_path("sequence", key))


def get_sequence_entity_paths(
    paths: "ProtosPaths",
    entity_name: str,
    *,
    extensions: Iterable[str] = ("fasta", "fa"),
) -> List[Path]:
    """Return candidate paths for a sequence entity."""

    base_name = sanitize_storage_name(entity_name, default="sequence")
    entity_dir = _sequence_dir(paths, "entity_fasta_dir")
    return [entity_dir / f"{base_name}.{ext}" for ext in extensions]


def get_sequence_entity_path(
    paths: "ProtosPaths",
    entity_name: str,
    *,
    extension: str = "fasta",
) -> Path:
    entity_dir = _sequence_dir(paths, "entity_fasta_dir")
    ext = extension.lstrip(".").lower()
    base_name = sanitize_storage_name(entity_name, default="sequence")
    if base_name.lower().endswith(f".{ext}"):
        base_name = base_name[: -(len(ext) + 1)]
    return entity_dir / f"{base_name}.{ext}"


def get_sequence_dataset_path(
    paths: "ProtosPaths",
    dataset_name: str,
    *,
    extension: str = "fasta",
) -> Path:
    dataset_dir = _sequence_dir(paths, "dataset_fasta_dir")
    ext = extension.lstrip(".").lower()
    base_name = sanitize_storage_name(dataset_name, default="dataset")
    if base_name.lower().endswith(f".{ext}"):
        base_name = base_name[: -(len(ext) + 1)]
    return dataset_dir / f"{base_name}.{ext}"


def get_sequence_dataset_paths(
    paths: "ProtosPaths",
    dataset_name: str,
    *,
    extensions: Iterable[str] = ("fasta", "fa"),
) -> List[Path]:
    dataset_dir = _sequence_dir(paths, "dataset_fasta_dir")
    base_name = sanitize_storage_name(dataset_name, default="dataset")

    candidates: List[Path] = []
    if "." in dataset_name:
        stem, ext = dataset_name.rsplit(".", 1)
        sanitized_stem = sanitize_storage_name(stem, default="dataset")
        candidates.append(dataset_dir / f"{sanitized_stem}.{ext}")

    candidates.append(dataset_dir / base_name)
    for ext in extensions:
        candidates.append(dataset_dir / f"{base_name}.{ext}")

    unique: List[Path] = []
    seen = set()
    for candidate in candidates:
        if candidate not in seen:
            seen.add(candidate)
            unique.append(candidate)
    return unique


def to_data_relative_path(paths: "ProtosPaths", target: Union[str, Path]) -> str:
    target_path = Path(target).resolve()
    data_root = Path(paths.data_root).resolve()
    return str(target_path.relative_to(data_root))


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


def reset_protos_data(
    data_root: Optional[str] = None,
    *,
    wipe: bool = False,
    reinstall_reference: bool = True,
    backup_registry: bool = True,
):
    """Reset the global Protos data directory and entity registry.

    Args:
        data_root: Optional override for the data root. Must match the
            configured root once processors have been initialized.
        wipe: Remove the existing data directory before recreating it.
        reinstall_reference: Reinstall bundled reference data after reset.
        backup_registry: Persist a timestamped copy of the registry before
            clearing it.
    """

    try:
        paths = get_protos_paths(data_root)
    except RuntimeError as exc:
        raise RuntimeError(
            "Cannot reset Protos data to a new location after initialization. "
            "Call protos.set_data_path() before instantiating processors."
        ) from exc

    from protos.io.core import get_registry

    registry = get_registry()
    registry.reset(backup=backup_registry)

    paths.reinitialize(wipe=wipe, reinstall_reference=reinstall_reference)

    registry.refresh()

    return paths
