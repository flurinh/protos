"""Shared helpers for managing Protos data directories."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Optional

from protos.io.paths import path_config


def reset_paths_singleton() -> None:
    path_config._paths_instance = None  # type: ignore[attr-defined]
    path_config.ProtosPaths._instance = None  # type: ignore[attr-defined]


def resolve_target_path(path_option: Optional[Path]) -> Path:
    if path_option is not None:
        return path_option.expanduser().resolve()
    default_root = Path(path_config.get_default_data_root()).expanduser().resolve()
    return default_root


def collect_summary(data_root: Path) -> Dict[str, object]:
    if not data_root.exists():
        return {
            "data_root": data_root,
            "processors": [],
            "directory_count": 0,
            "file_count": 0,
            "reference_installed": False,
        }

    processors = sorted(p.name for p in data_root.iterdir() if p.is_dir())
    directory_count = sum(1 for p in data_root.rglob("*") if p.is_dir())
    file_count = sum(1 for p in data_root.rglob("*") if p.is_file())
    reference_marker = data_root / ".protos_initialized"
    return {
        "data_root": data_root,
        "processors": processors,
        "directory_count": directory_count,
        "file_count": file_count,
        "reference_installed": reference_marker.exists(),
    }


def initialize_data_root(target_root: Path) -> Dict[str, object]:
    reset_paths_singleton()
    paths = path_config.get_protos_paths(str(target_root))
    paths.get_processor_path("structure")
    return collect_summary(Path(paths.data_root))


def remove_data_root(target_root: Path) -> None:
    if target_root.exists():
        shutil.rmtree(target_root)
        reset_paths_singleton()


def clear_data_root(target_root: Path, *, reinitialize: bool = True) -> Dict[str, object]:
    remove_data_root(target_root)
    if reinitialize:
        return initialize_data_root(target_root)
    return {
        "data_root": target_root,
        "processors": [],
        "directory_count": 0,
        "file_count": 0,
        "reference_installed": False,
    }
