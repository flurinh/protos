"""Ingest package providing loaders for converting raw data to processor formats."""

from __future__ import annotations

import importlib
from typing import Any

__all__ = [
    'alphafold_utils',
    'ccd_index_builder',
    'ccd_loader',
    'ccd_loader_unified',
    'chembl_loader',
    'download_structures',
    'enamine_loader',
    'gpcrdb_loader',
    'gpcrdb_loader_utils',
    'ligand_loader',
    'ligand_utils',
    'ncbi_loader',
    'ncbi_utils',
    'qm9_loader',
    'sequence_loader',
    'structure_loader',
    'uniprot_loader',
    'uniprot_utils',
]


def __getattr__(name: str) -> Any:
    if name in __all__:
        module = importlib.import_module(f'.{name}', __name__)
        globals()[name] = module
        return module
    raise AttributeError(name)


def __dir__():
    return sorted(__all__)
