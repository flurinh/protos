#!/usr/bin/env python3
"""TDD scaffold for zero-config dataset registration helpers."""
from __future__ import annotations

from pathlib import Path

import pytest

import sys
_PACKAGE_SRC = Path(__file__).resolve().parents[1] / 'protos' / 'src'
if _PACKAGE_SRC.exists() and str(_PACKAGE_SRC) not in sys.path:
    sys.path.insert(0, str(_PACKAGE_SRC))

import protos
from protos.datasets import (
    register_chembl_ligand_dataset,
    register_gpcr_property_dataset,
    register_gpcr_sequence_dataset,
    register_rhodopsin_graph_dataset,
    register_rhodopsin_structure_dataset,
)
from protos.processing.graph import GraphProcessor
from protos.processing.molecule import MoleculeProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor


@pytest.fixture()
def isolated_data_root(tmp_path) -> Path:
    protos.set_data_path(str(tmp_path))
    return tmp_path


def test_register_gpcr_sequence_dataset(isolated_data_root: Path) -> None:
    dataset_name = register_gpcr_sequence_dataset()
    processor = SequenceProcessor()
    assert dataset_name in processor.list_datasets()
    sequences = processor.load_dataset(dataset_name)
    assert len(sequences) > 0


def test_register_rhodopsin_structure_dataset(isolated_data_root: Path) -> None:
    dataset_name = register_rhodopsin_structure_dataset()
    processor = StructureProcessor()
    assert dataset_name in processor.list_datasets()
    structures = processor.dataset_manager.get_dataset_entities(dataset_name)
    assert len(structures) > 0


def test_register_chembl_ligand_dataset(isolated_data_root: Path) -> None:
    dataset_name = register_chembl_ligand_dataset()
    processor = MoleculeProcessor()
    datasets = processor.dataset_manager.list_datasets()
    assert dataset_name in datasets


def test_register_gpcr_property_dataset(isolated_data_root: Path) -> None:
    dataset_name = register_gpcr_property_dataset()
    processor = PropertyProcessor()
    datasets = processor.dataset_manager.list_datasets()
    assert dataset_name in datasets


def test_register_rhodopsin_graph_dataset(isolated_data_root: Path) -> None:
    dataset_name = register_rhodopsin_graph_dataset()
    processor = GraphProcessor()
    assert dataset_name in processor.list_graphs()
