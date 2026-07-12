from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import protos.io.core as core
import protos.io.paths.path_config as path_config
import protos.io.ingest.structure_loader as structure_loader_module
from protos.io.core.entity_registry import EntityRegistry
from protos.io.ingest.structure_loader import StructureLoader


@pytest.fixture(autouse=True)
def isolated_data_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setenv("PROTOS_DATA_ROOT", str(tmp_path))
    path_config._paths_instance = None
    path_config.ProtosPaths._instance = None
    core._registry_instance = None
    EntityRegistry._instance = None
    yield
    path_config._paths_instance = None
    path_config.ProtosPaths._instance = None
    core._registry_instance = None
    EntityRegistry._instance = None


def _parsed_structure() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "atom_id": [1],
            "group": ["ATOM"],
            "auth_chain_id": ["A"],
            "auth_seq_id": [1],
            "auth_comp_id": ["GLY"],
            "atom_name": ["CA"],
            "x": [1.0],
            "y": [2.0],
            "z": [3.0],
        }
    )


def test_download_and_register_repairs_raw_registration(monkeypatch: pytest.MonkeyPatch):
    loader = StructureLoader(name="test_structure_loader")
    raw_path = Path(loader.paths.data_root) / "structure" / "mmcif" / "1ubq.cif"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text("raw CIF fixture", encoding="utf-8")
    loader.entity_registry.register_entity(
        "1UBQ",
        "structure",
        str(raw_path),
        {"source": "rcsb", "existing_note": "preserve me"},
    )

    monkeypatch.setattr(
        loader,
        "fetch_entity",
        lambda *args, **kwargs: pytest.fail("existing raw CIF should be reused"),
    )
    monkeypatch.setattr(
        structure_loader_module,
        "load_structure_from_cif",
        lambda *args, **kwargs: _parsed_structure(),
    )

    result = loader.download_and_register(
        "1UBQ",
        metadata={"request_note": "new metadata"},
    )

    assert result == "1UBQ"
    info = loader.entity_registry.find_entity("1UBQ", "structure")
    assert info is not None
    assert Path(info.file_path).suffix == ".pkl"
    assert (Path(loader.paths.data_root) / info.file_path).is_file()
    assert info.metadata["source"] == "rcsb"
    assert info.metadata["source_file"] == str(raw_path)
    assert info.metadata["existing_note"] == "preserve me"
    assert info.metadata["request_note"] == "new metadata"
    loaded = loader._get_processor().load_entity("1UBQ")
    assert loaded is not None and not loaded.empty
    assert loaded.index.names == ["structure_id", "atom_id"]


def test_download_and_register_rolls_back_failed_canonicalization(
    monkeypatch: pytest.MonkeyPatch,
):
    loader = StructureLoader(name="test_structure_loader")
    raw_path = Path(loader.paths.data_root) / "structure" / "mmcif" / "bad.cif"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text("invalid CIF", encoding="utf-8")
    monkeypatch.setattr(loader, "fetch_entity", lambda *args, **kwargs: raw_path)

    def fail_parse(*args, **kwargs):
        raise ValueError("invalid atom_site data")

    monkeypatch.setattr(structure_loader_module, "load_structure_from_cif", fail_parse)

    assert loader.download_and_register("1BAD", name="broken") is None
    assert loader.entity_registry.find_entity("broken", "structure") is None
    assert not (loader._get_processor().path_pkl_dir / "broken.pkl").exists()
