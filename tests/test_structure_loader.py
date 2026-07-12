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
        {
            "source": "rcsb",
            "existing_note": "preserve me",
            "content_hash": "stale-existing-hash",
        },
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
    assert info.metadata["content_hash"] != "stale-existing-hash"
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


def test_import_local_structures_derives_safe_default_names(
    monkeypatch: pytest.MonkeyPatch,
):
    loader = StructureLoader(name="test_structure_loader")
    raw_path = Path(loader.paths.data_root) / "raw.cif"
    raw_path.write_text("raw CIF fixture", encoding="utf-8")
    monkeypatch.setattr(loader, "fetch_entity", lambda *args, **kwargs: raw_path)
    monkeypatch.setattr(
        structure_loader_module,
        "load_structure_from_cif",
        lambda *args, **kwargs: _parsed_structure(),
    )

    successful, failed = loader.import_local_structures(
        ["/outside/path/first.cif", "relative/path/second.mmcif.gz"]
    )

    assert successful == ["first", "second"]
    assert failed == []
    assert (loader._get_processor().path_pkl_dir / "first.pkl").is_file()
    assert (loader._get_processor().path_pkl_dir / "second.pkl").is_file()

    explicit, failed = loader.import_local_structures(
        ["/outside/path/third.cif"],
        names=["chosen-name"],
    )
    assert explicit == ["chosen-name"]
    assert failed == []


def test_download_and_register_replaces_stale_integrity_metadata(
    monkeypatch: pytest.MonkeyPatch,
):
    loader = StructureLoader(name="test_structure_loader")
    raw_path = Path(loader.paths.data_root) / "metadata.cif"
    raw_path.write_text("raw CIF fixture", encoding="utf-8")
    monkeypatch.setattr(loader, "fetch_entity", lambda *args, **kwargs: raw_path)
    monkeypatch.setattr(
        structure_loader_module,
        "load_structure_from_cif",
        lambda *args, **kwargs: _parsed_structure(),
    )

    result = loader.download_and_register(
        "1META",
        metadata={
            "content_hash": "stale-hash",
            "schema_version": "stale-schema",
            "saved_at": "stale-time",
            "atom_count": 999,
        },
    )

    assert result == "1META"
    info = loader.entity_registry.find_entity("1META", "structure")
    assert info is not None
    assert info.metadata["content_hash"] != "stale-hash"
    assert info.metadata["schema_version"] == "1.0"
    assert info.metadata["saved_at"] != "stale-time"
    assert info.metadata["atom_count"] == 1


def test_failure_after_save_restores_overwritten_pkl_and_raw_registration(
    monkeypatch: pytest.MonkeyPatch,
):
    loader = StructureLoader(name="test_structure_loader")
    processor = loader._get_processor()
    raw_path = Path(loader.paths.data_root) / "structure" / "mmcif" / "repair.cif"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text("raw CIF fixture", encoding="utf-8")
    previous_metadata = {"source": "local", "note": "original"}
    loader.entity_registry.register_entity(
        "repair",
        "structure",
        str(raw_path),
        previous_metadata,
    )
    pkl_path = processor.path_pkl_dir / "repair.pkl"
    original_pkl = b"original PKL bytes"
    pkl_path.write_bytes(original_pkl)
    monkeypatch.setattr(
        structure_loader_module,
        "load_structure_from_cif",
        lambda *args, **kwargs: _parsed_structure(),
    )
    monkeypatch.setattr(processor, "load_entity", lambda *args, **kwargs: None)

    assert (
        loader.download_and_register("repair.cif", name="repair", source="local")
        is None
    )
    assert pkl_path.read_bytes() == original_pkl
    info = loader.entity_registry.find_entity("repair", "structure")
    assert info is not None
    assert info.file_path == str(raw_path)
    assert info.metadata == previous_metadata
    assert not list(pkl_path.parent.glob(".repair.*.pkl.bak"))


def test_registry_restoration_failure_retains_restored_pkl_backup(
    monkeypatch: pytest.MonkeyPatch,
):
    loader = StructureLoader(name="test_structure_loader")
    processor = loader._get_processor()
    raw_path = Path(loader.paths.data_root) / "structure" / "mmcif" / "repair.cif"
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.write_text("raw CIF fixture", encoding="utf-8")
    loader.entity_registry.register_entity(
        "repair",
        "structure",
        str(raw_path),
        {"source": "local"},
    )
    pkl_path = processor.path_pkl_dir / "repair.pkl"
    original_pkl = b"original PKL bytes"
    pkl_path.write_bytes(original_pkl)
    monkeypatch.setattr(
        structure_loader_module,
        "load_structure_from_cif",
        lambda *args, **kwargs: _parsed_structure(),
    )
    monkeypatch.setattr(processor, "load_entity", lambda *args, **kwargs: None)
    register_entity = loader.entity_registry.register_entity

    def fail_raw_restoration(name, format_type, file_path, metadata=None):
        if file_path == str(raw_path):
            raise OSError("registry persistence unavailable")
        return register_entity(name, format_type, file_path, metadata)

    monkeypatch.setattr(loader.entity_registry, "register_entity", fail_raw_restoration)

    assert (
        loader.download_and_register("repair.cif", name="repair", source="local")
        is None
    )
    assert pkl_path.read_bytes() == original_pkl
    backups = list(pkl_path.parent.glob(".repair.*.pkl.bak"))
    assert len(backups) == 1
    assert backups[0].read_bytes() == original_pkl
