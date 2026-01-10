import os

from protos.io.ingest.structure_loader import StructureLoader


def test_structure_loader_accepts_source_aliases(monkeypatch, tmp_path):
    """Ensure user-friendly aliases route to the canonical sources."""

    monkeypatch.setenv("PROTOS_DATA_ROOT", str(tmp_path / "protos_data"))

    current_alias = {"value": None}
    calls = []

    def fake_fetch(self, pdb_id, **kwargs):
        alias = current_alias["value"] or "direct"
        calls.append((alias, pdb_id))
        path = tmp_path / f"{pdb_id}_{alias}.cif"
        path.write_text(alias)
        return path

    monkeypatch.setattr(StructureLoader, "_fetch_from_rcsb", fake_fetch)

    loader = StructureLoader(name="alias_test_loader")

    for alias in ("pdb", "cif", "mmcif"):
        current_alias["value"] = alias
        result = loader.fetch_entity("7Y89", source=alias)
        assert result.read_text() == alias

    assert calls == [("pdb", "7Y89"), ("cif", "7Y89"), ("mmcif", "7Y89")]
