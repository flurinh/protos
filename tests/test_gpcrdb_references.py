"""Integrity tests for the bundled GPCRdb/CGN/CAN reference tables."""

from __future__ import annotations

import json
import hashlib
import re
from pathlib import Path

import pandas as pd

from protos.processing.grn.grn_utils import validate_grn_string
from protos.io.paths.path_config import ProtosPaths


BUNDLE = Path(__file__).parents[1] / "src" / "protos" / "reference_data" / "grn"
MANIFEST = json.loads((BUNDLE / "manifest.json").read_text(encoding="utf-8"))
SUPPORTED_TABLES = {
    "can_arrestin_human.csv",
    "cgn_galpha_g1213_human.csv",
    "cgn_galpha_gio_human.csv",
    "cgn_galpha_gq11_human.csv",
    "cgn_galpha_gs_human.csv",
    "cgn_galpha_human.csv",
    "gpcrdb_class_a.csv",
    "gpcrdb_class_b1.csv",
    "gpcrdb_class_b2.csv",
    "gpcrdb_class_c.csv",
    "gpcrdb_class_d1.csv",
    "gpcrdb_class_f.csv",
    "gpcrdb_class_o1.csv",
    "gpcrdb_class_o2.csv",
    "gpcrdb_class_t2.csv",
    "gpcrdb_ref.csv",
    "gpcrdb_unclassified.csv",
    "type_I_opsins.csv",
    "type_II_opsins.csv",
}


def load_reference(name: str) -> pd.DataFrame:
    return pd.read_csv(BUNDLE / "reference" / name, index_col=0, dtype=str).fillna("-")


def test_signal_protein_dot_notation_is_valid() -> None:
    assert validate_grn_string("G.H5.26")[0]
    assert validate_grn_string("H.ha-hb.01")[0]
    assert validate_grn_string("N.S1.01")[0]
    assert validate_grn_string("C.s20c.12")[0]
    assert not validate_grn_string("G.H5.00")[0]


def test_human_cgn_table_and_family_splits_are_complete() -> None:
    all_human = load_reference("cgn_galpha_human.csv")
    splits = [
        load_reference("cgn_galpha_gs_human.csv"),
        load_reference("cgn_galpha_gio_human.csv"),
        load_reference("cgn_galpha_gq11_human.csv"),
        load_reference("cgn_galpha_g1213_human.csv"),
    ]
    assert all_human.shape == (16, 539)
    assert {row for split in splits for row in split.index} == set(all_human.index)
    assert sum(len(split) for split in splits) == len(all_human)
    assert all(validate_grn_string(column)[0] for column in all_human.columns)


def test_human_can_table_is_complete() -> None:
    table = load_reference("can_arrestin_human.csv")
    assert table.shape == (4, 437)
    assert set(table.index) == {
        "ARRB1_HUMAN",
        "ARRB2_HUMAN",
        "ARRC_HUMAN",
        "ARRS_HUMAN",
    }
    assert all(validate_grn_string(column)[0] for column in table.columns)


def test_receptor_class_tables_partition_the_export() -> None:
    provenance = json.loads((BUNDLE / "gpcrdb_provenance.json").read_text())
    counts = provenance["statistics"]["receptors"]["class_rows"]
    assert sum(counts.values()) == provenance["statistics"]["receptors"]["input_rows"] == 401
    assert counts["Class D1 (Ste2-like fungal pheromone)"] == 0
    assert provenance["source"]["runtime_api_used"] is False


def test_bundle_contains_only_supported_current_tables() -> None:
    present = {path.name for path in (BUNDLE / "reference").glob("*.csv")}
    assert present == SUPPORTED_TABLES


def test_gpcrdb_tables_exclude_upstream_loop_and_terminal_numbering() -> None:
    forbidden = re.compile(r"^(?:0|9)\.\d+$|^[1-9][1-9]\.\d+$")
    for path in (BUNDLE / "reference").glob("gpcrdb_*.csv"):
        table = load_reference(path.name)
        assert table.shape[1] == 424
        assert not any("x" in column or forbidden.fullmatch(column) for column in table.columns)


def test_type_ii_opsins_are_wild_type_only_and_have_protos_loop_policy() -> None:
    table = load_reference("type_II_opsins.csv")
    assert len(table) == 364
    assert "U57536.1" in table.index
    assert "NM_001014890" in table.index
    assert not any(re.fullmatch(r"O[0-9a-f]{9}", entity) for entity in table.index)
    assert not any(re.fullmatch(r"^[1-9][1-9]\.\d+$", column) for column in table.columns)
    provenance = json.loads((BUNDLE / "opsin_provenance.json").read_text())
    assert provenance["type_II"]["wt_rows"] == 364
    assert provenance["type_II"]["source_subset"].startswith("VPOD_1.3 WT")
    assert "Exact ungapped" in provenance["type_II"]["selection"]


def test_type_i_opsins_include_separate_tara_domains_and_pi_bulge() -> None:
    table = load_reference("type_I_opsins.csv")
    assert table.shape == (130, 1911)
    assert "7pl9" not in table.index
    assert {"TARA_A", "TARA_B"} <= set(table.index)
    assert table.loc["TARA_A", ["5.43", "5.44", "5.45", "5.451", "5.46"]].to_dict() == {
        "5.43": "-",
        "5.44": "Y187",
        "5.45": "A188",
        "5.451": "I189",
        "5.46": "F190",
    }
    assert table.loc["TARA_B", "5.451"] == "-"


def test_fresh_data_root_installs_and_annotates_without_network(
    tmp_path: Path, monkeypatch
) -> None:
    import protos.io.paths.path_config as path_config
    from protos.io.formats.fasta_utils import read_fasta
    from protos.processing.grn.grn_processor import GRNProcessor

    paths = ProtosPaths(str(tmp_path))
    paths.get_processor_path("grn")
    installed = tmp_path / "grn"
    assert (installed / "gpcrdb_provenance.json").is_file()
    assert (installed / "reference" / "gpcrdb_class_b2.csv").is_file()
    assert (installed / "reference" / "cgn_galpha_human.csv").is_file()
    assert (installed / "reference" / "can_arrestin_human.csv").is_file()
    assert not (installed / "reference" / "gpcr_a_core.csv").exists()
    assert not (installed / "reference" / "vpod1_2.csv").exists()

    monkeypatch.setattr(path_config, "_paths_instance", paths)
    monkeypatch.setattr(ProtosPaths, "_instance", paths)
    fixture = Path(__file__).parent / "fixtures" / "uniprot_grn_sequences.fasta"
    sequences = read_fasta(str(fixture))
    annotations, summary = GRNProcessor().annotate_sequences(
        {"P07550|ADRB2_HUMAN": sequences["P07550|ADRB2_HUMAN"]},
        reference_table="gpcrdb_class_a",
        protein_family="gpcrdb_class_a",
    )
    assert summary["per_sequence"]["P07550|ADRB2_HUMAN"]["status"] == "ok"
    assert (annotations.loc["P07550|ADRB2_HUMAN"] != "-").any()


def test_reference_install_repairs_damage_and_force_removes_unbundled_csvs(
    tmp_path: Path,
) -> None:
    paths = ProtosPaths(str(tmp_path))
    paths.get_processor_path("grn")
    reference = tmp_path / "grn" / "reference"
    damaged = reference / "gpcrdb_class_a.csv"
    damaged.write_text("damaged\n", encoding="utf-8")
    custom = reference / "custom.csv"
    custom.write_text("entity_name,custom.01\n", encoding="utf-8")

    paths._install_reference_data()  # pylint: disable=protected-access
    assert hashlib.sha256(damaged.read_bytes()).hexdigest() == MANIFEST["files"][damaged.name]
    assert custom.exists(), "automatic repair must preserve an unrelated custom table"

    paths._install_reference_data(force=True)  # pylint: disable=protected-access
    assert not custom.exists(), "explicit refresh is authoritative"
    marker = (tmp_path / ".protos_initialized").read_text(encoding="utf-8")
    assert f"bundle={MANIFEST['bundle_version']}" in marker


def test_reinitialize_honours_reference_install_switch(tmp_path: Path) -> None:
    paths = ProtosPaths(str(tmp_path))
    paths.get_processor_path("grn")
    reference = tmp_path / "grn" / "reference" / "gpcrdb_class_a.csv"
    reference.unlink()

    paths.reinitialize(reinstall_reference=False)
    assert not reference.exists()
    paths.reinitialize(reinstall_reference=True)
    assert reference.is_file()
