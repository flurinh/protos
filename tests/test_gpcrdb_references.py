"""Integrity tests for the bundled GPCRdb/CGN/CAN reference tables."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from protos.processing.grn.grn_utils import validate_grn_string
from protos.io.paths.path_config import ProtosPaths


BUNDLE = Path(__file__).parents[1] / "src" / "protos" / "reference_data" / "grn"


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
