"""Structure-level GRN projection and visualization tests."""

from __future__ import annotations

import pandas as pd
import logging

from protos.processing.grn import GRNProcessor
from protos.processing.structure import StructureProcessor
from protos.visualization.grn_structure_vis import plot_grn_ca_structure


def test_chain_residue_keys_include_modified_polymer_residues() -> None:
    atoms = pd.DataFrame(
        [
            {
                "atom_id": 1,
                "group": "ATOM",
                "atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": 11,
                "label_seq_id": 2,
                "res_name3l": "ALA",
                "res_name1l": "A",
                "insertion": "",
                "model_num": 1,
            },
            {
                "atom_id": 2,
                "group": "HETATM",
                "atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": 10,
                "label_seq_id": 1,
                "res_name3l": "MSE",
                "res_name1l": "M",
                "insertion": "",
                "model_num": 1,
            },
        ]
    )
    sequence, keys = StructureProcessor._extract_chain_residue_keys(atoms, "A")
    assert sequence == "MA"
    assert keys == [
        {"label_seq_id": 1, "auth_seq_id": 10, "insertion": ""},
        {"label_seq_id": 2, "auth_seq_id": 11, "insertion": ""},
    ]
    restricted_sequence, restricted_keys = StructureProcessor._extract_chain_residue_keys(
        atoms, "A", auth_residue_range=(11, 11)
    )
    assert restricted_sequence == "A"
    assert restricted_keys == [
        {"label_seq_id": 2, "auth_seq_id": 11, "insertion": ""}
    ]


def test_plotly_ca_visualization_has_grn_hovertext_and_chain_edges() -> None:
    structure = pd.DataFrame(
        [
            {
                "structure_id": "demo",
                "atom_id": 1,
                "atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": 10,
                "label_seq_id": 1,
                "res_name3l": "ALA",
                "insertion": "",
                "model_num": 1,
                "grn": "G.HN.01",
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
            },
            {
                "structure_id": "demo",
                "atom_id": 2,
                "atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": 11,
                "label_seq_id": 2,
                "res_name3l": "GLY",
                "insertion": "",
                "model_num": 1,
                "grn": "G.HN.02",
                "x": 3.8,
                "y": 0.0,
                "z": 0.0,
            },
            {
                "structure_id": "demo",
                "atom_id": 3,
                "atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": 20,
                "label_seq_id": 10,
                "res_name3l": "SER",
                "insertion": "",
                "model_num": 1,
                "grn": "",
                "x": 20.0,
                "y": 0.0,
                "z": 0.0,
            },
        ]
    )
    figure = plot_grn_ca_structure(structure)
    assert len(figure.data) == 2
    edge_trace, atom_trace = figure.data
    assert list(edge_trace.x) == [0.0, 3.8, None]
    assert "GRN: G.HN.01" in atom_trace.text[0]
    assert "GRN: unassigned" in atom_trace.text[2]


def test_direct_structure_annotation_maps_modified_residue_atoms(monkeypatch) -> None:
    structure = pd.DataFrame(
        [
            {
                "structure_id": "demo",
                "atom_id": 1,
                "group": "HETATM",
                "atom_name": "CA",
                "auth_chain_id": "A",
                "auth_seq_id": 10,
                "label_seq_id": 1,
                "res_name3l": "MSE",
                "res_name1l": "M",
                "insertion": "",
                "model_num": 1,
                "x": 0.0,
                "y": 0.0,
                "z": 0.0,
            },
            {
                "structure_id": "demo",
                "atom_id": 2,
                "group": "HETATM",
                "atom_name": "SE",
                "auth_chain_id": "A",
                "auth_seq_id": 10,
                "label_seq_id": 1,
                "res_name3l": "MSE",
                "res_name1l": "M",
                "insertion": "",
                "model_num": 1,
                "x": 1.0,
                "y": 0.0,
                "z": 0.0,
            },
            *[
                {
                    "structure_id": "demo",
                    "atom_id": index + 3,
                    "group": "ATOM",
                    "atom_name": "CA",
                    "auth_chain_id": "A",
                    "auth_seq_id": index + 11,
                    "label_seq_id": index + 2,
                    "res_name3l": "ALA",
                    "res_name1l": "A",
                    "insertion": "",
                    "model_num": 1,
                    "x": float(index + 2),
                    "y": 0.0,
                    "z": 0.0,
                }
                for index in range(9)
            ],
        ]
    )
    processor = object.__new__(StructureProcessor)
    processor.frames = {}
    processor._dirty = False
    processor.logger = logging.getLogger("test_structure_grn")
    processor.load_entity = lambda _structure_id: structure

    def fake_annotation(_self, sequences, **_kwargs):
        name = next(iter(sequences))
        columns = [f"S.{position:02d}" for position in range(1, 11)]
        values = [f"M1"] + [f"A{position}" for position in range(2, 11)]
        return (
            pd.DataFrame([values], index=[name], columns=columns),
            {"per_sequence": {name: {"status": "ok", "coverage": 1.0}}},
        )

    monkeypatch.setattr(GRNProcessor, "__init__", lambda self: None)
    monkeypatch.setattr(GRNProcessor, "annotate_sequences", fake_annotation)
    annotated, summary = processor.annotate_with_grn(
        "demo",
        reference_table="memory",
        protein_family="test",
        chains=["A"],
        save=False,
        return_summary=True,
    )
    result = annotated.reset_index()
    mse = result[result["label_seq_id"] == 1]
    assert set(mse["grn"]) == {"S.01"}
    assert summary["chains"]["A"]["status"] == "ok"
