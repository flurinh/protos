"""Behavioural tests for segment-based, gap-aware GRN annotation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from protos.analysis.sequence.alignment_engine import (
    SequenceAlignmentEngine,
    SequenceAlignmentResult,
)
from protos.io.formats.fasta_utils import read_fasta
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.assign_grns import assign_grns_to_sequences
from protos.processing.grn.grn_utils import (
    make_directional_flexible_grn,
    make_insertion_grn,
    normalize_grn_format,
    parse_directional_flexible_grn,
    parse_grn_str2float,
    sort_grns_str,
    validate_grn_string,
)


BUNDLE = Path(__file__).parents[1] / "src" / "protos" / "reference_data" / "grn"
FIXTURES = Path(__file__).parent / "fixtures"
MANIFEST = json.loads((BUNDLE / "manifest.json").read_text())
REFERENCE_TABLES = sorted(MANIFEST["files"])
FROZEN_UNIPROT = read_fasta(str(FIXTURES / "uniprot_grn_sequences.fasta"))

OFFLINE_UNIPROT_CASES = [
    ("P07550|ADRB2_HUMAN", "gpcrdb_class_a", "beta2-adrenoceptor-Human"),
    ("P43220|GLP1R_HUMAN", "gpcrdb_class_b1", "GLP-1-receptor-Human"),
    ("Q9Y653|ADGRG1_HUMAN", "gpcrdb_class_b2", "ADGRG1-Human"),
    ("P41180|CASR_HUMAN", "gpcrdb_class_c", "CaS-receptor-Human"),
    ("Q9UP38|FZD1_HUMAN", "gpcrdb_class_f", "FZD1-Human"),
    ("Q8TCB6|OR51E1_HUMAN", "gpcrdb_class_o1", "OR51E1-Human"),
    ("Q9P1Q5|OR1A1_HUMAN", "gpcrdb_class_o2", "Olfactory-receptor-1A1-Human"),
    ("Q9NYV8|TAS2R14_HUMAN", "gpcrdb_class_t2", "TAS2R14-Human"),
    ("P51810|GPR143_HUMAN", "gpcrdb_unclassified", "GPR143-Human"),
    ("P07550|ADRB2_HUMAN", "gpcrdb_ref", "beta2-adrenoceptor-Human"),
    ("P63092|GNAS2_HUMAN", "cgn_galpha_gs_human", "GNAS2_HUMAN"),
    ("P63096|GNAI1_HUMAN", "cgn_galpha_gio_human", "GNAI1_HUMAN"),
    ("P50148|GNAQ_HUMAN", "cgn_galpha_gq11_human", "GNAQ_HUMAN"),
    ("Q03113|GNA12_HUMAN", "cgn_galpha_g1213_human", "GNA12_HUMAN"),
    ("P63092|GNAS2_HUMAN", "cgn_galpha_human", "GNAS2_HUMAN"),
    ("P49407|ARRB1_HUMAN", "can_arrestin_human", "ARRB1_HUMAN"),
]


def projection(
    query: str,
    reference: str,
    labels: list[str],
    columns: list[str] | None = None,
    *,
    assign_insertions: bool = False,
    standard_regions: dict[str, tuple[str, str]] | None = None,
    strict_regions: dict[str, tuple[str, str]] | None = None,
    max_alignment_gap: int = 1,
):
    result = SequenceAlignmentResult(
        seq1_id="query",
        seq2_id="reference",
        score=0.0,
        alignment=[query, "".join("|" if q == r else "-" for q, r in zip(query, reference)), reference],
        method="test",
    )
    processor = object.__new__(GRNProcessor)
    return processor._project_alignment(
        "query",
        result,
        labels,
        columns or labels,
        sequence_grn_order=columns or labels,
        assign_unambiguous_insertions=assign_insertions,
        standard_regions=standard_regions,
        strict_regions=strict_regions,
        max_alignment_gap=max_alignment_gap,
    )


def test_manifest_covers_and_authenticates_every_bundled_reference_table() -> None:
    reference_paths = sorted((BUNDLE / "reference").glob("*.csv"))
    assert REFERENCE_TABLES == [path.name for path in reference_paths]
    for path in reference_paths:
        assert hashlib.sha256(path.read_bytes()).hexdigest() == MANIFEST["files"][path.name]


@pytest.mark.parametrize(
    "label",
    [
        "1.50",
        "12.003",
        "12x47",
        "1.411",
        "5x461",
        "n.09",
        "c.01",
        "TM1.50",
        "G.HN.03",
        "H.ha-hb.01",
        "A.H1.50",
        "arrestin.N.S1.01",
    ],
)
def test_segment_identifiers_can_be_numeric_or_arbitrary_strings(label: str) -> None:
    assert validate_grn_string(label)[0]


@pytest.mark.parametrize(
    "label",
    ["G..03", ".HN.03", "G.HN.", "G.HN.zero", "G.HN.00", "segment"],
)
def test_malformed_segment_labels_are_rejected(label: str) -> None:
    assert not validate_grn_string(label)[0]


def test_numeric_normalization_and_float_behaviour_are_unchanged() -> None:
    assert normalize_grn_format("1x50") == "1.50"
    assert normalize_grn_format("12x005") == "12.005"
    assert parse_grn_str2float("1.50") == 1.5
    assert parse_grn_str2float("1x50") == 1.5
    assert parse_grn_str2float("12x001") == 12.001
    assert sort_grns_str(["2.50", "12.001", "1.50"]) == [
        "1.50",
        "12.001",
        "2.50",
    ]
    processor = object.__new__(GRNProcessor)
    assert processor._normalize_columns(["1x50", "1x411", "12x001"]) == [
        "1.50",
        "1.411",
        "12.001",
    ]


def test_directional_flexible_coordinates_are_class_agnostic() -> None:
    assert parse_directional_flexible_grn("12.003") == (1, 2, 3)
    assert parse_directional_flexible_grn("21.003") == (2, 1, 3)
    assert parse_directional_flexible_grn("19.047") == (1, 9, 47)
    assert make_directional_flexible_grn(2, 1, 3) == "21.003"
    assert validate_grn_string("19.047") == (
        True,
        "Valid directional flexible-region GRN format",
    )


def test_insertion_coordinate_uses_final_position_for_any_segment_identifier() -> None:
    assert make_insertion_grn("1.41", 1) == "1.411"
    assert make_insertion_grn("G.HN.03", 1) == "G.HN.031"
    assert make_insertion_grn("arrestin.N.S1.01", 2) == "arrestin.N.S1.012"
    with pytest.raises(ValueError, match="indices 1 through 9"):
        make_insertion_grn("G.HN.03", 10)


def test_hierarchical_sort_uses_reference_segment_order_not_lexical_order() -> None:
    labels = ["Z.loop.02", "Z.loop.01", "A.helix.02", "A.helix.01"]
    assert sort_grns_str(labels) == [
        "Z.loop.01",
        "Z.loop.02",
        "A.helix.01",
        "A.helix.02",
    ]
    assert sort_grns_str(labels, segment_order=["A.helix", "Z.loop"]) == [
        "A.helix.01",
        "A.helix.02",
        "Z.loop.01",
        "Z.loop.02",
    ]


def test_reference_order_is_inferred_from_sequence_positions() -> None:
    table = pd.DataFrame(
        [["B2", "A1", "C3", "D4"]],
        index=["ref"],
        columns=["Z.segment.02", "A.first.01", "Z.segment.03", "loop.named.01"],
    )
    assert GRNProcessor._infer_reference_grn_order(table) == [
        "A.first.01",
        "Z.segment.02",
        "Z.segment.03",
        "loop.named.01",
    ]


def test_directional_flexible_columns_sort_between_intrinsically_ordered_segments() -> None:
    table = pd.DataFrame(
        [["A1", "D4", "C3", "B2", "E5", "F6"]],
        index=["ref"],
        columns=["G.HN.01", "21.001", "12.002", "12.001", "G.S1.01", "G.S1.02"],
    )
    assert GRNProcessor._infer_reference_grn_order(table) == [
        "G.HN.01",
        "12.001",
        "12.002",
        "21.001",
        "G.S1.01",
        "G.S1.02",
    ]


def test_empty_columns_inherit_segment_direction_and_gpcr_insertions_sort_locally() -> None:
    reverse_segment = pd.DataFrame(
        [["C3", "-", "A1"]],
        index=["ref"],
        columns=["named.01", "named.02", "named.03"],
    )
    assert GRNProcessor._infer_reference_grn_order(reverse_segment) == [
        "named.03",
        "named.02",
        "named.01",
    ]

    gpcr_insert = pd.DataFrame(
        [["A1", "B2", "C3"]],
        index=["ref"],
        columns=["1.41", "1.411", "1.42"],
    )
    assert GRNProcessor._infer_reference_grn_order(gpcr_insert) == [
        "1.41",
        "1.411",
        "1.42",
    ]


def test_reference_sequence_reconstruction_uses_cell_positions_not_columns() -> None:
    table = pd.DataFrame(
        [["C3", "A1", "B2"]],
        index=["ref"],
        columns=["segment.03", "segment.01", "segment.02"],
    )
    processor = object.__new__(GRNProcessor)
    sequences, mappings = processor._prepare_reference_sequences(table)
    assert sequences == {"ref": "ABC"}
    assert mappings == {"ref": ["segment.01", "segment.02", "segment.03"]}


def test_reference_sequence_reconstruction_rejects_duplicate_positions() -> None:
    table = pd.DataFrame(
        [["A1", "B1"]], index=["broken"], columns=["one.01", "two.01"]
    )
    with pytest.raises(ValueError, match="duplicate sequence positions"):
        GRNProcessor._prepare_reference_sequences(object.__new__(GRNProcessor), table)


def test_internal_insertion_is_detected_and_downstream_positions_are_corrected() -> None:
    row, coverage, assigned, indels = projection(
        "ACXDE",
        "AC-DE",
        ["S.01", "S.02", "S.03", "S.04"],
    )
    assert row.to_dict() == {"S.01": "A1", "S.02": "C2", "S.03": "D4", "S.04": "E5"}
    assert coverage == 1.0
    assert assigned == 4
    assert indels["insertion_residues"] == 1
    assert indels["insertions"][0] == {
        "after_grn": "S.02",
        "before_grn": "S.03",
        "query_start": 3,
        "query_end": 3,
        "residues": "X",
        "candidate_grns": [],
        "assignment": "unassigned",
        "generated_grns": [],
        "kind": "internal_insertion",
        "unassigned_reason": "insertion_annotation_disabled",
    }


def test_deletion_is_detected_and_deleted_grn_remains_empty() -> None:
    row, coverage, assigned, indels = projection(
        "AC-DE",
        "ACXDE",
        ["S.01", "S.02", "loop.01", "S.03", "S.04"],
    )
    assert row["loop.01"] == "-"
    assert row["S.03"] == "D3"
    assert row["S.04"] == "E4"
    assert coverage == pytest.approx(4 / 5)
    assert assigned == 4
    assert indels["deletion_residues"] == 1
    assert indels["deletions"] == [
        {
            "grns": ["loop.01"],
            "reference_residues": "X",
            "after_query_position": 2,
        }
    ]


def test_unambiguous_insertion_can_fill_an_existing_empty_reference_column() -> None:
    columns = ["S.01", "S.02", "loop.named.01", "S.03", "S.04"]
    row, coverage, assigned, indels = projection(
        "ACXDE",
        "AC-DE",
        ["S.01", "S.02", "S.03", "S.04"],
        columns,
        assign_insertions=True,
    )
    assert row["loop.named.01"] == "X3"
    assert row["S.03"] == "D4"
    assert coverage == 1.0
    assert assigned == 5
    assert indels["insertions"][0]["candidate_grns"] == ["loop.named.01"]
    assert indels["insertions"][0]["assignment"] == "assigned_exact_candidate_count"


def test_insertion_inside_string_segment_creates_new_decimal_columns() -> None:
    row, coverage, assigned, indels = projection(
        "ACXDE",
        "AC-DE",
        ["G.HN.01", "G.HN.02", "G.HN.03", "G.HN.04"],
        assign_insertions=True,
    )
    assert row["G.HN.021"] == "X3"
    assert row["G.HN.03"] == "D4"
    assert coverage == 1.0
    assert assigned == 5
    event = indels["insertions"][0]
    assert event["generated_grns"] == ["G.HN.021"]
    assert event["assignment"] == "assigned_generated_coordinates"


def test_insertion_after_existing_bulge_advances_without_collision() -> None:
    row, _, _, indels = projection(
        "ACXDE",
        "AC-DE",
        ["1.41", "1.411", "1.42", "1.43"],
        assign_insertions=True,
    )
    assert row["1.411"] == "C2"
    assert row["1.412"] == "X3"
    assert indels["insertions"][0]["generated_grns"] == ["1.412"]


def test_generated_flexible_columns_merge_independently_of_sequence_order() -> None:
    merged = ["12.001", "21.001"]
    GRNProcessor._merge_ordered_labels(
        merged,
        ["12.001", "12.002", "12.003", "21.002", "21.001"],
    )
    assert merged == ["12.001", "12.002", "12.003", "21.002", "21.001"]


def test_insertion_beyond_compact_coordinate_capacity_is_explicitly_unassigned() -> None:
    row, _, assigned, indels = projection(
        "ACXXXXXXXXXXDE",
        "AC----------DE",
        ["S.01", "S.02", "S.03", "S.04"],
        assign_insertions=True,
    )
    assert assigned == 4
    assert all(not column.startswith("S.02") or column == "S.02" for column in row.index)
    event = indels["insertions"][0]
    assert event["assignment"] == "unassigned"
    assert event["unassigned_reason"] == "coordinate_capacity_or_segment_mapping"


@pytest.mark.parametrize(
    ("residues", "expected"),
    [
        ("XY", [("12.001", "X3"), ("21.001", "Y4")]),
        (
            "XYZ",
            [("12.001", "X3"), ("12.002", "Y4"), ("21.001", "Z5")],
        ),
    ],
)
def test_inter_segment_residues_are_labelled_from_both_directions(
    residues: str,
    expected: list[tuple[str, str]],
) -> None:
    row, _, assigned, indels = projection(
        f"AC{residues}DE",
        f"AC{'-' * len(residues)}DE",
        ["G.HN.01", "G.HN.02", "G.S1.01", "G.S1.02"],
        assign_insertions=True,
    )
    for grn, value in expected:
        assert row[grn] == value
    assert assigned == 4 + len(residues)
    assert indels["insertions"][0]["generated_grns"] == [grn for grn, _ in expected]


def test_long_gap_inside_strict_region_is_compressed_across_standard_extent() -> None:
    standard = {"core": ("S.01", "S.04")}
    strict = {"core": ("S.02", "S.03")}
    row, coverage, assigned, indels = projection(
        "ACXYDE",
        "AC--DE",
        ["S.01", "S.02", "S.03", "S.04"],
        assign_insertions=True,
        standard_regions=standard,
        strict_regions=strict,
    )
    assert row.to_dict() == {
        "S.01": "A1",
        "S.02": "C2",
        "S.03": "X3",
        "S.04": "Y4",
    }
    assert coverage == 1.0
    assert assigned == 4
    event = indels["insertions"][0]
    assert event["kind"] == "alignment_gap_compressed"
    assert event["assignment"] == "compressed_to_standard_coordinates"
    assert event["corrections"] == [
        {"grn": "S.03", "from": "D5", "to": "X3"},
        {"grn": "S.04", "from": "E6", "to": "Y4"},
    ]
    assert indels["insertion_residues"] == 0


def test_long_gap_without_strict_anchors_remains_a_biological_insertion() -> None:
    row, _, assigned, indels = projection(
        "ACXYDE",
        "AC--DE",
        ["S.01", "S.02", "S.03", "S.04"],
        assign_insertions=True,
    )
    assert row["S.021"] == "X3"
    assert row["S.022"] == "Y4"
    assert row["S.03"] == "D5"
    assert assigned == 6
    assert indels["insertions"][0]["kind"] == "internal_insertion"


def test_single_residue_gap_is_not_compressed_even_inside_strict_region() -> None:
    regions = {"core": ("S.01", "S.04")}
    row, _, _, indels = projection(
        "ACXDE",
        "AC-DE",
        ["S.01", "S.02", "S.03", "S.04"],
        assign_insertions=True,
        standard_regions=regions,
        strict_regions=regions,
    )
    assert row["S.021"] == "X3"
    assert row["S.03"] == "D4"
    assert indels["insertions"][0]["kind"] == "internal_insertion"


def test_non_exact_candidate_set_is_preserved_while_new_insertion_column_is_created() -> None:
    columns = ["S.01", "S.02", "loop.01", "loop.02", "S.03", "S.04"]
    row, _, _, indels = projection(
        "ACXDE",
        "AC-DE",
        ["S.01", "S.02", "S.03", "S.04"],
        columns,
        assign_insertions=True,
    )
    assert row["loop.01"] == row["loop.02"] == "-"
    assert row["S.021"] == "X3"
    assert indels["insertions"][0]["candidate_grns"] == ["loop.01", "loop.02"]
    assert indels["insertions"][0]["assignment"] == "assigned_generated_coordinates"


def test_terminal_overhang_uses_protos_distance_coordinates() -> None:
    row, coverage, assigned, indels = projection(
        "XXACDEYY", "--ACDE--", ["S.01", "S.02", "S.03", "S.04"]
    )
    assert row[["n.2", "n.1", "c.1", "c.2"]].tolist() == ["X1", "X2", "Y7", "Y8"]
    assert assigned == 8
    assert coverage == 1.0
    assert indels["insertion_residues"] == 0
    assert indels["terminal_overhang_residues"] == 4
    assert [event["kind"] for event in indels["insertions"]] == [
        "n_terminal_overhang",
        "c_terminal_overhang",
    ]
    assert all(
        event["assignment"] == "assigned_terminal_coordinates"
        for event in indels["insertions"]
    )


def test_gap_parameters_are_actually_applied_to_biopython() -> None:
    from protos.processing.sequence.seq_alignment import get_end_gap_scores

    engine = SequenceAlignmentEngine(
        open_gap_score=-7.0, extend_gap_score=-0.25, end_gap_score=-1.5
    )
    assert engine.gap_parameters == {
        "open_gap_score": -7.0,
        "extend_gap_score": -0.25,
        "end_gap_score": -1.5,
    }
    internal_open = (
        engine._aligner.open_internal_gap_score
        if hasattr(engine._aligner, "open_internal_gap_score")
        else engine._aligner.internal_open_gap_score
    )
    internal_extend = (
        engine._aligner.extend_internal_gap_score
        if hasattr(engine._aligner, "extend_internal_gap_score")
        else engine._aligner.internal_extend_gap_score
    )
    assert internal_open == -7.0
    assert internal_extend == -0.25
    assert get_end_gap_scores(engine._aligner) == (-1.5, -1.5)


def test_legacy_expansion_entry_point_bypasses_numeric_math_for_string_segments() -> None:
    reference = pd.DataFrame(
        [["A1", "C2", "D3", "E4"]],
        index=["segment_reference"],
        columns=["G.HN.01", "G.HN.02", "G.HN.03", "G.HN.04"],
    )
    annotations, metadata = assign_grns_to_sequences(
        {"query": "ACDE"}, reference, protein_family="cgn", use_mmseqs=False
    )
    assert annotations.loc["query"].tolist() == ["A1", "C2", "D3", "E4"]
    assert metadata["query"]["expansion_method"] == "segment_projection"


def test_legacy_expansion_does_not_float_convert_gpcr_bulge_labels() -> None:
    reference = pd.DataFrame(
        [["A1", "C2", "D3"]],
        index=["bulge_reference"],
        columns=["1.41", "1.411", "1.42"],
    )
    annotations, metadata = assign_grns_to_sequences(
        {"query": "ACD"}, reference, protein_family="gpcr_a", use_mmseqs=False
    )
    assert annotations.loc["query"].tolist() == ["A1", "C2", "D3"]
    assert metadata["query"]["expansion_method"] == "segment_projection"


def public_processor_with_table(table: pd.DataFrame) -> GRNProcessor:
    processor = object.__new__(GRNProcessor)
    processor.load_reference_table = lambda _name: table
    return processor


def test_public_annotation_reports_invalid_and_empty_sequences_without_crashing() -> None:
    reference = pd.DataFrame(
        [["A1", "C2"]], index=["ref"], columns=["segment.01", "segment.02"]
    )
    annotations, summary = public_processor_with_table(reference).annotate_sequences(
        {"empty": "", "invalid": "AC?"},
        reference_table="memory",
        protein_family="arbitrary",
    )
    assert summary["per_sequence"]["empty"]["status"] == "empty_sequence"
    assert summary["per_sequence"]["invalid"]["status"] == "invalid_sequence"
    assert (annotations == "-").all(axis=None)


def test_public_annotation_threshold_rejects_result_without_losing_diagnostics() -> None:
    reference = pd.DataFrame(
        [["A1", "C2"]], index=["ref"], columns=["segment.01", "segment.02"]
    )
    annotations, summary = public_processor_with_table(reference).annotate_sequences(
        {"query": "AC"},
        reference_table="memory",
        protein_family="arbitrary",
        min_coverage=1.01,
    )
    info = summary["per_sequence"]["query"]
    assert info["status"] == "below_threshold"
    assert info["coverage"] == 1.0
    assert info["assigned_positions"] == 2
    assert summary["global"]["annotated"] == 0
    assert (annotations.loc["query"] == "-").all()


def test_public_annotation_adds_generated_insertion_column_in_sequence_order() -> None:
    reference = pd.DataFrame(
        [["A1", "C2", "D3", "E4"]],
        index=["ref"],
        columns=["G.HN.01", "G.HN.02", "G.HN.03", "G.HN.04"],
    )
    annotations, summary = public_processor_with_table(reference).annotate_sequences(
        {"query": "ACXDE"},
        reference_table="memory",
        protein_family="arbitrary",
    )
    assert annotations.columns.tolist() == [
        "G.HN.01",
        "G.HN.02",
        "G.HN.021",
        "G.HN.03",
        "G.HN.04",
    ]
    assert annotations.loc["query", "G.HN.021"] == "X3"
    assert summary["per_sequence"]["query"]["assigned_positions"] == 5


def test_public_annotation_rejects_unknown_search_and_empty_reference() -> None:
    reference = pd.DataFrame(
        [["A1"]], index=["ref"], columns=["segment.01"]
    )
    with pytest.raises(ValueError, match="Unsupported GRN reference search"):
        public_processor_with_table(reference).annotate_sequences(
            {"query": "A"},
            reference_table="memory",
            protein_family="arbitrary",
            search="magic",
        )
    with pytest.raises(ValueError, match="no usable reference sequences"):
        public_processor_with_table(reference.iloc[:0]).annotate_sequences(
            {"query": "A"},
            reference_table="memory",
            protein_family="arbitrary",
        )


@pytest.mark.parametrize("fixture_id,table_name,expected_reference", OFFLINE_UNIPROT_CASES)
def test_frozen_full_length_uniprot_member_annotates_each_biological_table(
    fixture_id: str,
    table_name: str,
    expected_reference: str,
) -> None:
    table = pd.read_csv(
        BUNDLE / "reference" / f"{table_name}.csv", index_col=0, dtype=str
    ).fillna("-")
    annotations, summary = public_processor_with_table(table).annotate_sequences(
        {fixture_id: FROZEN_UNIPROT[fixture_id]},
        reference_table=table_name,
        protein_family=table_name,
    )
    info = summary["per_sequence"][fixture_id]
    assert info["status"] == "ok"
    assert info["reference"] == expected_reference
    assert info["coverage"] == 1.0
    assert info["deletion_residues"] == 0
    assert (annotations.loc[fixture_id] != "-").any()


@pytest.mark.parametrize("table_name", REFERENCE_TABLES)
def test_a_real_member_round_trips_through_every_bundled_table(table_name: str) -> None:
    table = pd.read_csv(BUNDLE / "reference" / table_name, index_col=0, dtype=str).fillna("-")
    assert all(
        validate_grn_string(str(column))[0] for column in table.columns
    ), f"{table_name} contains a GRN rejected by the validator"
    if table.empty:
        assert table_name == "gpcrdb_class_d1.csv"
        return

    processor = object.__new__(GRNProcessor)
    sequences, mappings = processor._prepare_reference_sequences(table)
    assert sequences, f"{table_name} has no usable member sequences"
    member = max(sequences, key=lambda name: len(sequences[name]))
    engine = SequenceAlignmentEngine()
    best_ref, alignment, score = processor._find_best_reference(
        member, sequences[member], sequences, engine
    )
    assert alignment is not None
    assert best_ref is not None
    assert score > 0

    row, coverage, assigned, indels = processor._project_alignment(
        member,
        alignment,
        mappings[best_ref],
        [str(column) for column in table.columns],
        sequence_grn_order=processor._infer_reference_grn_order(table),
    )
    assert assigned == len(mappings[best_ref])
    assert coverage == 1.0
    assert indels["insertion_residues"] == 0
    assert indels["deletion_residues"] == 0
    assert sum(value != "-" for value in row) == assigned
