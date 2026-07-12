"""Live UniProt acquisition, registration, and GRN annotation tests.

Run with ``PROTOS_RUN_NETWORK_TESTS=1 pytest -m network``.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from protos.io.formats.fasta_utils import read_fasta


pytestmark = [
    pytest.mark.network,
    pytest.mark.skipif(
        os.environ.get("PROTOS_RUN_NETWORK_TESTS") != "1",
        reason="set PROTOS_RUN_NETWORK_TESTS=1 to query UniProt",
    ),
]


UNIPROT_CASES = [
    # GPCR top-level classes represented by the human export.
    ("P07550", "ADRB2_HUMAN", "gpcrdb_class_a", "beta2-adrenoceptor-Human"),
    ("P43220", "GLP1R_HUMAN", "gpcrdb_class_b1", "GLP-1-receptor-Human"),
    ("Q9Y653", "ADGRG1_HUMAN", "gpcrdb_class_b2", "ADGRG1-Human"),
    ("P41180", "CASR_HUMAN", "gpcrdb_class_c", "CaS-receptor-Human"),
    ("Q9UP38", "FZD1_HUMAN", "gpcrdb_class_f", "FZD1-Human"),
    ("Q8TCB6", "OR51E1_HUMAN", "gpcrdb_class_o1", "OR51E1-Human"),
    ("Q9P1Q5", "OR1A1_HUMAN", "gpcrdb_class_o2", "Olfactory-receptor-1A1-Human"),
    ("Q9NYV8", "TAS2R14_HUMAN", "gpcrdb_class_t2", "TAS2R14-Human"),
    ("P51810", "GPR143_HUMAN", "gpcrdb_unclassified", "GPR143-Human"),
    # Aggregate/compatibility GPCR tables.
    ("P07550", "ADRB2_ALL_HUMAN", "gpcrdb_ref", "beta2-adrenoceptor-Human"),
    ("P07550", "ADRB2_CORE_HUMAN", "gpcr_a_core", "beta2-adrenoceptor-Human"),
    # Every G-alpha family plus the aggregate table.
    ("P63092", "GNAS2_HUMAN", "cgn_galpha_gs_human", "GNAS2_HUMAN"),
    ("P63096", "GNAI1_HUMAN", "cgn_galpha_gio_human", "GNAI1_HUMAN"),
    ("P50148", "GNAQ_HUMAN", "cgn_galpha_gq11_human", "GNAQ_HUMAN"),
    ("Q03113", "GNA12_HUMAN", "cgn_galpha_g1213_human", "GNA12_HUMAN"),
    ("P63092", "GNAS2_ALL_HUMAN", "cgn_galpha_human", "GNAS2_HUMAN"),
    # Arrestin CAN.
    ("P49407", "ARRB1_HUMAN", "can_arrestin_human", "ARRB1_HUMAN"),
]

FROZEN_RECORDS = read_fasta(
    str(Path(__file__).parent / "fixtures" / "uniprot_grn_sequences.fasta")
)
FROZEN_BY_ACCESSION = {
    record_id.split("|", 1)[0]: sequence
    for record_id, sequence in FROZEN_RECORDS.items()
}


@pytest.fixture(scope="module")
def live_processors(tmp_path_factory):
    import protos
    import protos.io.paths.path_config as path_config
    from protos.io.ingest.sequence_loader import SequenceLoader
    from protos.processing.sequence import SequenceProcessor

    path_config._paths_instance = None
    path_config.ProtosPaths._instance = None
    data_root = tmp_path_factory.mktemp("protos-uniprot-grn")
    protos.set_data_path(str(data_root))
    sequence_processor = SequenceProcessor()
    loader = SequenceLoader(processor=sequence_processor)
    yield sequence_processor, loader
    path_config._paths_instance = None
    path_config.ProtosPaths._instance = None


@pytest.mark.parametrize("accession,entity,table,expected_reference", UNIPROT_CASES)
def test_uniprot_fetch_register_and_annotate_all_biological_tables(
    live_processors,
    accession: str,
    entity: str,
    table: str,
    expected_reference: str,
) -> None:
    sequence_processor, loader = live_processors
    saved = loader.download_and_register(
        accession,
        name=entity,
        materialize_entities=True,
        metadata={"test_accession": accession},
    )
    assert saved == entity
    sequence = sequence_processor.load_entity(entity)
    assert isinstance(sequence, str) and sequence
    assert sequence == FROZEN_BY_ACCESSION[accession]

    annotations, summary = sequence_processor.annotate_with_grn(
        sequences={entity: sequence},
        reference_table=table,
        protein_family=table,
        return_summary=True,
    )
    info = summary["per_sequence"][entity]
    assert info["status"] == "ok"
    assert info["reference"] == expected_reference
    assert info["coverage"] == 1.0
    assert info["deletion_residues"] == 0
    assert (annotations.loc[entity] != "-").any()


def test_fungal_class_d1_is_an_explicit_empty_table_exception(live_processors) -> None:
    sequence_processor, _ = live_processors
    with pytest.raises(ValueError, match="no usable reference sequences"):
        sequence_processor.annotate_with_grn(
            sequences={"STE2_YEAST": "M" * 431},
            reference_table="gpcrdb_class_d1",
            protein_family="gpcrdb_class_d1",
        )
