"""Smoke tests for ModelManager configurable adapters and providers.

Run with: pytest -q protos/test_models.py
"""

from __future__ import annotations

from pathlib import Path
import shutil
import sys
from typing import Dict

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager
from protos.io.formats.fasta_utils import read_fasta
from protos.models.model_specs import ArtifactSpec, ExecutionSpec, ModelCard
from protos.processing.sequence import SequenceProcessor
from protos.processing.graph.graph_processor import GraphProcessor
from protos.processing.structure import StructureProcessor
from protos.processing.property import PropertyProcessor


def ensure_data_root() -> Path:
    # Use the repository-managed Protos data directory
    data_root = Path(__file__).resolve().parent / "data"
    repo_data = Path(__file__).resolve().parent / "data"
    # If a top-level protos/data exists, prefer it
    top_level_data = Path(__file__).resolve().parents[1] / "data"
    if top_level_data.exists():
        data_root = top_level_data
    data_root.mkdir(parents=True, exist_ok=True)
    try:
        protos.set_data_path(str(data_root))
    except Exception:
        # Already initialized (possibly by earlier tests); leave as-is
        pass
    return data_root


def make_sequence_dataset(name: str = "test_sequences") -> str:
    seq_proc = SequenceProcessor()
    seqs: Dict[str, str] = {
        "SEQ_A": "MKTIIALSYIFCLVFADYKDDDDK",
        "SEQ_B": "MGHHHHHHSSGLVPRGSHMASMTGGQQMGRGSEF",
    }
    seq_proc.save_sequences(
        seqs,
        output_file=name,
        dataset_name=name,
        metadata={"kind": "test"},
        materialize_entities=True,
    )
    return name


def make_graph_entity(structure_id: str = "3sn6") -> str:
    # Ensure the structure PKL is auto-registered if present in cache
    sp = StructureProcessor()
    _ = sp.load_entity(structure_id)  # triggers auto-register when cache exists
    gp = GraphProcessor()
    entity = gp.generate_graph(structure_id, level="atom", cutoff=5.0)
    return entity


def test_configurable_external_prepares_job(tmp_path):
    ensure_data_root()
    dataset = make_sequence_dataset("mm_seq_ext")

    manager = ModelManager()
    card = ModelCard(
        name="external_test",
        version="0.1",
        description="Configurable external adapter test",
        execution=ExecutionSpec(mode="external_config", entrypoint="echo external_test", expected_config="json"),
        input_spec=[
            ArtifactSpec(name="sequence_dataset", kind="sequence", provider="sequence_dataset", format="fasta")
        ],
        output_spec=[],
    )
    manager.register_model(card, None)

    inv = manager.prepare(
        "external_test",
        inputs={"sequence_dataset": dataset},
        config={"foo": "bar"},
        metadata={"source": "test"},
    )
    assert inv.job is not None
    job = inv.job
    assert job.working_dir.exists()
    ctx = job.metadata.get("context", {})
    cfg = job.metadata.get("config_path")
    assert cfg is not None and Path(cfg).exists()
    packaged = job.metadata.get("packaged_inputs", {})
    assert "sequence_dataset" in packaged
    assert Path(packaged["sequence_dataset"]).exists()


def test_configurable_runtime_executes_and_ingests(tmp_path):
    ensure_data_root()
    dataset = make_sequence_dataset("mm_seq_rt")

    manager = ModelManager()
    card = ModelCard(
        name="dummy_runtime",
        version="0.0",
        description="Configurable runtime adapter test",
        execution=ExecutionSpec(mode="runtime", entrypoint="protos.models.dummy_runtime:run"),
        input_spec=[
            ArtifactSpec(name="sequence_dataset", kind="sequence", provider="sequence_dataset", format="fasta")
        ],
        output_spec=[
            ArtifactSpec(name="property_table", kind="property", provider="property_processor", format="csv", optional=True)
        ],
    )
    manager.register_model(card, None)

    inv = manager.prepare("dummy_runtime", inputs={"sequence_dataset": dataset})
    assert inv.runtime is not None
    # Work/outputs dirs present via metadata
    meta = inv.runtime.metadata
    assert meta.get("work_dir") and Path(meta["work_dir"]).exists()
    assert meta.get("outputs_dir") and Path(meta["outputs_dir"]).exists()

    # Ingest property table
    summary = manager.ingest_outputs(inv)
    tables = [e for e in summary.get("ingested", []) if e.get("type") == "property_table"]
    assert tables, "Expected a property table to be ingested"

    # Confirm table persisted
    table_name = tables[0]["name"]
    prop = PropertyProcessor()
    df = prop.load_property_table(table_name)
    assert not df.empty


def test_graph_entity_provider(tmp_path):
    ensure_data_root()
    graph_entity = make_graph_entity("3sn6")

    manager = ModelManager()
    card = ModelCard(
        name="graph_external",
        version="0.1",
        description="External with graph input",
        execution=ExecutionSpec(mode="external_config", entrypoint="echo graph_external", expected_config="yaml"),
        input_spec=[
            ArtifactSpec(name="pocket_graph", kind="graph", provider="graph_entity", format="pickle")
        ],
        output_spec=[],
    )
    manager.register_model(card, None)

    inv = manager.prepare("graph_external", inputs={"pocket_graph": graph_entity})
    assert inv.job is not None
    packaged = inv.job.metadata.get("packaged_inputs", {})
    assert "pocket_graph" in packaged
    assert Path(packaged["pocket_graph"]).exists()


def test_pocketdta_prepare_job_with_repo_dataset():
    ensure_data_root()
    dataset_dir = Path("protos/src/protos/models/PocketDTA/dataset")
    config_yaml = Path("protos/src/protos/models/PocketDTA/configs/default_config.yaml")
    if not dataset_dir.exists() or not config_yaml.exists():
        return
    manager = ModelManager()
    inv = manager.prepare(
        "pocketdta",
        inputs={"dataset_dir": str(dataset_dir), "config_yaml": str(config_yaml)},
        config={"task": "Davis"},
    )
    assert inv.job is not None
    meta = inv.job.metadata or {}
    assert meta.get("context")
    cfg_path = meta.get("context", {}).get("config_path")
    assert cfg_path and Path(cfg_path).exists()


def test_graphscoredta_prepare_job():
    ensure_data_root()
    repo_root = Path("protos/src/protos/models/GraphscoreDTA")
    graphs = repo_root / "test_set" / "all_in_one_graph_test2016.bin"
    labels = repo_root / "test_set" / "labels_test2016.csv"
    vina = repo_root / "test_set" / "Vina_terms2016.pkl"
    model = repo_root / "model" / "modelp.pth"
    if not (graphs.exists() and labels.exists() and vina.exists() and model.exists()):
        return

    manager = ModelManager()
    inv = manager.prepare(
        "graphscoredta",
        inputs={
            "graphs_bin": str(graphs),
            "labels_csv": str(labels),
            "vina_pkl": str(vina),
            "model_path": str(model),
        },
    )
    assert inv.job is not None
    meta = inv.job.metadata or {}
    assert meta.get("context")
    manifest = (meta.get("staged") or {})
    assert manifest.get("test_set") and manifest.get("model")
def test_boltz_yaml_with_ligand(tmp_path):
    ensure_data_root()
    # Load a real GPCR sequence from packaged Boltz demo FASTA
    fasta_path = Path("data/models/boltz2/ADRB2_HUMAN_wild_type/sequences.fasta")
    assert fasta_path.exists(), "Expected demo FASTA at data/models/boltz2/..."
    seqs = read_fasta(str(fasta_path))
    assert seqs, "No sequences found in demo FASTA"
    seq_proc = SequenceProcessor()
    dataset_name = "adrb2_boltz_demo"
    seq_proc.save_sequences(seqs, output_file=dataset_name, dataset_name=dataset_name, materialize_entities=True)

    manager = ModelManager()
    first_entity = next(iter(seqs.keys()))
    # Prepare Boltz job with a real ligand example (dopamine)
    inv = manager.prepare(
        "boltz2",
        inputs={"sequence_dataset": dataset_name, "entity": first_entity},
        config={
            "ligand": {"id": "DOP", "smiles": "CNCCC1=CC(=C(C=C1)O)O"},
        },
    )
    assert inv.job is not None
    yaml_path = None
    for b in inv.job.artifacts:
        if b.spec.name == "boltz_config":
            yaml_path = b.path
            break
    assert yaml_path is not None and Path(yaml_path).exists()
    text = Path(yaml_path).read_text(encoding="utf-8")
    assert "ligand:" in text and "affinity" in text


def test_ligandmpnn_prepare_job_with_real_cif():
    ensure_data_root()
    cif_path = Path("data/structure/mmcif/3sn6.cif")
    if not cif_path.exists():
        return
    try:
        import gemmi  # noqa: F401
    except Exception:
        return

    manager = ModelManager()
    inv = manager.prepare(
        "ligand_mpnn",
        inputs={"structure_pdb": str(cif_path)},
        config={"seed": 42, "batch_size": 1, "temperature": 0.1, "model_type": "protein_mpnn"},
    )
    assert inv.job is not None
    cmd = inv.job.command
    assert any("LigandMPNN" in part for part in cmd)
    assert "--pdb_path" in cmd and "--out_folder" in cmd


def test_pocket2mol_prepare_job_with_real_cif():
    ensure_data_root()
    cif_path = Path("data/structure/mmcif/3sn6.cif")
    if not cif_path.exists():
        return
    try:
        import gemmi  # noqa: F401
    except Exception:
        return

    manager = ModelManager()
    inv = manager.prepare(
        "pocket2mol",
        inputs={"structure_pdb": str(cif_path)},
        config={"bbox_size": 20.0},
    )
    assert inv.job is not None
    cmd = inv.job.command
    assert any("Pocket2Mol" in part for part in cmd)
    assert "--center" in cmd and "--bbox_size" in cmd


def test_unidock_prepare_job():
    ensure_data_root()
    # Use bundled Uni-Dock test assets if available
    receptor = Path("protos/src/protos/models/Uni-Dock/unidock/test/receptor/1a30_protein.pdbqt")
    ligand = Path("protos/src/protos/models/Uni-Dock/unidock/test/ligands/1a30_ligand.sdf")
    if not (receptor.exists() and ligand.exists()):
        return

    manager = ModelManager()
    inv = manager.prepare(
        "unidock",
        inputs={
            "receptor_pdb": str(receptor),
            "ligand_file": str(ligand),
        },
        config={
            "search_mode": "fast",
            "num_modes": 1,
            "scoring": "vina",
        },
    )
    assert inv.job is not None
    cmd = inv.job.command
    # Basic CLI shape
    assert cmd[0] == "unidock"
    joined = " ".join(cmd)
    assert "--receptor" in joined and "--gpu_batch" in joined and "--dir" in joined
    assert "--center_x" in joined and "--size_x" in joined and "--scoring" in joined
