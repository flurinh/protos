"""Generate a real Lambda graph and send it to the Docker runtime."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict

import importlib
import numpy as np
import pytest

SRC_DIR = Path(__file__).resolve().parent / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.processing.grn import GRNProcessor
from protos.processing.embedding.embedding_processor import EmbeddingProcessor

lambda_utils = importlib.import_module("protos.models.lambda.lmda.runtime_utils")
align_embeddings_to_grn = lambda_utils.align_embeddings_to_grn
build_grn_assignments = lambda_utils.build_grn_assignments
Predictor = importlib.import_module(
    "protos.models.lambda.lmda.inference.predictor"
).Predictor

RHO_SEQUENCE = (
    "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGF"
    "TTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGM"
)

DOCKER_IMAGE = os.environ.get("LAMBDA_DOCKER_IMAGE", "lambda:latest")
LAMBDA_SRC = Path(__file__).resolve().parent / "src/protos/models/lambda/lmda"


def _require_embedding_deps() -> None:
    missing = []
    try:
        import torch  # noqa: F401
    except Exception:
        missing.append("torch")
    try:
        import transformers  # noqa: F401
    except Exception:
        missing.append("transformers")
    if missing:
        raise RuntimeError("Missing dependencies: " + ", ".join(missing))


def ensure_sequence_dataset(tmp_path: Path) -> str:
    protos.set_data_path(str(tmp_path))
    seq_proc = SequenceProcessor()
    dataset_name = "docker_lambda_sequences"
    seq_proc.save_sequences(
        {"RHO": RHO_SEQUENCE},
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )
    return dataset_name


def build_real_graph(tmp_path: Path, dataset_name: str) -> Dict[str, object]:
    seq_proc = SequenceProcessor()
    grn_table_name = f"{dataset_name}_grn"
    seq_proc.annotate_with_grn(
        dataset_name=dataset_name,
        reference_table="vpod1_2",
        protein_family="gpcr_a",
        output_table=grn_table_name,
        allow_create=True,
        return_summary=False,
    )

    grn_proc = GRNProcessor()
    grn_table = grn_proc.load_table(grn_table_name)

    emb_proc = EmbeddingProcessor(model_name="esm2_t12_35m", batch_size=1)
    embeddings = emb_proc.embed_sequences(
        {"RHO": RHO_SEQUENCE},
        embedding_type="per_residue",
        register_entities=False,
    )
    embeddings_np = {
        key: tensor.detach().cpu().numpy() for key, tensor in embeddings.items()
    }

    assignments = build_grn_assignments(grn_table)
    grn_dict, aligned = align_embeddings_to_grn(assignments, embeddings_np)

    protein_id = next(iter(aligned))
    labels = grn_dict[protein_id]
    embedding_arr = aligned[protein_id]

    configs_dir = LAMBDA_SRC / "configs"
    predictor = Predictor(
        "007061",
        verbose=False,
        prefer_checkpoint=False,
        binding_config_path=str(configs_dir / "binding_domain2.json"),
        positional_encoding_path=str(configs_dir / "final_mapping7.csv"),
    )

    graph_data = predictor._prepare_graph_input(
        labels, embedding_arr, protein_id, "gpcr_a"
    )

    return {
        "protein_id": protein_id,
        "protein_family": "gpcr_a",
        "x": graph_data.x.detach().cpu().tolist(),
        "edge_index": graph_data.edge_index.detach().cpu().tolist(),
        "pos": (
            graph_data.pos.detach().cpu().tolist()
            if hasattr(graph_data, "pos")
            else None
        ),
        "grn_labels": labels,
    }


def send_to_container(message: Dict[str, object]) -> Dict[str, object]:
    docker_cmd = ["docker", "run", "--rm", "-i", DOCKER_IMAGE]
    proc = subprocess.run(
        docker_cmd,
        input=json.dumps(message).encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Docker run failed (code {proc.returncode}): {proc.stderr.decode().strip()}"
        )
    return json.loads(proc.stdout.decode("utf-8"))


def run_workflow(tmp_path: Path) -> Dict[str, object]:
    _require_embedding_deps()
    dataset_name = ensure_sequence_dataset(tmp_path)
    graph = build_real_graph(tmp_path, dataset_name)
    message = {
        "run_id": "docker_demo",
        "graph": graph,
    }
    response = send_to_container(message)
    return response


def main() -> int:
    tmp_dir = Path("./temp_docker_lambda")
    shutil.rmtree(tmp_dir, ignore_errors=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    try:
        _require_embedding_deps()
        response = run_workflow(tmp_dir)
    except Exception as exc:
        print(f"docker lambda workflow failed: {exc}")
        return 1
    print(json.dumps(response, indent=2))
    return 0


@pytest.mark.integration
def test_docker_lambda_workflow(tmp_path):
    if shutil.which("docker") is None:
        pytest.skip("Docker is not available")
    try:
        _require_embedding_deps()
    except RuntimeError as exc:
        pytest.skip(str(exc))
    try:
        response = run_workflow(tmp_path)
    except Exception as exc:
        pytest.skip(str(exc))
    assert response.get("status") == "received"
    assert response.get("graph_summary", {}).get("nodes", 0) > 0


if __name__ == "__main__":
    raise SystemExit(main())
