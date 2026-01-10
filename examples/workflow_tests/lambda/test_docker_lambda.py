#!/usr/bin/env python3
"""Lambda workflow using Docker container for model execution.

This script prepares inputs using protos ModelManager and runs the Lambda
model via the protos/lambda:latest Docker container.

Usage:
    python test_docker_lambda.py                    # Prepare job only
    python test_docker_lambda.py --run              # Prepare and run Docker
    python test_docker_lambda.py --run --gpu        # Run with GPU support
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor


RHO_SEQUENCE = (
    "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLA"
    "VADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFT"
    "WVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQES"
    "ATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTT"
    "ICCGKNPLGDDEASATVSKTETSQVAPA"
)


def ensure_data_root(base_dir: Optional[Path] = None) -> Path:
    data_root = base_dir or (REPO_ROOT / "data")
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def run_docker_lambda_workflow(
    data_root: Path,
    *,
    run_container: bool = False,
    use_gpu: bool = False,
    verbose: bool = True,
) -> dict:
    """Run Lambda prediction workflow via Docker.

    Args:
        data_root: Directory to store Protos artifacts
        run_container: If True, actually run the Docker container
        use_gpu: If True, run with GPU support (--gpus all)
        verbose: Print progress messages

    Returns:
        Dict with workflow results including job info and predictions
    """
    ensure_data_root(data_root)
    if verbose:
        print(f"[lambda-docker] data_root: {data_root}")

    manager = ModelManager(data_root=data_root)
    seq_proc = SequenceProcessor()

    dataset_name = "rhodopsin_docker_test"
    sequence_id = "OPSD_BOVIN"
    sequences = {sequence_id: RHO_SEQUENCE}

    # Step 1: Save sequences
    if verbose:
        print("[lambda-docker] saving sequence dataset")
    seq_proc.save_sequences(
        sequences,
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )

    # Step 2: Annotate with GRN
    grn_table_name = f"{dataset_name}_grn"
    if verbose:
        print("[lambda-docker] annotating sequences with GRN")
    seq_proc.annotate_with_grn(
        dataset_name=dataset_name,
        reference_table="vpod1_2",
        protein_family="gpcr_a",
        output_table=grn_table_name,
        allow_create=True,
        return_summary=False,
    )

    # Step 3: Generate embeddings locally (required for Lambda)
    embedding_dataset_name = f"{dataset_name}__ankh_large__per_residue"
    if verbose:
        print("[lambda-docker] generating Ankh Large embeddings locally")

    # Force CPU for Ankh Large (1.15B params) to avoid GPU OOM on small GPUs
    emb_proc = EmbeddingProcessor(model_name="ankh_large", device="cpu")
    emb_proc.embed_sequences(
        sequences,
        embedding_type="per_residue",
        save_dataset=embedding_dataset_name,
    )
    if verbose:
        print(f"[lambda-docker] saved embedding dataset: {embedding_dataset_name}")

    # Step 4: Prepare Lambda job with Docker mode
    if verbose:
        print("[lambda-docker] preparing Lambda job for Docker execution")

    invocation = manager.prepare(
        "lambda",
        inputs={
            "sequence_dataset": dataset_name,
            "grn_table": grn_table_name,
            "embedding_dataset": embedding_dataset_name,
            "protein_family": "gpcr_a",
        },
        config={
            "run_id": "007061",
            "batch_size": 1,
            "embedding_model": "ankh_large",
            "embedding_type": "per_residue",
            "use_docker": True,
            "use_gpu": use_gpu,
            "job_name": f"docker_lambda_{dataset_name}",
        },
    )

    job = invocation.job
    if job is None:
        raise RuntimeError("Lambda adapter did not produce a Docker job")

    result = {
        "job": job,
        "input_dir": Path(job.metadata.get("input_dir")),
        "output_dir": Path(job.metadata.get("output_dir")),
        "command": job.command,
        "predictions": None,
        "property_table": None,
    }

    if verbose:
        print("[lambda-docker] prepared Docker job:")
        print(f"  command: {' '.join(job.command)}")
        print(f"  input_dir: {result['input_dir']}")
        print(f"  output_dir: {result['output_dir']}")

    # Verify input files
    input_dir = result["input_dir"]
    assert (input_dir / "aligned_embeddings.npz").exists(), "Embeddings not created"
    assert (input_dir / "grn_table.csv").exists(), "GRN table not created"
    assert (input_dir / "config.json").exists(), "Config not created"

    if verbose:
        print("[lambda-docker] input files verified:")
        for f in input_dir.iterdir():
            print(f"    {f.name} ({f.stat().st_size} bytes)")

    if not run_container:
        if verbose:
            print("\n[lambda-docker] Job prepared. To run:")
            print(f"  {' '.join(job.command)}")
        return result

    # Step 5: Run Docker container
    if verbose:
        print("[lambda-docker] running Docker container...")

    proc = subprocess.run(
        job.command,
        cwd=job.working_dir,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        print(f"[lambda-docker] Docker failed (exit code {proc.returncode})")
        print(f"[lambda-docker] stdout: {proc.stdout}")
        print(f"[lambda-docker] stderr: {proc.stderr}")
        raise RuntimeError(f"Docker execution failed: {proc.stderr}")

    if verbose:
        print("[lambda-docker] Docker execution completed")

    # Step 6: Load predictions
    output_dir = result["output_dir"]
    predictions_path = output_dir / "predictions.csv"
    property_rows_path = output_dir / "property_rows.json"
    metadata_path = output_dir / "metadata.json"

    if predictions_path.exists():
        predictions = pd.read_csv(predictions_path)
        result["predictions"] = predictions
        if verbose:
            print(f"[lambda-docker] loaded predictions: {len(predictions)} rows")
            print(predictions.head())

    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if verbose:
            print(f"[lambda-docker] run metadata: {metadata}")

    # Step 7: Ingest property rows
    if property_rows_path.exists():
        property_rows = json.loads(property_rows_path.read_text())
        property_table = f"lambda_docker_{dataset_name}"

        prop_proc = PropertyProcessor()
        prop_proc.record_properties(
            property_table,
            property_rows,
            metadata={
                "model": "lambda",
                "execution": "docker",
                "run_id": "007061",
            },
            allow_create=True,
        )
        result["property_table"] = property_table
        if verbose:
            print(f"[lambda-docker] recorded property table: {property_table}")

    return result


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Lambda prediction via Docker container"
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=REPO_ROOT / "data",
        help="Directory to store Protos artifacts",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Actually run the Docker container",
    )
    parser.add_argument(
        "--gpu",
        action="store_true",
        help="Run with GPU support (--gpus all)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce console logging",
    )
    args = parser.parse_args(argv)

    try:
        result = run_docker_lambda_workflow(
            args.data_root,
            run_container=args.run,
            use_gpu=args.gpu,
            verbose=not args.quiet,
        )
    except Exception as exc:
        print(f"Lambda Docker workflow failed: {exc}")
        return 1

    if result["predictions"] is not None:
        print(f"\nSuccess! Predictions: {len(result['predictions'])} rows")
        if result["property_table"]:
            print(f"Property table: {result['property_table']}")
    else:
        print("\nJob prepared but not executed (use --run to execute)")

    return 0


# Pytest tests

@pytest.fixture
def lambda_docker_env(tmp_path):
    """Set up environment for Lambda Docker tests."""
    data_root = tmp_path / "protos_data"
    ensure_data_root(data_root)
    return data_root


def test_docker_lambda_job_preparation(lambda_docker_env):
    """Test that Lambda Docker job preparation works."""
    result = run_docker_lambda_workflow(
        lambda_docker_env,
        run_container=False,
        verbose=True,
    )

    job = result["job"]
    assert job is not None
    assert job.command
    assert "docker" in job.command[0]

    input_dir = result["input_dir"]
    assert input_dir.exists()
    assert (input_dir / "aligned_embeddings.npz").exists()
    assert (input_dir / "grn_table.csv").exists()
    assert (input_dir / "config.json").exists()


@pytest.mark.docker
@pytest.mark.slow
def test_docker_lambda_execution(lambda_docker_env):
    """Test actual Docker container execution."""
    # Check if Docker image exists
    proc = subprocess.run(
        ["docker", "images", "protos/lambda:latest", "-q"],
        capture_output=True,
        text=True,
    )
    if not proc.stdout.strip():
        pytest.skip("protos/lambda:latest Docker image not found")

    result = run_docker_lambda_workflow(
        lambda_docker_env,
        run_container=True,
        verbose=True,
    )

    assert result["predictions"] is not None
    assert isinstance(result["predictions"], pd.DataFrame)
    assert not result["predictions"].empty


if __name__ == "__main__":
    raise SystemExit(main())
