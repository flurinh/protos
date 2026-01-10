from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import pytest

# Ensure the Protos package resolves to the local source tree
REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager
from protos.models.model_specs import ModelInvocation
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor


RHO_SEQUENCE = (
    "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGF"
    "TTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGM"
    "QCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYI"
    "FTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"
)

def _require_dependencies() -> None:
    missing = []
    try:
        import torch  # noqa: F401
    except Exception:
        missing.append("torch")
    try:
        import torch_geometric  # noqa: F401
    except Exception:
        missing.append("torch_geometric")
    try:
        import transformers  # noqa: F401
    except Exception:
        missing.append("transformers")

    if missing:
        raise RuntimeError(
            "Missing optional dependencies: " + ", ".join(missing) +
            ". Install protos with the lambda/embedding extras or run pip install "
            "torch torch_geometric transformers"
        )


def ensure_data_root(base_dir: Path = None) -> Path:
    data_root = base_dir or (REPO_ROOT / "data")
    data_root.mkdir(parents=True, exist_ok=True)
    try:
        protos.set_data_path(str(data_root))
    except Exception:
        pass
    return data_root


def run_lambda_prediction(
    data_root: Path,
    *,
    verbose: bool = True,
    debug: bool = False,
    use_docker: bool = False,
    run_docker: bool = False,
) -> tuple[ModelInvocation, Optional[pd.DataFrame], Optional[str]]:
    """Run Lambda prediction workflow.

    Args:
        data_root: Directory to store Protos artifacts
        verbose: Print progress messages
        debug: Enable debug output
        use_docker: If True, prepare job for Docker execution instead of runtime
        run_docker: If True and use_docker=True, actually run the Docker container

    Returns:
        Tuple of (invocation, predictions_df, property_table_name)
    """
    ensure_data_root(data_root)
    if verbose:
        print(f"[lambda] data_root: {data_root}")

    manager = ModelManager(data_root=data_root)
    seq_proc = SequenceProcessor()

    dataset_name = "rhodopsin_bovine_uniprot"
    sequence_id = "OPSD_BOVIN"
    sequences = {sequence_id: RHO_SEQUENCE}

    if verbose:
        print("[lambda] saving sequence dataset")
    seq_proc.save_sequences(
        sequences,
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )

    grn_table_name = f"{dataset_name}_grn"
    if verbose:
        print("[lambda] annotating sequences with GRN")
    seq_proc.annotate_with_grn(
        dataset_name=dataset_name,
        reference_table="vpod1_2",
        protein_family="gpcr_a",
        output_table=grn_table_name,
        allow_create=True,
        return_summary=False,
    )

    invocation = manager.prepare(
        "lambda",
        inputs={
            "sequence_dataset": dataset_name,
            "grn_table": grn_table_name,
            "protein_family": "gpcr_a",
        },
        config={
            "run_id": "007061",
            "property_table": "lambda_rhodopsin_predictions",
            "batch_size": 1,
            "embedding_model": "ankh_large",
            "embedding_type": "per_residue",
            "ingest_embeddings": False,
            "debug": debug,
            "use_docker": use_docker,
            "job_name": f"lambda_rhodopsin_{dataset_name}",
        },
    )

    # Handle Docker execution mode
    if use_docker and invocation.job is not None:
        job = invocation.job
        if verbose:
            print("[lambda] prepared Docker job:")
            print("  command:", " ".join(job.command))
            print("  working_dir:", job.working_dir)
            print("  input_dir:", job.metadata.get("input_dir"))
            print("  output_dir:", job.metadata.get("output_dir"))

        if run_docker:
            if verbose:
                print("[lambda] running Docker container...")

            result = subprocess.run(
                job.command,
                cwd=job.working_dir,
                capture_output=True,
                text=True,
            )

            if result.returncode != 0:
                print(f"[lambda] Docker failed with code {result.returncode}")
                print(f"[lambda] stderr: {result.stderr}")
                raise RuntimeError(f"Docker execution failed: {result.stderr}")

            if verbose:
                print("[lambda] Docker execution completed")

            # Load predictions from Docker output
            output_dir = Path(job.metadata.get("output_dir"))
            predictions_path = output_dir / "predictions.csv"
            property_rows_path = output_dir / "property_rows.json"

            if predictions_path.exists():
                predictions = pd.read_csv(predictions_path)
                if verbose:
                    print(f"[lambda] loaded predictions: {len(predictions)} rows")

                # Ingest property rows if available
                if property_rows_path.exists():
                    property_rows = json.loads(property_rows_path.read_text())
                    property_table = f"lambda_rhodopsin_docker_{dataset_name}"

                    prop_proc = PropertyProcessor()
                    prop_proc.record_properties(
                        property_table,
                        property_rows,
                        metadata={
                            "model": "lambda",
                            "execution": "docker",
                            "run_id": job.metadata.get("run_config", {}).get("run_id"),
                        },
                        allow_create=True,
                    )
                    if verbose:
                        print(f"[lambda] recorded property table: {property_table}")

                    return invocation, predictions, property_table

                return invocation, predictions, None

        return invocation, None, None

    # Handle runtime execution mode (original behavior)
    runtime = invocation.runtime
    if runtime is not None:
        predictions = runtime.outputs["predictions"]
        property_table = runtime.metadata["property_table"]

        if verbose:
            print(f"[lambda] predictions rows: {len(predictions)}")

        prop_proc = PropertyProcessor()
        table_df = prop_proc.load_table(property_table)
        if verbose:
            print(
                f"[lambda] property table '{property_table}' rows: {len(table_df)}"
            )
        return invocation, predictions, property_table

    if verbose:
        job = invocation.job
        if job:
            print("[lambda] prepared external job:")
            print("  command:", " ".join(job.command))
            print("  working_dir:", job.working_dir)

    return invocation, None, None




def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run Lambda prediction on rhodopsin")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=REPO_ROOT / "data",
        help="Directory to store Protos artifacts",
    )
    parser.add_argument("--quiet", action="store_true", help="Reduce console logging")
    parser.add_argument(
        "--docker",
        action="store_true",
        help="Prepare job for Docker execution instead of runtime",
    )
    parser.add_argument(
        "--run-docker",
        action="store_true",
        help="Actually run the Docker container (implies --docker)",
    )
    args = parser.parse_args(argv)

    use_docker = args.docker or args.run_docker
    run_docker = args.run_docker

    # Only require local dependencies if not using Docker
    if not use_docker:
        try:
            _require_dependencies()
        except RuntimeError as exc:
            print(exc)
            print("Hint: Use --docker to prepare a job for Docker execution")
            return 1

    debug = not args.quiet
    try:
        invocation, predictions, property_table = run_lambda_prediction(
            args.data_root,
            verbose=not args.quiet,
            debug=debug,
            use_docker=use_docker,
            run_docker=run_docker,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"Lambda workflow failed: {exc}")
        return 1

    if predictions is not None and property_table:
        prop_proc = PropertyProcessor()
        property_path = prop_proc.tables_dir / f"{property_table}.csv"
        print(f"Prediction rows: {len(predictions)}")
        print(f"Property table saved to: {property_path}")
    elif predictions is not None:
        print(f"Prediction rows: {len(predictions)}")
        print("(Property table not created)")
    elif invocation.job is not None:
        job = invocation.job
        print("Prepared Lambda Docker job:")
        print("  command:", " ".join(job.command))
        print("  working_dir:", job.working_dir)
        print("  input_dir:", job.metadata.get("input_dir"))
        print("  output_dir:", job.metadata.get("output_dir"))
        print("  artifacts:")
        for bundle in job.artifacts:
            print("   -", bundle.spec.name, bundle.spec.kind, bundle.path)
        if not run_docker:
            print("\nTo run the job, execute:")
            print(f"  {' '.join(job.command)}")
    elif invocation.runtime is not None:
        print("Lambda runtime executed but produced no predictions")
    else:
        print("Lambda invocation produced no runtime or job metadata")
        return 1

    return 0


@pytest.mark.integration
def test_lambda_workflow(tmp_path):
    """Test Lambda workflow with local runtime (requires torch/torch_geometric)."""
    try:
        _require_dependencies()
    except RuntimeError as exc:
        pytest.skip(str(exc))

    data_root = tmp_path / "protos_data"
    try:
        invocation, predictions, property_table = run_lambda_prediction(
            data_root, debug=True
        )
    except RuntimeError as exc:
        pytest.skip(str(exc))

    if invocation.runtime is not None and predictions is not None and property_table:
        assert isinstance(predictions, pd.DataFrame)
        assert not predictions.empty

        prop_proc = PropertyProcessor()
        table_df = prop_proc.load_table(property_table)
        assert not table_df.empty
        assert any(col.endswith("007061") for col in table_df.columns)
    else:
        job = invocation.job
        assert job is not None, "Lambda submission did not include a job"
        assert job.command, "Prepared job missing command"
        assert job.working_dir is not None
        assert job.working_dir.exists(), "Job working dir missing"
        assert job.artifacts, "Lambda submission should package artifacts"


@pytest.mark.integration
@pytest.mark.docker
def test_lambda_workflow_docker(tmp_path):
    """Test Lambda workflow with Docker execution."""
    data_root = tmp_path / "protos_data"

    # First test: prepare job only (no Docker run)
    invocation, predictions, property_table = run_lambda_prediction(
        data_root,
        verbose=True,
        debug=True,
        use_docker=True,
        run_docker=False,
    )

    job = invocation.job
    assert job is not None, "Docker mode should produce a job"
    assert job.command, "Prepared job missing command"
    assert "docker" in job.command[0], "Command should use docker"
    assert job.working_dir is not None
    assert job.working_dir.exists(), "Job working dir missing"
    assert job.artifacts, "Lambda submission should package artifacts"

    # Verify input files were created
    input_dir = Path(job.metadata.get("input_dir"))
    assert (input_dir / "aligned_embeddings.npz").exists(), "Embeddings file missing"
    assert (input_dir / "grn_table.csv").exists(), "GRN table file missing"
    assert (input_dir / "config.json").exists(), "Config file missing"


@pytest.mark.integration
@pytest.mark.docker
@pytest.mark.slow
def test_lambda_workflow_docker_run(tmp_path):
    """Test Lambda workflow with actual Docker container execution."""
    import subprocess

    # Check if Docker is available
    try:
        result = subprocess.run(
            ["docker", "images", "protos/lambda:latest", "-q"],
            capture_output=True,
            text=True,
        )
        if not result.stdout.strip():
            pytest.skip("protos/lambda:latest Docker image not found")
    except FileNotFoundError:
        pytest.skip("Docker not available")

    data_root = tmp_path / "protos_data"

    invocation, predictions, property_table = run_lambda_prediction(
        data_root,
        verbose=True,
        debug=True,
        use_docker=True,
        run_docker=True,
    )

    assert predictions is not None, "Docker execution should produce predictions"
    assert isinstance(predictions, pd.DataFrame)
    assert not predictions.empty

    if property_table:
        prop_proc = PropertyProcessor()
        table_df = prop_proc.load_table(property_table)
        assert not table_df.empty


if __name__ == "__main__":
    raise SystemExit(main())
