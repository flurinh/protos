from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
import pytest

import protos
from protos.io.paths import get_protos_paths
from protos.io.formats.fasta_utils import read_fasta
from protos.models.model_manager import ModelManager
from protos.models.model_specs import ModelInvocation
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor


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
            "Missing optional dependencies: "
            + ", ".join(missing)
            + ". Install protos with the lambda/embedding extras or run pip install "
            "torch torch_geometric transformers"
        )


def _set_repo_data_root() -> Path:
    """Point Protos to this repo's managed data directory (protos/data).

    This relies entirely on ProtosPaths; all processor paths are derived from it.
    """
    data_root = Path(__file__).resolve().parent / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def _load_fish2_sequences() -> Tuple[str, dict[str, str]]:
    """Load sequences from data/input/fish2.fasta using ProtosPaths only.

    Returns:
        Tuple of (dataset_name, sequences_dict)
    """
    paths = get_protos_paths()
    fasta_path = Path(paths.data_root) / "input" / "fish2.fasta"
    if not fasta_path.exists():
        raise FileNotFoundError(
            f"Input FASTA not found: {fasta_path}. Place fish2.fasta under data/input/."
        )
    sequences = read_fasta(str(fasta_path))
    if not sequences:
        raise ValueError("fish2.fasta contained no sequences")
    return "fish2", sequences


def run_fish2_lambda_prediction(
    *, verbose: bool = True, debug: bool = False
) -> tuple[ModelInvocation, Optional[pd.DataFrame], Optional[str]]:
    _set_repo_data_root()
    if verbose:
        print(f"[fish2] data_root: {get_protos_paths().data_root}")

    dataset_name, sequences = _load_fish2_sequences()

    # Register sequences without removing the original input file
    seq_proc = SequenceProcessor()
    if verbose:
        print("[fish2] registering sequence dataset from input/fish2.fasta")
    seq_proc.save_sequences(
        sequences,
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )

    # Annotate with GRN
    grn_table_name = f"{dataset_name}_grn"
    if verbose:
        print("[fish2] annotating sequences with GRN")
    seq_proc.annotate_with_grn(
        dataset_name=dataset_name,
        reference_table="vpod1_2",
        protein_family="gpcr_a",
        output_table=grn_table_name,
        allow_create=True,
        return_summary=False,
    )

    # Prepare and run Lambda prediction
    manager = ModelManager()
    invocation = manager.prepare(
        "lambda",
        inputs={
            "sequence_dataset": dataset_name,
            "grn_table": grn_table_name,
            "protein_family": "gpcr_a",
        },
        config={
            "run_id": "007061",
            "property_table": "lambda_fish2_predictions",
            "batch_size": 1,
            "embedding_model": "ankh_large",
            "embedding_type": "per_residue",
            "ingest_embeddings": False,
            "debug": debug,
        },
    )

    runtime = invocation.runtime
    if runtime is not None:
        predictions = runtime.outputs["predictions"]
        property_table = runtime.metadata["property_table"]

        if verbose:
            print(f"[fish2] predictions rows: {len(predictions)}")

        prop_proc = PropertyProcessor()
        table_df = prop_proc.load_table(property_table)
        if verbose:
            print(
                f"[fish2] property table '{property_table}' rows: {len(table_df)}"
            )
        return invocation, predictions, property_table

    if verbose and invocation.job is not None:
        print("[fish2] prepared Lambda submission (no runtime)")
        print("  command:", " ".join(invocation.job.command))
        print("  working_dir:", invocation.job.working_dir)

    return invocation, None, None


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run Lambda prediction on fish2 dataset")
    parser.add_argument("--quiet", action="store_true", help="Reduce console logging")
    args = parser.parse_args(argv)

    try:
        _require_dependencies()
    except RuntimeError as exc:
        print(exc)
        return 1

    try:
        invocation, predictions, property_table = run_fish2_lambda_prediction(
            verbose=not args.quiet, debug=not args.quiet
        )
    except Exception as exc:  # noqa: BLE001
        print(f"Fish2 workflow failed: {exc}")
        return 1

    if invocation.runtime is not None and predictions is not None and property_table:
        prop_proc = PropertyProcessor()
        property_path = prop_proc.tables_dir / f"{property_table}.csv"
        print(f"Prediction rows: {len(predictions)}")
        print(f"Property table saved to: {property_path}")
    else:
        job = invocation.job
        if job is None:
            print("Lambda invocation produced neither runtime results nor job metadata")
            return 1
        print("Prepared Lambda submission (fish2)")
        print("  command:", " ".join(job.command))
        print("  working_dir:", job.working_dir)
        print("  artifacts:")
        for bundle in job.artifacts:
            print("   -", bundle.spec.name, bundle.spec.kind, bundle.path)
    return 0


@pytest.mark.integration
def test_fish2_workflow():
    try:
        _require_dependencies()
    except RuntimeError as exc:
        pytest.skip(str(exc))

    # Use repo data root; ensure input file remains intact
    try:
        invocation, predictions, property_table = run_fish2_lambda_prediction(
            debug=True
        )
    except RuntimeError as exc:
        pytest.skip(str(exc))

    if invocation.runtime is not None and predictions is not None and property_table:
        assert isinstance(predictions, pd.DataFrame)
        assert not predictions.empty

        prop_proc = PropertyProcessor()
        table_df = prop_proc.load_table(property_table)
        assert not table_df.empty
        assert "prediction_run_id" in table_df.columns
        assert any(str(v).endswith("007061") for v in table_df["prediction_run_id"]) 
    else:
        job = invocation.job
        assert job is not None
        assert job.command
        assert job.working_dir is not None
        assert job.working_dir.exists()
        assert job.artifacts


if __name__ == "__main__":
    raise SystemExit(main())
