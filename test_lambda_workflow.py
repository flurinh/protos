from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd
import pytest

# Ensure the Protos package resolves to the local source tree
SRC_DIR = Path(__file__).resolve().parent / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager
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


def ensure_data_root(base_dir: Path) -> Path:
    # Always prefer protos/data per repository guidelines
    protos_data = Path(__file__).resolve().parent / "data"
    data_root = protos_data
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
) -> tuple[pd.DataFrame, str]:
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
        },
    )

    runtime = invocation.runtime
    if runtime is None:
        raise RuntimeError("Lambda invocation did not return runtime results")

    predictions = runtime.outputs["predictions"]
    property_table = runtime.metadata["property_table"]

    if verbose:
        print(f"[lambda] predictions rows: {len(predictions)}")

    prop_proc = PropertyProcessor()
    table_df = prop_proc.load_table(property_table)
    if verbose:
        print(f"[lambda] property table '{property_table}' rows: {len(table_df)}")

    return predictions, property_table




def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run Lambda prediction on rhodopsin")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
        help="Directory to store Protos artifacts",
    )
    parser.add_argument("--quiet", action="store_true", help="Reduce console logging")
    args = parser.parse_args(argv)

    try:
        _require_dependencies()
    except RuntimeError as exc:
        print(exc)
        return 1

    debug = not args.quiet
    try:
        predictions, property_table = run_lambda_prediction(
            args.data_root,
            verbose=not args.quiet,
            debug=debug,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"Lambda workflow failed: {exc}")
        return 1

    prop_proc = PropertyProcessor()
    property_path = prop_proc.tables_dir / f"{property_table}.csv"

    print(f"Prediction rows: {len(predictions)}")
    print(f"Property table saved to: {property_path}")
    return 0


@pytest.mark.integration
def test_lambda_workflow(tmp_path):
    try:
        _require_dependencies()
    except RuntimeError as exc:
        pytest.skip(str(exc))

    data_root = tmp_path / "protos_data"
    try:
        predictions, property_table = run_lambda_prediction(data_root, debug=True)
    except RuntimeError as exc:
        pytest.skip(str(exc))

    assert isinstance(predictions, pd.DataFrame)
    assert not predictions.empty

    prop_proc = PropertyProcessor()
    table_df = prop_proc.load_table(property_table)
    assert not table_df.empty
    assert any(col.endswith("007061") for col in table_df.columns)


if __name__ == "__main__":
    raise SystemExit(main())
