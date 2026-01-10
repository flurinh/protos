#!/usr/bin/env python3
"""Opsin Lambda prediction workflow using Docker container.

This script:
1. Fetches opsin sequences from UniProt (Rh1, Lw, Sw, etc.)
2. Registers them as a protos dataset
3. Annotates with GRN (generic residue numbering)
4. Generates Ankh Large embeddings locally (CPU)
5. Runs Lambda prediction via Docker container
6. Ingests predictions into property table

Usage:
    python test_docker_opsin.py                    # Prepare job only
    python test_docker_opsin.py --run              # Prepare and run Docker
    python test_docker_opsin.py --run --gpu        # Run with GPU support
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.io.ingest.utils.uniprot_utils import get_uniprot
from protos.models.model_manager import ModelManager
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.sequence import SequenceProcessor


# Curated list of opsin UniProt IDs for testing
# These are well-characterized opsins with known spectral properties
OPSIN_UNIPROT_IDS = {
    # Rhodopsins (Rh1) - ~500nm
    "OPSD_BOVIN": "P02699",   # Bovine rhodopsin - reference structure
    "OPSD_HUMAN": "P08100",   # Human rhodopsin
    "OPSD_MOUSE": "P15409",   # Mouse rhodopsin

    # Long-wavelength sensitive (LWS/Lw) - ~560nm
    "OPSR_HUMAN": "P04000",   # Human red opsin
    "OPSR_MOUSE": "P51489",   # Mouse red opsin (M-cone)

    # Middle-wavelength sensitive (MWS) - ~530nm
    "OPSG_HUMAN": "P04001",   # Human green opsin

    # Short-wavelength sensitive (SWS1/Sw) - ~420nm
    "OPSB_HUMAN": "P03999",   # Human blue opsin
    "OPSB_MOUSE": "P51491",   # Mouse blue opsin (S-cone)

    # Other vertebrate opsins
    "OPN4_HUMAN": "Q9UHM6",   # Melanopsin (non-visual)
    "RGR_HUMAN": "P47804",    # Retinal G protein-coupled receptor
}


def ensure_data_root(base_dir: Optional[Path] = None) -> Path:
    """Set up protos data root directory."""
    data_root = base_dir or (REPO_ROOT / "data")
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def fetch_opsin_sequences(
    uniprot_ids: Dict[str, str],
    verbose: bool = True,
) -> Dict[str, str]:
    """Fetch opsin sequences from UniProt.

    Args:
        uniprot_ids: Dict mapping entry names to UniProt accessions
        verbose: Print progress messages

    Returns:
        Dict mapping entry names to sequences
    """
    sequences = {}

    for entry_name, accession in uniprot_ids.items():
        if verbose:
            print(f"  Fetching {entry_name} ({accession})...")

        try:
            df = get_uniprot(accession, batchsize=1, reviewed=True)

            if not df.empty and "Sequence" in df.columns:
                seq = df.iloc[0]["Sequence"]
                sequences[entry_name] = seq
                if verbose:
                    print(f"    OK ({len(seq)} aa)")
            else:
                # Try without reviewed filter
                df = get_uniprot(accession, batchsize=1, reviewed=False)
                if not df.empty and "Sequence" in df.columns:
                    seq = df.iloc[0]["Sequence"]
                    sequences[entry_name] = seq
                    if verbose:
                        print(f"    OK ({len(seq)} aa, unreviewed)")
                else:
                    if verbose:
                        print(f"    SKIP: no sequence found")
        except Exception as exc:
            if verbose:
                print(f"    ERROR: {exc}")

    return sequences


def run_opsin_docker_workflow(
    data_root: Path,
    *,
    run_container: bool = False,
    use_gpu: bool = False,
    verbose: bool = True,
    max_sequences: Optional[int] = None,
) -> dict:
    """Run Lambda prediction on opsins via Docker.

    Args:
        data_root: Directory to store Protos artifacts
        run_container: If True, actually run the Docker container
        use_gpu: If True, run with GPU support (--gpus all)
        verbose: Print progress messages
        max_sequences: Limit number of sequences (for testing)

    Returns:
        Dict with workflow results including job info and predictions
    """
    ensure_data_root(data_root)
    if verbose:
        print(f"[opsin-docker] data_root: {data_root}")

    # Select sequences
    opsin_ids = dict(list(OPSIN_UNIPROT_IDS.items())[:max_sequences]) if max_sequences else OPSIN_UNIPROT_IDS

    # Step 1: Fetch sequences from UniProt
    if verbose:
        print(f"\n[opsin-docker] Step 1: Fetching {len(opsin_ids)} opsin sequences from UniProt")

    sequences = fetch_opsin_sequences(opsin_ids, verbose=verbose)

    if not sequences:
        raise RuntimeError("No sequences fetched from UniProt")

    if verbose:
        print(f"  Retrieved {len(sequences)} sequences")

    # Step 2: Register sequences
    dataset_name = "opsin_docker_test"
    if verbose:
        print(f"\n[opsin-docker] Step 2: Registering sequence dataset '{dataset_name}'")

    seq_proc = SequenceProcessor()
    seq_proc.save_sequences(
        sequences,
        output_file=dataset_name,
        dataset_name=dataset_name,
        materialize_entities=True,
    )

    # Step 3: Annotate with GRN
    grn_table_name = f"{dataset_name}_grn"
    if verbose:
        print(f"\n[opsin-docker] Step 3: Annotating with GRN (vpod1_2 reference)")

    seq_proc.annotate_with_grn(
        dataset_name=dataset_name,
        reference_table="vpod1_2",
        protein_family="gpcr_a",
        output_table=grn_table_name,
        allow_create=True,
        return_summary=False,
    )

    # Step 4: Generate embeddings locally
    embedding_dataset_name = f"{dataset_name}__ankh_large__per_residue"
    if verbose:
        print(f"\n[opsin-docker] Step 4: Generating Ankh Large embeddings (CPU)")

    # Force CPU for Ankh Large (1.15B params) to avoid GPU OOM on small GPUs
    emb_proc = EmbeddingProcessor(model_name="ankh_large", device="cpu")
    emb_proc.embed_sequences(
        sequences,
        embedding_type="per_residue",
        save_dataset=embedding_dataset_name,
    )
    if verbose:
        print(f"  Saved embedding dataset: {embedding_dataset_name}")

    # Step 5: Prepare Lambda job with Docker mode
    if verbose:
        print(f"\n[opsin-docker] Step 5: Preparing Lambda Docker job")

    manager = ModelManager(data_root=data_root)

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
        "sequences": sequences,
        "predictions": None,
        "property_table": None,
    }

    if verbose:
        print(f"  Docker command: {' '.join(job.command[:6])}...")
        print(f"  Input dir: {result['input_dir']}")
        print(f"  Output dir: {result['output_dir']}")

    # Verify input files
    input_dir = result["input_dir"]
    assert (input_dir / "aligned_embeddings.npz").exists(), "Embeddings not created"
    assert (input_dir / "grn_table.csv").exists(), "GRN table not created"
    assert (input_dir / "config.json").exists(), "Config not created"

    if verbose:
        print(f"  Input files:")
        for f in input_dir.iterdir():
            print(f"    {f.name} ({f.stat().st_size:,} bytes)")

    if not run_container:
        if verbose:
            print(f"\n[opsin-docker] Job prepared. To run:")
            print(f"  {' '.join(job.command)}")
        return result

    # Step 6: Run Docker container
    if verbose:
        print(f"\n[opsin-docker] Step 6: Running Docker container...")

    proc = subprocess.run(
        job.command,
        cwd=job.working_dir,
        capture_output=True,
        text=True,
    )

    if proc.returncode != 0:
        print(f"[opsin-docker] Docker failed (exit code {proc.returncode})")
        print(f"[opsin-docker] stdout: {proc.stdout}")
        print(f"[opsin-docker] stderr: {proc.stderr}")
        raise RuntimeError(f"Docker execution failed: {proc.stderr}")

    if verbose:
        print("[opsin-docker] Docker execution completed")

    # Step 7: Load predictions
    output_dir = result["output_dir"]
    predictions_path = output_dir / "predictions.csv"
    property_rows_path = output_dir / "property_rows.json"
    metadata_path = output_dir / "metadata.json"

    if predictions_path.exists():
        predictions = pd.read_csv(predictions_path)
        result["predictions"] = predictions
        if verbose:
            print(f"\n[opsin-docker] Predictions: {len(predictions)} sequences")
            print(predictions.to_string())

    if metadata_path.exists():
        metadata = json.loads(metadata_path.read_text())
        if verbose:
            print(f"\n[opsin-docker] Run metadata: {metadata}")

    # Step 8: Ingest property rows
    if property_rows_path.exists():
        property_rows = json.loads(property_rows_path.read_text())
        property_table = f"lambda_opsin_{dataset_name}"

        prop_proc = PropertyProcessor()
        prop_proc.record_properties(
            property_table,
            property_rows,
            metadata={
                "model": "lambda",
                "execution": "docker",
                "run_id": "007061",
                "sequences": list(sequences.keys()),
            },
            allow_create=True,
        )
        result["property_table"] = property_table
        if verbose:
            print(f"[opsin-docker] Recorded property table: {property_table}")

    return result


def main(argv: Optional[list[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Lambda prediction on opsins via Docker container"
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
        "--max-sequences",
        type=int,
        default=None,
        help="Limit number of sequences (for quick testing)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce console logging",
    )
    args = parser.parse_args(argv)

    try:
        result = run_opsin_docker_workflow(
            args.data_root,
            run_container=args.run,
            use_gpu=args.gpu,
            verbose=not args.quiet,
            max_sequences=args.max_sequences,
        )
    except Exception as exc:
        print(f"Opsin Docker workflow failed: {exc}")
        import traceback
        traceback.print_exc()
        return 1

    if result["predictions"] is not None:
        print(f"\nSuccess! Predictions: {len(result['predictions'])} sequences")
        if result["property_table"]:
            print(f"Property table: {result['property_table']}")

        # Show summary
        preds = result["predictions"]
        print("\nSpectral predictions (Lambda max in nm):")
        for _, row in preds.iterrows():
            name = row["protein_id"]
            lmax = row.get("Lmax_11_cis_pred_007061", row.get("lmax_11_cis_pred_007061", "N/A"))
            print(f"  {name}: {lmax:.1f} nm" if isinstance(lmax, float) else f"  {name}: {lmax}")
    else:
        print("\nJob prepared but not executed (use --run to execute)")

    return 0


# Pytest tests

@pytest.fixture
def opsin_docker_env(tmp_path):
    """Set up environment for opsin Docker tests."""
    data_root = tmp_path / "protos_data"
    ensure_data_root(data_root)
    return data_root


def test_opsin_docker_job_preparation(opsin_docker_env):
    """Test that opsin Docker job preparation works."""
    result = run_opsin_docker_workflow(
        opsin_docker_env,
        run_container=False,
        verbose=True,
        max_sequences=2,  # Quick test with just 2 sequences
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
def test_opsin_docker_execution(opsin_docker_env):
    """Test actual Docker container execution."""
    # Check if Docker image exists
    proc = subprocess.run(
        ["docker", "images", "protos/lambda:latest", "-q"],
        capture_output=True,
        text=True,
    )
    if not proc.stdout.strip():
        pytest.skip("protos/lambda:latest Docker image not found")

    result = run_opsin_docker_workflow(
        opsin_docker_env,
        run_container=True,
        verbose=True,
        max_sequences=3,
    )

    assert result["predictions"] is not None
    assert isinstance(result["predictions"], pd.DataFrame)
    assert not result["predictions"].empty


if __name__ == "__main__":
    raise SystemExit(main())
