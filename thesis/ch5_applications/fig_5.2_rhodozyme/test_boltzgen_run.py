#!/usr/bin/env python3
"""Test run for BoltzGen rhodozyme design.

Submits a single BoltzGen job to verify the pipeline works before
running the full grid search.
"""
from __future__ import annotations

import sys
import time
import yaml
from pathlib import Path
from datetime import datetime

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"
CONFIGS_DIR = OUTPUT_DIR / "boltzgen_configs"

import protos
from protos.models import ModelManager, JobStatus


def main() -> int:
    """Run a single BoltzGen test job."""
    print("=" * 70)
    print("BOLTZGEN TEST RUN - RHODOZYME DESIGN")
    print("=" * 70)

    # Initialize Protos - must pass data_root explicitly to ModelManager
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))
    print(f"Data root: {data_root}")

    # Use the first trypsin config for testing
    test_config_path = CONFIGS_DIR / "trypsin_rot000_000_000.yaml"

    if not test_config_path.exists():
        print(f"ERROR: Config file not found: {test_config_path}")
        return 1

    print(f"\nTest config: {test_config_path.name}")

    # Load config
    with open(test_config_path) as f:
        config = yaml.safe_load(f)

    print(f"  Entities: {len(config.get('entities', []))}")
    print(f"  Num designs: {config.get('num_designs', 1)}")

    # Initialize ModelManager with explicit data_root
    mm = ModelManager(data_root=data_root, use_api=False)
    print(f"\nExecutor: {type(mm._executor).__name__}")
    print(f"GPU enabled: {mm._executor.use_gpu}")
    print(f"Data root: {mm.paths.data_root}")

    # Prepare the job
    print("\n[1] Preparing job...")
    try:
        invocation = mm.prepare(
            "boltzgen",
            config=config,
            metadata={"job_name": "test_trypsin_rot000"}
        )
        print(f"  Working dir: {invocation.job.working_dir}")
        print(f"  Command: {' '.join(invocation.job.command[:5])}...")
    except Exception as e:
        print(f"ERROR preparing job: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Submit the job
    print("\n[2] Submitting job...")
    try:
        job_state = mm.submit_job(invocation, persistent=True)
        print(f"  Job ID: {job_state.job_id}")
        print(f"  Status: {job_state.status}")
    except Exception as e:
        print(f"ERROR submitting job: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Wait for completion with progress updates
    print("\n[3] Waiting for completion...")
    start_time = time.time()
    last_status = None

    while True:
        state = mm.job_status(job_state.job_id)

        if state.status != last_status:
            elapsed = time.time() - start_time
            print(f"  [{elapsed:.0f}s] Status: {state.status.value}")
            last_status = state.status

        if state.status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED):
            break

        time.sleep(10)  # Poll every 10 seconds

    # Show results
    elapsed = time.time() - start_time
    print(f"\n[4] Job completed in {elapsed:.1f}s")
    print(f"  Final status: {state.status.value}")

    if state.error:
        print(f"  Error: {state.error}")

    # Get result details
    result = mm.job_result(job_state.job_id)
    if result:
        print(f"\n  Exit code: {result.exit_code}")
        print(f"  Duration: {result.duration_seconds:.1f}s")
        print(f"  Output files: {len(result.output_files)}")

        if result.output_files:
            print("\n  Generated files:")
            for f in result.output_files[:10]:  # Show first 10
                print(f"    - {f}")
            if len(result.output_files) > 10:
                print(f"    ... and {len(result.output_files) - 10} more")

        # Show stderr if there were issues
        if result.exit_code != 0 and result.stderr:
            print(f"\n  Stderr (last 500 chars):")
            print(f"    {result.stderr[-500:]}")

    # Ingest outputs to register designs in data/structure/mmcif/
    if state.status == JobStatus.COMPLETED:
        print("\n[5] Ingesting outputs...")
        summary = mm.ingest_job_outputs(job_state.job_id, cleanup=False)
        print(f"  Ingested: {len(summary.get('ingested', []))} items")
        for item in summary.get('ingested', []):
            item_type = item.get('type', 'unknown')
            item_name = item.get('name', item.get('path', '?'))
            print(f"    - {item_type}: {item_name}")

        # Show where structures were saved
        mmcif_dir = data_root / "structure" / "mmcif"
        print(f"\n  Structures saved to: {mmcif_dir}")

    # Summary
    print("\n" + "=" * 70)
    if state.status == JobStatus.COMPLETED:
        print("TEST RUN SUCCESSFUL")
        print(f"Output directory: {invocation.job.working_dir}")
        print(f"Designs registered in: {data_root}/structure/mmcif/")
        return 0
    else:
        print("TEST RUN FAILED")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
