#!/usr/bin/env python3
"""Run BoltzGen with corrected configs (no duplicate ligands)."""
from __future__ import annotations

import sys
import time
import yaml
from pathlib import Path

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.models import ModelManager

# Directories
OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"
CONFIG_DIR = OUTPUT_DIR / "boltzgen_configs_corrected"
DATA_ROOT = REPO_ROOT / "data"


def main() -> int:
    print("=" * 70)
    print("RUNNING BOLTZGEN WITH CORRECTED CONFIGS")
    print("(No duplicate retinal, substrate as separate chain)")
    print("=" * 70)

    # Initialize
    protos.set_data_path(str(DATA_ROOT))

    # Use Docker executor
    mm = ModelManager(data_root=DATA_ROOT, use_api=False)

    # Get config files
    config_files = sorted(CONFIG_DIR.glob("*.yaml"))
    print(f"\nFound {len(config_files)} configs to run:")
    for cf in config_files:
        print(f"  - {cf.name}")

    # Run each config
    results = []

    for i, config_path in enumerate(config_files):
        print(f"\n{'='*70}")
        print(f"[{i+1}/{len(config_files)}] Running: {config_path.name}")
        print("=" * 70)

        # Load config
        with open(config_path) as f:
            config = yaml.safe_load(f)

        print(f"  Entities: {len(config['entities'])}")
        print(f"  Num designs: {config['num_designs']}")

        # Show what's being designed
        file_entity = config['entities'][0]['file']
        print(f"  Structure: {Path(file_entity['path']).name}")
        print(f"  Designable: {file_entity['design'][0]['chain']['res_index']}")
        print(f"  Fixed: {file_entity['not_design'][0]['chain']['res_index']}")

        # Prepare and submit job
        try:
            print("\n  Preparing and submitting job...")
            state = mm.prepare_and_submit("boltzgen", config=config, persistent=True)
            job_id = state.job_id
            print(f"  Job submitted: {job_id}")

            # Wait for completion
            print("  Waiting for completion...")
            start_time = time.time()

            while state.status not in ["completed", "failed", "error"]:
                time.sleep(10)
                state = mm.job_status(job_id)
                elapsed = time.time() - start_time
                print(f"    Status: {state.status} ({elapsed:.0f}s)")

            elapsed = time.time() - start_time
            print(f"\n  Completed in {elapsed:.1f}s with status: {state.status}")

            if state.status == "completed":
                # Ingest outputs
                print("  Ingesting outputs...")
                summary = mm.ingest_job_outputs(job_id, cleanup=False)
                print(f"  Ingested {len(summary.get('structures', []))} structures")

                results.append({
                    "config": config_path.name,
                    "job_id": job_id,
                    "status": "success",
                    "structures": summary.get("structures", []),
                    "elapsed": elapsed,
                })
            else:
                error_msg = getattr(state, 'error', 'Unknown error')
                print(f"  ERROR: {error_msg}")
                results.append({
                    "config": config_path.name,
                    "job_id": job_id,
                    "status": "failed",
                    "error": error_msg,
                })

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            results.append({
                "config": config_path.name,
                "status": "error",
                "error": str(e),
            })

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    for r in results:
        status_icon = "✓" if r["status"] == "success" else "✗"
        print(f"\n{status_icon} {r['config']}")
        if r["status"] == "success":
            print(f"  Job: {r['job_id']}")
            print(f"  Structures: {len(r['structures'])}")
            print(f"  Time: {r['elapsed']:.1f}s")
        else:
            print(f"  Error: {r.get('error', 'Unknown')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
