#!/usr/bin/env python3
"""Full grid search for BoltzGen rhodozyme designs.

Runs all 144 configurations (48 rotations × 3 enzymes) and aggregates results.

Usage:
    python run_boltzgen_gridsearch.py                    # Run all
    python run_boltzgen_gridsearch.py --enzyme trypsin   # Run only trypsin
    python run_boltzgen_gridsearch.py --max-jobs 10      # Run first 10
    python run_boltzgen_gridsearch.py --resume           # Resume from last run
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import yaml
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field, asdict

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

OUTPUT_DIR = THESIS_DIR / "outputs" / "boltzgen_rhodozyme"
CONFIGS_DIR = OUTPUT_DIR / "boltzgen_configs"
RESULTS_DIR = OUTPUT_DIR / "results"

import protos
from protos.models import ModelManager, JobStatus


@dataclass
class JobTracker:
    """Track job progress and results."""
    config_name: str
    enzyme: str
    rotation: str
    config_path: str
    status: str = "pending"
    job_id: Optional[str] = None
    start_time: Optional[str] = None
    end_time: Optional[str] = None
    duration_seconds: float = 0.0
    exit_code: Optional[int] = None
    num_designs: int = 0
    output_dir: Optional[str] = None
    error: Optional[str] = None
    metrics: Dict[str, Any] = field(default_factory=dict)


def load_progress(progress_file: Path) -> Dict[str, JobTracker]:
    """Load progress from previous run."""
    if not progress_file.exists():
        return {}

    with open(progress_file) as f:
        data = json.load(f)

    trackers = {}
    for name, info in data.get("jobs", {}).items():
        trackers[name] = JobTracker(**info)
    return trackers


def save_progress(progress_file: Path, trackers: Dict[str, JobTracker], summary: Dict[str, Any]):
    """Save progress to file."""
    data = {
        "summary": summary,
        "jobs": {name: asdict(t) for name, t in trackers.items()}
    }
    with open(progress_file, "w") as f:
        json.dump(data, f, indent=2)


def collect_metrics(output_dir: Path) -> Dict[str, Any]:
    """Collect key metrics from BoltzGen output."""
    metrics = {}

    # Load final metrics CSV
    metrics_file = output_dir / "predictions" / "final_ranked_designs" / "final_designs_metrics_30.csv"
    if metrics_file.exists():
        import pandas as pd
        df = pd.read_csv(metrics_file)
        if len(df) > 0:
            # Get best design metrics
            best = df.iloc[0]
            metrics["best_rank"] = int(best.get("final_rank", 1))
            metrics["best_iptm"] = float(best.get("design_to_target_iptm", 0))
            metrics["best_ptm"] = float(best.get("design_ptm", 0))
            metrics["best_rmsd"] = float(best.get("filter_rmsd", 999))
            metrics["num_designs"] = len(df)

            # Average metrics
            metrics["avg_iptm"] = float(df["design_to_target_iptm"].mean())
            metrics["avg_ptm"] = float(df["design_ptm"].mean())
            metrics["avg_rmsd"] = float(df["filter_rmsd"].mean())

    return metrics


def run_job(mm: ModelManager, config_name: str, config: Dict, tracker: JobTracker) -> JobTracker:
    """Run a single BoltzGen job."""
    tracker.start_time = datetime.now().isoformat()
    tracker.status = "running"

    try:
        # Prepare
        invocation = mm.prepare(
            "boltzgen",
            config=config,
            metadata={"job_name": config_name}
        )

        # Submit
        job_state = mm.submit_job(invocation, persistent=True)
        tracker.job_id = job_state.job_id
        tracker.output_dir = str(invocation.job.working_dir)

        # Wait for completion
        while True:
            state = mm.job_status(job_state.job_id)
            if state.status in (JobStatus.COMPLETED, JobStatus.FAILED, JobStatus.CANCELLED):
                break
            time.sleep(10)

        # Get result
        result = mm.job_result(job_state.job_id)
        tracker.end_time = datetime.now().isoformat()
        tracker.exit_code = result.exit_code if result else -1
        tracker.duration_seconds = result.duration_seconds if result else 0

        if state.status == JobStatus.COMPLETED:
            tracker.status = "completed"
            # Collect metrics
            tracker.metrics = collect_metrics(Path(tracker.output_dir))
            tracker.num_designs = tracker.metrics.get("num_designs", 0)

            # Ingest outputs to register designs in data/structure/mmcif/
            try:
                summary = mm.ingest_job_outputs(job_state.job_id, cleanup=False)
                tracker.metrics["ingested_count"] = len(summary.get("ingested", []))
            except Exception as e:
                tracker.metrics["ingestion_error"] = str(e)
        else:
            tracker.status = "failed"
            tracker.error = state.error

    except Exception as e:
        tracker.status = "failed"
        tracker.error = str(e)
        tracker.end_time = datetime.now().isoformat()

    return tracker


def main() -> int:
    parser = argparse.ArgumentParser(description="Run BoltzGen grid search")
    parser.add_argument("--enzyme", choices=["trypsin", "subtilisin", "papain"],
                        help="Only run specific enzyme")
    parser.add_argument("--max-jobs", type=int, default=None,
                        help="Maximum number of jobs to run")
    parser.add_argument("--resume", action="store_true",
                        help="Resume from previous run")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be run without executing")
    args = parser.parse_args()

    print("=" * 70)
    print("BOLTZGEN GRID SEARCH - RHODOZYME DESIGN")
    print("=" * 70)

    # Initialize - must pass data_root explicitly to ModelManager
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Data root: {data_root}")

    progress_file = RESULTS_DIR / "gridsearch_progress.json"

    # Load configs
    config_files = sorted(CONFIGS_DIR.glob("*.yaml"))
    print(f"\nTotal configurations: {len(config_files)}")

    # Filter by enzyme
    if args.enzyme:
        config_files = [f for f in config_files if f.stem.startswith(args.enzyme)]
        print(f"Filtered to {args.enzyme}: {len(config_files)}")

    # Limit jobs
    if args.max_jobs:
        config_files = config_files[:args.max_jobs]
        print(f"Limited to: {len(config_files)}")

    # Load progress
    if args.resume:
        trackers = load_progress(progress_file)
        print(f"Resuming: {sum(1 for t in trackers.values() if t.status == 'completed')} completed")
    else:
        trackers = {}

    # Build job list
    jobs_to_run = []
    for config_file in config_files:
        name = config_file.stem

        # Skip if already completed
        if name in trackers and trackers[name].status == "completed":
            continue

        # Parse enzyme and rotation from name (e.g., "trypsin_rot000_000_000")
        parts = name.split("_")
        enzyme = parts[0]
        rotation = "_".join(parts[1:])

        # Load config
        with open(config_file) as f:
            config = yaml.safe_load(f)

        # Create tracker
        tracker = JobTracker(
            config_name=name,
            enzyme=enzyme,
            rotation=rotation,
            config_path=str(config_file),
        )
        trackers[name] = tracker
        jobs_to_run.append((name, config, tracker))

    print(f"Jobs to run: {len(jobs_to_run)}")

    if args.dry_run:
        print("\nDry run - jobs that would be executed:")
        for name, _, _ in jobs_to_run[:10]:
            print(f"  {name}")
        if len(jobs_to_run) > 10:
            print(f"  ... and {len(jobs_to_run) - 10} more")
        return 0

    if not jobs_to_run:
        print("\nNo jobs to run - all completed!")
        return 0

    # Initialize ModelManager with explicit data_root
    mm = ModelManager(data_root=data_root, use_api=False)
    print(f"\nExecutor: {type(mm._executor).__name__}")
    print(f"GPU enabled: {mm._executor.use_gpu}")
    print(f"Models dir: {mm.paths.data_root}/models")

    # Estimate time
    est_time_per_job = 150  # seconds based on test run
    est_total = len(jobs_to_run) * est_time_per_job / 3600
    print(f"\nEstimated time: {est_total:.1f} hours ({len(jobs_to_run)} jobs × ~{est_time_per_job}s)")

    # Run jobs
    print("\n" + "-" * 70)
    print("RUNNING JOBS")
    print("-" * 70)

    start_time = time.time()
    completed = 0
    failed = 0

    for i, (name, config, tracker) in enumerate(jobs_to_run):
        print(f"\n[{i+1}/{len(jobs_to_run)}] {name}")

        tracker = run_job(mm, name, config, tracker)
        trackers[name] = tracker

        if tracker.status == "completed":
            completed += 1
            metrics = tracker.metrics
            print(f"  ✓ Completed in {tracker.duration_seconds:.0f}s")
            print(f"    Designs: {tracker.num_designs}")
            if metrics:
                print(f"    Best: RMSD={metrics.get('best_rmsd', 'N/A'):.2f}Å, "
                      f"iPTM={metrics.get('best_iptm', 'N/A'):.3f}, "
                      f"pTM={metrics.get('best_ptm', 'N/A'):.3f}")
        else:
            failed += 1
            print(f"  ✗ Failed: {tracker.error}")

        # Save progress after each job
        elapsed = time.time() - start_time
        remaining = (len(jobs_to_run) - i - 1) * (elapsed / (i + 1))
        summary = {
            "started": datetime.fromtimestamp(start_time).isoformat(),
            "elapsed_hours": elapsed / 3600,
            "remaining_hours": remaining / 3600,
            "total_jobs": len(trackers),
            "completed": completed,
            "failed": failed,
            "pending": len(jobs_to_run) - i - 1,
        }
        save_progress(progress_file, trackers, summary)

    # Final summary
    total_time = time.time() - start_time
    print("\n" + "=" * 70)
    print("GRID SEARCH COMPLETE")
    print("=" * 70)
    print(f"Total time: {total_time/3600:.2f} hours")
    print(f"Completed: {completed}")
    print(f"Failed: {failed}")
    print(f"Results: {progress_file}")

    # Aggregate best results
    if completed > 0:
        print("\nTOP 10 DESIGNS (by iPTM):")
        results = [(name, t) for name, t in trackers.items()
                   if t.status == "completed" and t.metrics]
        results.sort(key=lambda x: x[1].metrics.get("best_iptm", 0), reverse=True)

        for i, (name, t) in enumerate(results[:10]):
            m = t.metrics
            print(f"  {i+1}. {name}: iPTM={m.get('best_iptm', 0):.3f}, "
                  f"pTM={m.get('best_ptm', 0):.3f}, RMSD={m.get('best_rmsd', 999):.2f}Å")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
