#!/usr/bin/env python3
"""Test script for BoltzGen Docker job execution.

This script demonstrates the full workflow:
1. Use ModelManager to prepare a BoltzGen job
2. Submit to Docker executor
3. Monitor job status
4. Retrieve results

Usage:
    python scripts/test_boltzgen_docker.py [--build] [--run]
"""

import argparse
import subprocess
import sys
from pathlib import Path

# Add protos to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from protos.models import ModelManager, JobStatus


def build_boltzgen_image(download_weights: bool = False) -> bool:
    """Build the BoltzGen Docker image."""
    boltzgen_dir = Path(__file__).parent.parent / "src/protos/models/boltzgen"

    print(f"Building BoltzGen Docker image from {boltzgen_dir}")

    cmd = [
        "docker", "build",
        "-t", "protos/boltzgen:latest",
        "--build-arg", f"DOWNLOAD_WEIGHTS={'true' if download_weights else 'false'}",
        str(boltzgen_dir),
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(boltzgen_dir))
    return result.returncode == 0


def check_image_exists() -> bool:
    """Check if BoltzGen image exists."""
    result = subprocess.run(
        ["docker", "images", "-q", "protos/boltzgen:latest"],
        capture_output=True,
        text=True,
    )
    return bool(result.stdout.strip())


def prepare_test_job(mm: ModelManager) -> dict:
    """Prepare a simple BoltzGen test job."""

    # Simple de novo protein design config
    config = {
        "job_name": "test_boltzgen_denovo",
        "entities": [
            {
                "protein": {
                    "id": "A",
                    "sequence": "40..60",  # Design a small protein (40-60 residues)
                }
            }
        ],
        # Fast settings for testing
        "recycling_steps": 1,
        "sampling_steps": 10,
        "diffusion_samples": 2,
        "num_designs": 2,  # Small number for testing
        # devices defaults to all available GPUs
    }

    print("\n=== Preparing BoltzGen Job ===")
    print(f"Config: {config}")

    invocation = mm.prepare("boltzgen", config=config)

    print(f"\nJob prepared:")
    print(f"  Model: {invocation.model}")
    print(f"  Is external: {invocation.is_external()}")

    if invocation.job:
        print(f"  Working dir: {invocation.job.working_dir}")
        print(f"  Command: {' '.join(invocation.job.command)}")
        print(f"  Output dir: {invocation.job.metadata.get('output_dir')}")

    return invocation


def run_test_job(mm: ModelManager, invocation) -> None:
    """Submit and monitor the test job."""

    print("\n=== Submitting Job ===")
    state = mm.submit_job(invocation)

    print(f"Job submitted: {state.job_id}")
    print(f"Status: {state.status.value}")
    print(f"Executor: {state.executor}")

    print("\n=== Monitoring Job ===")
    print("(Press Ctrl+C to stop monitoring)")

    try:
        final_state = mm.wait_for_job(
            state.job_id,
            timeout_seconds=3600,  # 1 hour
            poll_interval=5.0,
        )

        print(f"\nJob completed!")
        print(f"Final status: {final_state.status.value}")

        if final_state.result:
            print(f"Exit code: {final_state.result.exit_code}")
            print(f"Duration: {final_state.result.duration_seconds:.1f}s")
            print(f"Output files: {len(final_state.result.output_files)}")

            for f in final_state.result.output_files[:10]:
                print(f"  - {f}")

            if final_state.result.stderr:
                print(f"\nStderr (last 500 chars):")
                print(final_state.result.stderr[-500:])

        if final_state.error:
            print(f"\nError: {final_state.error}")

    except KeyboardInterrupt:
        print("\n\nMonitoring stopped. Job continues in background.")
        print(f"Check status with: mm.job_status('{state.job_id}')")


def list_jobs(mm: ModelManager) -> None:
    """List all jobs."""
    print("\n=== Job History ===")
    jobs = mm.list_jobs()

    if not jobs:
        print("No jobs found.")
        return

    for job in jobs[:10]:
        print(f"  {job.job_id}: {job.status.value} ({job.model})")


def main():
    parser = argparse.ArgumentParser(description="Test BoltzGen Docker workflow")
    parser.add_argument("--build", action="store_true", help="Build Docker image")
    parser.add_argument("--build-with-weights", action="store_true",
                       help="Build Docker image with model weights")
    parser.add_argument("--run", action="store_true", help="Run test job")
    parser.add_argument("--prepare-only", action="store_true",
                       help="Only prepare job, don't submit")
    parser.add_argument("--list", action="store_true", help="List existing jobs")
    parser.add_argument("--status", type=str, help="Check status of job ID")
    args = parser.parse_args()

    # Build image if requested
    if args.build or args.build_with_weights:
        success = build_boltzgen_image(download_weights=args.build_with_weights)
        if not success:
            print("Failed to build Docker image")
            sys.exit(1)
        print("Docker image built successfully!")
        if not args.run:
            return

    # Check if image exists
    if not check_image_exists() and (args.run or args.prepare_only):
        print("BoltzGen Docker image not found!")
        print("Run with --build to build the image first:")
        print("  python scripts/test_boltzgen_docker.py --build")
        sys.exit(1)

    # Initialize ModelManager
    mm = ModelManager()

    # List jobs
    if args.list:
        list_jobs(mm)
        return

    # Check specific job status
    if args.status:
        try:
            state = mm.job_status(args.status)
            print(f"\nJob: {state.job_id}")
            print(f"Model: {state.model}")
            print(f"Status: {state.status.value}")
            print(f"Created: {state.created_at}")
            if state.completed_at:
                print(f"Completed: {state.completed_at}")
            if state.error:
                print(f"Error: {state.error}")
            if state.result:
                print(f"Exit code: {state.result.exit_code}")
                print(f"Output files: {len(state.result.output_files)}")
        except KeyError:
            print(f"Job not found: {args.status}")
        return

    # Prepare job
    if args.run or args.prepare_only:
        invocation = prepare_test_job(mm)

        if args.prepare_only:
            print("\nJob prepared but not submitted (--prepare-only)")
            return

        # Run job
        run_test_job(mm, invocation)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
