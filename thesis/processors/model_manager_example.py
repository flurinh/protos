#!/usr/bin/env python3
"""ModelManager Example: Structure prediction job preparation.

Demonstrates ProtOS capabilities:
- SequenceProcessor: Load and manage sequence data
- ModelManager: Prepare jobs for external compute resources
- Zero-configuration job packaging (YAML, FASTA, metadata)

Question: "How do I prepare Boltz2 structure predictions for cluster submission?"

This example shows how ProtOS wraps processor data into compute-ready jobs
that can be submitted to cloud, cluster, or containerized environments.
The jobs are fully self-contained with all inputs and configurations.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Setup paths
THESIS_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = THESIS_DIR.parent
if (REPO_ROOT / "src").exists():
    sys.path.insert(0, str(REPO_ROOT / "src"))

import protos
from protos.processing.sequence import SequenceProcessor
from protos.io.ingest.sequence_loader import SequenceLoader
from protos.models.model_manager import ModelManager


# =============================================================================
# Configuration
# =============================================================================
OUTPUT_DIR = THESIS_DIR / "outputs" / "model_manager"
FIGURES_DIR = THESIS_DIR / "processors" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

DATASET_NAME = "rhodopsin_prediction_targets"

# Target sequences for structure prediction
PREDICTION_TARGETS = {
    "RHO_BOVIN": "P02699",    # Bovine rhodopsin (has crystal structure for validation)
    "RHO_HUMAN": "P08100",    # Human rhodopsin
    "OPN4_HUMAN": "Q9UHM6",   # Melanopsin (no crystal structure)
}

# Mutation variants to predict
MUTATIONS = [
    {"position": 296, "wild_type": "K", "mutant": "A", "name": "K296A"},  # Retinal binding site
    {"position": 113, "wild_type": "E", "mutant": "Q", "name": "E113Q"},  # Counterion
]


def main() -> int:
    """Run the ModelManager example."""
    print("=" * 70)
    print("MODEL MANAGER EXAMPLE")
    print("Preparing Boltz2 Structure Prediction Jobs")
    print("=" * 70)

    # Initialize ProtOS
    data_root = REPO_ROOT / "data"
    protos.set_data_path(str(data_root))

    # -------------------------------------------------------------------------
    # Step 1: Load sequences using ProtOS SequenceProcessor
    # -------------------------------------------------------------------------
    print("\n[1] Loading target sequences...")
    seq_proc = SequenceProcessor()
    loader = SequenceLoader(processor=seq_proc)

    sequences = {}
    for name, uniprot_id in PREDICTION_TARGETS.items():
        loader.download_and_register(
            f"uniprot:{uniprot_id}",
            name=name,
            materialize_entities=True,
        )
        seq = seq_proc.load_entity(name)
        if seq:
            sequences[name] = seq
            print(f"  {name}: {len(seq)} aa")

    # Register as dataset (required for ModelManager)
    seq_proc.save_sequences(
        sequences,
        output_file=DATASET_NAME,
        dataset_name=DATASET_NAME,
    )
    print(f"  Dataset registered: {DATASET_NAME}")

    # -------------------------------------------------------------------------
    # Step 2: Initialize ModelManager
    # -------------------------------------------------------------------------
    print("\n[2] Initializing ModelManager...")
    manager = ModelManager(data_root=data_root)

    # List available models
    print("  Available models:")
    for model_name, card in manager.cards.items():
        print(f"    - {model_name}: {card.description}")

    # -------------------------------------------------------------------------
    # Step 3: Prepare wild-type prediction jobs
    # -------------------------------------------------------------------------
    print("\n[3] Preparing wild-type structure prediction jobs...")

    wild_type_jobs = []
    for entity_name in sequences.keys():
        try:
            # ModelManager wraps sequence data into Boltz2 job
            invocation = manager.prepare_input(
                "boltz2",
                entity_name=entity_name,
                dataset_name=DATASET_NAME,
                config={
                    "use_msa_server": True,
                    "recycling_steps": 3,
                    "sampling_steps": 200,
                },
                metadata={"variant": "wild_type"},
            )

            if invocation.job:
                job = invocation.job
                wild_type_jobs.append({
                    "entity": entity_name,
                    "working_dir": str(job.working_dir),
                    "command": " ".join(job.command),
                    "artifacts": [a.spec.name for a in job.artifacts],
                })

                print(f"\n  {entity_name}:")
                print(f"    Working dir: {job.working_dir}")
                print(f"    Command: {' '.join(job.command[:4])}...")
                print(f"    Artifacts: {len(job.artifacts)}")

                # Show generated files
                if job.working_dir.exists():
                    for f in job.working_dir.iterdir():
                        print(f"      - {f.name}")

        except Exception as e:
            print(f"  {entity_name}: Failed - {e}")

    # -------------------------------------------------------------------------
    # Step 4: Prepare mutant prediction jobs
    # -------------------------------------------------------------------------
    print("\n[4] Preparing mutant structure prediction jobs...")

    mutant_jobs = []
    reference = "RHO_BOVIN"

    for mutation in MUTATIONS:
        try:
            invocation = manager.prepare_input(
                "boltz2",
                entity_name=reference,
                dataset_name=DATASET_NAME,
                config={
                    "mutations": [mutation],
                    "use_msa_server": True,
                    "recycling_steps": 3,
                },
                metadata={"variant": mutation["name"]},
            )

            if invocation.job:
                job = invocation.job
                mutant_jobs.append({
                    "entity": reference,
                    "mutation": mutation["name"],
                    "working_dir": str(job.working_dir),
                    "command": " ".join(job.command),
                })

                print(f"\n  {reference} {mutation['name']}:")
                print(f"    Working dir: {job.working_dir}")
                print(f"    Mutation: {mutation['wild_type']}{mutation['position']}{mutation['mutant']}")

        except Exception as e:
            print(f"  {mutation['name']}: Failed - {e}")

    # -------------------------------------------------------------------------
    # Step 5: Prepare batch submission manifest
    # -------------------------------------------------------------------------
    print("\n[5] Creating batch submission manifest...")

    all_jobs = wild_type_jobs + mutant_jobs
    manifest = {
        "batch_name": "rhodopsin_structure_predictions",
        "model": "boltz2",
        "total_jobs": len(all_jobs),
        "dataset": DATASET_NAME,
        "jobs": all_jobs,
        "execution": {
            "container": "ghcr.io/jwohlwend/boltz:0.4.0",
            "gpu_required": True,
            "estimated_time_per_job": "30-60 min",
            "recommended_resources": {
                "gpu": "A100 or V100",
                "memory": "32GB",
                "storage": "50GB per job",
            },
        },
        "submission_examples": {
            "local": "boltz predict config.yaml --out_dir predictions",
            "slurm": "sbatch --gres=gpu:1 --mem=32G run_boltz.sh",
            "docker": "docker run --gpus all -v $PWD:/data boltz predict /data/config.yaml",
        },
    }

    manifest_path = OUTPUT_DIR / "batch_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  Saved: {manifest_path}")

    # -------------------------------------------------------------------------
    # Step 6: Show job directory structure
    # -------------------------------------------------------------------------
    print("\n[6] Job directory structure (example)...")

    if wild_type_jobs:
        example_dir = Path(wild_type_jobs[0]["working_dir"])
        if example_dir.exists():
            print(f"\n  {example_dir.name}/")
            for f in sorted(example_dir.iterdir()):
                size = f.stat().st_size if f.is_file() else 0
                print(f"    {f.name:25s} ({size:,} bytes)")

            # Show config.yaml content
            config_file = example_dir / "config.yaml"
            if config_file.exists():
                print(f"\n  config.yaml preview:")
                with open(config_file) as f:
                    lines = f.readlines()[:15]
                    for line in lines:
                        print(f"    {line.rstrip()}")
                    if len(lines) == 15:
                        print("    ...")

    # -------------------------------------------------------------------------
    # Step 7: Create summary visualization
    # -------------------------------------------------------------------------
    print("\n[7] Creating job summary visualization...")

    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("Jobs by Type", "Job Workflow"),
        specs=[[{"type": "pie"}, {"type": "sankey"}]],
    )

    # Pie chart of job types
    fig.add_trace(go.Pie(
        labels=["Wild-type", "Mutants"],
        values=[len(wild_type_jobs), len(mutant_jobs)],
        marker_colors=["#1f77b4", "#d62728"],
        hole=0.4,
    ), row=1, col=1)

    # Sankey diagram of workflow
    fig.add_trace(go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            label=["Sequences", "ModelManager", "Boltz2 Jobs", "Cluster/Cloud"],
            color=["#2ecc71", "#3498db", "#9b59b6", "#e74c3c"],
        ),
        link=dict(
            source=[0, 1, 2],
            target=[1, 2, 3],
            value=[len(sequences), len(all_jobs), len(all_jobs)],
            color=["rgba(46,204,113,0.4)", "rgba(52,152,219,0.4)", "rgba(155,89,182,0.4)"],
        ),
    ), row=1, col=2)

    fig.update_layout(
        title="Boltz2 Job Preparation Summary",
        height=400,
        width=900,
        paper_bgcolor="white",
    )
    fig.write_image(str(FIGURES_DIR / "model_manager_jobs.png"), scale=2)
    print(f"  Saved: {FIGURES_DIR / 'model_manager_jobs.png'}")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("COMPLETE - Jobs Ready for External Submission")
    print("=" * 70)
    print(f"\nSummary:")
    print(f"  Sequences loaded: {len(sequences)}")
    print(f"  Wild-type jobs: {len(wild_type_jobs)}")
    print(f"  Mutant jobs: {len(mutant_jobs)}")
    print(f"  Total jobs: {len(all_jobs)}")
    print(f"\nJob artifacts location: {data_root / 'models' / 'boltz2'}")
    print(f"Batch manifest: {manifest_path}")
    print(f"\nNext steps:")
    print(f"  1. Transfer job directories to compute cluster")
    print(f"  2. Run: boltz predict config.yaml --out_dir predictions")
    print(f"  3. Or use container: docker run --gpus all boltz ...")
    print(f"\nOutputs: {OUTPUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
