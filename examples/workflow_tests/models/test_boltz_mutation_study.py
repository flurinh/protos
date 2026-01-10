#!/usr/bin/env python3
"""Exercise Boltz mutation batch preparation end-to-end."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager
from protos.processing.sequence import SequenceProcessor


def ensure_data_root(base_dir: Path) -> tuple[ModelManager, SequenceProcessor]:
    base_dir.mkdir(parents=True, exist_ok=True)
    manager = ModelManager(data_root=base_dir)
    manager.paths.reinitialize(wipe=True)
    protos.set_data_path(str(base_dir))
    return manager, SequenceProcessor()


def prepare_dataset(seq_proc: SequenceProcessor, name: str, sequences: Dict[str, str]) -> None:
    seq_proc.save_sequences(
        sequences,
        output_file=name,
        dataset_name=name,
        materialize_entities=True,
    )


def run_mutation_batch(
    manager: ModelManager,
    dataset_name: str,
    mutations: List[Dict[str, object]],
) -> List[Path]:
    invocations = manager.prepare_boltz_mutations(dataset_name, mutations)
    config_paths: List[Path] = []

    for invocation, expected_name in zip(invocations, ("SEQ_ALPHA", "SEQ_BETA")):
        job = invocation.job
        if job is None:
            raise AssertionError("Boltz invocation missing job payload")

        config_bundle = next(
            bundle for bundle in job.artifacts if bundle.spec.name == "boltz_config"
        )
        fasta_bundle = next(
            bundle for bundle in job.artifacts if bundle.spec.name == "boltz_sequences"
        )

        yaml_data = yaml.safe_load(config_bundle.path.read_text())
        assert "sequences" in yaml_data

        fasta_content = fasta_bundle.path.read_text()
        assert any(
            line.startswith(f">{expected_name}") for line in fasta_content.splitlines()
        )

        for mutation in invocation.metadata["mutation_entry"]["mutations"]:
            position = mutation["position"]
            mutant = mutation["mutant"]
            seq_entry = yaml_data["sequences"][0]["protein"]["sequence"]
            assert seq_entry[position - 1] == mutant

        config_paths.append(config_bundle.path)

    return config_paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate Boltz configs for mutation screens")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=REPO_ROOT / "data",
        help="Directory to store generated configs",
    )
    args = parser.parse_args()

    manager, seq_proc = ensure_data_root(args.data_root)

    sequences = {
        "SEQ_ALPHA": "MKTAYIAKQRQISFVKSHFSRQDILDLWQ"[:15],
        "SEQ_BETA": "GAVLVTGIVLDSGDGVTHVVPIYEGY"[:15],
    }
    dataset_name = "test_mutations"
    prepare_dataset(seq_proc, dataset_name, sequences)

    mutations = [
        {
            "entity": "SEQ_ALPHA",
            "mutations": [
                {"position": 3, "mutant": "L", "original": "T", "name": "T3L"}
            ],
        },
        {
            "entity": "SEQ_BETA",
            "mutations": [
                {"position": 5, "mutant": "R", "original": "V", "name": "V5R"}
            ],
            "config": {"default_sequence_type": "protein"},
        },
    ]

    config_paths = run_mutation_batch(manager, dataset_name, mutations)

    print("Generated Boltz configs:")
    for path in config_paths:
        print(f"  - {path}")

    return 0


@pytest.fixture
def boltz_env(tmp_path):
    base = tmp_path / "protos_data"
    manager, seq_proc = ensure_data_root(base)
    return manager, seq_proc


def test_prepare_boltz_mutations(boltz_env):
    manager, seq_proc = boltz_env

    sequences = {
        "SEQ_ALPHA": "MKTAYIAKQRQISFVKSHFSRQDILDLWQ"[:15],
        "SEQ_BETA": "GAVLVTGIVLDSGDGVTHVVPIYEGY"[:15],
    }

    dataset_name = "test_mutations"
    prepare_dataset(seq_proc, dataset_name, sequences)

    mutations = [
        {
            "entity": "SEQ_ALPHA",
            "mutations": [
                {"position": 3, "mutant": "L", "original": "T", "name": "T3L"}
            ],
        },
        {
            "entity": "SEQ_BETA",
            "mutations": [
                {"position": 5, "mutant": "R", "original": "V", "name": "V5R"}
            ],
            "config": {"default_sequence_type": "protein"},
        },
    ]

    config_paths = run_mutation_batch(manager, dataset_name, mutations)
    assert len(config_paths) == 2


if __name__ == "__main__":
    raise SystemExit(main())
