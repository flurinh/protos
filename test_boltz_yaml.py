#!/usr/bin/env python3
"""Verify Boltz YAML generation matches the expected schema."""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Tuple

import pytest
import yaml

SRC_DIR = Path(__file__).resolve().parent / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models import model_manager as mm


def ensure_data_root(base_dir: Path) -> Tuple[mm.ModelManager, mm.BoltzAdapter]:
    base_dir.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(base_dir))
    manager = mm.ModelManager(data_root=base_dir)
    manager.paths.reinitialize(wipe=True)
    adapter = manager.adapters["boltz2"]
    return manager, adapter


def ensure_global_data_root() -> Path:
    tmp_root = Path.cwd() / "data"
    tmp_root.mkdir(parents=True, exist_ok=True)
    return tmp_root


def generate_yaml(samples: OrderedDict[str, str], config: dict | None = None) -> OrderedDict:
    _, adapter = ensure_data_root(ensure_global_data_root())
    config = config or {}
    return adapter._generate_yaml(samples, config)


def _dump_yaml(data: OrderedDict) -> str:
    return yaml.dump(
        data,
        Dumper=mm.BoltzYamlDumper,
        default_flow_style=False,
        sort_keys=False,
        indent=2,
    )


def run_default_structure_check() -> None:
    samples = OrderedDict([
        ("SEQ_ALPHA", "AAA"),
        ("SEQ_BETA", "BBB"),
    ])
    yaml_data = generate_yaml(samples)
    assert list(yaml_data.keys()) == ["sequences"]
    assert yaml_data["sequences"] == [
        {"protein": OrderedDict([("id", ["A1"]), ("sequence", "AAA")])},
        {"protein": OrderedDict([("id", ["B1"]), ("sequence", "BBB")])},
    ]


def run_override_check() -> None:
    samples = OrderedDict([("SEQ_ALPHA", "AAA")])
    config = {
        "sequence_overrides": {
            "SEQ_ALPHA": {"id": ["A3"], "fields": {"chain": "A"}},
        },
        "extra_sequences": [
            {"ligand": OrderedDict([("id", ["L1"]), ("ccd", "EKY")])},
        ],
        "constraints": [
            {"pocket": OrderedDict([("binder", "L1"), ("contacts", [["A3", 42]])])},
        ],
        "properties": [{"affinity": OrderedDict([("binder", "L1")])}],
    }
    yaml_data = generate_yaml(samples, config)
    assert list(yaml_data.keys()) == ["sequences", "constraints", "properties"]
    seq_entry = yaml_data["sequences"][0]["protein"]
    assert seq_entry["id"] == ["A3"]
    assert seq_entry["chain"] == "A"
    ligand_entry = yaml_data["sequences"][1]["ligand"]
    assert ligand_entry["id"] == ["L1"]
    assert ligand_entry["ccd"] == "EKY"
    contacts = yaml_data["constraints"][0]["pocket"]["contacts"]
    assert contacts == [["A3", 42]]


def run_reference_template_check() -> None:
    reference_path = (
        Path(__file__).resolve().parent / "src/protos/reference_data/models/boltz/example.yaml"
    )
    reference_text = reference_path.read_text()
    reference_data = yaml.safe_load(reference_text)
    yaml_data = generate_yaml(
        OrderedDict(),
        config={
            "sequences": reference_data["sequences"],
            "constraints": reference_data.get("constraints"),
            "properties": reference_data.get("properties"),
        },
    )
    generated_text = _dump_yaml(yaml_data)
    assert yaml.safe_load(generated_text) == reference_data
    assert generated_text.strip().startswith("sequences:")


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Boltz YAML serialization")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path(__file__).resolve().parent / "data",
        help="Directory storing temporary Protos artifacts",
    )
    args = parser.parse_args()
    ensure_data_root(args.data_root)
    run_default_structure_check()
    run_override_check()
    run_reference_template_check()
    print("Boltz YAML checks passed")
    return 0


@pytest.fixture
def boltz_env(tmp_path):
    return ensure_data_root(tmp_path / "protos_data")


def test_boltz_yaml_default_structure(boltz_env):
    run_default_structure_check()


def test_boltz_yaml_supports_overrides_and_extras(boltz_env):
    run_override_check()


def test_boltz_yaml_matches_reference_template(boltz_env):
    run_reference_template_check()


if __name__ == "__main__":
    raise SystemExit(main())
