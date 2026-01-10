#!/usr/bin/env python3
"""Verify BoltzGen YAML generation matches the expected schema."""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models import model_manager as mm


def ensure_data_root(base_dir: Path) -> mm.ModelManager:
    base_dir.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(base_dir))
    manager = mm.ModelManager(data_root=base_dir)
    manager.paths.reinitialize(wipe=True)
    return manager


def ensure_global_data_root() -> Path:
    tmp_root = Path.cwd() / "data"
    tmp_root.mkdir(parents=True, exist_ok=True)
    return tmp_root


def get_adapter() -> mm.BoltzGenAdapter:
    manager = ensure_data_root(ensure_global_data_root())
    return manager.adapters["boltzgen"]


def _dump_yaml(data: OrderedDict) -> str:
    return yaml.dump(
        data,
        Dumper=mm.BoltzYamlDumper,
        default_flow_style=False,
        sort_keys=False,
        indent=2,
    )


def run_designed_protein_check() -> None:
    """Test designed protein with length range specification."""
    adapter = get_adapter()
    config = {
        "designed_proteins": [
            {"id": "B", "min_length": 80, "max_length": 140},
        ],
    }
    yaml_data = adapter._generate_yaml(config)

    assert "entities" in yaml_data
    entities = yaml_data["entities"]
    assert len(entities) == 1

    protein = entities[0]["protein"]
    assert protein["id"] == "B"
    assert protein["sequence"] == "80..140"


def run_file_entity_check() -> None:
    """Test file entity with binding types and structure groups."""
    adapter = get_adapter()
    config = {
        "target_structure": {
            "path": "target.cif",
            "include": [
                {"chain": {"id": "A"}},
            ],
            "binding_types": {
                "A": "binder",
            },
            "structure_groups": {
                "A": "pocket",
            },
        },
    }
    yaml_data = adapter._generate_yaml(config)

    assert "entities" in yaml_data
    entities = yaml_data["entities"]
    assert len(entities) == 1

    file_entity = entities[0]["file"]
    assert file_entity["path"] == "target.cif"
    assert file_entity["include"] == [{"chain": {"id": "A"}}]
    assert file_entity["binding_types"] == {"A": "binder"}
    assert file_entity["structure_groups"] == {"A": "pocket"}


def run_combined_design_check() -> None:
    """Test combined design with protein, target structure, and ligand."""
    adapter = get_adapter()
    config = {
        "designed_proteins": [
            {"id": "B", "min_length": 50, "max_length": 100},
        ],
        "target_structure": {
            "path": "receptor.cif",
            "include": [{"chain": {"id": "A"}}],
            "binding_types": {"A": "target"},
        },
        "ligands": [
            {"id": "LIG", "ccd": "ATP"},
        ],
    }
    yaml_data = adapter._generate_yaml(config)

    assert "entities" in yaml_data
    entities = yaml_data["entities"]
    assert len(entities) == 3

    # Check protein
    protein = entities[0]["protein"]
    assert protein["id"] == "B"
    assert protein["sequence"] == "50..100"

    # Check file
    file_entity = entities[1]["file"]
    assert file_entity["path"] == "receptor.cif"

    # Check ligand
    ligand = entities[2]["ligand"]
    assert ligand["id"] == "LIG"
    assert ligand["ccd"] == "ATP"


def run_constraints_check() -> None:
    """Test constraints section for bond specifications."""
    adapter = get_adapter()
    config = {
        "designed_proteins": [
            {"id": "A", "sequence": "3..5C6C3"},
        ],
        "constraints": [
            {
                "bond": {
                    "atom1": ["A", 1, "SG"],
                    "atom2": ["A", 8, "SG"],
                }
            },
        ],
    }
    yaml_data = adapter._generate_yaml(config)

    assert "entities" in yaml_data
    assert "constraints" in yaml_data

    # Check designed protein with explicit sequence pattern
    protein = yaml_data["entities"][0]["protein"]
    assert protein["sequence"] == "3..5C6C3"

    # Check constraints
    constraints = yaml_data["constraints"]
    assert len(constraints) == 1
    bond = constraints[0]["bond"]
    assert bond["atom1"] == ["A", 1, "SG"]
    assert bond["atom2"] == ["A", 8, "SG"]


def run_explicit_entities_check() -> None:
    """Test direct entities specification (passthrough mode)."""
    adapter = get_adapter()
    explicit_entities = [
        {"protein": {"id": "B", "sequence": "80..140"}},
        {
            "file": {
                "path": "target.cif",
                "include": [{"chain": {"id": "A"}}],
                "binding_types": {"A": "binder"},
                "structure_groups": {"A": "pocket"},
                "design": {"A": [[10, 20], [30, 40]]},
                "secondary_structure": {"A": "10..20:H,30..40:E"},
            }
        },
    ]
    config = {"entities": explicit_entities}
    yaml_data = adapter._generate_yaml(config)

    assert yaml_data["entities"] == explicit_entities


def run_convenience_method_check() -> None:
    """Test the generate_design_config convenience method."""
    adapter = get_adapter()
    yaml_data = adapter.generate_design_config(
        designed_protein_id="X",
        min_length=60,
        max_length=120,
        target_structure="protein.cif",
        target_chain="A",
        ligand_ccd="FAD",
    )

    assert "entities" in yaml_data
    entities = yaml_data["entities"]

    # Check designed protein
    protein = entities[0]["protein"]
    assert protein["id"] == "X"
    assert protein["sequence"] == "60..120"

    # Check target structure
    file_entity = entities[1]["file"]
    assert file_entity["path"] == "protein.cif"
    assert file_entity["include"] == [{"chain": {"id": "A"}}]

    # Check ligand
    ligand = entities[2]["ligand"]
    assert ligand["id"] == "LIG"
    assert ligand["ccd"] == "FAD"


def run_yaml_serialization_check() -> None:
    """Verify YAML output is properly formatted."""
    adapter = get_adapter()
    config = {
        "designed_proteins": [{"id": "B", "min_length": 80, "max_length": 140}],
        "target_structure": {"path": "target.cif"},
    }
    yaml_data = adapter._generate_yaml(config)
    yaml_text = _dump_yaml(yaml_data)

    # Should start with entities
    assert yaml_text.strip().startswith("entities:")

    # Should be parseable
    parsed = yaml.safe_load(yaml_text)
    assert "entities" in parsed


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate BoltzGen YAML serialization")
    parser.add_argument(
        "--data-root",
        type=Path,
        default=REPO_ROOT / "data",
        help="Directory storing temporary Protos artifacts",
    )
    args = parser.parse_args()
    ensure_data_root(args.data_root)

    run_designed_protein_check()
    run_file_entity_check()
    run_combined_design_check()
    run_constraints_check()
    run_explicit_entities_check()
    run_convenience_method_check()
    run_yaml_serialization_check()

    print("BoltzGen YAML checks passed")
    return 0


@pytest.fixture
def boltzgen_env(tmp_path):
    return ensure_data_root(tmp_path / "protos_data")


def test_boltzgen_designed_protein(boltzgen_env):
    run_designed_protein_check()


def test_boltzgen_file_entity(boltzgen_env):
    run_file_entity_check()


def test_boltzgen_combined_design(boltzgen_env):
    run_combined_design_check()


def test_boltzgen_constraints(boltzgen_env):
    run_constraints_check()


def test_boltzgen_explicit_entities(boltzgen_env):
    run_explicit_entities_check()


def test_boltzgen_convenience_method(boltzgen_env):
    run_convenience_method_check()


def test_boltzgen_yaml_serialization(boltzgen_env):
    run_yaml_serialization_check()


if __name__ == "__main__":
    raise SystemExit(main())
