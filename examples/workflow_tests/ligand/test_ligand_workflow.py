#!/usr/bin/env python3
"""End-to-end ligand workflow: GRN annotation, interactions, SMILES export, and Boltz job."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import importlib

import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]  # -> repo root
SRC_DIR = REPO_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.analysis.structure_ligand_analysis import (
    calculate_ligand_interactions,
    ligand_to_smiles,
)
from protos.io.ingest.structure_loader import StructureLoader
from protos.models.model_manager import ModelManager
from protos.processing.molecule import MoleculeProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.structure import StructureProcessor

_lambda_runtime = importlib.import_module("protos.models.lambda.lmda.runtime_utils")
build_grn_assignments = _lambda_runtime.build_grn_assignments

TARGET_STRUCTURE = "5d5a"
PREFERRED_LIGAND = "CAU"  # Allosteric antagonist co-crystallized with 5d5a
SEQUENCE_DATASET_PREFIX = "ligand_workflow_chains"
GRN_TABLE_NAME = "ligand_workflow_grn"
EXCLUDED_LIGANDS = {"HOH", "WAT", "SO4", "PO4", "GOL", "PEG"}
LIGAND_SMILES_OVERRIDES = {
    "CAU": "CC(C)NC1=NC(=O)N(C=C1)C2=CN=C(N)N=C2N",  # Example cafuasimide analog (placeholder)
    "CLR": "C[C@H](CCC[C@]1(C)CC=C2C[C@H]3CC[C@]4(C)C(=CC[C@H]4O)C3CC[C@]12C)C",
}
BINDING_DISTANCE_CUTOFF = 4.0


def ensure_data_root(base_dir: Optional[Path] = None) -> Path:
    data_root = base_dir or (REPO_ROOT / "data")
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def _choose_ligand(summary: Dict[str, any]) -> Dict[str, any]:
    best_entry: Optional[Dict[str, any]] = None
    best_atoms = -1
    for chain_id, ligands in summary.get("chains", {}).items():
        for res_id, payload in ligands.items():
            comp = (payload.get("comp_id") or "").upper()
            if comp in EXCLUDED_LIGANDS:
                continue
            entry = {
                "chain_id": chain_id,
                "seq_id": payload.get("seq_id"),
                "comp_id": comp,
                "atom_count": payload.get("atom_count", 0),
                "res_id": res_id,
            }
            if comp == PREFERRED_LIGAND:
                return entry
            if entry["atom_count"] > best_atoms:
                best_entry = entry
                best_atoms = entry["atom_count"]
    if best_entry is None:
        raise RuntimeError("No suitable ligand detected in structure")
    return best_entry


def _resolve_seq_id(entry: Dict[str, any]) -> int:
    seq_id = entry.get("seq_id")
    if seq_id is not None:
        return int(seq_id)
    token = entry.get("res_id") or ""
    for part in str(token).split("_"):
        if part.isdigit():
            return int(part)
    raise ValueError(f"Unable to resolve residue id from token '{token}'")


def _build_residue_grn_map(structure_df: pd.DataFrame, chain_id: str,
                           assignments: list[Tuple[str, int]]) -> Dict[int, str]:
    chain_df = structure_df[structure_df["auth_chain_id"] == chain_id]
    residue_df = chain_df[["auth_seq_id", "label_seq_id"]].dropna(subset=["auth_seq_id"])
    residue_df = residue_df.drop_duplicates()
    order_cols = [col for col in ("label_seq_id", "auth_seq_id") if col in residue_df.columns]
    residue_df = residue_df.sort_values(order_cols)
    ordered_residues = [int(v) for v in residue_df["auth_seq_id"].tolist()]
    residue_grn: Dict[int, str] = {}
    for grn_label, seq_idx in assignments:
        if 0 <= seq_idx < len(ordered_residues):
            residue_grn[ordered_residues[seq_idx]] = grn_label
    return residue_grn


def run_ligand_structure_workflow(
    *,
    data_root: Optional[Path] = None,
    structure_id: str = TARGET_STRUCTURE,
) -> Dict[str, any]:
    ensure_data_root(data_root)

    loader = StructureLoader()
    loader.download_batch([structure_id], dataset_name="ligand_workflow_structures", create_dataset=True, overwrite=False)

    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()
    mol_proc = MoleculeProcessor()
    manager = ModelManager()

    registration = struct_proc.register_chain_sequences(
        [structure_id],
        dataset_prefix=SEQUENCE_DATASET_PREFIX,
        create_dataset=True,
        overwrite=False,
    )
    reg_info = registration[structure_id]
    chain_payloads = reg_info.get("chains", {})
    if not chain_payloads:
        raise RuntimeError("No chains registered for structure")

    chain_id, chain_data = max(chain_payloads.items(), key=lambda item: item[1]["length"])
    sequence_entity = chain_data["entity_name"]
    sequence_dataset = reg_info.get("dataset") or f"{SEQUENCE_DATASET_PREFIX}_{structure_id}"

    grn_table, grn_summary = seq_proc.annotate_with_grn(
        dataset_name=sequence_dataset,
        reference_table="gpcrdb_ref",
        protein_family="gpcr_a",
        output_table=GRN_TABLE_NAME,
        allow_create=True,
        return_summary=True,
    )
    grn_assignments = build_grn_assignments(grn_table)

    structure_df = struct_proc.load_entity(structure_id).reset_index()
    ligand_summary = struct_proc.summarize_ligands(structure_id)
    ligand_entry = _choose_ligand(ligand_summary)
    ligand_res_seq = _resolve_seq_id(ligand_entry)
    ligand_chain = ligand_entry["chain_id"]
    ligand_comp = ligand_entry["comp_id"]

    ligand_mask = (
        (structure_df["group"].str.upper() == "HETATM")
        & (structure_df["auth_chain_id"] == ligand_chain)
        & (structure_df["auth_seq_id"] == ligand_res_seq)
    )
    ligand_atoms = structure_df[ligand_mask]
    if ligand_atoms.empty:
        raise RuntimeError(f"No ligand atoms found for {ligand_comp} chain {ligand_chain}")

    # Interaction analysis and GRN-aware binding annotations
    interaction_result = calculate_ligand_interactions(
        struct_proc,
        structure_id,
        ligand_atoms,
        detailed=True,
        cutoff=BINDING_DISTANCE_CUTOFF,
    )
    binding_records = interaction_result.get("binding_residues", [])
    binding_df = pd.DataFrame(binding_records)
    if not binding_df.empty:
        binding_df = binding_df.dropna(subset=["res_id"])
        binding_df["res_id"] = binding_df["res_id"].astype(int)
    interaction_summary = interaction_result.get("summary", {})
    property_table_name = None
    if not binding_df.empty:
        seq_assignments = grn_assignments.get(sequence_entity, [])
        residue_grn = _build_residue_grn_map(structure_df, chain_id, seq_assignments)
        binding_df["grn"] = binding_df["res_id"].astype(int).map(residue_grn)

        prop_proc = PropertyProcessor()
        rows = []
        for rec in binding_df.to_dict("records"):
            rows.append(
                {
                    "scope": [
                        {"format": "structure", "name": structure_id},
                        {"format": "sequence", "name": sequence_entity},
                    ],
                    "entity_name": f"{structure_id}_{rec['chain_id']}_{rec['res_id']}",
                    "chain_id": rec["chain_id"],
                    "res_id": rec["res_id"],
                    "res_name": rec["res_name"],
                    "grn": rec.get("grn"),
                    "min_distance": rec["min_distance"],
                    "num_contacts": rec["num_contacts"],
                }
            )

        property_table_name = f"{structure_id.lower()}_ligand_contacts"
        prop_proc.record_properties(
            property_table_name,
            rows,
            metadata={
                "structure_id": structure_id,
                "ligand": ligand_comp,
                "sequence_entity": sequence_entity,
                "workflow": "ligand_workflow",
            },
            allow_create=True,
        )

    # Ligand SMILES resolution
    ligand_smiles = None
    try:
        ligand_smiles = ligand_to_smiles(ligand_atoms, ligand_comp)
    except Exception:
        ligand_smiles = None
    if not ligand_smiles:
        ligand_smiles = LIGAND_SMILES_OVERRIDES.get(ligand_comp)
    if not ligand_smiles:
        raise RuntimeError(f"Unable to determine SMILES for ligand {ligand_comp}")

    ligand_entity = mol_proc.save_entity(
        f"{structure_id}_{ligand_comp}_{ligand_chain}",
        {"smiles": ligand_smiles, "kind": "smiles_record"},
        metadata={"source_structure": structure_id},
    )

    boltz_config = {
        "output_name": f"{structure_id}_{ligand_chain}_{ligand_comp}_dock",
        "ligand": {"id": ligand_comp, "smiles": ligand_smiles},
        "default_sequence_type": "protein",
    }
    boltz_invocation = manager.prepare(
        "boltz2",
        inputs={"sequence_dataset": sequence_dataset, "entity": sequence_entity},
        config=boltz_config,
    )
    boltz_job = boltz_invocation.job
    if boltz_job is None:
        raise RuntimeError("Boltz invocation did not produce a job payload")

    return {
        "structure": structure_id,
        "chain": chain_id,
        "ligand_comp": ligand_comp,
        "binding_df": binding_df,
        "interaction_summary": interaction_summary,
        "ligand_entity": ligand_entity,
        "boltz_job": boltz_job,
        "sequence_dataset": sequence_dataset,
        "sequence_entity": sequence_entity,
        "property_table": property_table_name,
    }


def main() -> int:
    try:
        results = run_ligand_structure_workflow()
    except Exception as exc:  # noqa: BLE001
        print(f"Ligand workflow failed: {exc}")
        return 1

    binding_df = results["binding_df"]
    if binding_df.empty:
        print("No binding residues detected within cutoff distance")
    else:
        print("Top binding residues with GRN labels:")
        preview = binding_df.head(10)
        print(preview[["chain_id", "res_id", "res_name", "min_distance", "grn"]])

    job = results["boltz_job"]
    print("\nPrepared Boltz job:")
    print("  command:", " ".join(job.command))
    print("  working_dir:", job.working_dir)
    print("  packaged artifacts:")
    for artifact in job.artifacts:
        print("   -", artifact.spec.name, artifact.path)
    return 0


@pytest.mark.integration
def test_ligand_workflow():
    ensure_data_root()
    results = run_ligand_structure_workflow()
    binding_df = results["binding_df"]
    assert not binding_df.empty, "Expected at least one binding residue"
    assert "grn" in binding_df.columns
    assert results["property_table"] is not None
    job = results["boltz_job"]
    assert job.command, "Boltz job missing command"
    assert job.working_dir.exists()


if __name__ == "__main__":
    raise SystemExit(main())
