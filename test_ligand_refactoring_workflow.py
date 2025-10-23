#!/usr/bin/env python3
"""Ligand ingest/export workflow using real downloads via chembl_loader.

All path management via ProtosPaths; no hardcoded examples.
"""

from __future__ import annotations

import sys
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.processing.structure import StructureProcessor
from protos.processing.molecule import MoleculeProcessor
from protos.io.ingest.ligand_loader import LigandLoader
from protos.io.ingest.structure_loader import StructureLoader
from protos.io.export.ligand_exporter import LigandExporter
from protos.io.ingest import chembl_loader
from protos.io.formats.structure_schema import STRUCT_COLUMN_DTYPE
from protos.models.model_manager import ModelManager

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover - skip if rdkit missing
    Chem = None


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def _fetch_top3_smiles_for_protein(protein_id: str) -> dict[str, str]:
    """Use ChEMBL to fetch 3 ligands for a protein (no hardcoded examples)."""
    results = chembl_loader.query_protein_ligands(
        protein_id,
        activity_types=None,
        min_pchembl=5.0,
        max_value_nm=10000,
        limit=3,
    )
    smiles_map: dict[str, str] = {}
    for idx, rec in enumerate(results, start=1):
        name = rec.get('chembl_id') or f"ligand_{idx}"
        smi = rec.get('smiles')
        if smi:
            smiles_map[name] = smi
    return smiles_map


def main() -> None:
    if Chem is None:
        print("RDKit is not available; skipping ligand refactoring check.")
        return

    data_root = ensure_data_root()
    print(f"Data root initialised at {data_root}")

    structure_proc = StructureProcessor()
    molecule_proc = MoleculeProcessor()
    ligand_loader = LigandLoader(ligand_processor=molecule_proc)
    structure_loader = StructureLoader()

    # Get protein ID from environment/CLI; fall back to a sane default for examples
    protein_id = (sys.argv[1] if len(sys.argv) > 1 else None) or os.environ.get("PROTOS_PROTEIN_ID")
    if not protein_id:
        protein_id = os.environ.get("PROTOS_DEFAULT_PROTEIN_ID", "P24941")  # EGFR_HUMAN as default example
        print(f"Using default protein ID for example run: {protein_id}")

    smiles_map = _fetch_top3_smiles_for_protein(protein_id)
    if not smiles_map:
        print(f"No ligands returned for protein '{protein_id}'.")
        return

    dataset_name = f"chembl_{protein_id}_ligands"
    dataset_id, entities = ligand_loader.import_smiles(smiles_map, dataset_name=dataset_name)
    print(f"Registered dataset: {dataset_id} -> {entities}")
    molecule_dataset_entities = molecule_proc.dataset_manager.get_dataset_entities(f"{dataset_id}_molecules")
    print(f"Molecule dataset entities: {molecule_dataset_entities}")

    if len(entities) != 3:
        raise AssertionError(f"Expected 3 structure entities, got {len(entities)}")
    if len(molecule_dataset_entities) != 3:
        raise AssertionError(
            f"Expected 3 molecule descriptors, got {len(molecule_dataset_entities)}"
        )

    for entity in entities:
        df = structure_proc.load_entity(entity)
        assert df is not None
        df = df.reset_index()
        for column in STRUCT_COLUMN_DTYPE:
            if column in ('phi', 'psi', 'omega', 'grn', 'res_name1l'):
                continue
            if column not in df.columns:
                raise AssertionError(f"Missing column {column} in structure frame")

        pkl_path = structure_proc.path_pkl_dir / f"{entity}.pkl"
        if not pkl_path.exists():
            raise AssertionError(f"Canonical cache missing for {entity}: {pkl_path}")
    print("Canonical frame validation succeeded.")

    # Quick ligand summaries for each entity
    for ent in entities:
        summary = structure_proc.summarize_ligands(ent)
        print(f"Ligand summary for {ent}: groups={summary['total_groups']}")

    # Export dataset to a multi-mol SDF (one record per ligand)
    ds_export = structure_proc.export_dataset(
        dataset_id,
        format="sdf",
        overwrite=True,
        ligand_group_by='res_id',
        min_atoms=1,
    )
    ds_sdf = Path(ds_export.get("dataset_file"))
    if not ds_sdf.exists():
        raise AssertionError(f"Dataset SDF export missing: {ds_sdf}")
    print(f"Dataset SDF export: {ds_sdf}")

    # Re-import the exported SDF to verify round-trip
    re_ds_name = f"{dataset_id}_reimport"
    re_dataset_id, re_entities = ligand_loader.import_sdf(str(ds_sdf), dataset_name=re_ds_name)
    print(f"Re-imported dataset: {re_dataset_id} -> count={len(re_entities)}")
    if len(re_entities) == 0:
        raise AssertionError("Re-import produced no entities")

    # Protein (PDB) with ligands; export ligands only to SDF (use default example if not provided)
    # Default to a GPCR example if none provided
    pdb_id = os.environ.get("PROTOS_PDB_ID") or os.environ.get("PROTOS_DEFAULT_PDB_ID", "5D5A")
    reg_name = structure_loader.download_and_register(pdb_id)
    if reg_name:
        # Print ligand information (compact) and choose main ligand by atom count
        summary = structure_proc.summarize_ligands(reg_name)
        print(f"Ligand info for {reg_name}: total_groups={summary['total_groups']}")
        # Flatten groups
        candidates = []
        for chain, groups in summary.get('chains', {}).items():
            for resid, info in groups.items():
                comp = info.get('comp_id')
                atomc = info.get('atom_count', 0)
                candidates.append((resid, comp, atomc, chain))
        # Show a compact overview
        preview = ", ".join([f"{rid}({cmp},{cnt})" for rid, cmp, cnt, _ in candidates[:10]])
        print(f"Ligand groups (first 10): {preview}")

        # Exclude common solvents/ions when picking the main ligand
        exclude_comp = {
            'HOH', 'WAT', 'DOD', 'SO4', 'PO4', 'GOL', 'PEG', 'MPD', 'ACT', 'BME',
            'NA', 'K', 'CL', 'CA', 'MG', 'ZN', 'MN', 'CU', 'FE',
        }
        filtered = [(rid, cmp, cnt) for rid, cmp, cnt, _ in candidates if (cmp or '').upper() not in exclude_comp]
        choice = None
        if filtered:
            choice = max(filtered, key=lambda x: x[2])
        elif candidates:
            choice = max([(rid, cmp, cnt) for rid, cmp, cnt, _ in candidates], key=lambda x: x[2])

        if not choice:
            print(f"No ligand groups found to export for {reg_name}.")
        else:
            chosen_resid, chosen_comp, chosen_atoms = choice
            print(f"Exporting main ligand: res_id={chosen_resid} comp_id={chosen_comp} atoms={chosen_atoms}")
            lig_sdf = structure_proc.export_entity(
                reg_name,
                format="sdf",
                overwrite=True,
                include_res_ids=[chosen_resid],
            )
            if not Path(lig_sdf).exists():
                raise AssertionError(f"Ligand SDF export missing for PDB {pdb_id}: {lig_sdf}")
            print(f"Exported main ligand for {pdb_id} to {lig_sdf}")

            # Prepare Uni-Dock job using zero-config paths (no hardcoded examples)
            try:
                receptor_pdb = structure_proc.export_entity(reg_name, format="pdb", overwrite=True)
                manager = ModelManager()
                inv = manager.prepare(
                    "unidock",
                    inputs={
                        "receptor_pdb": str(receptor_pdb),
                        "ligand_file": str(lig_sdf),
                    },
                    config={
                        # Follow Uni-Dock docs: valid search_mode and scoring
                        "search_mode": "fast",
                        "num_modes": 1,
                        "scoring": "vina",
                    },
                )
                if inv.job is None:
                    raise AssertionError("Uni-Dock prepare() did not produce a job")
                print("Uni-Dock command:")
                print(" ".join(inv.job.command))
            except Exception as e:  # pragma: no cover - environment dependent
                print(f"Uni-Dock prepare failed: {e}")


if __name__ == "__main__":
    main()
