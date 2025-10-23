"""Central place for registering reference datasets with the Protos registry.

Functions defined here should encapsulate the zero-configuration data-management
pattern: callers supply only the target data root via ``protos.set_data_path``
and then invoke one of these helpers to ensure the dataset is present. The
helpers return the canonical dataset name registered with the appropriate
processor.

Implementation is currently pending; the stubs raise ``NotImplementedError`` so
that tests capturing the desired behaviour can be written first.
"""
from __future__ import annotations

import csv
import json
import pickle
import shutil
from importlib import resources
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from protos.io.formats.fasta_utils import read_fasta
from protos.io.ingest.sequence_loader import SequenceLoader
from protos.io.ingest.structure_loader import StructureLoader
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor
from protos.processing.molecule import MoleculeProcessor
from protos.processing.property import PropertyProcessor
from protos.processing.graph import GraphProcessor


def _entity_exists(processor, name: str) -> bool:
    return processor.entity_registry.find_entity(name, processor.processor_type) is not None


def register_gpcr_sequence_dataset(dataset_name: str = "gpcr_agonist_antagonist_sequences") -> str:
    """Ensure the GPCR agonist/antagonist sequence dataset is registered."""

    processor = SequenceProcessor()
    if processor.dataset_manager.dataset_exists(dataset_name):
        return dataset_name

    fasta_pkg = "protos.data.sequence.fasta"
    fasta_resource = "gpcr_agonist_antagonist_sequences.fasta"

    with resources.as_file(resources.files(fasta_pkg).joinpath(fasta_resource)) as fasta_path:
        sequences = read_fasta(str(fasta_path))

    loader = SequenceLoader(processor=processor)
    records = [
        {
            "name": name,
            "sequence": sequence,
            "metadata": {
                "source": "package:protos.data.sequence.fasta",
                "resource": fasta_resource,
            },
        }
        for name, sequence in sequences.items()
    ]

    loader.register_sequence_records(
        records,
        dataset_name=dataset_name,
        dataset_metadata={
            "source": "protos.datasets.register_gpcr_sequence_dataset",
            "resource": fasta_resource,
        },
        overwrite=False,
    )

    return dataset_name


def register_rhodopsin_structure_dataset(dataset_name: str = "rhodopsin_states") -> str:
    """Ensure the rhodopsin state structure dataset is registered."""

    processor = StructureProcessor()
    existing = set(processor.dataset_manager.get_dataset_entities(dataset_name)) if processor.dataset_manager.dataset_exists(dataset_name) else set()

    pdb_ids = ["1U19", "1F88", "3PQR", "3PXO", "2I37", "6CMO", "4ZWJ"]

    loader = StructureLoader(processor=processor)
    mmcif_pkg = "protos.data.structure.mmcif"

    for pdb_id in pdb_ids:
        if pdb_id in existing or _entity_exists(processor, pdb_id):
            continue

        resource_name = f"{pdb_id.lower()}.cif"
        try:
            resource_path = resources.files(mmcif_pkg).joinpath(resource_name)
        except FileNotFoundError as exc:  # pragma: no cover - safety
            raise FileNotFoundError(f"Packaged mmCIF resource missing: {resource_name}") from exc

        with resources.as_file(resource_path) as src_path:
            temp_dir = Path(processor.path_temp_dir)
            temp_dir.mkdir(parents=True, exist_ok=True)
            temp_target = temp_dir / resource_name
            shutil.copyfile(src_path, temp_target)
            try:
                loader.download_and_register(
                    str(temp_target),
                    name=pdb_id,
                    source="local",
                    metadata={
                        "source": "package:protos.data.structure.mmcif",
                        "resource": resource_name,
                    },
                )
            finally:
                temp_target.unlink(missing_ok=True)

    available_entities = [pdb for pdb in pdb_ids if _entity_exists(processor, pdb)]
    if not available_entities:
        raise RuntimeError(
            "Failed to register any rhodopsin structures; ensure packaged resources are available."
        )

    processor.create_dataset(
        dataset_name,
        available_entities,
        {
            "source": "protos.datasets.register_rhodopsin_structure_dataset",
            "description": "Rhodopsin state reference structures packaged with Protos",
        },
    )

    return dataset_name


def register_chembl_ligand_dataset(dataset_name: str = "P24941_chembl_ligands") -> str:
    """Ensure the ChEMBL ligand dataset for GPCR benchmarking is registered."""

    processor = MoleculeProcessor()
    if processor.dataset_manager.dataset_exists(dataset_name):
        return dataset_name

    dataset_pkg = "protos.data.ligand.datasets"
    dataset_resource = "P24941_chembl_ligands.json"

    with resources.as_file(resources.files(dataset_pkg).joinpath(dataset_resource)) as dataset_path:
        dataset_info = json.loads(Path(dataset_path).read_text(encoding="utf-8"))

    smiles_list: List[str] = dataset_info.get("entities", [])
    metadata = dataset_info.get("metadata", {}).copy()
    metadata.setdefault("source", "protos.datasets.register_chembl_ligand_dataset")

    for smiles in smiles_list:
        if _entity_exists(processor, smiles):
            continue
        processor.save_entity(
            smiles,
            {"smiles": smiles},
            metadata={
                "source": "package:protos.data.ligand.datasets",
                "dataset": dataset_resource,
            },
        )

    processor.create_dataset(dataset_name, smiles_list, metadata)

    return dataset_name


def register_gpcr_property_dataset(dataset_name: str = "gpcr_ligand_binding_analysis") -> str:
    """Ensure the GPCR ligand binding property dataset is registered."""

    processor = PropertyProcessor()
    if processor.dataset_manager.dataset_exists(dataset_name):
        return dataset_name

    dataset_pkg = "protos.data.property.datasets"
    dataset_resource = "gpcr_ligand_binding_analysis.json"

    with resources.as_file(resources.files(dataset_pkg).joinpath(dataset_resource)) as dataset_path:
        dataset_info = json.loads(Path(dataset_path).read_text(encoding="utf-8"))

    table_rel = dataset_info.get("table_file")
    if not table_rel:
        raise ValueError("Property dataset metadata missing 'table_file'")

    table_name = Path(table_rel).name
    tables_pkg = "protos.data.property.tables"

    with resources.as_file(resources.files(tables_pkg).joinpath(table_name)) as table_path:
        with open(table_path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)

    property_rows: List[Dict[str, Any]] = []
    for row in rows:
        entity_id = row.pop("entity_id")
        structure_id = entity_id.split("_")[0]
        property_rows.append(
            {
                **row,
                "entity_name": entity_id,
                "scope": [
                    {
                        "format": "structure",
                        "name": structure_id,
                    }
                ],
            }
        )

    processor.record_properties(
        dataset_name,
        property_rows,
        metadata={
            "source": "protos.datasets.register_gpcr_property_dataset",
            **dataset_info.get("metadata", {}),
        },
        allow_create=True,
    )

    return dataset_name


def register_rhodopsin_graph_dataset(dataset_name: str = "rhodopsin_states_residue_graphs") -> str:
    """Ensure the rhodopsin residue-level graph dataset is registered."""

    processor = GraphProcessor()
    if processor.dataset_manager.dataset_exists(dataset_name):
        return dataset_name

    dataset_pkg = "protos.data.graph.datasets"
    dataset_resource = "rhodopsin_states_residue_graphs.json"

    with resources.as_file(resources.files(dataset_pkg).joinpath(dataset_resource)) as dataset_path:
        dataset_info = json.loads(Path(dataset_path).read_text(encoding="utf-8"))

    graph_names: List[str] = dataset_info.get("entities", [])
    metadata = dataset_info.get("metadata", {}).copy()
    metadata.setdefault("source", "protos.datasets.register_rhodopsin_graph_dataset")

    graphs_pkg = "protos.data.graph.graphs"
    for graph_name in graph_names:
        if _entity_exists(processor, graph_name):
            continue
        resource_name = f"{graph_name}.pkl"
        with resources.as_file(resources.files(graphs_pkg).joinpath(resource_name)) as graph_path:
            with open(graph_path, "rb") as handle:
                graph_payload = pickle.load(handle)
        processor.save_entity(
            graph_name,
            graph_payload,
            metadata=metadata,
        )

    processor.create_dataset(dataset_name, graph_names, metadata)

    return dataset_name
