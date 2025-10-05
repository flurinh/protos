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

from typing import Optional


def register_gpcr_sequence_dataset(dataset_name: str = "gpcr_agonist_antagonist_sequences") -> str:
    """Ensure the GPCR agonist/antagonist sequence dataset is registered.

    Args:
        dataset_name: Optional override for the dataset identifier.

    Returns:
        The dataset name once registration succeeds.
    """

    raise NotImplementedError("register_gpcr_sequence_dataset is not implemented yet")


def register_rhodopsin_structure_dataset(dataset_name: str = "rhodopsin_states") -> str:
    """Ensure the rhodopsin state structure dataset is registered."""

    raise NotImplementedError("register_rhodopsin_structure_dataset is not implemented yet")


def register_chembl_ligand_dataset(dataset_name: str = "P24941_chembl_ligands") -> str:
    """Ensure the ChEMBL ligand dataset for GPCR benchmarking is registered."""

    raise NotImplementedError("register_chembl_ligand_dataset is not implemented yet")


def register_gpcr_property_dataset(dataset_name: str = "gpcr_ligand_binding_analysis") -> str:
    """Ensure the GPCR ligand binding property dataset is registered."""

    raise NotImplementedError("register_gpcr_property_dataset is not implemented yet")


def register_rhodopsin_graph_dataset(dataset_name: str = "rhodopsin_states_residue_graphs") -> str:
    """Ensure the rhodopsin residue-level graph dataset is registered."""

    raise NotImplementedError("register_rhodopsin_graph_dataset is not implemented yet")
