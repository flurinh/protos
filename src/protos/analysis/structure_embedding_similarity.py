"""Structure embedding similarity helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from protos.analysis.structure import alignment as structure_alignment

try:  # Optional torch dependency
    import torch
    import torch.nn.functional as F

    HAS_TORCH = True
except ImportError:  # pragma: no cover - environment without torch
    torch = None  # type: ignore
    F = None  # type: ignore
    HAS_TORCH = False


@dataclass
class ChainSelection:
    structure_id: str
    chain_id: str
    sequence_id: str


def _build_chain_dataframe(struct_proc, structure_id: str, chain_id: str) -> pd.DataFrame:
    df = struct_proc.load_entity(structure_id)
    if df is None:
        raise ValueError(f"Structure {structure_id} not available")

    reset = df.reset_index()
    chain_df = reset[(reset["auth_chain_id"] == chain_id) & (reset["atom_name"] == "CA")].copy()
    chain_df = chain_df.sort_values(
        ["model_num", "auth_seq_id", "insertion"], na_position="last"
    ).reset_index(drop=True)
    chain_df["seq_index"] = range(len(chain_df))
    return chain_df


def _compute_cosine_similarity(
    ref_embeddings: "torch.Tensor",
    target_embeddings: "torch.Tensor",
    ref_index: int,
    target_index: int,
) -> float:
    if not HAS_TORCH:
        raise RuntimeError("PyTorch is required for embedding similarity computation")
    ref_vec = ref_embeddings[ref_index]
    tgt_vec = target_embeddings[target_index]
    sim = F.cosine_similarity(ref_vec.unsqueeze(0), tgt_vec.unsqueeze(0), dim=1)
    return float(sim.item())


def compute_structure_embedding_similarity(
    structure_proc,
    embedding_proc,
    selection: Iterable[ChainSelection],
    *,
    reference_structure: str,
    reference_chain: str,
    embedding_dataset: str,
    window_size: int = 8,
    max_gap: int = 30,
    property_table_name: Optional[str] = None,
    property_metadata: Optional[Dict[str, Any]] = None,
    record_property_table: bool = False,
) -> Dict[str, Any]:
    """Compute per-residue embedding similarity against a reference chain."""

    if not HAS_TORCH:
        raise RuntimeError("PyTorch + torch.nn.functional required for similarity computation")

    selection_map: Dict[str, ChainSelection] = {
        sel.structure_id: sel for sel in selection
    }

    if reference_structure not in selection_map:
        raise ValueError("Reference structure must be included in selection")

    embeddings = embedding_proc.load_embeddings(embedding_dataset)

    reference_sel = selection_map[reference_structure]
    ref_sequence_id = reference_sel.sequence_id or reference_sel.structure_id
    if ref_sequence_id not in embeddings:
        raise ValueError(f"Reference sequence '{ref_sequence_id}' not found in embeddings")

    ref_embeddings = embeddings[ref_sequence_id]
    ref_chain = reference_sel.chain_id or reference_chain
    ref_df = _build_chain_dataframe(structure_proc, reference_structure, ref_chain)

    ref_plot = {
        "structure_id": reference_structure,
        "chain": ref_chain,
        "sequence": ref_sequence_id,
        "coordinates": ref_df[["x", "y", "z"]].to_dict(orient="list"),
        "values": [[] for _ in range(len(ref_df))],
    }

    records: List[Dict[str, Any]] = []
    rmsd_map: Dict[str, float] = {}
    plot_points: Dict[str, Dict[str, Any]] = {reference_structure: ref_plot}

    for structure_id, selection_entry in selection_map.items():
        if structure_id == reference_structure:
            continue

        sequence_id = selection_entry.sequence_id
        if sequence_id not in embeddings:
            raise ValueError(f"Sequence '{sequence_id}' not present in embedding dataset")

        target_embeddings = embeddings[sequence_id]
        chain_id = selection_entry.chain_id

        target_df = _build_chain_dataframe(structure_proc, structure_id, chain_id)

        ref_coords = ref_df[["x", "y", "z"]]
        target_coords = target_df[["x", "y", "z"]]

        aligned_coords, _rotation, _translation, path, rmsd = structure_alignment.align_structures(
            ref_coords,
            target_coords,
            window_size=window_size,
            max_gap=max_gap,
        )
        rmsd_map[structure_id] = float(rmsd)

        aligned_df = target_df.copy()
        aligned_df[["x", "y", "z"]] = aligned_coords

        plot_payload = {
            "structure_id": structure_id,
            "chain": chain_id,
            "sequence": sequence_id,
            "coordinates": aligned_df[["x", "y", "z"]].to_dict(orient="list"),
            "values": [[] for _ in range(len(aligned_df))],
        }
        plot_points[structure_id] = plot_payload

        ref_indices, target_indices = path
        for idx_r, idx_t in zip(ref_indices, target_indices):
            ref_row = ref_df.iloc[int(idx_r)]
            tgt_row = aligned_df.iloc[int(idx_t)]

            ref_idx = int(ref_row["seq_index"])
            tgt_idx = int(tgt_row["seq_index"])
            similarity = _compute_cosine_similarity(ref_embeddings, target_embeddings, ref_idx, tgt_idx)

            records.append(
                {
                    "reference_structure": reference_structure,
                    "reference_chain": ref_chain,
                    "reference_sequence": ref_sequence_id,
                    "reference_auth_seq_id": int(ref_row["auth_seq_id"]),
                    "reference_residue": ref_row["res_name"],
                    "target_structure": structure_id,
                    "target_chain": chain_id,
                    "target_sequence": sequence_id,
                    "target_auth_seq_id": int(tgt_row["auth_seq_id"]),
                    "target_residue": tgt_row["res_name"],
                    "cosine_similarity": similarity,
                }
            )

            ref_plot["values"][ref_idx].append(similarity)
            plot_payload["values"][tgt_idx].append(similarity)

    df = pd.DataFrame.from_records(records)

    per_structure_summary = (
        df.groupby("target_structure")["cosine_similarity"]
        .agg(["mean", "min", "max", "count"])
        .reset_index()
        .rename(columns={
            "mean": "mean_cosine",
            "min": "min_cosine",
            "max": "max_cosine",
            "count": "aligned_pairs",
        })
    )

    for payload in plot_points.values():
        payload["similarity"] = [
            sum(vals) / len(vals) if vals else 0.0
            for vals in payload["values"]
        ]
        payload.pop("values", None)

    property_table_recorded: Optional[str] = None
    if record_property_table and property_table_name:
        from protos.processing.property import PropertyProcessor

        prop_proc = PropertyProcessor()
        rows: List[Dict[str, Any]] = []

        for row in records:
            scope = [
                {"format": "structure", "name": row["reference_structure"]},
                {"format": "sequence", "name": row["reference_sequence"]},
                {"format": "structure", "name": row["target_structure"]},
                {"format": "sequence", "name": row["target_sequence"]},
            ]
            rows.append(
                {
                    "scope": scope,
                    "entity_name": f"{row['reference_structure']}->{row['target_structure']}:{row['reference_auth_seq_id']}:{row['target_auth_seq_id']}",
                    "cosine_similarity": row["cosine_similarity"],
                    "reference_residue": row["reference_residue"],
                    "target_residue": row["target_residue"],
                    "reference_auth_seq_id": row["reference_auth_seq_id"],
                    "target_auth_seq_id": row["target_auth_seq_id"],
                }
            )

        metadata = property_metadata.copy() if property_metadata else {}
        metadata.update(
            {
                "reference_structure": reference_structure,
                "reference_chain": ref_chain,
                "embedding_dataset": embedding_dataset,
            }
        )

        prop_proc.record_properties(
            property_table_name,
            rows,
            metadata=metadata,
            allow_create=True,
        )
        property_table_recorded = property_table_name

    return {
        "records": df,
        "summary": per_structure_summary,
        "rmsd": rmsd_map,
        "plot_points": plot_points,
        "property_table": property_table_recorded,
    }


__all__ = [
    "ChainSelection",
    "compute_structure_embedding_similarity",
]
