#!/usr/bin/env python3
"""Compute per-residue embedding similarities across aligned GPCR structures."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

try:
    import torch
    import torch.nn.functional as F
except ImportError:  # pragma: no cover - optional dependency
    torch = None  # type: ignore[assignment]
    F = None  # type: ignore[assignment]

import plotly.graph_objects as go

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos

from protos.processing.structure import StructureProcessor
from protos.processing.sequence import SequenceProcessor
from protos.processing.embedding import EmbeddingProcessor
from protos.processing.property import PropertyProcessor
from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine
from protos.analysis.structure import alignment as structure_alignment
from protos.io.paths import get_protos_paths

STRUCTURE_IDS = ["3sn6", "5d5a", "6b73"]
REFERENCE_STRUCTURE = "5d5a"
REFERENCE_CHAIN = "A"
REFERENCE_SEQUENCE_ID = f"{REFERENCE_STRUCTURE}_chain_{REFERENCE_CHAIN}"
SEQUENCE_THRESHOLD = 0.35
EMBEDDING_MODEL = "esm2_t12_35m"
EMBEDDING_TYPE = "per_residue"
SEQUENCE_DATASET = "gpcr_gpcr_chains"
EMBEDDING_DATASET = f"{SEQUENCE_DATASET}__{EMBEDDING_MODEL}__{EMBEDDING_TYPE}"
PROPERTY_TABLE = "gpcr_structure_embedding_similarity"


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def ensure_structures(struct_proc: StructureProcessor, structure_ids: Iterable[str]) -> None:
    try:
        existing = set(struct_proc.list_entities())
    except FileNotFoundError:
        existing = set()

    missing = [sid for sid in structure_ids if sid not in existing]
    if not missing:
        return

    from protos.io.ingest.structure_loader import StructureLoader

    loader = StructureLoader()
    loader.download_batch(missing, dataset_name="gpcr_structures")


def ensure_gpcr_dataset(struct_proc: StructureProcessor, structure_ids: Iterable[str]) -> None:
    dataset_name = "gpcr_structures"
    if not struct_proc.dataset_manager.dataset_exists(dataset_name):
        struct_proc.create_dataset(dataset_name, list(structure_ids), metadata={"source": "script"})


def ensure_chain_sequences(struct_proc: StructureProcessor, structure_ids: Iterable[str]) -> None:
    relations = struct_proc.list_dataset_related_sequences("gpcr_structures")
    missing = [sid for sid in structure_ids if not relations.get(sid)]
    if not missing:
        return
    struct_proc.register_chain_sequences(
        missing,
        dataset_prefix="gpcr_chain_dataset",
        create_dataset=True,
        overwrite=False,
    )


def materialize_chain_sequences(
    struct_proc: StructureProcessor,
    seq_proc: SequenceProcessor,
    structure_ids: Iterable[str],
) -> None:
    collected = struct_proc.collect_chain_sequences(structure_ids)
    for chains in collected.values():
        for payload in chains.values():
            seq_name = payload["entity_name"]
            sequence = payload["sequence"]
            existing = seq_proc.load_entity(seq_name)
            if isinstance(existing, str):
                continue
            seq_proc.save_entity(seq_name, sequence)


def load_structure_chain_map(struct_proc: StructureProcessor, structure_ids: Iterable[str]) -> Dict[str, List[Dict[str, object]]]:
    mapping: Dict[str, List[Dict[str, object]]] = {}
    relations = struct_proc.list_related_sequences(structure_ids, include_unloaded=True)
    for struct_id, entries in relations.items():
        chain_list: List[Dict[str, object]] = []
        for entry in entries:
            meta = entry.get("metadata", {})
            chain_list.append(
                {
                    "sequence_name": entry["name"],
                    "chain_id": meta.get("chain_id"),
                    "sequence_length": meta.get("sequence_length"),
                }
            )
        mapping[struct_id] = chain_list
    return mapping


def classify_gpcr_chains(
    struct_map: Dict[str, List[Dict[str, object]]],
    seq_proc: SequenceProcessor,
    reference_sequence_id: str,
    threshold: float,
) -> Dict[str, Dict[str, object]]:
    classifier = SequenceAlignmentEngine()
    reference_sequence = seq_proc.load_entity(reference_sequence_id)
    if not isinstance(reference_sequence, str):
        raise ValueError(f"Reference sequence '{reference_sequence_id}' not available")

    selection: Dict[str, Dict[str, object]] = {}
    for structure_id, chains in struct_map.items():
        best: Optional[Dict[str, object]] = None
        for chain in chains:
            seq_name = chain["sequence_name"]
            try:
                sequence = seq_proc.load_entity(seq_name)
            except FileNotFoundError:
                continue
            if not isinstance(sequence, str):
                continue
            result = classifier.align_pairwise(
                seq_name,
                sequence,
                reference_sequence_id,
                reference_sequence,
            )
            norm_score = result.score / max(len(sequence), 1)
            payload = {
                "sequence_name": seq_name,
                "chain_id": chain.get("chain_id"),
                "sequence": sequence,
                "score": norm_score,
                "alignment": result.alignment,
            }
            if norm_score >= threshold and (best is None or norm_score > best["score"]):
                best = payload
        if best is not None:
            selection[structure_id] = best
    return selection


def ensure_sequence_dataset(seq_proc: SequenceProcessor, sequences: Dict[str, str]) -> None:
    if not sequences:
        raise ValueError("No sequences selected for GPCR chains")

    for name, sequence in sequences.items():
        existing = seq_proc.load_entity(name)
        if isinstance(existing, str):
            continue
        seq_proc.save_entity(name, sequence)

    if seq_proc.dataset_manager.dataset_exists(SEQUENCE_DATASET):
        existing = seq_proc.load_dataset(SEQUENCE_DATASET)
        if set(existing.keys()) == set(sequences.keys()):
            return
    seq_proc.save_sequences(
        sequences,
        output_file=SEQUENCE_DATASET,
        dataset_name=SEQUENCE_DATASET,
        metadata={"kind": "gpcr_gpcr_chains"},
        materialize_entities=False,
    )


def ensure_embeddings(
    emb_proc: EmbeddingProcessor,
    sequences: Dict[str, str],
    *,
    embedding_type: str,
) -> Dict[str, torch.Tensor]:
    if torch is None or F is None:
        raise RuntimeError("PyTorch is required for embedding similarity; install torch and transformers")
    if not emb_proc.dependencies_available:
        raise RuntimeError("Embedding dependencies missing; install torch + transformers")

    if not emb_proc.dataset_manager.dataset_exists(EMBEDDING_DATASET):
        emb_proc.embed_sequences(
            sequences,
            embedding_type=embedding_type,
            save_dataset=EMBEDDING_DATASET,
            register_entities=True,
        )
    return emb_proc.load_embeddings(EMBEDDING_DATASET)


def load_embedding_entity_map(emb_proc: EmbeddingProcessor) -> Dict[str, str]:
    info = emb_proc.dataset_manager.load_dataset(EMBEDDING_DATASET)
    entities = info.get("entities", [])
    metadata = info.get("metadata", {})
    sequence_ids = metadata.get("sequence_ids", [])
    return {seq_id: entity for seq_id, entity in zip(sequence_ids, entities)}


def build_chain_dataframe(
    struct_proc: StructureProcessor,
    structure_id: str,
    chain_id: str,
) -> pd.DataFrame:
    df = struct_proc.load_entity(structure_id)
    if df is None:
        raise ValueError(f"Structure {structure_id} not available")
    reset = df.reset_index()
    chain_df = reset[(reset["auth_chain_id"] == chain_id) & (reset["atom_name"] == "CA")].copy()
    chain_df = chain_df.sort_values(["model_num", "auth_seq_id", "insertion"], na_position="last").reset_index(drop=True)
    chain_df["seq_index"] = range(len(chain_df))
    return chain_df


def compute_cosine_similarity(
    ref_embeddings: torch.Tensor,
    target_embeddings: torch.Tensor,
    ref_index: int,
    target_index: int,
) -> float:
    if torch is None or F is None:
        raise RuntimeError("PyTorch not available for cosine similarity computation")
    ref_vec = ref_embeddings[ref_index]
    tgt_vec = target_embeddings[target_index]
    sim = F.cosine_similarity(ref_vec.unsqueeze(0), tgt_vec.unsqueeze(0), dim=1)
    return float(sim.item())


def structure_alignment_records(
    struct_proc: StructureProcessor,
    embeddings: Dict[str, torch.Tensor],
    selection: Dict[str, Dict[str, object]],
    embedding_entities: Dict[str, str],
) -> Tuple[pd.DataFrame, Dict[str, float], Dict[str, Dict[str, object]]]:
    records: List[Dict[str, object]] = []
    rmsd_map: Dict[str, float] = {}
    plot_points: Dict[str, Dict[str, object]] = {}

    reference_payload = selection.get(REFERENCE_STRUCTURE)
    if reference_payload is None:
        raise ValueError(f"Reference structure {REFERENCE_STRUCTURE} not classified")

    ref_chain_id = reference_payload.get("chain_id") or REFERENCE_CHAIN
    ref_sequence_id = reference_payload["sequence_name"]
    ref_embeddings = embeddings[ref_sequence_id]
    ref_df = build_chain_dataframe(struct_proc, REFERENCE_STRUCTURE, ref_chain_id)

    ref_plot = {
        "df": ref_df[["x", "y", "z"]].copy(),
        "values": [[] for _ in range(len(ref_df))],
        "chain": ref_chain_id,
        "sequence": ref_sequence_id,
    }
    plot_points[REFERENCE_STRUCTURE] = ref_plot

    for structure_id, payload in selection.items():
        if structure_id == REFERENCE_STRUCTURE:
            continue
        chain_id = payload.get("chain_id")
        sequence_id = payload["sequence_name"]
        target_embeddings = embeddings.get(sequence_id)
        if target_embeddings is None:
            continue

        target_df = build_chain_dataframe(struct_proc, structure_id, chain_id)

        ref_coords = ref_df[["x", "y", "z"]]
        target_coords = target_df[["x", "y", "z"]]
        aligned_coords, _rotation, _translation, path, rmsd = structure_alignment.align_structures(
            ref_coords,
            target_coords,
            window_size=8,
            max_gap=30,
        )
        rmsd_map[structure_id] = rmsd

        aligned_df = target_df.copy()
        aligned_df[["x", "y", "z"]] = aligned_coords
        target_plot = {
            "df": aligned_df[["x", "y", "z"]].copy(),
            "values": [[] for _ in range(len(aligned_df))],
            "chain": chain_id,
            "sequence": sequence_id,
        }
        plot_points[structure_id] = target_plot

        ref_indices, target_indices = path
        for idx_r, idx_t in zip(ref_indices, target_indices):
            ref_row = ref_df.iloc[int(idx_r)]
            tgt_row = aligned_df.iloc[int(idx_t)]
            ref_idx = int(ref_row["seq_index"])
            tgt_idx = int(tgt_row["seq_index"])
            similarity = compute_cosine_similarity(ref_embeddings, target_embeddings, ref_idx, tgt_idx)
            records.append(
                {
                    "reference_structure": REFERENCE_STRUCTURE,
                    "reference_chain": ref_chain_id,
                    "reference_sequence": ref_sequence_id,
                    "reference_auth_seq_id": int(ref_row["auth_seq_id"]),
                    "reference_residue": ref_row["res_name"],
                    "reference_embedding": embedding_entities.get(ref_sequence_id),
                    "target_structure": structure_id,
                    "target_chain": chain_id,
                    "target_sequence": sequence_id,
                    "target_auth_seq_id": int(tgt_row["auth_seq_id"]),
                    "target_residue": tgt_row["res_name"],
                    "target_embedding": embedding_entities.get(sequence_id),
                    "cosine_similarity": similarity,
                    "alignment_index": len(records),
                }
            )
            ref_plot["values"][ref_idx].append(similarity)
            target_plot["values"][tgt_idx].append(similarity)

    df = pd.DataFrame.from_records(records)

    for structure_id, payload in plot_points.items():
        averages = [sum(vals) / len(vals) if vals else 0.0 for vals in payload["values"]]
        payload["similarity"] = averages

    return df, rmsd_map, plot_points


def register_similarity_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        raise ValueError("No similarity records computed")

    def make_scope(row: pd.Series) -> List[Dict[str, str]]:
        scope_items: List[Dict[str, str]] = []
        embed_ref = row["reference_embedding"]
        embed_tgt = row["target_embedding"]
        if embed_ref:
            scope_items.append({"format": "embedding", "name": embed_ref})
        scope_items.append({"format": "sequence", "name": row["reference_sequence"]})
        scope_items.append({"format": "structure", "name": row["reference_structure"]})
        if embed_tgt:
            scope_items.append({"format": "embedding", "name": embed_tgt})
        scope_items.append({"format": "sequence", "name": row["target_sequence"]})
        scope_items.append({"format": "structure", "name": row["target_structure"]})
        return scope_items

    payload = df.copy()
    payload["scope"] = payload.apply(make_scope, axis=1)
    payload["entity_name"] = payload.apply(
        lambda row: f"{row['reference_structure']}_{row['target_structure']}_{row['reference_auth_seq_id']}_{row['target_auth_seq_id']}",
        axis=1,
    )

    prop_proc = PropertyProcessor()
    recorded = prop_proc.record_properties(
        PROPERTY_TABLE,
        payload,
        metadata={"models": EMBEDDING_MODEL, "embedding_type": EMBEDDING_TYPE},
    )
    return recorded


def summarize(df: pd.DataFrame, rmsd_map: Dict[str, float]) -> None:
    if df.empty:
        print("No similarity data available")
        return
    summary = (
        df.groupby(["target_structure"])["cosine_similarity"].agg(["mean", "min", "max", "count"]).reset_index()
    )
    print("\nPer-structure cosine similarity summary:")
    print(summary.to_string(index=False))
    print("\nAlignment RMSD (Å):")
    for structure_id, rmsd in rmsd_map.items():
        print(f"  {structure_id}: {rmsd:.3f}")


def build_similarity_figure(plot_points: Dict[str, Dict[str, object]]) -> go.Figure:
    fig = go.Figure()
    colorscale = [[0.0, "red"], [1.0, "blue"]]
    show_colorbar = True

    for structure_id, payload in plot_points.items():
        coords = payload["df"]
        similarity = payload["similarity"]
        chain_id = payload.get("chain")
        sequence_id = payload.get("sequence")
        text = [
            f"{structure_id} chain {chain_id} | seq idx {idx} | {sequence_id} | sim {value:.3f}"
            for idx, value in enumerate(similarity)
        ]
        marker = dict(
            size=4,
            color=similarity,
            colorscale=colorscale,
            cmin=0.0,
            cmax=1.0,
            showscale=show_colorbar,
            colorbar=dict(title="Cosine similarity", len=0.5) if show_colorbar else None,
        )
        fig.add_trace(
            go.Scatter3d(
                x=coords["x"],
                y=coords["y"],
                z=coords["z"],
                mode="markers",
                marker=marker,
                name=f"{structure_id} chain {chain_id}",
                text=text,
                hoverinfo="text",
            )
        )
        show_colorbar = False

    fig.update_layout(
        title="Aligned GPCR chains colored by cosine similarity",
        scene=dict(aspectmode="data"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1.0),
    )
    return fig


def save_similarity_figure(fig: go.Figure, filename: str = "gpcr_embedding_similarity.html") -> Path:
    paths = get_protos_paths()
    embedding_root = Path(paths.get_processor_path("embedding"))
    vis_dir = embedding_root / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)
    output_path = vis_dir / filename
    fig.write_html(str(output_path))
    return output_path


def main() -> None:
    ensure_data_root()

    struct_proc = StructureProcessor()
    seq_proc = SequenceProcessor()

    ensure_structures(struct_proc, STRUCTURE_IDS)
    ensure_gpcr_dataset(struct_proc, STRUCTURE_IDS)
    ensure_chain_sequences(struct_proc, STRUCTURE_IDS)
    materialize_chain_sequences(struct_proc, seq_proc, STRUCTURE_IDS)

    struct_map = load_structure_chain_map(struct_proc, STRUCTURE_IDS)
    selection = classify_gpcr_chains(struct_map, seq_proc, REFERENCE_SEQUENCE_ID, SEQUENCE_THRESHOLD)
    if len(selection) < 2:
        raise RuntimeError("Insufficient GPCR chain selections to compute similarities")

    sequences = {payload["sequence_name"]: payload["sequence"] for payload in selection.values()}
    ensure_sequence_dataset(seq_proc, sequences)

    emb_proc = EmbeddingProcessor(model_name=EMBEDDING_MODEL)
    embeddings = ensure_embeddings(emb_proc, sequences, embedding_type=EMBEDDING_TYPE)
    embedding_entities = load_embedding_entity_map(emb_proc)

    similarity_df, rmsd_map, plot_points = structure_alignment_records(
        struct_proc,
        embeddings,
        selection,
        embedding_entities,
    )
    if similarity_df.empty:
        raise RuntimeError("No aligned residues found between structures")

    recorded = register_similarity_table(similarity_df)
    summarize(similarity_df, rmsd_map)
    print(f"\nRegistered property table '{PROPERTY_TABLE}' with {len(recorded)} rows")
    artifact = PropertyProcessor().tables_dir / f"{PROPERTY_TABLE}.csv"
    print(f"Property table written to: {artifact}")

    figure = build_similarity_figure(plot_points)
    figure_path = save_similarity_figure(figure)
    print(f"Visualization saved to: {figure_path}")


if __name__ == "__main__":
    main()
