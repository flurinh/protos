#!/usr/bin/env python3
"""Align two GPCR structures on retinal (RET) and visualize the result."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import plotly.graph_objects as go

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos

from protos.processing.structure import StructureProcessor
from protos.io.ingest.structure_loader import StructureLoader
from protos.analysis.structure import alignment as structure_alignment
from protos.io.paths import get_protos_paths

PDB_IDS = ["4FBZ", "1F88"]
REFERENCE_ID = "4FBZ"
TARGET_ID = "1F88"
LIGAND_RESNAME = "RET"
OUTPUT_HTML = "ret_alignment.html"


def ensure_data_root() -> Path:
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def download_structures(structure_ids: list[str]) -> None:
    loader = StructureLoader()
    try:
        loader.download_batch(structure_ids, dataset_name="ret_alignment_structures")
    except Exception as exc:  # noqa: BLE001
        print(f"Warning: issue downloading structures: {exc}")


def load_structure_frames(struct_proc: StructureProcessor, structure_ids: list[str]) -> Dict[str, pd.DataFrame]:
    frames: Dict[str, pd.DataFrame] = {}
    for pdb_id in structure_ids:
        df = struct_proc.load_entity(pdb_id)
        if df is None:
            raise RuntimeError(f"Failed to load structure {pdb_id}")
        frames[pdb_id] = df.reset_index()
    return frames


def select_ligand_subset(df, resname: str):
    ligand = df[(df["group"] == "HETATM") & (df["res_name"] == resname)]
    if ligand.empty:
        raise RuntimeError(f"No {resname} ligand found in structure")
    # choose the first residue occurrence (based on res_id)
    first_residue = ligand["res_id"].iloc[0]
    subset = ligand[ligand["res_id"] == first_residue].copy()
    if subset.empty:
        raise RuntimeError(f"Failed to isolate residue {first_residue}")
    return subset


def align_on_ligand(reference_df, target_df):
    ref_subset = select_ligand_subset(reference_df, LIGAND_RESNAME)
    tgt_subset = select_ligand_subset(target_df, LIGAND_RESNAME)

    ref_coords = ref_subset[["x", "y", "z"]]
    tgt_coords = tgt_subset[["x", "y", "z"]]

    aligned_coords_df, rotation, translation, path, rmsd = structure_alignment.align_structures(
        ref_coords,
        tgt_coords,
        window_size=8,
        max_gap=30,
    )
    return {
        "aligned_subset": aligned_coords_df,
        "rotation": rotation,
        "translation": translation,
        "path": path,
        "rmsd": rmsd,
        "reference_subset": ref_subset,
        "target_subset": tgt_subset,
    }


def apply_transform(df: pd.DataFrame, rotation, translation) -> pd.DataFrame:
    rotation = np.asarray(rotation, dtype=float)
    translation = np.asarray(translation, dtype=float).reshape(1, 3)
    coords = df[["x", "y", "z"]].to_numpy(dtype=float)
    rotated = coords @ rotation
    translated = rotated + translation
    result = df.copy()
    result[["x", "y", "z"]] = translated
    return result


def build_plot(reference_df, reference_ligand, aligned_df, aligned_ligand, output_path: Path) -> None:
    fig = go.Figure()

    def add_structure(df, name, color):
        ca = df[df["atom_name"] == "CA"]
        fig.add_trace(
            go.Scatter3d(
                x=ca["x"],
                y=ca["y"],
                z=ca["z"],
                mode="markers",
                marker=dict(size=3, color=color, opacity=0.4),
                name=name,
                hoverinfo="skip",
            )
        )

    def add_ligand(df, name, color):
        fig.add_trace(
            go.Scatter3d(
                x=df["x"],
                y=df["y"],
                z=df["z"],
                mode="markers",
                marker=dict(size=6, color=color),
                name=name,
                text=df["atom_name"],
            )
        )

    add_structure(reference_df, f"{REFERENCE_ID} (CA)", "#1f77b4")
    add_ligand(reference_ligand, f"{REFERENCE_ID} RET", "#ff7f0e")

    add_structure(aligned_df, f"{TARGET_ID} aligned (CA)", "#2ca02c")
    add_ligand(aligned_ligand, f"{TARGET_ID} aligned RET", "#d62728")

    fig.update_layout(
        title=f"Alignment of {TARGET_ID} onto {REFERENCE_ID} using {LIGAND_RESNAME}",
        scene=dict(aspectmode="data"),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1.0),
    )

    fig.write_html(str(output_path))
    print(f"Visualization saved to {output_path}")


def main() -> None:
    data_root = ensure_data_root()
    print(f"Data root configured at {data_root}")

    download_structures(PDB_IDS)

    struct_proc = StructureProcessor()
    frames = load_structure_frames(struct_proc, PDB_IDS)

    reference_df = frames[REFERENCE_ID]
    target_df = frames[TARGET_ID]

    alignment_result = align_on_ligand(reference_df, target_df)
    rotation = alignment_result["rotation"]
    translation = alignment_result["translation"]
    rmsd = alignment_result["rmsd"]
    print(f"Alignment RMSD (RET): {rmsd:.3f} Å")

    aligned_target_df = apply_transform(target_df, rotation, translation)
    aligned_ligand_df = apply_transform(alignment_result["target_subset"], rotation, translation)

    paths = get_protos_paths()
    vis_dir = Path(paths.get_processor_path("structure")) / "alignments"
    vis_dir.mkdir(parents=True, exist_ok=True)
    output_path = vis_dir / OUTPUT_HTML

    build_plot(
        reference_df,
        alignment_result["reference_subset"],
        aligned_target_df,
        aligned_ligand_df,
        output_path,
    )


if __name__ == "__main__":
    main()
