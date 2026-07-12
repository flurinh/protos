"""Interactive GRN-aware C-alpha structure visualization."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import qualitative

__all__ = ["plot_grn_ca_structure", "write_grn_ca_html"]


def _ca_residues(
    structure: pd.DataFrame,
    chains: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Select one ordered C-alpha record per polymer residue."""

    df = structure.reset_index() if isinstance(structure.index, pd.MultiIndex) else structure.copy()
    atom_col = "atom_name" if "atom_name" in df.columns else "res_atom_name"
    ca = df[df[atom_col].astype(str).str.upper() == "CA"].copy()
    if chains is not None:
        allowed = {str(chain) for chain in chains}
        ca = ca[ca["auth_chain_id"].astype(str).isin(allowed)]
    if ca.empty:
        raise ValueError("Structure contains no C-alpha atoms for the selected chains")

    if "model_num" in ca.columns and ca["model_num"].notna().any():
        ca = ca[ca["model_num"] == ca["model_num"].dropna().min()]

    ordered = []
    for _, chain in ca.groupby("auth_chain_id", sort=False):
        if "label_seq_id" in chain.columns and chain["label_seq_id"].notna().any():
            chain = chain[chain["label_seq_id"].notna()].sort_values(
                ["label_seq_id", "atom_id"]
            )
            chain = chain.drop_duplicates("label_seq_id", keep="first")
        else:
            chain["_insertion_order"] = chain.get(
                "insertion", pd.Series("", index=chain.index)
            ).fillna("").astype(str)
            chain = chain.sort_values(["auth_seq_id", "_insertion_order", "atom_id"])
            chain = chain.drop_duplicates(
                ["auth_seq_id", "_insertion_order"], keep="first"
            )
        ordered.append(chain)
    return pd.concat(ordered, ignore_index=True)


def plot_grn_ca_structure(
    structure: pd.DataFrame,
    *,
    structure_id: Optional[str] = None,
    chains: Optional[Iterable[str]] = None,
    max_ca_distance: Optional[float] = 4.5,
    marker_size: float = 4.0,
    edge_width: float = 3.0,
    show_unassigned: bool = True,
) -> go.Figure:
    """Plot C-alpha atoms with GRNs in hover text and sequential CA edges.

    Edges connect consecutive polymer residues within a chain.  By default an
    edge is omitted when its C-alpha distance exceeds 4.5 Å, preventing a line
    from bridging unresolved structure gaps.  Set ``max_ca_distance=None`` to
    connect every consecutive observed C-alpha pair.
    """

    ca = _ca_residues(structure, chains)
    if "grn" not in ca.columns:
        ca["grn"] = ""
    ca["grn"] = ca["grn"].fillna("").astype(str)
    if not show_unassigned:
        ca = ca[ca["grn"] != ""]
    if ca.empty:
        raise ValueError("No C-alpha atoms remain after GRN filtering")

    if structure_id is None and "structure_id" in ca.columns:
        ids = ca["structure_id"].dropna().astype(str).unique()
        structure_id = ids[0] if len(ids) == 1 else None

    figure = go.Figure()
    palette = qualitative.Plotly + qualitative.Safe + qualitative.Dark24
    chain_ids = list(dict.fromkeys(ca["auth_chain_id"].astype(str)))

    for chain_index, chain_id in enumerate(chain_ids):
        chain = ca[ca["auth_chain_id"].astype(str) == chain_id].copy()
        color = palette[chain_index % len(palette)]
        coords = chain[["x", "y", "z"]].to_numpy(dtype=float)

        edge_x, edge_y, edge_z = [], [], []
        for left, right in zip(coords[:-1], coords[1:]):
            distance = float(np.linalg.norm(right - left))
            if max_ca_distance is not None and distance > max_ca_distance:
                continue
            edge_x.extend([left[0], right[0], None])
            edge_y.extend([left[1], right[1], None])
            edge_z.extend([left[2], right[2], None])
        if edge_x:
            figure.add_trace(
                go.Scatter3d(
                    x=edge_x,
                    y=edge_y,
                    z=edge_z,
                    mode="lines",
                    line={"color": color, "width": edge_width},
                    opacity=0.55,
                    hoverinfo="skip",
                    name=f"Chain {chain_id} edges",
                    legendgroup=chain_id,
                    showlegend=False,
                )
            )

        residue_name = chain.get(
            "res_name3l", chain.get("auth_comp_id", pd.Series("", index=chain.index))
        ).fillna("").astype(str)
        insertion = chain.get(
            "insertion", pd.Series("", index=chain.index)
        ).fillna("").astype(str)
        auth_position = chain["auth_seq_id"].astype(str) + insertion
        grn_text = chain["grn"].where(chain["grn"] != "", "unassigned")
        hovertext = [
            f"GRN: {grn}<br>Chain: {chain_id}<br>Residue: {name} {position}"
            for grn, name, position in zip(grn_text, residue_name, auth_position)
        ]
        marker_colors = [color if grn else "#B8B8B8" for grn in chain["grn"]]
        figure.add_trace(
            go.Scatter3d(
                x=coords[:, 0],
                y=coords[:, 1],
                z=coords[:, 2],
                mode="markers",
                marker={"size": marker_size, "color": marker_colors, "opacity": 0.9},
                text=hovertext,
                hoverinfo="text",
                name=f"Chain {chain_id}",
                legendgroup=chain_id,
            )
        )

    title = "GRN C-alpha structure"
    if structure_id:
        title = f"{structure_id}: GRN-annotated C-alpha trace"
    figure.update_layout(
        title=title,
        scene={
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "aspectmode": "data",
            "bgcolor": "rgba(0,0,0,0)",
        },
        paper_bgcolor="white",
        margin={"l": 0, "r": 0, "t": 50, "b": 0},
        legend={"itemsizing": "constant"},
    )
    return figure


def write_grn_ca_html(
    structure: pd.DataFrame,
    output_path: Union[str, Path],
    **plot_options,
) -> Path:
    """Write a self-contained interactive GRN C-alpha Plotly document."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure = plot_grn_ca_structure(structure, **plot_options)
    figure.write_html(output, include_plotlyjs=True, full_html=True)
    return output
