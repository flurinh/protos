#!/usr/bin/env python
"""Fetch representative signalling structures, assign GRNs, and plot C-alpha traces."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import plotly.io as pio


STRUCTURES = {
    "gpcr": {
        "pdb": "2RH1",
        "chains": ["A"],
        "annotations": [("A", "gpcrdb_class_a", "gpcrdb_class_a", (29, 342))],
    },
    "arrestin": {
        "pdb": "4JQI",
        "chains": ["A"],
        "annotations": [("A", "can_arrestin_human", "can_arrestin_human", None)],
    },
    "gprotein_trimer": {
        "pdb": "1GOT",
        "chains": ["A", "B", "G"],
        "annotations": [("A", "cgn_galpha_gio_human", "cgn_galpha_gio_human", None)],
    },
    "gpcr_gs_complex": {
        "pdb": "3SN6",
        "chains": ["R", "A", "B", "G"],
        "annotations": [
            ("R", "gpcrdb_class_a", "gpcrdb_class_a", (29, 341)),
            ("A", "cgn_galpha_gs_human", "cgn_galpha_gs_human", None),
        ],
    },
    "gpcr_arrestin_complex": {
        "pdb": "5W0P",
        "chains": ["A"],
        "annotations": [
            ("A", "gpcrdb_class_a", "gpcrdb_class_a", (1, 343)),
            ("A", "can_arrestin_human", "can_arrestin_human", (2012, 2362)),
        ],
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("data/grn_structure_demo"),
        help="Isolated ProtOS data root used for downloaded structures",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/visualizations/grn_structures"),
        help="Directory for individual Plotly documents and the gallery",
    )
    args = parser.parse_args()

    os.environ["PROTOS_DATA_ROOT"] = str(args.data_root.resolve())

    from protos.io.ingest.structure_loader import StructureLoader
    from protos.processing.structure import StructureProcessor
    from protos.visualization.grn_structure_vis import plot_grn_ca_structure

    args.output_dir.mkdir(parents=True, exist_ok=True)
    processor = StructureProcessor(name="grn_structure_demo")
    loader = StructureLoader(processor=processor)
    figures = []
    summary = {}

    for label, specification in STRUCTURES.items():
        pdb_id = specification["pdb"]
        structure_id = pdb_id.lower()
        registered = loader.download_and_register(pdb_id, name=structure_id)
        if registered is None:
            raise RuntimeError(f"RCSB download failed for {pdb_id}")

        annotation_runs = []
        for chain, reference_table, protein_family, residue_range in specification["annotations"]:
            _, annotation_summary = processor.annotate_with_grn(
                structure_id,
                reference_table=reference_table,
                protein_family=protein_family,
                chains=[chain],
                residue_ranges={chain: residue_range} if residue_range else None,
                save=True,
                return_summary=True,
            )
            info = annotation_summary["chains"].get(chain, {})
            annotation_runs.append(
                {
                    "chain": chain,
                    "reference_table": reference_table,
                    "auth_residue_range": residue_range,
                    "status": info.get("status"),
                    "reference": info.get("reference"),
                    "coverage": info.get("coverage"),
                    "assigned_positions": info.get("assigned_positions"),
                    "insertion_residues": info.get("insertion_residues"),
                    "deletion_residues": info.get("deletion_residues"),
                }
            )

        structure = processor.load_entity(structure_id)
        figure = plot_grn_ca_structure(
            structure,
            structure_id=f"{pdb_id} — {label.replace('_', ' ')}",
            chains=specification["chains"],
        )
        output = args.output_dir / f"{label}_{structure_id}.html"
        figure.write_html(output, include_plotlyjs=True, full_html=True)
        figures.append((label, pdb_id, figure, output))

        frame = structure.reset_index()
        ca = frame[
            frame["atom_name"].astype(str).str.upper().eq("CA")
            & frame["auth_chain_id"].astype(str).isin(specification["chains"])
        ]
        if "label_seq_id" in ca.columns and ca["label_seq_id"].notna().any():
            ca = ca[ca["label_seq_id"].notna()]
        assigned = ca["grn"].fillna("").astype(str).ne("")
        summary[label] = {
            "pdb_id": pdb_id,
            "chains": specification["chains"],
            "ca_atoms": int(len(ca)),
            "grn_assigned_ca_atoms": int(assigned.sum()),
            "annotations": annotation_runs,
            "output": str(output.resolve()),
        }

    sections = []
    for index, (label, pdb_id, figure, _) in enumerate(figures):
        plot = pio.to_html(
            figure,
            full_html=False,
            include_plotlyjs=True if index == 0 else False,
        )
        sections.append(f"<h2>{pdb_id} — {label.replace('_', ' ')}</h2>{plot}")
    gallery = args.output_dir / "grn_structure_gallery.html"
    gallery.write_text(
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>ProtOS GRN structure gallery</title>"
        "<style>body{font-family:sans-serif;margin:2rem}h2{margin-top:3rem}</style>"
        "</head><body><h1>ProtOS GRN-annotated signalling structures</h1>"
        + "".join(sections)
        + "</body></html>",
        encoding="utf-8",
    )
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(gallery.resolve())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
