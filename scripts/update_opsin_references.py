#!/usr/bin/env python3
"""Normalize the bundled type-I and WT-only type-II opsin GRN tables."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


DIRECTIONAL_LOOP = re.compile(r"^[1-9][1-9]\.\d+$")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_fasta_sequences(path: Path) -> set[str]:
    sequences: list[str] = []
    current: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if current:
                sequences.append("".join(current).replace("-", ""))
            current = []
        elif line:
            current.append(line)
    if current:
        sequences.append("".join(current).replace("-", ""))
    if len(sequences) != len(set(sequences)):
        raise ValueError("WT FASTA contains duplicate protein sequences")
    return set(sequences)


def row_sequence(row: pd.Series) -> str:
    return "".join(str(value)[0] for value in row if value != "-")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--type-i", type=Path, required=True)
    parser.add_argument("--type-ii-all", type=Path, required=True)
    parser.add_argument("--type-ii-wt-fasta", type=Path, required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    type_i = pd.read_csv(args.type_i, index_col=0, dtype=str).fillna("-")
    type_ii = pd.read_csv(args.type_ii_all, index_col=0, dtype=str).fillna("-")
    wt_sequences = read_fasta_sequences(args.type_ii_wt_fasta)

    matched = type_ii.apply(row_sequence, axis=1).isin(wt_sequences)
    wt = type_ii.loc[matched]
    if len(wt) != len(wt_sequences):
        raise ValueError(
            f"Expected {len(wt_sequences)} WT sequences, matched {len(wt)} table rows"
        )
    excluded = [column for column in wt.columns if DIRECTIONAL_LOOP.fullmatch(column)]
    wt = wt.drop(columns=excluded)

    type_i_path = output / "type_I_opsins.csv"
    type_ii_path = output / "type_II_opsins.csv"
    type_i.to_csv(type_i_path)
    wt.to_csv(type_ii_path)
    provenance = {
        "schema_version": 1,
        "type_I": {
            "policy": "Retain the curated microbial type-I opsin table unchanged.",
            "rows": len(type_i),
            "columns": len(type_i.columns),
        },
        "type_II": {
            "source_repository": "https://github.com/VisualPhysiologyDB/visual-physiology-opsin-db",
            "source_revision": args.source_revision,
            "source_subset": "VPOD_1.2 WT heterologous aligned FASTA",
            "source_license": "GPL-3.0",
            "wt_fasta_sha256": sha256(args.type_ii_wt_fasta),
            "selection": "Exact ungapped amino-acid sequence match to the upstream WT subset.",
            "input_rows": len(type_ii),
            "wt_rows": len(wt),
            "excluded_directional_loop_columns": excluded,
        },
        "files": {
            type_i_path.name: sha256(type_i_path),
            type_ii_path.name: sha256(type_ii_path),
        },
    }
    args.provenance.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
