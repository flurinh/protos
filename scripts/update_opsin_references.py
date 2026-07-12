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


def read_fasta_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    entity_id: str | None = None
    current: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if entity_id is not None:
                sequences[entity_id] = "".join(current).replace("-", "").upper()
            entity_id = line[1:].split()[0]
            if not entity_id or entity_id in sequences:
                raise ValueError(f"Duplicate or empty FASTA entity ID: {entity_id!r}")
            current = []
        elif line:
            current.append(line)
    if entity_id is not None:
        sequences[entity_id] = "".join(current).replace("-", "").upper()
    if not sequences:
        raise ValueError("WT FASTA contains no sequences")
    if any(not sequence for sequence in sequences.values()):
        raise ValueError("WT FASTA contains an empty protein sequence")
    if len(sequences) != len(set(sequences.values())):
        raise ValueError("WT FASTA contains duplicate protein sequences")
    return sequences


def read_wt_metadata(path: Path) -> dict[str, str]:
    """Map exact WT protein sequences to their unique source accessions."""

    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"Accession", "Protein"}
    if not required <= set(table.columns):
        raise ValueError(f"WT metadata is missing columns: {sorted(required - set(table))}")
    if table["Accession"].duplicated().any() or table["Protein"].duplicated().any():
        raise ValueError("WT metadata accessions and protein sequences must be unique")
    if table["Accession"].str.strip().eq("").any() or table["Protein"].eq("").any():
        raise ValueError("WT metadata contains an empty accession or protein sequence")
    if "Mutations" in table and table["Mutations"].str.strip().ne("").any():
        raise ValueError("WT metadata unexpectedly contains mutant entries")
    return dict(zip(table["Protein"].str.upper(), table["Accession"]))


def row_sequence(row: pd.Series) -> str:
    return "".join(str(value)[0] for value in row if value != "-")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--type-i", type=Path, required=True)
    parser.add_argument("--type-ii-all", type=Path, required=True)
    parser.add_argument("--type-ii-wt-fasta", type=Path, required=True)
    parser.add_argument("--type-ii-wt-metadata", type=Path)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--source-revision", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    type_i = pd.read_csv(args.type_i, index_col=0, dtype=str).fillna("-")
    type_ii = pd.read_csv(args.type_ii_all, index_col=0, dtype=str).fillna("-")
    wt_fasta = read_fasta_sequences(args.type_ii_wt_fasta)
    if args.type_ii_wt_metadata:
        entity_by_sequence = read_wt_metadata(args.type_ii_wt_metadata)
        if set(entity_by_sequence) != set(wt_fasta.values()):
            raise ValueError("WT FASTA and metadata contain different protein sequences")
    else:
        entity_by_sequence = {sequence: entity for entity, sequence in wt_fasta.items()}

    source_sequences = type_ii.apply(row_sequence, axis=1)
    matched = source_sequences.isin(entity_by_sequence)
    wt = type_ii.loc[matched]
    if len(wt) != len(entity_by_sequence):
        raise ValueError(
            f"Expected {len(entity_by_sequence)} WT sequences, matched {len(wt)} table rows"
        )
    wt.index = [entity_by_sequence[sequence] for sequence in source_sequences[matched]]
    wt.index.name = "entity_name"
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
            "source_subset": f"{args.source_version} WT heterologous aligned FASTA",
            "source_license": "GPL-3.0",
            "wt_fasta_sha256": sha256(args.type_ii_wt_fasta),
            "selection": (
                "Exact ungapped amino-acid sequence match to the upstream WT subset; "
                "opaque GRN row hashes replaced by unique source accessions."
            ),
            "wt_metadata_sha256": (
                sha256(args.type_ii_wt_metadata) if args.type_ii_wt_metadata else None
            ),
            "source_table_sha256": sha256(args.type_ii_all),
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
