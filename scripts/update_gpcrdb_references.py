#!/usr/bin/env python3
"""Build ProtOS GRN tables from a pinned gpcrdb_data checkout.

The script deliberately has no network code.  GPCR tables need an explicit
residue export because gpcrdb_data stores the inputs to the Protwis build
(anchors and numbering maps), not the built receptor-residue table itself.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


CLASS_FILES = {
    "Class A (Rhodopsin)": "gpcrdb_class_a.csv",
    "Class B1 (Secretin)": "gpcrdb_class_b1.csv",
    "Class B2 (Adhesion)": "gpcrdb_class_b2.csv",
    "Class C (Glutamate)": "gpcrdb_class_c.csv",
    "Class D1 (Ste2-like fungal pheromone)": "gpcrdb_class_d1.csv",
    "Class F (Frizzled)": "gpcrdb_class_f.csv",
    "Class O1 (fish-like odorant)": "gpcrdb_class_o1.csv",
    "Class O2 (tetrapod specific odorant)": "gpcrdb_class_o2.csv",
    "Class T2 (Taste 2)": "gpcrdb_class_t2.csv",
    "Unclassified": "gpcrdb_unclassified.csv",
}

GPROTEIN_FAMILIES = {
    "gs": {"GNAS2", "GNAL"},
    "gio": {"GNAI1", "GNAI2", "GNAI3", "GNAT1", "GNAT2", "GNAT3", "GNAO", "GNAZ"},
    "gq11": {"GNAQ", "GNA11", "GNA14", "GNA15"},
    "g1213": {"GNA12", "GNA13"},
}

ARRESTINS = {
    "P49407": "ARRB1_HUMAN",
    "P32121": "ARRB2_HUMAN",
    "P36575": "ARRC_HUMAN",
    "P10523": "ARRS_HUMAN",
}


def source_revision(root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_signal_label(label: str) -> str:
    """Use ProtOS dot notation and two digits for a signal-protein position."""
    label = label.strip().replace("(", "").replace(")", "")
    parts = label.split(".")
    if len(parts) != 3 or not parts[2].isdigit():
        raise ValueError(f"Unexpected signal-protein GRN: {label!r}")
    return f"{parts[0]}.{parts[1]}.{int(parts[2]):02d}"


def normalized_name(value: str) -> str:
    value = html.unescape(re.sub(r"<[^>]+>", "", value)).lower()
    value = re.sub(r"olfactory[- ]receptor", "or", value)
    replacements = {
        "α": "alpha",
        "β": "beta",
        "δ": "delta",
        "κ": "kappa",
        "μ": "mu",
        "–": "-",
    }
    for old, new in replacements.items():
        value = value.replace(old, new)
    value = re.sub(r"(?:-?human|-?receptor)", "", value)
    return re.sub(r"[^a-z0-9]", "", value)


def receptor_name_index(hierarchy_path: Path) -> dict[str, set[str]]:
    """Map every display name/alias/gene symbol to its top-level class."""
    index: dict[str, set[str]] = defaultdict(set)
    current_class: str | None = None
    for line in hierarchy_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(" ") and " | " in line:
            current_class = line.split(" | ", 1)[0].strip()
            continue
        stripped = line.strip()
        if not stripped.startswith('"gpcr"') or current_class is None:
            continue
        fields = next(csv.reader([stripped]))
        candidates = [fields[4]]
        if len(fields) > 8:
            candidates.extend(fields[8].split("|"))
        if len(fields) > 10:
            candidates.append(fields[10])
        for candidate in candidates:
            key = normalized_name(candidate)
            if key:
                index[key].add(current_class)
    return index


def build_receptor_tables(
    export_path: Path, hierarchy_path: Path, output_dir: Path
) -> tuple[list[Path], dict[str, object]]:
    table = pd.read_csv(export_path, index_col=0, dtype=str).fillna("-")
    table.index.name = "entity_name"
    index = receptor_name_index(hierarchy_path)
    assignments: dict[str, list[str]] = defaultdict(list)
    unmatched: list[str] = []
    ambiguous: dict[str, list[str]] = {}

    for entity in table.index:
        classes = index.get(normalized_name(str(entity)), set())
        if len(classes) == 1:
            assignments[next(iter(classes))].append(str(entity))
        elif not classes:
            unmatched.append(str(entity))
        else:
            ambiguous[str(entity)] = sorted(classes)

    if unmatched or ambiguous:
        raise ValueError(
            "Receptor export could not be classified exactly: "
            f"unmatched={unmatched}, ambiguous={ambiguous}"
        )

    written: list[Path] = []
    counts: dict[str, int] = {}
    for class_name, filename in CLASS_FILES.items():
        rows = assignments.get(class_name, [])
        destination = output_dir / filename
        table.loc[rows].to_csv(destination)
        written.append(destination)
        counts[class_name] = len(rows)

    return written, {
        "input_rows": len(table),
        "class_rows": counts,
        "empty_classes": [name for name, count in counts.items() if count == 0],
    }


def grn_table_from_residues(
    rows: pd.DataFrame, *, column_order: list[str] | None = None
) -> pd.DataFrame:
    rows = rows.copy()
    rows["label"] = rows["CGN"].map(canonical_signal_label)
    rows["cell"] = rows["Residue"].astype(str) + rows["Position"].astype(int).astype(str)
    order = column_order or (
        rows[["label", "sortColumn"]]
        .drop_duplicates()
        .sort_values("sortColumn")["label"]
        .tolist()
    )
    table = rows.pivot(index="Uniprot_ID", columns="label", values="cell")
    table = table.reindex(columns=order).fillna("-")
    table.index.name = "entity_name"
    return table


def build_cgn_tables(
    source: Path, lookup_source: Path, output_dir: Path
) -> tuple[list[Path], dict[str, object]]:
    residues = pd.read_csv(source, sep="\t", dtype=str).fillna("")
    residues["sortColumn"] = residues["sortColumn"].astype(int)
    lookup = pd.read_csv(lookup_source, dtype=str).fillna("")
    cgn_by_sort = {
        int(row["Unnamed: 0"]): row["CGN_new"] for _, row in lookup.iterrows()
    }
    all_cgns = [
        canonical_signal_label(cgn_by_sort[position])
        for position in sorted(cgn_by_sort)
    ]
    missing_sort_columns = sorted(set(residues["sortColumn"]) - set(cgn_by_sort))
    if missing_sort_columns:
        raise ValueError(f"CGN lookup is missing sort columns: {missing_sort_columns}")
    # PDB_UNIPROT_ENSEMBLE_ALL.txt contains a legacy CGN column.  Protwis itself
    # refreshes it from CGN_lookup.csv via sortColumn; do the same here.
    residues["CGN"] = residues["sortColumn"].map(cgn_by_sort)
    residues["label"] = residues["CGN"].map(canonical_signal_label)
    human = residues[residues["Uniprot_ID"].str.endswith("_HUMAN")].copy()
    table = grn_table_from_residues(human, column_order=all_cgns)
    written = [output_dir / "cgn_galpha_human.csv"]
    table.to_csv(written[0])

    family_counts: dict[str, int] = {}
    for family, genes in GPROTEIN_FAMILIES.items():
        selected = [name for name in table.index if name.split("_", 1)[0] in genes]
        destination = output_dir / f"cgn_galpha_{family}_human.csv"
        table.loc[selected].to_csv(destination)
        written.append(destination)
        family_counts[family] = len(selected)

    return written, {
        "human_proteins": len(table),
        "generic_positions": len(table.columns),
        "family_rows": family_counts,
    }


def build_can_table(source: Path, output_dir: Path) -> tuple[list[Path], dict[str, object]]:
    alignment = pd.read_csv(source, dtype=str).fillna("-")
    labels = [canonical_signal_label(label) for label in alignment["CAN_id"]]
    rows: dict[str, list[str]] = {}
    for accession, entity in ARRESTINS.items():
        sequence_position = 0
        values: list[str] = []
        for residue in alignment[accession]:
            residue = str(residue)
            if residue == "-":
                values.append("-")
            else:
                sequence_position += 1
                values.append(f"{residue}{sequence_position}")
        rows[entity] = values
    table = pd.DataFrame.from_dict(rows, orient="index", columns=labels)
    table.index.name = "entity_name"
    destination = output_dir / "can_arrestin_human.csv"
    table.to_csv(destination)
    return [destination], {
        "human_proteins": len(table),
        "generic_positions": len(table.columns),
    }


def validate_tables(paths: Iterable[Path]) -> None:
    cell_pattern = re.compile(r"^[A-Z][1-9][0-9]*$")
    for path in paths:
        table = pd.read_csv(path, index_col=0, dtype=str).fillna("-")
        if table.index.has_duplicates or table.columns.has_duplicates:
            raise ValueError(f"Duplicate rows or GRNs in {path.name}")
        for entity, row in table.iterrows():
            positions: list[int] = []
            for value in row:
                if value == "-":
                    continue
                if not cell_pattern.fullmatch(value):
                    raise ValueError(f"Malformed cell {value!r} in {path.name}/{entity}")
                positions.append(int(value[1:]))
            if len(positions) != len(set(positions)):
                raise ValueError(f"Duplicate sequence positions in {path.name}/{entity}")
            if path.name.startswith(("cgn_", "can_")) and positions != sorted(positions):
                raise ValueError(f"Non-monotonic sequence positions in {path.name}/{entity}")


def update_bundle_manifest(path: Path, output_dir: Path, provenance_path: Path) -> None:
    manifest = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}
    manifest.update(
        {
            "schema_version": 2,
            "bundle_version": "2026.07.12",
            "distribution": "bundled",
            "files": {
                item.name: sha256(item) for item in sorted(output_dir.glob("*.csv"))
            },
            "provenance": {
                "gpcrdb": provenance_path.name,
                "policy": "No runtime API or network fetch; regenerate from pinned raw sources.",
            },
        }
    )
    manifest.setdefault("aliases", {})
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gpcrdb-data", type=Path, required=True)
    parser.add_argument("--receptor-export", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--provenance", type=Path)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()

    root = args.gpcrdb_data.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    receptor_paths, receptor_stats = build_receptor_tables(
        args.receptor_export,
        root / "protein_data" / "proteins_and_families.txt",
        output,
    )
    cgn_paths, cgn_stats = build_cgn_tables(
        root / "g_protein_data" / "PDB_UNIPROT_ENSEMBLE_ALL.txt",
        root / "g_protein_data" / "CGN_lookup.csv",
        output,
    )
    can_paths, can_stats = build_can_table(
        root / "arrestin_data" / "CAN_aln.csv", output
    )
    generated = receptor_paths + cgn_paths + can_paths
    validate_tables(generated)

    provenance_path = args.provenance or output.parent / "gpcrdb_provenance.json"
    provenance = {
        "schema_version": 1,
        "source": {
            "repository": "https://github.com/protwis/gpcrdb_data",
            "commit": source_revision(root),
            "license": "CC-BY-4.0",
            "runtime_api_used": False,
        },
        "inputs": {
            "receptors": {
                "generated_export": str(args.receptor_export),
                "generated_export_sha256": sha256(args.receptor_export),
                "classification": "protein_data/proteins_and_families.txt",
                "numbering_schemes": "residue_data/generic_numbers/",
                "curated_anchors": "residue_data/reference_positions/",
            },
            "gprotein": [
                "g_protein_data/PDB_UNIPROT_ENSEMBLE_ALL.txt",
                "g_protein_data/CGN_lookup.csv",
            ],
            "arrestin": "arrestin_data/CAN_aln.csv",
        },
        "notation": {
            "gpcr": "existing ProtOS numeric dot notation",
            "cgn": "G.<segment>.<two-digit-position>",
            "can": "<domain>.<segment>.<two-digit-position>",
        },
        "statistics": {
            "receptors": receptor_stats,
            "cgn": cgn_stats,
            "can": can_stats,
        },
        "files": {path.name: sha256(path) for path in sorted(generated)},
    }
    provenance_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    if args.manifest:
        update_bundle_manifest(args.manifest, output, provenance_path)
    print(f"generated {len(generated)} tables from gpcrdb_data {provenance['source']['commit'][:12]}")
    print(f"provenance: {provenance_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
