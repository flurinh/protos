"""High-level helpers for assigning Generic Residue Numbers (GRNs)."""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from protos.analysis.sequence.alignment_utils import init_biopython_aligner, perform_pairwise_alignment
from protos.io.formats.fasta_utils import read_fasta
from protos.processing.grn.grn_table_utils import expand_annotation, init_row_from_alignment
from protos.processing.grn.grn_utils import (
    GRNConfigManager,
    GRN_GAP_SYMBOL,
    GRN_UNKNOWN_SYMBOL,
    get_seq,
    init_grn_intervals,
    sort_grns_str,
    validate_grn_string,
)
from protos.processing.sequence.seq_alignment import mmseqs2_align2

logger = logging.getLogger(__name__)

_GAP_TOKENS = {GRN_GAP_SYMBOL, GRN_UNKNOWN_SYMBOL, '-', 'X', None}


def _sanitize_sequence(sequence: str) -> str:
    """Return an uppercase sequence without whitespace or gap characters."""

    return sequence.replace('\n', '').replace('-', '').replace(' ', '').upper()


def _ensure_directory(path: str) -> str:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return str(directory)


def _select_mmseqs_hits(hits: pd.DataFrame) -> Dict[str, Dict[str, float]]:
    if hits is None or hits.empty:
        return {}

    df = hits.copy()
    for column in ("e_value", "bit_score", "sequence_identity"):
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    df = df.sort_values(
        by=["query_id", "e_value", "bit_score"],
        ascending=[True, True, False],
        na_position="last",
    )

    best_rows = df.groupby("query_id", as_index=False).first()
    summary: Dict[str, Dict[str, float]] = {}
    for row in best_rows.itertuples(index=False):
        summary[row.query_id] = {
            "target_id": row.target_id,
            "bit_score": getattr(row, "bit_score", None),
            "e_value": getattr(row, "e_value", None),
            "sequence_identity": getattr(row, "sequence_identity", None),
        }

    return summary


def _build_reference_sequences(reference_table: pd.DataFrame) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    for ref_id in reference_table.index:
        seq = get_seq(ref_id, reference_table)
        if not seq:
            continue
        sanitized = _sanitize_sequence(seq)
        if sanitized:
            sequences[ref_id] = sanitized
    return sequences


def _determine_strict_grns(
    *,
    strict_grns: Optional[Sequence[str]] = None,
    protein_family: Optional[str] = None,
    config_manager: Optional[GRNConfigManager] = None,
    reference_table: Optional[pd.DataFrame] = None,
) -> List[str]:
    if strict_grns:
        return list(strict_grns)

    family = protein_family or (reference_table.attrs.get("protein_family") if reference_table is not None else None)
    manager = config_manager or GRNConfigManager()
    if family:
        config = manager.get_config(family, strict=True)
        if config:
            return init_grn_intervals(config)

    if reference_table is None:
        return []

    valid_columns = []
    for column in reference_table.columns:
        is_valid, _ = validate_grn_string(str(column))
        if is_valid:
            valid_columns.append(str(column))
    return valid_columns


def _construct_seed_row(
    reference_row: pd.Series,
    alignment_lines: Sequence[str],
) -> pd.Series:
    ref_dict = {
        grn: val
        for grn, val in reference_row.items()
        if val not in _GAP_TOKENS
    }
    if not ref_dict:
        return pd.Series(dtype=object)

    seq_pos2grn = {idx + 1: grn for idx, grn in enumerate(ref_dict.keys())}
    seed_row = init_row_from_alignment(list(alignment_lines), seq_pos2grn)
    return seed_row.replace({None: '-'})


def _expand_row(
    seed_row: pd.Series,
    sequence: str,
    alignment_lines: Sequence[str],
    *,
    protein_family: str,
    verbose: int = 0,
) -> Tuple[List[str], List[str]]:
    grn_list, rn_list, _missing = expand_annotation(
        seed_row,
        sequence,
        list(alignment_lines),
        protein_family=protein_family,
        max_alignment_gap=1,
        verbose=verbose,
    )
    return grn_list, rn_list


def assign_grns_to_sequences(
    query_sequences: Dict[str, str],
    reference_table: pd.DataFrame,
    *,
    protein_family: Optional[str] = None,
    strict_grns: Optional[Sequence[str]] = None,
    use_mmseqs: bool = True,
    aligner=None,
    temp_folder: Optional[str] = None,
    verbose: int = 0,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, object]]]:
    """Assign GRNs to a set of query sequences using a reference table.

    Returns the resulting annotation table and per-sequence alignment metadata.
    """

    if not query_sequences:
        return pd.DataFrame(columns=reference_table.columns), {}

    if reference_table.empty:
        raise ValueError("Reference table is empty; cannot assign GRNs")

    sanitized_queries = {
        seq_id: _sanitize_sequence(seq)
        for seq_id, seq in query_sequences.items()
        if _sanitize_sequence(seq)
    }

    if not sanitized_queries:
        raise ValueError("No valid sequences supplied for GRN assignment")

    reference_sequences = _build_reference_sequences(reference_table)
    if not reference_sequences:
        raise ValueError("Reference table does not contain usable sequences")

    aligner = aligner or init_biopython_aligner()
    strict_positions = _determine_strict_grns(
        strict_grns=strict_grns,
        protein_family=protein_family,
        reference_table=reference_table,
    )

    mmseqs_summary: Dict[str, Dict[str, float]] = {}
    if use_mmseqs:
        temp_dir = _ensure_directory(temp_folder or os.path.join("temp", "grn_assign"))
        hits = mmseqs2_align2(sanitized_queries, reference_sequences, temp_folder=temp_dir)
        mmseqs_summary = _select_mmseqs_hits(hits)

    ref_valid_columns = [
        str(col)
        for col in reference_table.columns
        if validate_grn_string(str(col))[0]
    ]
    total_reference_positions = len(ref_valid_columns)

    metadata: Dict[str, Dict[str, object]] = {}
    rows: Dict[str, pd.Series] = {}

    for seq_id, sequence in sanitized_queries.items():
        preferred_ref_id: Optional[str] = None
        mmseqs_meta = mmseqs_summary.get(seq_id)
        if mmseqs_meta:
            preferred_ref_id = mmseqs_meta.get("target_id")

        best_result = None
        best_ref_id = None

        if preferred_ref_id and preferred_ref_id in reference_sequences:
            best_result = perform_pairwise_alignment(
                sequence,
                reference_sequences[preferred_ref_id],
                aligner,
                seq_id,
                preferred_ref_id,
            )
            best_ref_id = preferred_ref_id

        if best_result is None:
            best_score = float("-inf")
            for ref_id, ref_sequence in reference_sequences.items():
                result = perform_pairwise_alignment(sequence, ref_sequence, aligner, seq_id, ref_id)
                if result.score > best_score:
                    best_score = result.score
                    best_result = result
                    best_ref_id = ref_id

        if best_result is None or best_ref_id is None:
            logger.warning("No alignment produced for %s; skipping", seq_id)
            continue

        if best_ref_id not in reference_table.index:
            logger.warning(
                "Reference '%s' selected for %s is not present in table; skipping",
                best_ref_id,
                seq_id,
            )
            continue

        alignment_lines = list(best_result.alignment_lines)
        seed_row = _construct_seed_row(reference_table.loc[best_ref_id], alignment_lines)

        if seed_row.empty:
            logger.warning("Seed row for %s is empty; skipping", seq_id)
            continue

        if strict_positions:
            available = [grn for grn in strict_positions if grn in seed_row.index]
            filtered_seed = seed_row.loc[available] if available else seed_row
        else:
            filtered_seed = seed_row

        family = protein_family or 'gpcr_a'
        try:
            grn_list, rn_list = _expand_row(
                filtered_seed,
                sequence,
                alignment_lines,
                protein_family=family,
                verbose=verbose,
            )
        except Exception as exc:  # pragma: no cover - defensive fallback
            logger.warning("GRN expansion failed for %s (%s); using seed mapping", seq_id, exc)
            grn_list = [grn for grn, value in filtered_seed.items() if value not in _GAP_TOKENS]
            rn_list = [filtered_seed[grn] for grn in grn_list]

        if not grn_list:
            logger.warning("No GRNs could be assigned for %s", seq_id)
            continue

        # Remove duplicate GRNs (keep first occurrence)
        seen_grns = set()
        unique_grn_list = []
        unique_rn_list = []
        for grn, rn in zip(grn_list, rn_list):
            if grn not in seen_grns:
                seen_grns.add(grn)
                unique_grn_list.append(grn)
                unique_rn_list.append(rn)

        annotation_row = pd.Series(unique_rn_list, index=unique_grn_list, dtype=object)
        rows[seq_id] = annotation_row

        assigned_valid = [
            grn
            for grn in annotation_row.index
            if validate_grn_string(str(grn))[0]
            and annotation_row[grn] not in _GAP_TOKENS
        ]
        meta_entry: Dict[str, object] = {
            "reference_id": best_ref_id,
            "alignment_score": float(best_result.score),
            "assigned_positions": int(len(assigned_valid)),
            "reference_positions": total_reference_positions,
            "coverage_fraction": (
                float(len(assigned_valid)) / total_reference_positions
                if total_reference_positions > 0
                else 0.0
            ),
        }

        if mmseqs_meta:
            meta_entry.update(mmseqs_meta)
            meta_entry.setdefault("method", "mmseqs+biopython")
        else:
            meta_entry["method"] = "biopython"

        metadata[seq_id] = meta_entry

    if not rows:
        return pd.DataFrame(columns=reference_table.columns), metadata

    annotation_table = pd.DataFrame(rows).T
    annotation_table.index.name = "entity_name"
    annotation_table = annotation_table.fillna('-')

    # Remove duplicate columns (keep first occurrence)
    annotation_table = annotation_table.loc[:, ~annotation_table.columns.duplicated()]

    # Collect all columns (from reference and annotation table)
    all_columns = set(str(col) for col in reference_table.columns)
    all_columns.update(str(col) for col in annotation_table.columns)

    # Filter to valid GRN columns and sort them
    valid_grn_columns = [col for col in all_columns if validate_grn_string(col)[0]]
    sorted_columns = sort_grns_str(valid_grn_columns)

    # Add any non-GRN columns at the end (shouldn't happen normally)
    non_grn_columns = [col for col in all_columns if not validate_grn_string(col)[0]]
    ordered_columns = sorted_columns + sorted(non_grn_columns)

    # Add missing columns with '-' fill value
    for col in ordered_columns:
        if col not in annotation_table.columns:
            annotation_table[col] = '-'

    # Reorder columns in sorted GRN order
    annotation_table = annotation_table[ordered_columns]
    annotation_table = annotation_table.astype(str)

    return annotation_table, metadata


def _parse_sequences(path: str) -> Dict[str, str]:
    sequences = read_fasta(path)
    return {seq_id: seq for seq_id, seq in sequences.items() if seq}


def _load_reference_table(ref_name: str, base_path: Path) -> pd.DataFrame:
    path = base_path / f"{ref_name}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Reference table not found: {path}")
    table = pd.read_csv(path, index_col=0)
    table.index.name = "entity_name"
    return table


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Assign GRNs to a FASTA dataset")
    parser.add_argument("--protein-family", "-p", required=True, help="Protein family key (e.g. gpcr_a)")
    parser.add_argument("--dataset", "-d", required=True, help="Path to query FASTA file")
    parser.add_argument("--reference", "-r", required=True, help="Name of reference table (CSV under grn/reference)")
    parser.add_argument("--output", "-o", help="Output CSV path (defaults to dataset name + _grn.csv)")
    parser.add_argument("--data-root", default="data", help="Project data root containing grn/reference")
    parser.add_argument("--no-mmseqs", action="store_true", help="Disable MMseqs acceleration")

    args = parser.parse_args(argv)

    data_root = Path(args.data_root)
    ref_table = _load_reference_table(args.reference, data_root / "grn" / "reference")
    queries = _parse_sequences(args.dataset)

    annotation_table, metadata = assign_grns_to_sequences(
        queries,
        ref_table,
        protein_family=args.protein_family,
        use_mmseqs=not args.no_mmseqs,
    )

    if annotation_table.empty:
        logger.error("No annotations were produced")
        return

    output_path = args.output
    if not output_path:
        output_path = Path(args.dataset).with_suffix("_grn.csv")

    annotation_table.to_csv(output_path)
    logger.info("Saved GRN annotations to %s", output_path)
    logger.debug("Alignment summary: %s", metadata)


if __name__ == "__main__":  # pragma: no cover - CLI helper
    logging.basicConfig(level=logging.INFO)
    main()
