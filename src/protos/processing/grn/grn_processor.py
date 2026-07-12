"""GRNProcessor: manages GRN reference tables and sequence annotations."""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine
from protos.processing.grn.grn_utils import (
    make_directional_flexible_grn,
    make_insertion_grn,
    normalize_grn_format,
    parse_directional_flexible_grn,
    validate_grn_string,
)

GRN_INDEX = "entity_name"
GRN_COLUMNS = "grn_columns"
ARTIFACT_PATH = "artifact_path"


class GRNProcessor(BaseProcessor):
    """Stores reference GRN tables and produces new annotations."""

    processor_type = "grn"

    def __init__(self, name: str = "grn_processor") -> None:
        super().__init__(name=name)
        self._table_cache: Dict[str, pd.DataFrame] = {}

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------
    @property
    def tables_dir(self) -> Path:
        return Path(self.get_subdirectory_path("table_dir"))

    @property
    def reference_dir(self) -> Path:
        return Path(self.get_subdirectory_path("reference_dir"))

    @property
    def configs_dir(self) -> Path:
        return Path(self.get_subdirectory_path("configs_dir"))

    def _table_path(self, table_name: str) -> Path:
        return self.tables_dir / f"{table_name}.csv"

    def _relative_path(self, path: Path) -> str:
        return str(path.relative_to(self.paths.data_root))

    # ------------------------------------------------------------------
    # GRN parsing and mapping utilities
    # ------------------------------------------------------------------
    @staticmethod
    def parse_grn_value(value: Any) -> Optional[Tuple[str, int]]:
        """Parse a GRN table cell like 'M123' into (residue, 123).

        Returns None for missing or malformed entries ('-').
        """
        if not isinstance(value, str) or not value or value == '-':
            return None
        head = value[0]
        tail = value[1:]
        try:
            pos = int(tail)
        except (TypeError, ValueError):
            return None
        return head, pos

    @staticmethod
    def build_grn_to_seq_index(grn_table: pd.DataFrame, *, sequence_id: str) -> Dict[str, int]:
        """Construct mapping from GRN label to 1-based sequence indices for one sequence.

        Expects GRN table to have sequence ids as index and GRN labels as columns,
        with cell values formatted like 'A123' (amino acid + sequence index).
        """
        if sequence_id not in grn_table.index:
            return {}
        row = grn_table.loc[sequence_id]
        mapping: Dict[str, int] = {}
        for label, value in row.items():
            parsed = GRNProcessor.parse_grn_value(value)
            if parsed is None:
                continue
            _, seq_idx = parsed
            mapping[str(label)] = int(seq_idx)
        return mapping

    # ------------------------------------------------------------------
    # Public API: recording / loading tables
    # ------------------------------------------------------------------
    def record_table(
        self,
        table_name: str,
        table: pd.DataFrame,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        per_entity_metadata: Optional[Dict[str, Dict[str, Any]]] = None,
        allow_create: bool = False,
        link_entities: bool = True,
    ) -> pd.DataFrame:
        """Persist a GRN table and register dataset-level relationships."""

        if table.index.name != GRN_INDEX:
            table = table.copy()
            table.index.name = GRN_INDEX

        normalized_columns = self._normalize_columns(table.columns)
        table = table.copy()
        table.columns = normalized_columns
        table = table.fillna('-')
        self._write_table(table_name, table)

        artifact_rel = self._relative_path(self._table_path(table_name))
        dataset_metadata = {
            ARTIFACT_PATH: artifact_rel,
            GRN_COLUMNS: list(table.columns),
            "row_count": int(len(table)),
        }
        if metadata:
            dataset_metadata.update(metadata)

        self._register_table_entity(table_name, artifact_rel, dataset_metadata)

        if self.dataset_manager.dataset_exists(table_name):
            self.dataset_manager.update_metadata(table_name, dataset_metadata)
            self.dataset_manager.add_to_dataset(table_name, [table_name])
        else:
            self.create_dataset(table_name, [table_name], dataset_metadata)

        if link_entities:
            for entity_name in table.index:
                entity_meta = per_entity_metadata.get(entity_name) if per_entity_metadata else None
                self._link_table_to_entity(
                    table_name,
                    entity_name,
                    metadata=entity_meta,
                    allow_create=allow_create,
                )

        return table

    def load_table(self, table_name: str) -> pd.DataFrame:
        return self._load_table(table_name).copy()

    def list_tables(self) -> List[str]:
        return self.dataset_manager.list_datasets()

    # ------------------------------------------------------------------
    # Reference utilities
    # ------------------------------------------------------------------
    def load_reference_table(self, reference_name: str) -> pd.DataFrame:
        path = self.reference_dir / f"{reference_name}.csv"
        if not path.exists():
            raise FileNotFoundError(f"Reference table not found: {path}")
        df = pd.read_csv(path, index_col=0)
        df.index.name = GRN_INDEX
        df.columns = self._normalize_columns(df.columns)
        return df

    # ------------------------------------------------------------------
    # Annotation helpers
    # ------------------------------------------------------------------
    def annotate_sequences(
        self,
        sequences: Dict[str, str],
        *,
        reference_table: str,
        protein_family: str,
        search: str = "pairwise",
        open_gap_score: float = -10.0,
        extend_gap_score: float = -0.05,
        end_gap_score: float = 0.0,
        min_normalized_score: Optional[float] = None,
        min_coverage: Optional[float] = None,
        assign_unambiguous_insertions: bool = True,
        standard_regions: Optional[Dict[str, Tuple[str, str]]] = None,
        strict_regions: Optional[Dict[str, Tuple[str, str]]] = None,
        max_alignment_gap: int = 1,
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Annotate sequences by projecting a gap-aware reference alignment.

        Insertions and deletions are reported in the per-sequence summary.
        Insertions inside a segment create ``position + .001`` columns, while
        residues between segments receive directional flexible-region labels.
        Terminal overhangs receive ``n.<distance>`` and ``c.<distance>`` labels.
        Deletion columns remain ``-``.  Alignment gaps larger than
        ``max_alignment_gap`` are compressed only when both anchors are inside
        the same configured strict non-flexible region.
        """

        if not sequences:
            raise ValueError("No sequences provided for GRN annotation")
        if max_alignment_gap < 0:
            raise ValueError("max_alignment_gap must be zero or greater")

        reference_df = self.load_reference_table(reference_table)
        ref_sequences, ref_grn_map = self._prepare_reference_sequences(reference_df)
        if not ref_sequences:
            raise ValueError(
                f"Reference table '{reference_table}' has no usable reference sequences"
            )

        if search != "pairwise":
            raise ValueError(f"Unsupported GRN reference search method: {search!r}")

        aligner = SequenceAlignmentEngine(
            open_gap_score=open_gap_score,
            extend_gap_score=extend_gap_score,
            end_gap_score=end_gap_score,
        )
        rows: List[pd.Series] = []
        per_sequence: Dict[str, Dict[str, Any]] = {}
        column_order = reference_df.columns.tolist()
        sequence_grn_order = self._infer_reference_grn_order(reference_df)
        configured_standard, configured_strict = self._load_region_configuration(
            protein_family
        )
        standard_regions = (
            configured_standard if standard_regions is None else standard_regions
        )
        strict_regions = configured_strict if strict_regions is None else strict_regions
        generated_after: Dict[str, List[str]] = {}
        generated_n_terminal: List[str] = []
        generated_c_terminal: List[str] = []

        for seq_name, raw_sequence in sequences.items():
            try:
                sequence = self._sanitize_sequence(raw_sequence)
            except ValueError as exc:
                per_sequence[seq_name] = {
                    "reference": None,
                    "normalized_score": None,
                    "coverage": 0.0,
                    "status": "invalid_sequence",
                    "error": str(exc),
                }
                rows.append(pd.Series({col: '-' for col in column_order}, name=seq_name))
                continue
            if not sequence:
                per_sequence[seq_name] = {
                    "reference": None,
                    "normalized_score": None,
                    "coverage": 0.0,
                    "status": "empty_sequence",
                }
                rows.append(pd.Series({col: '-' for col in column_order}, name=seq_name))
                continue

            best_ref, best_alignment, best_score = self._find_best_reference(
                seq_name,
                sequence,
                ref_sequences,
                aligner,
            )

            if best_alignment is None or best_ref is None:
                per_sequence[seq_name] = {
                    "reference": None,
                    "normalized_score": None,
                    "coverage": 0.0,
                    "status": "alignment_failed",
                }
                rows.append(pd.Series({col: '-' for col in column_order}, name=seq_name))
                continue

            series, coverage, assigned, indels = self._project_alignment(
                seq_name,
                best_alignment,
                ref_grn_map[best_ref],
                column_order,
                sequence_grn_order=sequence_grn_order,
                assign_unambiguous_insertions=assign_unambiguous_insertions,
                standard_regions=standard_regions,
                strict_regions=strict_regions,
                max_alignment_gap=max_alignment_gap,
            )

            for event in indels["insertions"]:
                if event.get("kind") == "n_terminal_overhang":
                    self._merge_ordered_labels(
                        generated_n_terminal, event.get("generated_grns", [])
                    )
                    continue
                if event.get("kind") == "c_terminal_overhang":
                    self._merge_ordered_labels(
                        generated_c_terminal, event.get("generated_grns", [])
                    )
                    continue
                anchor = event.get("after_grn")
                if anchor is None:
                    continue
                bucket = generated_after.setdefault(anchor, [])
                self._merge_ordered_labels(bucket, event.get("generated_grns", []))

            series.name = seq_name
            rows.append(series)
            accepted = (
                (min_normalized_score is None or best_score >= min_normalized_score)
                and (min_coverage is None or coverage >= min_coverage)
            )
            if not accepted:
                series = pd.Series({col: '-' for col in column_order}, name=seq_name)
                rows[-1] = series

            per_sequence[seq_name] = {
                "reference": best_ref,
                "normalized_score": best_score,
                "coverage": coverage,
                "assigned_positions": assigned,
                "insertions": indels["insertions"],
                "deletions": indels["deletions"],
                "insertion_residues": indels["insertion_residues"],
                "deletion_residues": indels["deletion_residues"],
                "terminal_overhang_residues": indels["terminal_overhang_residues"],
                "status": (
                    "below_threshold"
                    if not accepted
                    else "ok" if coverage > 0 else "no_overlap"
                ),
            }

        annotations = pd.DataFrame(rows)
        annotations.index.name = GRN_INDEX
        annotations = annotations.fillna('-')
        output_columns: List[str] = list(generated_n_terminal)
        for column in column_order:
            output_columns.append(column)
            output_columns.extend(generated_after.get(column, []))
        output_columns.extend(
            column
            for column in annotations.columns
            if column not in output_columns and column not in generated_c_terminal
        )
        output_columns.extend(generated_c_terminal)
        annotations = annotations.reindex(columns=output_columns, fill_value='-')

        summary = {
            "per_sequence": per_sequence,
            "global": {
                "annotated": sum(
                    1 for info in per_sequence.values() if info["status"] == "ok"
                ),
                "total": len(per_sequence),
                "reference_table": reference_table,
                "protein_family": protein_family,
                "search": search,
                "gap_parameters": dict(aligner.gap_parameters),
                "min_normalized_score": min_normalized_score,
                "min_coverage": min_coverage,
                "assign_unambiguous_insertions": assign_unambiguous_insertions,
                "standard_regions": standard_regions,
                "strict_regions": strict_regions,
                "max_alignment_gap": max_alignment_gap,
            },
        }

        return annotations, summary

    # ------------------------------------------------------------------
    # BaseProcessor abstract methods
    # ------------------------------------------------------------------
    def load_entity(self, name: str) -> Optional[Dict[str, Any]]:
        tables = self.dataset_manager.list_datasets()
        for table_name in tables:
            table = self._load_table(table_name)
            if name in table.index:
                return table.loc[name].to_dict()
        return None

    def save_entity(
        self,
        name: str,
        data: Dict[str, Any],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        table_name = (metadata or {}).get("table")
        if not table_name:
            raise ValueError("Metadata must include table name for GRN rows")

        table = self._load_table(table_name)
        series = pd.Series(data, name=name)
        table.loc[name] = series
        self.record_table(table_name, table, metadata=metadata)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _register_table_entity(
        self,
        table_name: str,
        artifact_rel: str,
        metadata: Dict[str, Any],
    ) -> None:
        entity_metadata = dict(metadata)
        entity_metadata.setdefault("table", table_name)
        self.entity_registry.register_entity(
            name=table_name,
            format_type=self.processor_type,
            file_path=artifact_rel,
            metadata=entity_metadata,
        )

    def _link_table_to_entity(
        self,
        table_name: str,
        entity_name: str,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        allow_create: bool,
    ) -> None:
        if not self._ensure_entity_exists(entity_name, default_format="sequence", allow_create=allow_create):
            return

        relationship_metadata = {"table": table_name}
        related_structures = []
        if metadata:
            relationship_metadata.update(
                {k: v for k, v in metadata.items() if k != "related_structures"}
            )
            related_structures = metadata.get("related_structures", [])

        try:
            self.entity_registry.add_relationship(
                source_name=table_name,
                target_name=entity_name,
                rel_type="annotated_by",
                metadata=relationship_metadata,
            )
        except ValueError:
            pass

        for related in related_structures:
            if isinstance(related, dict):
                struct_name = related.get("name")
                chain_id = related.get("chain_id")
            else:
                struct_name = related
                chain_id = None

            if not struct_name:
                continue

            if not self._ensure_entity_exists(struct_name, default_format="structure", allow_create=allow_create):
                continue

            structure_metadata = {
                "table": table_name,
                "sequence": entity_name,
            }
            if chain_id is not None:
                structure_metadata["chain_id"] = chain_id

            try:
                self.entity_registry.add_relationship(
                    source_name=table_name,
                    target_name=struct_name,
                    rel_type="annotated_by",
                    metadata=structure_metadata,
                )
            except ValueError:
                continue

    def _load_table(self, table_name: str) -> pd.DataFrame:
        if table_name in self._table_cache:
            return self._table_cache[table_name]

        path = self._table_path(table_name)
        if not path.exists():
            df = pd.DataFrame()
            self._table_cache[table_name] = df
            return df

        df = pd.read_csv(path, index_col=0)
        df.index.name = GRN_INDEX
        df = df.fillna('-')
        self._table_cache[table_name] = df
        return df

    def _write_table(self, table_name: str, table: pd.DataFrame) -> None:
        path = self._table_path(table_name)
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path)
        self._table_cache[table_name] = table

    def _normalize_columns(self, columns: Iterable[Any]) -> List[str]:
        normalized = []
        for col in columns:
            col_str = str(col)
            # Canonical storage is dot notation even when a valid legacy x
            # label is supplied by an older table.
            if "x" in col_str:
                col_str = normalize_grn_format(col_str)
            is_valid, _ = validate_grn_string(col_str)
            if not is_valid:
                col_str = normalize_grn_format(col_str)
            normalized.append(col_str)
        return normalized

    def _prepare_reference_sequences(
        self,
        reference_df: pd.DataFrame,
    ) -> Tuple[Dict[str, str], Dict[str, List[str]]]:
        sequences: Dict[str, str] = {}
        residue_maps: Dict[str, List[str]] = {}

        for ref_name, row in reference_df.iterrows():
            mapped_residues: List[Tuple[int, str, str]] = []
            for grn, value in row.items():
                parsed = self.parse_grn_value(value)
                if parsed is not None:
                    residue, sequence_position = parsed
                    mapped_residues.append((sequence_position, residue, str(grn)))
            mapped_residues.sort(key=lambda item: item[0])
            if not mapped_residues:
                continue
            positions = [item[0] for item in mapped_residues]
            if len(positions) != len(set(positions)):
                raise ValueError(f"Reference row '{ref_name}' has duplicate sequence positions")
            sequences[ref_name] = ''.join(item[1] for item in mapped_residues)
            residue_maps[ref_name] = [item[2] for item in mapped_residues]

        return sequences, residue_maps

    @staticmethod
    def _split_segment_position(label: str) -> Optional[Tuple[str, float]]:
        """Split ``<segment>.<position>`` without assuming a numeric segment."""

        if "." not in label:
            return None
        if parse_directional_flexible_grn(label) is not None:
            return None
        segment, position = label.rsplit(".", 1)
        if not segment or not position.isdigit():
            return None
        if len(position) > 2:
            # Structure-corrected insertion: 1.411 or G.HN.031 sorts between
            # its two-digit anchor and the following standard coordinate.
            suffix = position[2:]
            order_value = int(position[:2]) + int(suffix) / (10 ** len(suffix))
        else:
            order_value = float(int(position))
        return segment, order_value

    @classmethod
    def _infer_reference_grn_order(cls, reference_df: pd.DataFrame) -> List[str]:
        """Infer sequence order from residue positions, then table appearance.

        This handles arbitrary segment identifiers and repairs tables whose
        physical column order is not sequence order.  Empty GRN columns inherit
        the observed order and direction of their segment.
        """

        columns = [str(column) for column in reference_df.columns]
        observations: Dict[str, List[int]] = {column: [] for column in columns}
        segment_labels: Dict[str, List[Tuple[str, float, int]]] = {}

        for column_index, column in enumerate(columns):
            parsed_label = cls._split_segment_position(column)
            for value in reference_df.iloc[:, column_index]:
                parsed_value = cls.parse_grn_value(value)
                if parsed_value is not None:
                    observations[column].append(parsed_value[1])
            if parsed_label is None:
                continue
            segment, position = parsed_label
            segment_labels.setdefault(segment, []).append((column, position, column_index))

        def median(values: Sequence[float]) -> Optional[float]:
            if not values:
                return None
            ordered = sorted(values)
            middle = len(ordered) // 2
            if len(ordered) % 2:
                return float(ordered[middle])
            return (ordered[middle - 1] + ordered[middle]) / 2.0

        segment_rank: Dict[str, Tuple[float, int]] = {}
        segment_direction: Dict[str, int] = {}
        for segment, labels in segment_labels.items():
            observed_pairs = [
                (position, median(observations[label]))
                for label, position, _ in labels
                if observations[label]
            ]
            observed_medians = [value for _, value in observed_pairs if value is not None]
            first_column = min(item[2] for item in labels)
            segment_rank[segment] = (
                median(observed_medians)
                if observed_medians
                else float("inf"),
                first_column,
            )
            if len(observed_pairs) >= 2:
                first = min(observed_pairs, key=lambda item: item[0])
                last = max(observed_pairs, key=lambda item: item[0])
                segment_direction[segment] = 1 if last[1] >= first[1] else -1
            else:
                segment_direction[segment] = 1

        ordered_segments = sorted(segment_labels, key=lambda segment: segment_rank[segment])
        segment_index = {segment: index for index, segment in enumerate(ordered_segments)}
        ordinal_index = {
            int(segment): index
            for segment, index in segment_index.items()
            if segment.isdigit() and 1 <= int(segment) <= 9
        }
        if not ordinal_index:
            ordinal_index = {
                ordinal: index
                for ordinal, index in enumerate(range(len(ordered_segments)), start=1)
            }

        def label_key(item: Tuple[str, int]):
            label, original_index = item
            flexible = parse_directional_flexible_grn(label)
            if flexible is not None:
                nearer, farther, distance = flexible
                if nearer in ordinal_index and farther in ordinal_index:
                    nearer_index = ordinal_index[nearer]
                    farther_index = ordinal_index[farther]
                    left_index = min(nearer_index, farther_index)
                    from_left = nearer_index == left_index
                    return (
                        2 * left_index + 1,
                        0 if from_left else 1,
                        distance if from_left else -distance,
                        original_index,
                    )
            parsed_label = cls._split_segment_position(label)
            if parsed_label is None:
                observed = median(observations[label])
                return (2 * len(segment_index), 0, observed or original_index, original_index)
            segment, position = parsed_label
            position_key = segment_direction[segment] * position
            return 2 * segment_index[segment], 0, position_key, original_index

        return [
            label
            for _, label in sorted(
                enumerate(columns), key=lambda item: label_key((item[1], item[0]))
            )
        ]

    def _load_region_configuration(
        self,
        protein_family: str,
    ) -> Tuple[Dict[str, Tuple[str, str]], Dict[str, Tuple[str, str]]]:
        """Load optional standard/strict non-flexible extents for a family."""

        family_aliases = {
            "gpcrdb_class_a": "gpcr_a",
        }
        config_family = family_aliases.get(protein_family, protein_family)
        try:
            config_path = self.configs_dir / "config.json"
            if not config_path.exists():
                return {}, {}
            family = json.loads(config_path.read_text()).get(config_family, {})
        except (AttributeError, OSError, TypeError, ValueError):
            return {}, {}

        def intervals(name: str) -> Dict[str, Tuple[str, str]]:
            result: Dict[str, Tuple[str, str]] = {}
            for segment, bounds in family.get(name, {}).items():
                if isinstance(bounds, list) and len(bounds) == 2:
                    result[str(segment)] = (str(bounds[0]), str(bounds[1]))
            return result

        return intervals("standard"), intervals("strict")

    @staticmethod
    def _merge_ordered_labels(existing: List[str], ordered_new: List[str]) -> None:
        """Merge a per-event label order without scrambling earlier events."""

        for index, label in enumerate(ordered_new):
            if label in existing:
                continue
            previous = next(
                (candidate for candidate in reversed(ordered_new[:index]) if candidate in existing),
                None,
            )
            following = next(
                (candidate for candidate in ordered_new[index + 1:] if candidate in existing),
                None,
            )
            if previous is not None:
                existing.insert(existing.index(previous) + 1, label)
            elif following is not None:
                existing.insert(existing.index(following), label)
            else:
                existing.append(label)

    def _find_best_reference(
        self,
        seq_name: str,
        sequence: str,
        reference_sequences: Dict[str, str],
        aligner: SequenceAlignmentEngine,
    ) -> Tuple[Optional[str], Optional[Any], float]:
        best_ref = None
        best_alignment = None
        best_score = float("-inf")

        # Exact matches are common in reference integrity tests and avoid an
        # unnecessary all-vs-all scan without changing selection semantics.
        for ref_name, ref_sequence in reference_sequences.items():
            if sequence == ref_sequence:
                result = aligner.align_pairwise(seq_name, sequence, ref_name, ref_sequence)
                norm_score = result.score / max(len(sequence), 1)
                return ref_name, result, norm_score

        for ref_name, ref_sequence in reference_sequences.items():
            if not ref_sequence:
                continue
            result = aligner.align_pairwise(seq_name, sequence, ref_name, ref_sequence)
            norm_score = result.score / max(len(sequence), len(ref_sequence), 1)
            if norm_score > best_score:
                best_score = norm_score
                best_ref = ref_name
                best_alignment = result
        return best_ref, best_alignment, best_score

    def _project_alignment(
        self,
        seq_name: str,
        alignment_result,
        reference_grns: List[str],
        column_order: List[str],
        *,
        sequence_grn_order: Optional[List[str]] = None,
        assign_unambiguous_insertions: bool = True,
        standard_regions: Optional[Dict[str, Tuple[str, str]]] = None,
        strict_regions: Optional[Dict[str, Tuple[str, str]]] = None,
        max_alignment_gap: int = 1,
    ) -> Tuple[pd.Series, float, int, Dict[str, Any]]:
        query_alignment, _, reference_alignment = alignment_result.alignment

        row_data = {col: '-' for col in column_order}
        reference_index = 0
        query_index = 0  # Track position in query sequence
        assigned = 0
        insertion_events: List[Dict[str, Any]] = []
        deletion_events: List[Dict[str, Any]] = []
        active_insertion: Optional[Dict[str, Any]] = None
        active_deletion: Optional[Dict[str, Any]] = None
        previous_grn: Optional[str] = None

        total_reference = len(reference_grns)

        for q_char, r_char in zip(query_alignment, reference_alignment):
            if r_char != '-':
                if reference_index < total_reference:
                    current_grn = reference_grns[reference_index]
                else:
                    current_grn = None
            else:
                current_grn = None

            if q_char != '-':
                query_index += 1  # Increment for non-gap positions
                
            if q_char != '-' and r_char != '-' and current_grn is not None:
                # Store amino acid with position (e.g., 'M1', 'E2', etc.)
                row_data[current_grn] = f"{q_char}{query_index}"
                assigned += 1

            if q_char != '-' and r_char == '-':
                if active_insertion is None:
                    active_insertion = {
                        "after_grn": previous_grn,
                        "before_grn": None,
                        "query_start": query_index,
                        "query_end": query_index,
                        "residues": q_char,
                    }
                    insertion_events.append(active_insertion)
                else:
                    active_insertion["query_end"] = query_index
                    active_insertion["residues"] += q_char
            elif active_insertion is not None and current_grn is not None:
                active_insertion["before_grn"] = current_grn
                active_insertion = None

            if q_char == '-' and r_char != '-' and current_grn is not None:
                if active_deletion is None:
                    active_deletion = {
                        "grns": [current_grn],
                        "reference_residues": r_char,
                        "after_query_position": query_index,
                    }
                    deletion_events.append(active_deletion)
                else:
                    active_deletion["grns"].append(current_grn)
                    active_deletion["reference_residues"] += r_char
            else:
                active_deletion = None

            if r_char != '-':
                if current_grn is not None:
                    previous_grn = current_grn
                reference_index += 1

        ordered_grns = sequence_grn_order or column_order
        order_index = {grn: index for index, grn in enumerate(ordered_grns)}
        query_sequence = query_alignment.replace('-', '')
        standard_regions = standard_regions or {}
        strict_regions = strict_regions or {}
        ordered_segments = list(
            dict.fromkeys(
                parsed[0]
                for grn in ordered_grns
                if (parsed := self._split_segment_position(grn)) is not None
            )
        )
        insertion_assignments = 0
        terminal_assignments = 0
        for event in insertion_events:
            after = event["after_grn"]
            before = event["before_grn"]
            candidates: List[str] = []
            if after in order_index and before in order_index:
                left = order_index[after] + 1
                right = order_index[before]
                if left <= right:
                    candidates = [
                        grn
                        for grn in ordered_grns[left:right]
                        if row_data.get(grn, "-") == "-"
                    ]
            event["candidate_grns"] = candidates
            event["assignment"] = "unassigned"
            event["generated_grns"] = []
            if event["after_grn"] is None and event["before_grn"] is not None:
                event["kind"] = "n_terminal_overhang"
            elif event["after_grn"] is not None and event["before_grn"] is None:
                event["kind"] = "c_terminal_overhang"
            elif event["after_grn"] is None and event["before_grn"] is None:
                event["kind"] = "unaligned_query"
            else:
                event["kind"] = "internal_insertion"

            if (
                event["kind"] == "internal_insertion"
                and len(event["residues"]) > max_alignment_gap
                and self._anchors_share_region(
                    after,
                    before,
                    strict_regions,
                    order_index,
                )
            ):
                corrections = self._compress_alignment_gap(
                    row_data,
                    event,
                    query_sequence,
                    ordered_grns,
                    order_index,
                    standard_regions,
                    strict_regions,
                )
                if corrections:
                    event["kind"] = "alignment_gap_compressed"
                    event["assignment"] = "compressed_to_standard_coordinates"
                    event["corrections"] = corrections
                    continue

            if (
                event["kind"] == "internal_insertion"
                and assign_unambiguous_insertions
                and len(candidates) == len(event["residues"])
            ):
                for offset, (grn, residue) in enumerate(zip(candidates, event["residues"])):
                    row_data[grn] = f"{residue}{event['query_start'] + offset}"
                    assigned += 1
                    insertion_assignments += 1
                event["assignment"] = "assigned_exact_candidate_count"
                continue

            if event["kind"] == "n_terminal_overhang":
                length = len(event["residues"])
                generated = [f"n.{distance}" for distance in range(length, 0, -1)]
                for offset, (grn, residue) in enumerate(zip(generated, event["residues"])):
                    row_data[grn] = f"{residue}{event['query_start'] + offset}"
                event["generated_grns"] = generated
                event["assignment"] = "assigned_terminal_coordinates"
                assigned += length
                terminal_assignments += length
                continue
            if event["kind"] == "c_terminal_overhang":
                length = len(event["residues"])
                generated = [f"c.{distance}" for distance in range(1, length + 1)]
                for offset, (grn, residue) in enumerate(zip(generated, event["residues"])):
                    row_data[grn] = f"{residue}{event['query_start'] + offset}"
                event["generated_grns"] = generated
                event["assignment"] = "assigned_terminal_coordinates"
                assigned += length
                terminal_assignments += length
                continue
            if event["kind"] != "internal_insertion":
                event["unassigned_reason"] = "unaligned_run"
                continue
            if not assign_unambiguous_insertions:
                event["unassigned_reason"] = "insertion_annotation_disabled"
                continue

            after_parsed = self._split_segment_position(after) if after else None
            before_parsed = self._split_segment_position(before) if before else None
            generated: List[str] = []
            if after_parsed and before_parsed and after_parsed[0] == before_parsed[0]:
                segment, position = after.rsplit('.', 1)
                anchor = f"{segment}.{position[:2]}"
                for insertion_index in range(1, 10):
                    try:
                        candidate = make_insertion_grn(anchor, insertion_index)
                    except ValueError:
                        generated = []
                        break
                    if row_data.get(candidate, '-') == '-':
                        generated.append(candidate)
                    if len(generated) == len(event["residues"]):
                        break
            elif after_parsed and before_parsed:
                near_n = self._segment_ordinal(after_parsed[0], ordered_segments)
                near_c = self._segment_ordinal(before_parsed[0], ordered_segments)
                if near_n is not None and near_c is not None:
                    n_count = (len(event["residues"]) + 1) // 2
                    c_count = len(event["residues"]) - n_count
                    generated.extend(
                        make_directional_flexible_grn(near_n, near_c, distance)
                        for distance in range(1, n_count + 1)
                    )
                    generated.extend(
                        make_directional_flexible_grn(near_c, near_n, distance)
                        for distance in range(c_count, 0, -1)
                    )

            if len(generated) == len(event["residues"]):
                if all(row_data.get(grn, '-') == '-' for grn in generated):
                    for offset, (grn, residue) in enumerate(zip(generated, event["residues"])):
                        row_data[grn] = f"{residue}{event['query_start'] + offset}"
                    assigned += len(generated)
                    insertion_assignments += len(generated)
                    event["generated_grns"] = generated
                    event["assignment"] = "assigned_generated_coordinates"
                else:
                    event["unassigned_reason"] = "generated_coordinate_collision"
            else:
                event["unassigned_reason"] = "coordinate_capacity_or_segment_mapping"

        projectable_positions = total_reference + insertion_assignments + terminal_assignments
        coverage = assigned / projectable_positions if projectable_positions else 0.0
        series = pd.Series(row_data, name=seq_name)
        indels = {
            "insertions": insertion_events,
            "deletions": deletion_events,
            "insertion_residues": sum(
                len(event["residues"])
                for event in insertion_events
                if event["kind"] == "internal_insertion"
            ),
            "terminal_overhang_residues": sum(
                len(event["residues"])
                for event in insertion_events
                if event["kind"] in {"n_terminal_overhang", "c_terminal_overhang"}
            ),
            "deletion_residues": sum(len(event["grns"]) for event in deletion_events),
        }
        return series, coverage, assigned, indels

    @staticmethod
    def _segment_ordinal(segment: str, ordered_segments: List[str]) -> Optional[int]:
        """Map a real segment identifier to its compact 1-based sequence index."""

        if segment.isdigit() and 1 <= int(segment) <= 9:
            return int(segment)
        try:
            ordinal = ordered_segments.index(segment) + 1
        except ValueError:
            return None
        return ordinal if ordinal <= 9 else None

    @staticmethod
    def _anchors_share_region(
        after: Optional[str],
        before: Optional[str],
        regions: Dict[str, Tuple[str, str]],
        order_index: Dict[str, int],
    ) -> Optional[str]:
        """Return the configured region containing both alignment anchors."""

        if after not in order_index or before not in order_index:
            return None
        left, right = sorted((order_index[after], order_index[before]))
        for name, (start, end) in regions.items():
            if start not in order_index or end not in order_index:
                continue
            region_left, region_right = sorted((order_index[start], order_index[end]))
            if region_left <= left <= right <= region_right:
                return name
        return None

    @classmethod
    def _compress_alignment_gap(
        cls,
        row_data: Dict[str, str],
        event: Dict[str, Any],
        query_sequence: str,
        ordered_grns: List[str],
        order_index: Dict[str, int],
        standard_regions: Dict[str, Tuple[str, str]],
        strict_regions: Dict[str, Tuple[str, str]],
    ) -> List[Dict[str, str]]:
        """Shift projected residues left across a suspect long alignment gap."""

        region = cls._anchors_share_region(
            event.get("after_grn"),
            event.get("before_grn"),
            strict_regions,
            order_index,
        )
        if region is None:
            return []
        bounds = standard_regions.get(region) or strict_regions.get(region)
        if (
            not bounds
            or bounds[0] not in order_index
            or bounds[1] not in order_index
            or event["before_grn"] not in order_index
        ):
            return []
        start_index = order_index[event["before_grn"]]
        region_left, end_index = sorted((order_index[bounds[0]], order_index[bounds[1]]))
        if not region_left <= start_index <= end_index:
            return []

        shift = len(event["residues"])
        after_position = event["query_start"] - 1
        corrections: List[Dict[str, str]] = []
        replacements: List[Tuple[str, str]] = []
        for grn in ordered_grns[start_index:end_index + 1]:
            old_value = row_data.get(grn, '-')
            parsed = cls.parse_grn_value(old_value)
            if parsed is None:
                continue
            _, old_position = parsed
            new_position = old_position - shift
            if new_position <= after_position or not 1 <= new_position <= len(query_sequence):
                return []
            new_value = f"{query_sequence[new_position - 1]}{new_position}"
            replacements.append((grn, new_value))
            corrections.append({"grn": grn, "from": old_value, "to": new_value})
        for grn, new_value in replacements:
            row_data[grn] = new_value
        return corrections

    @staticmethod
    def _sanitize_sequence(sequence: Any) -> str:
        if sequence is None:
            return ""
        if not isinstance(sequence, str):
            raise ValueError(f"Sequence must be a string, got {type(sequence).__name__}")
        cleaned = re.sub(r"[\s.\-]", "", sequence).upper()
        invalid = sorted(set(cleaned) - set("ACDEFGHIKLMNPQRSTVWYBJZX*"))
        if invalid:
            raise ValueError(f"Unsupported residue symbols: {''.join(invalid)}")
        return cleaned

    def _ensure_entity_exists(self, name: str, *, default_format: str, allow_create: bool) -> bool:
        entity_info = self.entity_registry.find_entity(name, default_format)
        if not entity_info:
            entity_info = self.entity_registry.find_entity(name)
        if entity_info:
            return True
        if allow_create:
            self.entity_registry.register_entity(
                name=name,
                format_type=default_format,
                file_path="",
                metadata={"placeholder": True},
            )
            return True
        return False
