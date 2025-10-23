"""GRNProcessor: manages GRN reference tables and sequence annotations."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.analysis.sequence.alignment_engine import SequenceAlignmentEngine
from protos.processing.grn.grn_utils import normalize_grn_format, validate_grn_string

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
    ) -> Tuple[pd.DataFrame, Dict[str, Any]]:
        """Annotate sequences with GRN positions using a reference alignment."""

        if not sequences:
            raise ValueError("No sequences provided for GRN annotation")

        reference_df = self.load_reference_table(reference_table)
        ref_sequences, ref_grn_map = self._prepare_reference_sequences(reference_df)

        aligner = SequenceAlignmentEngine()
        rows: List[pd.Series] = []
        per_sequence: Dict[str, Dict[str, Any]] = {}
        column_order = reference_df.columns.tolist()

        for seq_name, sequence in sequences.items():
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

            series, coverage, assigned = self._project_alignment(
                seq_name,
                best_alignment,
                ref_grn_map[best_ref],
                column_order,
            )

            series.name = seq_name
            rows.append(series)
            per_sequence[seq_name] = {
                "reference": best_ref,
                "normalized_score": best_score,
                "coverage": coverage,
                "assigned_positions": assigned,
                "status": "ok" if coverage > 0 else "no_overlap",
            }

        annotations = pd.DataFrame(rows)
        annotations.index.name = GRN_INDEX
        annotations = annotations.fillna('-')

        summary = {
            "per_sequence": per_sequence,
            "global": {
                "annotated": sum(1 for info in per_sequence.values() if info["coverage"] > 0),
                "total": len(per_sequence),
                "reference_table": reference_table,
                "protein_family": protein_family,
                "search": search,
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
            residues: List[str] = []
            grns: List[str] = []
            for grn, value in row.items():
                if isinstance(value, str) and value and value != '-':
                    residues.append(value[0])
                    grns.append(grn)
            sequences[ref_name] = ''.join(residues)
            residue_maps[ref_name] = grns

        return sequences, residue_maps

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
    ) -> Tuple[pd.Series, float, int]:
        query_alignment, _, reference_alignment = alignment_result.alignment

        row_data = {col: '-' for col in column_order}
        reference_index = 0
        query_index = 0  # Track position in query sequence
        assigned = 0

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

            if r_char != '-':
                reference_index += 1

        coverage = assigned / total_reference if total_reference else 0.0
        series = pd.Series(row_data, name=seq_name)
        return series, coverage, assigned

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
