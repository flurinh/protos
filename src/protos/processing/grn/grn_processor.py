"""Registry-aware GRNProcessor for tabular GRN annotations."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.processing.grn.grn_utils import (
    normalize_grn_format,
    validate_grn_string,
)

GRN_INDEX = "entity_name"
SCOPE_COLUMN = "scope"
GRN_COLUMNS = "grn_columns"
ARTIFACT_PATH = "artifact_path"


class GRNProcessor(BaseProcessor):
    """Stores and retrieves GRN tables aligned with the entity registry."""

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
    # Recording & loading tables
    # ------------------------------------------------------------------
    def record_table(
        self,
        table_name: str,
        table: pd.DataFrame,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        allow_create: bool = False,
    ) -> pd.DataFrame:
        """Persist a GRN table and register relationships."""

        if table.index.name != GRN_INDEX:
            table = table.copy()
            table.index.name = GRN_INDEX

        normalized_columns = self._normalize_columns(table.columns)
        table.columns = normalized_columns
        self._write_table(table_name, table)

        artifact_rel = self._relative_path(self._table_path(table_name))
        dataset_metadata = {
            ARTIFACT_PATH: artifact_rel,
            GRN_COLUMNS: list(table.columns),
            "row_count": len(table),
        }
        if metadata:
            dataset_metadata.update(metadata)

        if self.dataset_manager.dataset_exists(table_name):
            self.dataset_manager.update_metadata(table_name, dataset_metadata)
        else:
            self.create_dataset(table_name, [], dataset_metadata)

        self._register_table_entity(table_name, artifact_rel, dataset_metadata)

        for row_index, entity_name in enumerate(table.index):
            self._register_row(
                table_name,
                entity_name,
                row_index,
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
    # Relationships & lookups
    # ------------------------------------------------------------------
    def get_annotations(
        self,
        entity_name: str,
        *,
        table_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Return GRN annotations for a given entity."""

        if table_name:
            table = self._load_table(table_name)
            if table.empty:
                return table
            if entity_name in table.index:
                return table.loc[[entity_name]].copy()
            return pd.DataFrame(columns=table.columns)

        relationships = self.entity_registry.get_relationships(
            entity_name,
            rel_type="annotated_by",
            direction="incoming",
        )

        rows: List[pd.DataFrame] = []
        for rel in relationships:
            metadata = rel.get("metadata", {})
            table = metadata.get("table")
            row_index = metadata.get("row_index")
            if not table:
                continue
            table_df = self._load_table(table)
            if row_index is not None and entity_name in table_df.index:
                rows.append(table_df.loc[[entity_name]])
            elif entity_name in table_df.index:
                rows.append(table_df.loc[[entity_name]])

        if rows:
            return pd.concat(rows, axis=0)
        return pd.DataFrame(columns=self._load_table_columns())

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
    def _register_row(
        self,
        table_name: str,
        entity_name: str,
        row_index: int,
        *,
        allow_create: bool,
    ) -> None:
        entity_info = self.entity_registry.find_entity(entity_name, "sequence")
        if not entity_info:
            entity_info = self.entity_registry.find_entity(entity_name)
        if not entity_info:
            if not allow_create:
                raise ValueError(f"Entity '{entity_name}' not found in registry")
            self.entity_registry.register_entity(
                name=entity_name,
                format_type=self.processor_type,
                file_path="",
                metadata={"placeholder": True},
            )

        self.entity_registry.add_relationship(
            source_name=table_name,
            target_name=entity_name,
            rel_type="annotated_by",
            metadata={"table": table_name, "row_index": int(row_index)},
        )

    def _register_table_entity(
        self,
        table_name: str,
        artifact_rel: str,
        metadata: Dict[str, Any],
    ) -> None:
        existing = self.entity_registry.find_entity(table_name, self.processor_type)
        entity_metadata = dict(metadata)
        entity_metadata.setdefault("table", table_name)
        self.entity_registry.register_entity(
            name=table_name,
            format_type=self.processor_type,
            file_path=artifact_rel,
            metadata=entity_metadata,
        )

    def _load_table(self, table_name: str) -> pd.DataFrame:
        if table_name in self._table_cache:
            return self._table_cache[table_name]

        path = self._table_path(table_name)
        if not path.exists():
            df = pd.DataFrame(columns=[])
            self._table_cache[table_name] = df
            return df

        df = pd.read_csv(path, index_col=0)
        df.index.name = GRN_INDEX
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

    def _load_table_columns(self) -> List[str]:
        if self._table_cache:
            return next(iter(self._table_cache.values())).columns.tolist()
        return []
