"""PropertyProcessor: registry-aware tabular property storage.

This processor treats property tables as datasets (CSV files) and
records each row as a property entry that annotates an existing entity.
"""

from __future__ import annotations

import uuid
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

import pandas as pd
import json

from protos.io.core.base_processor import BaseProcessor


PROPERTY_ENTRY_ID = "property_entry_id"
ENTITY_NAME = "entity_name"
SCOPE_COLUMN = "scope"


def _normalize_scope(value: Any) -> List[Dict[str, str]]:
    if isinstance(value, str):
        value = json.loads(value)
    if not isinstance(value, Iterable):
        raise ValueError("scope must be a list of {format, name} mappings")
    normalized: List[Dict[str, str]] = []
    for item in value:
        if not isinstance(item, dict):
            raise ValueError("Each scope entry must be a dict with 'format' and 'name'")
        fmt = item.get("format")
        name = item.get("name")
        if not name:
            raise ValueError("Scope entry missing 'name'")
        normalized.append({"format": fmt, "name": name})
    if not normalized:
        raise ValueError("Scope list cannot be empty")
    return normalized


class PropertyProcessor(BaseProcessor):
    """Manage property tables and relationships to existing entities."""

    processor_type = "property"

    def __init__(self, name: str = "property_processor") -> None:
        super().__init__(name=name)
        self._table_cache: Dict[str, pd.DataFrame] = {}

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------
    @property
    def tables_dir(self) -> Path:
        return Path(self.get_subdirectory_path("tables_dir"))

    @property
    def datasets_dir(self) -> Path:
        return Path(self.get_subdirectory_path("datasets_dir"))

    def _table_path(self, table_name: str) -> Path:
        return self.tables_dir / f"{table_name}.csv"

    def _relative_path(self, path: Path) -> str:
        return str(path.relative_to(self.paths.data_root))

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def record_properties(
        self,
        table_name: str,
        rows: Sequence[Dict[str, Any]] | pd.DataFrame,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        allow_create: bool = False,
        materialize_entries: bool = False,
    ) -> pd.DataFrame:
        """Insert rows into a property table and register relationships.

        When ``materialize_entries`` is ``False`` (default) the registry receives a
        single dataset-level entity that annotates each referenced structure or
        sequence. This keeps the registry lean even for residue-level tables.
        Set ``materialize_entries=True`` to retain the legacy behaviour of
        registering every row as its own entity.
        """

        if isinstance(rows, pd.DataFrame):
            new_df = rows.copy()
        else:
            new_df = pd.DataFrame(list(rows))

        if new_df.empty:
            raise ValueError("No property rows provided")
        if SCOPE_COLUMN not in new_df.columns:
            raise ValueError("Each property row must include 'scope'")

        table = self._load_table(table_name)
        existing_ids = set(
            table[PROPERTY_ENTRY_ID].dropna().astype(str)
        ) if not table.empty and PROPERTY_ENTRY_ID in table.columns else set()

        entry_ids: List[str] = []
        scopes: List[List[Dict[str, str]]] = []
        unique_targets: Dict[str, set[str]] = {}

        for idx, row in new_df.iterrows():
            scope_items = _normalize_scope(row[SCOPE_COLUMN])
            scopes.append(scope_items)

            entity_name = row.get(ENTITY_NAME) or scope_items[-1]["name"]
            new_df.at[idx, ENTITY_NAME] = entity_name

            for scope_item in scope_items:
                fmt = scope_item.get("format")
                name = scope_item["name"]
                entity_info = self.entity_registry.find_entity(name, fmt)
                if not entity_info:
                    if not allow_create:
                        raise ValueError(
                            f"Entity '{name}' (format={fmt}) not found; set allow_create=True to create placeholder"
                        )
                    self.entity_registry.register_entity(
                        name=name,
                        format_type=fmt or self.processor_type,
                        file_path="",
                        metadata={"placeholder": True},
                    )
                fmt_key = fmt or "unknown"
                unique_targets.setdefault(fmt_key, set()).add(name)

            if materialize_entries:
                entry_id = row.get(PROPERTY_ENTRY_ID)
                if not entry_id or str(entry_id).strip() == "":
                    entry_id = self._generate_entry_id(table_name)
                entry_id = str(entry_id)
                while entry_id in existing_ids or entry_id in entry_ids:
                    entry_id = self._generate_entry_id(table_name)

                entry_ids.append(entry_id)
                new_df.at[idx, PROPERTY_ENTRY_ID] = entry_id

        if not materialize_entries:
            new_df[PROPERTY_ENTRY_ID] = pd.Series([None] * len(new_df))

        new_df[SCOPE_COLUMN] = scopes
        updated = pd.concat([table, new_df], ignore_index=True)
        if materialize_entries and PROPERTY_ENTRY_ID in updated.columns:
            updated[PROPERTY_ENTRY_ID] = updated[PROPERTY_ENTRY_ID].astype(str)
        self._write_table(table_name, updated)

        artifact_rel = self._relative_path(self._table_path(table_name))

        registered_entry_ids: List[str] = []
        if materialize_entries:
            for entry_id, scope_items in zip(entry_ids, scopes):
                row_index = updated.index[updated[PROPERTY_ENTRY_ID] == entry_id][0]
                self.entity_registry.register_entity(
                    name=entry_id,
                    format_type=self.processor_type,
                    file_path=artifact_rel,
                    metadata={"table": table_name, "row_index": int(row_index)},
                )
                for scope_index, scope_item in enumerate(scope_items):
                    self.entity_registry.add_relationship(
                        source_name=entry_id,
                        target_name=scope_item["name"],
                        rel_type="annotated_by",
                        metadata={
                            "table": table_name,
                            "row_index": int(row_index),
                            "scope_index": scope_index,
                            "scope_format": scope_item.get("format"),
                        },
                    )
                registered_entry_ids.append(entry_id)

        dataset_metadata = {
            "artifact_path": artifact_rel,
            "row_count": len(updated),
            "columns": updated.columns.tolist(),
            "materialize_entries": materialize_entries,
        }

        index_rel_path: Optional[str] = None
        if not materialize_entries:
            index_path = self.datasets_dir / f"{table_name}__index.json"
            index_data: Dict[str, Dict[str, List[int]]] = {}
            if index_path.exists():
                with open(index_path, "r", encoding="utf-8") as handle:
                    index_data = json.load(handle)

            start_idx = len(table)
            for offset, scope_items in enumerate(scopes):
                row_index = start_idx + offset
                for scope_item in scope_items:
                    fmt = scope_item.get("format") or "unknown"
                    name = scope_item["name"]
                    entries = index_data.setdefault(fmt, {}).setdefault(name, [])
                    entries.append(int(row_index))

            for fmt_key, fmt_mapping in index_data.items():
                for entity_name, positions in fmt_mapping.items():
                    fmt_mapping[entity_name] = sorted(set(positions))

            index_path.parent.mkdir(parents=True, exist_ok=True)
            with open(index_path, "w", encoding="utf-8") as handle:
                json.dump(index_data, handle, indent=2)

            index_rel_path = str(index_path.relative_to(self.paths.data_root))
            dataset_metadata["index_artifact"] = index_rel_path

        if metadata:
            dataset_metadata.update(metadata)

        dataset_entity_name = table_name
        entity_metadata = {
            "table": table_name,
            "artifact_path": artifact_rel,
            "materialize_entries": materialize_entries,
        }
        if index_rel_path:
            entity_metadata["index_artifact"] = index_rel_path

        self.entity_registry.register_entity(
            name=dataset_entity_name,
            format_type=self.processor_type,
            file_path=artifact_rel,
            metadata=entity_metadata,
        )

        if not materialize_entries:
            for fmt_key, names in unique_targets.items():
                for target_name in names:
                    rel_metadata = {"table": table_name, "scope_format": fmt_key}
                    try:
                        self.entity_registry.add_relationship(
                            source_name=dataset_entity_name,
                            target_name=target_name,
                            rel_type="annotated_by",
                            metadata=rel_metadata,
                        )
                    except ValueError:
                        continue

        if self.dataset_manager.dataset_exists(table_name):
            if materialize_entries and registered_entry_ids:
                self.dataset_manager.add_to_dataset(table_name, registered_entry_ids)
            self.dataset_manager.update_metadata(table_name, dataset_metadata)
        else:
            entities_for_dataset = registered_entry_ids if materialize_entries else [dataset_entity_name]
            self.create_dataset(table_name, entities_for_dataset, dataset_metadata)

        return updated

    def load_table(self, table_name: str) -> pd.DataFrame:
        return self._load_table(table_name).copy()

    def get_properties(
        self,
        entity_name: str,
        *,
        table_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """Return property rows associated with an entity."""

        if table_name:
            table = self._load_table(table_name)
            if table.empty:
                return table
            return table[table[ENTITY_NAME] == entity_name].copy()

        relationships = self.entity_registry.get_relationships(
            entity_name,
            rel_type="annotated_by",
            direction="incoming",
        )
        rows: List[pd.DataFrame] = []
        for rel in relationships:
            rel_meta = rel.get("metadata", {})
            table = rel_meta.get("table")
            row_index = rel_meta.get("row_index")
            if table is None:
                continue
            table_df = self._load_table(table)
            if row_index is not None and 0 <= row_index < len(table_df):
                rows.append(table_df.iloc[[row_index]])
            else:
                subset = table_df[table_df[ENTITY_NAME] == entity_name]
                if not subset.empty:
                    rows.append(subset)

        if rows:
            return pd.concat(rows, ignore_index=True)
        return pd.DataFrame(columns=[PROPERTY_ENTRY_ID, ENTITY_NAME])

    def load_dataset_rows(
        self,
        table_name: str,
        entity_name: Optional[str] = None,
        *,
        format_type: Optional[str] = None,
    ) -> pd.DataFrame:
        """Load property rows for an optional target entity using the dataset index.

        Args:
            table_name: Property table name.
            entity_name: Optional entity name to filter rows by scope membership.
            format_type: Optional scope format filter (e.g. "structure").

        Returns:
            A DataFrame of matching rows. When ``entity_name`` is ``None`` the
            full table is returned.
        """

        table = self._load_table(table_name).copy()
        if entity_name is None:
            return table

        index_path = self.datasets_dir / f"{table_name}__index.json"
        candidate_rows: List[int] = []
        if index_path.exists():
            with open(index_path, "r", encoding="utf-8") as handle:
                index_data = json.load(handle)

            if format_type:
                candidate_rows.extend(index_data.get(format_type, {}).get(entity_name, []))
            else:
                for fmt_map in index_data.values():
                    candidate_rows.extend(fmt_map.get(entity_name, []))

        if candidate_rows:
            unique_idx = sorted(set(candidate_rows))
            return table.iloc[unique_idx].reset_index(drop=True)

        # Fallback: scan table scopes
        matched_idx: List[int] = []
        for idx, scopes in table[SCOPE_COLUMN].items():
            for scope in scopes:
                fmt = scope.get("format")
                name = scope.get("name")
                if name == entity_name and (format_type is None or fmt == format_type):
                    matched_idx.append(idx)
                    break

        if not matched_idx:
            return table.iloc[0:0]
        return table.iloc[matched_idx].reset_index(drop=True)

    def list_tables(self) -> List[str]:
        return self.dataset_manager.list_datasets()

    def list_datasets(self) -> List[str]:
        """Expose dataset listings via the dataset manager."""
        return super().list_datasets()

    # ------------------------------------------------------------------
    # Abstract method implementations
    # ------------------------------------------------------------------
    def load_entity(self, name: str) -> Optional[Dict[str, Any]]:
        entry_info = self.entity_registry.find_entity(name, self.processor_type)
        if not entry_info:
            return None
        metadata = entry_info.metadata or {}
        table_name = metadata.get("table")
        row_index = metadata.get("row_index")
        if table_name is None:
            return None
        table = self._load_table(table_name)
        if row_index is not None and 0 <= row_index < len(table):
            return table.iloc[int(row_index)].to_dict()
        match = table[table[PROPERTY_ENTRY_ID] == name]
        if not match.empty:
            return match.iloc[0].to_dict()
        return None

    def save_entity(
        self,
        name: str,
        data: Dict[str, Any],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        if not isinstance(data, dict):
            raise ValueError("PropertyProcessor.save_entity expects a dict")
        if ENTITY_NAME not in data:
            raise ValueError("Property entity data must include 'entity_name'")
        table_name = (metadata or {}).get("table", "default_properties")
        row = dict(data)
        row[PROPERTY_ENTRY_ID] = name
        self.record_properties(
            table_name,
            [row],
            metadata=metadata,
            allow_create=metadata.get("allow_create", False) if metadata else False,
        )

    # ------------------------------------------------------------------
    # Compatibility helpers
    # ------------------------------------------------------------------
    def create_property_table(
        self,
        table_name: str,
        data: pd.DataFrame | Sequence[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
        allow_create: bool = False,
    ) -> pd.DataFrame:
        table_file = self._table_path(table_name)
        if table_file.exists():
            table_file.unlink()
        self._table_cache.pop(table_name, None)
        return self.record_properties(
            table_name,
            data,
            metadata=metadata,
            allow_create=allow_create,
        )

    def load_property_table(self, table_name: str) -> pd.DataFrame:
        return self.load_table(table_name)

    def save_property_table(
        self,
        table_name: str,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> None:
        table = self._load_table(table_name)
        self._write_table(table_name, table)
        artifact_rel = self._relative_path(self._table_path(table_name))
        dataset_metadata = {
            "artifact_path": artifact_rel,
            "row_count": len(table),
            "columns": table.columns.tolist(),
        }
        if metadata:
            dataset_metadata.update(metadata)
        if self.dataset_manager.dataset_exists(table_name):
            self.dataset_manager.update_metadata(table_name, dataset_metadata)
        else:
            self.create_dataset(table_name, [], dataset_metadata)

    def add_property_column(
        self,
        table_name: str,
        property_name: str,
        values: Any,
    ) -> pd.DataFrame:
        table = self._load_table(table_name)
        if isinstance(values, pd.Series):
            table[property_name] = values
        elif isinstance(values, dict):
            table[property_name] = pd.Series(values)
        else:
            table[property_name] = values
        self._write_table(table_name, table)
        self.save_property_table(table_name)
        return table

    def filter_by_property(
        self,
        table_name: str,
        property_name: str,
        condition,
    ) -> pd.DataFrame:
        table = self._load_table(table_name)
        if property_name not in table.columns:
            raise ValueError(f"Property '{property_name}' not found in table '{table_name}'")
        return table[table[property_name].apply(condition)]

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _write_table(self, table_name: str, df: pd.DataFrame) -> None:
        storage_df = df.copy()
        if SCOPE_COLUMN in storage_df.columns:
            storage_df[SCOPE_COLUMN] = storage_df[SCOPE_COLUMN].apply(json.dumps)
        storage_df.to_csv(self._table_path(table_name), index=False)
        cache_df = df.copy()
        if SCOPE_COLUMN in cache_df.columns:
            cache_df[SCOPE_COLUMN] = cache_df[SCOPE_COLUMN].apply(_normalize_scope)
        self._table_cache[table_name] = cache_df

    # ------------------------------------------------------------------
    # Internal utilities
    # ------------------------------------------------------------------
    def _load_table(self, table_name: str) -> pd.DataFrame:
        if table_name in self._table_cache:
            return self._table_cache[table_name]

        table_path = self._table_path(table_name)
        if not table_path.exists():
            df = pd.DataFrame(columns=[PROPERTY_ENTRY_ID, ENTITY_NAME, SCOPE_COLUMN])
            self._table_cache[table_name] = df
            return df

        df = pd.read_csv(table_path)
        if SCOPE_COLUMN in df.columns:
            df[SCOPE_COLUMN] = df[SCOPE_COLUMN].apply(json.loads)
        if PROPERTY_ENTRY_ID not in df.columns:
            df.insert(0, PROPERTY_ENTRY_ID, pd.Series(dtype=str))
        if ENTITY_NAME not in df.columns:
            df.insert(1, ENTITY_NAME, pd.Series(dtype=str))
        self._table_cache[table_name] = df
        return df

    def _load_table_columns(self) -> List[str]:
        df = self._table_cache.get(next(iter(self._table_cache), ""))
        if df is not None:
            return df.columns.tolist()
        return [PROPERTY_ENTRY_ID, ENTITY_NAME, SCOPE_COLUMN]

    def _generate_entry_id(self, table_name: str) -> str:
        return f"{table_name}#{uuid.uuid4().hex[:8]}"
