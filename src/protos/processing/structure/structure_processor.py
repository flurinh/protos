from __future__ import annotations

import os
import time
import json
import warnings
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
import requests

from protos.core.base_processor import BaseProcessor
from protos.io import cif_utils
from protos.processing.structure.struct_utils import (
    load_structure as load_structure_util,
    STRUCT_COLUMN_DTYPE,
)

# ---------- helpers (tiny, local) ----------

def _safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-._" else "_" for c in s)

def _is_multiindexed(df: pd.DataFrame) -> bool:
    return isinstance(df.index, pd.MultiIndex) and list(df.index.names) == ["structure_id", "atom_id"]

class StructureProcessor(BaseProcessor):
    """
    IO + Storage only:
      - Single IO entry points (load_entity / save_entity / download_structure)
      - Cache priority: PKL cache > CIF file > Registry
      - Every DataFrame indexed as MultiIndex (structure_id, atom_id)
      - Datasets load to a single stacked MultiIndex DataFrame
      - Per-structure PKL caching (datasets are not materialized as PKL by default)
    """

    # ---------- paths ----------

    @property
    def path_structure_dir(self) -> Path:
        return self.get_subdirectory_path("structure_dir")

    @property
    def path_dataset_dir(self) -> Path:
        return self.get_subdirectory_path("dataset_dir")

    @property
    def path_cache_dir(self) -> Path:
        return self.get_subdirectory_path("cache_dir")

    @property
    def temp_dir(self) -> Path:
        return self.get_subdirectory_path("temp_dir")

    # ---------- init / storage state ----------

    def __init__(
        self,
        name: str = "structure_processor",
        paths=None,
        *,
        preload_ids: bool = False,
        limit: Optional[int] = None,
        remove_hetatm: bool = False,
        allow_exception: bool = False,
    ):
        super().__init__(name=name, paths=paths)

        self.limit = limit
        self.remove_hetatm = remove_hetatm
        self.allow_exception = allow_exception

        self.structure_ids: List[str] = []
        self.dfl: List[pd.DataFrame] = []   # each item is a single-structure DF with MultiIndex
        self.data: Optional[pd.DataFrame] = None  # stacked DF with MultiIndex (structure_id, atom_id)

        if preload_ids:
            self.structure_ids = self.get_available_structures()
            if self.limit is not None:
                self.structure_ids = self.structure_ids[: self.limit]

    # ---------- public IO (single entry) ----------

    def load_entity(self, structure_id: str) -> Optional[pd.DataFrame]:
        """
        Load structure as a MultiIndex DF. Priority: cache > CIF > registry.
        Also updates in-memory storage (replaces or inserts).
        """
        sid = str(structure_id)
        df = (
            self._load_from_cache(sid)
            or self._load_from_cif(sid)
            or self._load_from_registry(sid)
        )
        if df is None:
            return None

        df = self._ensure_schema(df, sid, drop_hetatm=self.remove_hetatm)
        df = self._ensure_multiindex(df, sid)

        self._update_storage(sid, df)
        return df

    def save_entity(
        self,
        structure_id: str,
        df: pd.DataFrame,
        *,
        format: str = "pkl",
        metadata: Optional[Dict] = None,
    ) -> str:
        """
        Save per-structure DF and register entity.
        format: 'pkl' | 'cif' | 'both'
        """
        sid = str(structure_id)
        df = self._ensure_schema(df, sid, drop_hetatm=False)
        df = self._ensure_multiindex(df, sid)

        saved: List[Tuple[str, Path]] = []

        if format in ("pkl", "both"):
            p = self._save_cache(sid, df)
            saved.append(("pkl", p))

        if format in ("cif", "both"):
            p = self._save_cif(sid, df)
            saved.append(("cif", p))

        # register primary file (first saved)
        primary_fmt, primary_path = saved[0]
        rel = primary_path.relative_to(self.paths.data_root)

        md = {
            **(metadata or {}),
            "formats_available": [fmt for fmt, _ in saved],
            "atom_count": df.shape[0],
            "chain_ids": sorted(df.reset_index().auth_chain_id.dropna().unique().tolist())
            if "auth_chain_id" in df.columns
            else [],
        }

        self.entity_registry.register_entity(
            name=sid,
            format_type=self.processor_type,
            file_path=str(rel),
            metadata=md,
        )

        # maintain memory mirror
        self._update_storage(sid, df)
        return sid

    def download_structure(
        self,
        structure_id: str,
        *,
        source: str = "rcsb",
        save_to_cache: bool = True,
        metadata: Optional[Dict] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Download, parse, index, save+register (pkl if requested + cif), and update storage.
        """
        sid = str(structure_id)
        try:
            tmp = self._download_to_temp(sid, source=source, metadata=metadata)
            df = self._parse_structure_file(tmp)
        finally:
            try:
                if "tmp" in locals() and tmp and Path(tmp).exists():
                    Path(tmp).unlink(missing_ok=True)
            except Exception:
                pass

        if df is None or df.empty:
            return None

        df = self._ensure_schema(df, sid, drop_hetatm=self.remove_hetatm)
        df = self._ensure_multiindex(df, sid)

        # always persist CIF; PKL optional
        fmt = "both" if save_to_cache else "cif"
        self.save_entity(
            sid,
            df,
            format=fmt,
            metadata={
                "source": source,
                "downloaded_at": datetime.utcnow().isoformat(),
                **(metadata or {}),
            },
        )
        return df

    # ---------- dataset IO ----------

    def load_dataset(self, dataset_name: str) -> Dict[str, pd.DataFrame]:
        """
        Load all entities listed in the dataset into memory.
        Returns {structure_id: per-structure DF}, and sets self.data as stacked DF (MultiIndex).
        """
        ds = self.get_dataset(dataset_name)
        if not ds:
            return {}
        ids: List[str] = ds.get("entities") or ds.get("pdb_ids") or ds.get("content") or []
        if self.limit is not None:
            ids = ids[: self.limit]

        loaded: Dict[str, pd.DataFrame] = {}
        for sid in ids:
            df = self.load_entity(sid)
            if df is not None and not df.empty:
                loaded[sid] = df

        # stacked already maintained in self.data by _update_storage
        return loaded

    def save_dataset(self, dataset_name: str, structure_ids: Optional[List[str]] = None, metadata: Optional[Dict] = None) -> str:
        """
        Datasets are logical; per our rule, saving a dataset does not write a dataset PKL.
        It persists the membership via the BaseProcessor dataset manager.
        """
        if structure_ids is None:
            structure_ids = list(self.structure_ids)
        # ensure registration exists for each id
        for sid in structure_ids:
            if not self.entity_exists(sid):
                self.logger.warning(f"Structure '{sid}' is not registered; saving anyway.")
        return super().create_dataset(dataset_name, structure_ids, metadata or {})

    # ---------- legacy name shims (kept short) ----------

    @property
    def pdb_ids(self) -> List[str]:
        warnings.warn("pdb_ids is deprecated, use structure_ids", DeprecationWarning, stacklevel=2)
        return self.structure_ids

    def load_structure(self, identifier: str, **kwargs) -> Optional[pd.DataFrame]:
        warnings.warn("load_structure is deprecated, use load_entity", DeprecationWarning, stacklevel=2)
        return self.load_entity(identifier)

    def save_structure(self, name: str, structure_df: pd.DataFrame, **kwargs) -> str:
        warnings.warn("save_structure is deprecated, use save_entity", DeprecationWarning, stacklevel=2)
        fmt = kwargs.pop("format", "pkl")
        return self.save_entity(name, structure_df, format=fmt, metadata=kwargs.get("metadata"))

    # ---------- private: cache / cif / registry ----------

    def _load_from_cache(self, structure_id: str) -> Optional[pd.DataFrame]:
        p = self.path_cache_dir / f"{_safe_name(structure_id)}.pkl"
        if not p.exists():
            return None
        try:
            df = pd.read_pickle(p)
        except Exception as e:
            self.logger.warning(f"Failed to read cache for {structure_id}: {e}")
            return None

        # auto-register if missing
        if not self.entity_exists(structure_id):
            self.entity_registry.register_entity(
                name=structure_id,
                format_type=self.processor_type,
                file_path=str(p.relative_to(self.paths.data_root)),
                metadata={"cached": True, "auto_discovered": True},
            )
        return df

    def _load_from_cif(self, structure_id: str) -> Optional[pd.DataFrame]:
        p = self.path_structure_dir / f"{_safe_name(structure_id)}.cif"
        if not p.exists():
            return None
        try:
            df = self._parse_cif_file(p)
        except Exception as e:
            self.logger.error(f"Failed to parse CIF for {structure_id}: {e}")
            return None

        # cache once parsed
        try:
            self._save_cache(structure_id, df)
        except Exception as e:
            self.logger.warning(f"Failed to write cache for {structure_id}: {e}")

        # auto-register if missing
        if not self.entity_exists(structure_id):
            self.entity_registry.register_entity(
                name=structure_id,
                format_type=self.processor_type,
                file_path=str(p.relative_to(self.paths.data_root)),
                metadata={"auto_discovered": True},
            )
        return df

    def _load_from_registry(self, structure_id: str) -> Optional[pd.DataFrame]:
        info = self.entity_registry.find_entity(structure_id, self.processor_type)
        if not info:
            return None
        file_path = Path(info.file_path)
        if not file_path.is_absolute():
            file_path = Path(self.paths.data_root) / file_path
        if not file_path.exists():
            return None

        if file_path.suffix.lower() == ".pkl":
            try:
                return pd.read_pickle(file_path)
            except Exception as e:
                self.logger.error(f"Failed to read registry PKL for {structure_id}: {e}")
                return None
        else:
            try:
                return self._parse_structure_file(file_path)
            except Exception as e:
                self.logger.error(f"Failed to parse registry file for {structure_id}: {e}")
                return None

    # ---------- private: save helpers ----------

    def _save_cache(self, structure_id: str, df: pd.DataFrame) -> Path:
        p = self.path_cache_dir / f"{_safe_name(structure_id)}.pkl"
        df.to_pickle(p)
        return p

    def _save_cif(self, structure_id: str, df: pd.DataFrame) -> Path:
        p = self.path_structure_dir / f"{_safe_name(structure_id)}.cif"
        cif_utils.dataframe_to_cif(df.reset_index(), str(p))
        return p

    # ---------- private: parsing ----------

    def _parse_cif_file(self, path: Path) -> pd.DataFrame:
        # struct_utils expects (pdb_id, folder); use filename stem & parent
        return load_structure_util(path.stem, folder=str(path.parent))

    def _parse_structure_file(self, path: Path) -> pd.DataFrame:
        ext = path.suffix.lower()
        if ext in (".cif", ".mmcif"):
            return self._parse_cif_file(path)
        elif ext == ".pkl":
            return pd.read_pickle(path)
        else:
            raise ValueError(f"Unsupported structure format: {ext}")

    def _download_to_temp(self, structure_id: str, *, source: str, metadata: Optional[Dict]) -> Path:
        if source == "rcsb":
            url = f"https://files.rcsb.org/download/{structure_id}.cif"
        elif source == "alphafold":
            url = f"https://alphafold.ebi.ac.uk/files/AF-{structure_id}-F1-model_v4.cif"
        elif source == "custom_url":
            if not metadata or "url" not in metadata:
                raise ValueError("custom_url requires metadata['url']")
            url = metadata["url"]
        else:
            raise ValueError(f"Unknown source '{source}'")

        r = requests.get(url, timeout=60)
        r.raise_for_status()
        tmp = self.temp_dir / f"{_safe_name(structure_id)}_{int(time.time())}.cif"
        tmp.write_bytes(r.content)
        return tmp

    # ---------- private: schema / typing / index ----------

    def _ensure_schema(self, df: pd.DataFrame, structure_id: str, *, drop_hetatm: bool) -> pd.DataFrame:
        # add structure_id column if needed
        if "structure_id" not in df.columns:
            if "pdb_id" in df.columns:
                df = df.copy()
                df["structure_id"] = df["pdb_id"]
            else:
                df = df.copy()
                df["structure_id"] = structure_id

        # standardize atom id column name if needed (expect 'atom_id' in struct_utils)
        if "atom_id" not in df.columns:
            raise ValueError("Structure DataFrame must contain 'atom_id' column")

        # optional HETATM removal if the parsed schema has 'group'
        if drop_hetatm and "group" in df.columns:
            df = df[df["group"] == "ATOM"].copy()

        # dtypes
        for col, dtype in STRUCT_COLUMN_DTYPE.items():
            if col in df.columns:
                try:
                    df[col] = df[col].astype(dtype)
                except Exception:
                    pass

        # coords sanity
        for c in ("x", "y", "z"):
            if c in df.columns:
                df[c] = pd.to_numeric(df[c], errors="coerce")

        return df

    def _ensure_multiindex(self, df: pd.DataFrame, structure_id: str) -> pd.DataFrame:
        if _is_multiindexed(df):
            return df

        if "structure_id" not in df.columns:
            df = df.copy()
            df["structure_id"] = structure_id

        if "atom_id" not in df.columns:
            raise ValueError("Structure DataFrame must contain 'atom_id' to build MultiIndex")

        df = df.copy()
        df.sort_values(["structure_id", "atom_id"], inplace=True)
        df.set_index(["structure_id", "atom_id"], inplace=True)
        df.index = df.index.set_names(["structure_id", "atom_id"])
        return df

    # ---------- private: storage updates ----------

    def _update_storage(self, structure_id: str, df: pd.DataFrame) -> None:
        if not _is_multiindexed(df):
            raise ValueError("Internal storage requires MultiIndex (structure_id, atom_id)")

        if structure_id in self.structure_ids:
            pos = self.structure_ids.index(structure_id)
            self.dfl[pos] = df
        else:
            self.structure_ids.append(structure_id)
            self.dfl.append(df)

        self._rebuild_stacked()

    def _rebuild_stacked(self) -> None:
        if not self.dfl:
            self.data = None
            return
        # concat preserves MultiIndex; verify disjoint structure_id values
        self.data = pd.concat(self.dfl, axis=0).sort_index()

    # ---------- discovery / utilities ----------

    def get_available_structures(self) -> List[str]:
        if not self.path_structure_dir.exists():
            return []
        ids = [p.stem for p in self.path_structure_dir.glob("*.cif")]
        ids += [p.stem for p in self.path_cache_dir.glob("*.pkl")]
        return sorted(list(dict.fromkeys(ids)))  # stable unique

    # ---------- bulk helpers ----------

    def load_structures(self, structure_ids: List[str]) -> pd.DataFrame:
        warnings.warn("load_structures is deprecated; prefer load_dataset or loop load_entity", DeprecationWarning, stacklevel=2)
        for sid in structure_ids:
            self.load_entity(sid)
        return self.data if self.data is not None else pd.DataFrame()

    # ---------- persistence for stacked data (optional convenience) ----------

    def save_data(self, dataset_id: str, data: Optional[pd.DataFrame] = None, *, file_format: str = "pkl") -> str:
        """
        Persist the current stacked DataFrame to structure_dataset/. Optional convenience.
        Does not change per-structure caching policy.
        """
        if data is None:
            if self.data is None:
                raise ValueError("No data to save")
            data = self.data

        if not _is_multiindexed(data):
            raise ValueError("save_data expects a stacked MultiIndex DataFrame")

        out = self.path_dataset_dir / f"{dataset_id}.{file_format}"
        if file_format == "pkl":
            data.to_pickle(out)
        elif file_format == "csv":
            data.reset_index().to_csv(out, index=False)
        else:
            raise ValueError(f"Unsupported format: {file_format}")

        self.logger.info(f"Saved dataset '{dataset_id}' to {out}")
        return str(out)

    def load_data(self, dataset_id: str, *, file_format: str = "pkl") -> Optional[pd.DataFrame]:
        """
        Load a previously persisted stacked dataset from structure_dataset/.
        This is orthogonal to load_dataset() which resolves membership and per-structure IO.
        """
        p = self.path_dataset_dir / f"{dataset_id}.{file_format}"
        if not p.exists():
            self.logger.warning(f"Dataset file not found: {p}")
            return None

        if file_format == "pkl":
            df = pd.read_pickle(p)
        elif file_format == "csv":
            df = pd.read_csv(p)
        else:
            raise ValueError(f"Unsupported format: {file_format}")

        if "structure_id" in df.columns and "atom_id" in df.columns:
            df = self._ensure_multiindex(df, "<unknown>")
        elif not _is_multiindexed(df):
            raise ValueError("Loaded dataset must contain structure_id/atom_id columns or MultiIndex")

        # refresh in-memory mirror from stacked df
        self._ingest_stacked(df)
        return self.data

    def _ingest_stacked(self, stacked: pd.DataFrame) -> None:
        self.structure_ids = []
        self.dfl = []
        self.data = None

        for sid, sdf in stacked.groupby(level=0, sort=False):
            self.structure_ids.append(str(sid))
            self.dfl.append(sdf.sort_index())

        self._rebuild_stacked()
