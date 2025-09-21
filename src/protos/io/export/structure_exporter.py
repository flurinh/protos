"""Structure exporter implementation.

Handles filesystem exports for StructureProcessor instances while keeping the
processor itself path-agnostic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from protos.io.core.base_exporter import BaseExporter
from protos.io.formats.cif_utils import write_cif_file


class StructureExporter(BaseExporter):
    """Exporter that serializes structures to CIF files."""

    def export_entity(
        self,
        name: str,
        out_path: Path,
        format: Optional[str] = None,
        overwrite: bool = True,
        chains: Optional[list[str]] = None,
    ) -> Path:
        export_format = format or "cif"
        if export_format != "cif":
            raise ValueError("StructureExporter supports only 'cif' exports")

        df = self._get_structure(name)
        if chains:
            df = self._filter_by_chains(df, chains)

        out_path = Path(out_path)
        if out_path.exists() and not overwrite:
            raise FileExistsError(
                f"File {out_path} already exists. Use overwrite=True to replace it."
            )

        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_cif_file(str(out_path), df, force_overwrite=overwrite)
        self.logger.info("Exported structure '%s' to %s", name, out_path)
        return out_path

    def _get_structure(self, name: str) -> pd.DataFrame:
        df = self.processor.frames.get(name)
        if df is None:
            df = self.processor.load_entity(name)
        if df is None:
            raise ValueError(f"Structure '{name}' not found")
        if not isinstance(df, pd.DataFrame):
            raise TypeError("StructureProcessor must return DataFrame instances")
        return df.reset_index()

    def _filter_by_chains(self, df: pd.DataFrame, chains: list[str]) -> pd.DataFrame:
        chains_set = set(chains)
        mask = df['auth_chain_id'].isin(chains_set)
        if 'label_chain_id' in df.columns:
            mask |= df['label_chain_id'].isin(chains_set)
        filtered = df.loc[mask]
        if filtered.empty:
            self.logger.warning("No atoms found for chains %s; exporting full structure", chains)
            return df
        return filtered
