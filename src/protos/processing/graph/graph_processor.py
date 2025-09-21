
"""GraphProcessor: manages structure-derived graph entities."""

from __future__ import annotations

import json
import math
import pickle
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from protos.io.core.base_processor import BaseProcessor
from protos.processing.structure import StructureProcessor

try:  # Optional PyTorch / PyG support
    import torch
except ImportError:  # pragma: no cover - optional dependency
    torch = None

try:  # pragma: no cover - optional dependency
    from torch_geometric.data import Data as PyGData
except ImportError:  # pragma: no cover - optional dependency
    PyGData = None


_ELEMENT_TO_Z: Dict[str, int] = {
    'H': 1,
    'C': 6,
    'N': 7,
    'O': 8,
    'P': 15,
    'S': 16,
    'SE': 34,
}


class GraphProcessor(BaseProcessor):
    """Processor for PyTorch Geometric compatible graphs derived from structures."""

    processor_type = "graph"

    def __init__(
        self,
        name: str = "graph_processor",
        *,
        default_cutoff: float = 5.0,
        default_level: str = "atom",
    ) -> None:
        super().__init__(name=name)

        self.graphs_dir = self.get_subdirectory_path('graphs_dir')
        self.datasets_dir = self.get_subdirectory_path('datasets_dir')
        self.graphs_dir.mkdir(parents=True, exist_ok=True)
        self.datasets_dir.mkdir(parents=True, exist_ok=True)

        self.default_cutoff = float(default_cutoff)
        self.default_level = default_level

        self._check_dependencies()

    # ------------------------------------------------------------------
    # Dependency helpers
    # ------------------------------------------------------------------
    def _check_dependencies(self) -> None:
        if torch is None or PyGData is None:
            self.logger.warning(
                "PyTorch Geometric is not available. Graphs will be saved as serialised "
                "NumPy dictionaries. Install with: pip install torch torch_geometric"
            )

    @property
    def dependencies_available(self) -> bool:
        return torch is not None and PyGData is not None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def generate_graph(
        self,
        structure_id: str,
        *,
        structure: Optional[pd.DataFrame] = None,
        structure_processor: Optional[StructureProcessor] = None,
        level: str = "atom",
        cutoff: Optional[float] = None,
        include_hydrogens: bool = False,
        dataset_name: Optional[str] = None,
    ) -> str:
        """Generate and persist a graph for a single structure."""

        cutoff = float(cutoff or self.default_cutoff)
        level = level or self.default_level

        if structure is None:
            struct_proc = structure_processor or StructureProcessor()
            structure = struct_proc.load_entity(structure_id)
        if structure is None:
            raise ValueError(f"Structure '{structure_id}' could not be loaded")

        graph_payload = self._structure_to_graph(
            structure_id,
            structure,
            level=level,
            cutoff=cutoff,
            include_hydrogens=include_hydrogens,
        )

        entity_name = self._persist_graph(
            structure_id,
            graph_payload,
            dataset_name=dataset_name,
        )

        return entity_name

    def generate_graphs_from_dataset(
        self,
        structure_dataset: str,
        *,
        structure_processor: Optional[StructureProcessor] = None,
        dataset_name: Optional[str] = None,
        level: str = "atom",
        cutoff: Optional[float] = None,
        include_hydrogens: bool = False,
    ) -> Tuple[str, List[str]]:
        """Generate graphs for every structure in a dataset."""

        cutoff = float(cutoff or self.default_cutoff)
        level = level or self.default_level
        struct_proc = structure_processor or StructureProcessor()

        structures = struct_proc.load_dataset(structure_dataset)
        if isinstance(structures, pd.DataFrame):
            structures = {structure_dataset: structures}
        if not isinstance(structures, dict) or len(structures) == 0:
            raise ValueError(f"Structure dataset '{structure_dataset}' is empty")

        graph_names: List[str] = []
        for structure_id, structure_df in structures.items():
            try:
                graph_name = self.generate_graph(
                    structure_id,
                    structure=structure_df,
                    structure_processor=struct_proc,
                    level=level,
                    cutoff=cutoff,
                    include_hydrogens=include_hydrogens,
                    dataset_name=None,  # postpone dataset registration until the end
                )
                graph_names.append(graph_name)
            except Exception as exc:
                self.logger.warning(
                    "Failed to generate graph for %s: %s", structure_id, exc
                )

        if not graph_names:
            raise RuntimeError("No graphs generated from the structure dataset")

        graph_dataset_name = dataset_name or self._sanitize_filename(
            f"{structure_dataset}__{level}__{cutoff:.1f}A"
        )

        metadata = {
            "source_structure_dataset": structure_dataset,
            "level": level,
            "cutoff": cutoff,
            "include_hydrogens": include_hydrogens,
            "entity_count": len(graph_names),
        }

        if self.dataset_manager.dataset_exists(graph_dataset_name):
            self.dataset_manager.delete_dataset(graph_dataset_name)

        self.create_dataset(graph_dataset_name, graph_names, metadata)

        return graph_dataset_name, graph_names

    # ------------------------------------------------------------------
    # BaseProcessor abstract implementations
    # ------------------------------------------------------------------
    def load_entity(self, name: str) -> Optional[Dict[str, Any]]:
        try:
            payload = self.load_graph(name)
            return payload
        except Exception:
            return None

    def save_entity(
        self,
        name: str,
        data: Dict[str, Any],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        graph_meta = dict(data.get("graph_metadata", {}))
        if metadata:
            graph_meta.update(metadata)

        structure_id = graph_meta.get("structure_id", name)
        dataset_name = graph_meta.get("dataset")
        entity_name = self._sanitize_filename(name)

        return self._persist_graph(
            structure_id,
            data,
            dataset_name=dataset_name,
            entity_name=entity_name,
        )

    def load_graph(self, entity_name: str) -> Dict[str, Any]:
        """Load a persisted graph entity."""

        entity_info = self.entity_registry.find_entity(entity_name, self.processor_type)
        if entity_info is None:
            raise ValueError(f"Graph entity '{entity_name}' not found")

        file_path = Path(self.paths.data_root) / entity_info.file_path
        if not file_path.exists():
            raise FileNotFoundError(f"Graph file not found: {file_path}")

        with open(file_path, "rb") as handle:
            return pickle.load(handle)

    def list_graphs(self, dataset: Optional[str] = None) -> List[str]:
        """List registered graph entities."""

        if dataset:
            return self.dataset_manager.get_dataset_entities(dataset)
        return self.entity_registry.list_entities(self.processor_type)

    def to_pyg(self, graph_dict: Dict[str, Any]) -> PyGData:
        """Convert a stored graph dictionary to a PyG Data object."""

        if not self.dependencies_available:
            raise RuntimeError("PyTorch Geometric is not installed")

        pos = torch.tensor(graph_dict["node_positions"], dtype=torch.float32)
        edge_index = torch.tensor(graph_dict["edge_index"], dtype=torch.long)
        edge_weight = torch.tensor(
            graph_dict.get("edge_weight"), dtype=torch.float32
        ).unsqueeze(1)

        x_data = graph_dict.get("node_features")
        if x_data is not None:
            x = torch.tensor(x_data, dtype=torch.float32)
        else:
            x = None

        pyg_data = PyGData(
            x=x,
            edge_index=edge_index,
            edge_attr=edge_weight,
            pos=pos,
        )
        return pyg_data

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _structure_to_graph(
        self,
        structure_id: str,
        structure_df: pd.DataFrame,
        *,
        level: str,
        cutoff: float,
        include_hydrogens: bool,
    ) -> Dict[str, Any]:
        """Convert a structure DataFrame into a graph payload."""

        df = structure_df.copy()

        coord_cols = [col for col in ("x_coord", "y_coord", "z_coord") if col in df.columns]
        if len(coord_cols) != 3:
            coord_cols = [col for col in ("x", "y", "z") if col in df.columns]
        if len(coord_cols) != 3:
            raise ValueError("Structure dataframe missing coordinate columns")

        element_col = None
        for candidate in ("type_symbol", "element", "atom_name"):
            if candidate in df.columns:
                element_col = candidate
                break
        if element_col is None:
            element_col = "element"
            df[element_col] = "C"  # default fallback

        if not include_hydrogens:
            df = df[df[element_col].str.upper() != "H"]

        if df.empty:
            raise ValueError("Structure has no atoms after filtering")

        if level == "residue":
            group_cols = [
                col
                for col in ("auth_asym_id", "auth_seq_id", "auth_comp_id")
                if col in df.columns
            ] or ["auth_seq_id"]
            grouped = df.groupby(group_cols, dropna=False)
            pos = grouped[coord_cols].mean().to_numpy(dtype=np.float32)

            elements = grouped[element_col].agg(
                lambda series: Counter(series).most_common(1)[0][0]
            ).tolist()

            node_metadata = [
                {
                    "chain": key[0] if isinstance(key, tuple) and len(key) > 0 else None,
                    "residue_number": key[1] if isinstance(key, tuple) and len(key) > 1 else None,
                    "residue_name": key[2] if isinstance(key, tuple) and len(key) > 2 else None,
                }
                for key in grouped.indices.keys()
            ]
        else:
            pos = df[coord_cols].to_numpy(dtype=np.float32)
            elements = df[element_col].tolist()
            node_metadata = [
                {
                    "atom_name": df.iloc[i].get("atom_name"),
                    "residue": df.iloc[i].get("auth_comp_id"),
                    "chain": df.iloc[i].get("auth_asym_id"),
                    "residue_number": df.iloc[i].get("auth_seq_id"),
                }
                for i in range(len(df))
            ]

        node_count = len(pos)
        if node_count < 1:
            raise ValueError("No nodes available for graph construction")

        atomic_numbers = np.array(
            [self._atomic_number(element) for element in elements], dtype=np.float32
        ).reshape(-1, 1)

        edge_index, edge_weight = self._build_edges(pos, cutoff=cutoff)

        graph_metadata = {
            "structure_id": structure_id,
            "level": level,
            "cutoff": cutoff,
            "include_hydrogens": include_hydrogens,
            "node_count": node_count,
            "edge_count": int(edge_index.shape[1]) if edge_index.size else 0,
        }

        payload = {
            "node_positions": pos,
            "node_features": atomic_numbers,
            "node_metadata": node_metadata,
            "edge_index": edge_index,
            "edge_weight": edge_weight,
            "graph_metadata": graph_metadata,
        }

        return payload

    def _persist_graph(
        self,
        structure_id: str,
        graph_payload: Dict[str, Any],
        *,
        dataset_name: Optional[str],
        entity_name: Optional[str] = None,
    ) -> str:
        """Persist a graph dictionary and register it in the entity registry."""

        graph_meta = graph_payload["graph_metadata"]
        if entity_name is None:
            name_parts = [
                structure_id,
                graph_meta["level"],
                f"{graph_meta['cutoff']:.1f}A",
            ]
            if dataset_name:
                name_parts.append(dataset_name)
            entity_name = self._sanitize_filename("__".join(map(str, name_parts)))
        else:
            entity_name = self._sanitize_filename(entity_name)

        file_path = self.graphs_dir / f"{entity_name}.pkl"
        with open(file_path, "wb") as handle:
            pickle.dump(graph_payload, handle)

        relative_path = str(file_path.relative_to(self.paths.data_root))

        metadata = {
            **graph_meta,
            "dependencies_available": self.dependencies_available,
        }
        if dataset_name:
            metadata["dataset"] = dataset_name

        self.entity_registry.register_entity(
            name=entity_name,
            format_type=self.processor_type,
            file_path=relative_path,
            metadata=metadata,
        )

        try:
            self.entity_registry.add_relationship(
                source_name=entity_name,
                target_name=structure_id,
                rel_type="derived_from",
                metadata={
                    "level": graph_meta["level"],
                    "cutoff": graph_meta["cutoff"],
                },
            )
        except ValueError:
            self.logger.warning(
                "Structure '%s' is not registered; relationship skipped", structure_id
            )

        return entity_name

    def _build_edges(
        self,
        positions: np.ndarray,
        *,
        cutoff: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute edges based on Euclidean distance cutoff."""

        if positions.shape[0] <= 1:
            return np.empty((2, 0), dtype=np.int64), np.empty((0,), dtype=np.float32)

        diff = positions[:, None, :] - positions[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)

        mask = (dist <= cutoff) & (dist > 0)
        src, dst = np.where(mask)
        edge_index = np.vstack([src, dst]).astype(np.int64)
        edge_weight = dist[src, dst].astype(np.float32)

        return edge_index, edge_weight

    def _atomic_number(self, element: Any) -> float:
        if isinstance(element, str):
            key = element.strip().upper()
            return float(_ELEMENT_TO_Z.get(key, 0))
        return 0.0
