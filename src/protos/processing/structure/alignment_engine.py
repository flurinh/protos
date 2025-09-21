"""Structure alignment engine that wraps available alignment algorithms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from protos.analysis.structure import alignment as alignment_utils


@dataclass
class AlignmentResult:
    """Container for an individual alignment outcome."""

    structure_id: str
    aligned_id: str
    rotation: np.ndarray
    translation: np.ndarray
    rmsd: float
    algorithm: str
    alignment_path: Optional[Tuple] = None
    error: Optional[str] = None

    @property
    def success(self) -> bool:
        return self.error is None


class StructureAlignmentEngine:
    """High-level facade for structure alignment operations."""

    def __init__(self, processor: "StructureProcessor") -> None:
        # Import locally to avoid circular dependency at module import time
        from protos.processing.structure.structure_processor import StructureProcessor

        if not isinstance(processor, StructureProcessor):  # pragma: no cover - safety
            raise TypeError("StructureAlignmentEngine requires a StructureProcessor instance")

        self._processor = processor
        self._logger = processor.logger

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def align(
        self,
        structure_ids: Iterable[str],
        reference_id: str,
        *,
        method: str = "cealign",
        atom_selection: str = "CA",
        chain_id: Optional[str] = None,
        apply_transform: bool = True,
        cealign_window: int = 8,
        cealign_max_gap: int = 30,
    ) -> Dict[str, AlignmentResult]:
        """Align a collection of structures to a common reference."""

        ordered_ids: List[str] = [reference_id]
        ordered_ids.extend(sid for sid in structure_ids if sid != reference_id)

        reference_df = self._load_structure(reference_id)
        if reference_df is None:
            raise ValueError(f"Reference structure {reference_id} not found")

        algorithm = method.lower()
        results: Dict[str, AlignmentResult] = {}

        for structure_id in ordered_ids:
            if structure_id == reference_id:
                results[structure_id] = AlignmentResult(
                    structure_id=reference_id,
                    aligned_id=reference_id,
                    rotation=np.eye(3),
                    translation=np.zeros(3),
                    rmsd=0.0,
                    algorithm=algorithm,
                )
                continue

            mobile_df = self._load_structure(structure_id)
            if mobile_df is None:
                self._logger.warning("Structure %s not found, skipping alignment", structure_id)
                continue

            try:
                result = self._run_alignment(
                    reference_df,
                    mobile_df,
                    structure_id=structure_id,
                    algorithm=algorithm,
                    atom_selection=atom_selection,
                    chain_id=chain_id,
                    apply_transform=apply_transform,
                    cealign_window=cealign_window,
                    cealign_max_gap=cealign_max_gap,
                    target_id=structure_id,
                )
            except Exception as exc:  # noqa: BLE001 - propagate as informational result
                self._logger.error("Failed to align %s: %s", structure_id, exc)
                result = AlignmentResult(
                    structure_id=structure_id,
                    aligned_id=structure_id,
                    rotation=np.eye(3),
                    translation=np.zeros(3),
                    rmsd=float("nan"),
                    algorithm=algorithm,
                    error=str(exc),
                )

            results[structure_id] = result

        return results

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _load_structure(self, structure_id: str) -> Optional[pd.DataFrame]:
        if structure_id in self._processor.frames:
            return self._processor.frames[structure_id]
        return self._processor.load_entity(structure_id)

    def _run_alignment(
        self,
        reference_df: pd.DataFrame,
        mobile_df: pd.DataFrame,
        *,
        structure_id: str,
        algorithm: str,
        atom_selection: str,
        chain_id: Optional[str],
        apply_transform: bool,
        cealign_window: int,
        cealign_max_gap: int,
        target_id: str,
    ) -> AlignmentResult:
        alg = algorithm
        rotation: Optional[np.ndarray] = None
        translation: Optional[np.ndarray] = None
        rmsd: Optional[float] = None
        alignment_path: Optional[Tuple] = None

        if alg not in {"cealign", "kabsch", "simple"}:
            raise ValueError("method must be one of {'cealign', 'kabsch', 'simple'}")

        if alg == "cealign":
            try:
                ref_coords = self._select_alignment_coordinates(
                    reference_df, atom_selection=atom_selection, chain_id=chain_id
                )
                mob_coords = self._select_alignment_coordinates(
                    mobile_df, atom_selection=atom_selection, chain_id=chain_id
                )
                (
                    _,
                    rotation,
                    translation,
                    alignment_path,
                    rmsd,
                ) = alignment_utils.align_structures(
                    ref_coords,
                    mob_coords,
                    window_size=cealign_window,
                    max_gap=cealign_max_gap,
                )
            except ImportError as exc:
                self._logger.warning(
                    "CEalign unavailable (%s). Falling back to simple alignment.", exc
                )
                alg = "simple"
            else:
                rotation = np.asarray(rotation, dtype=float)
                translation = np.asarray(translation, dtype=float).reshape(-1)
                if translation.size != 3:
                    raise ValueError("CEalign translation vector must have length 3")

                aligned_id = structure_id
                if apply_transform:
                    aligned_id = self._store_transformed_frame(
                        structure_id,
                        rotation=rotation,
                        translation=translation,
                        save_as=target_id,
                    )

                return AlignmentResult(
                    structure_id=structure_id,
                    aligned_id=aligned_id,
                    rotation=rotation,
                    translation=translation,
                    rmsd=rmsd if rmsd is not None else float("nan"),
                    algorithm=alg,
                    alignment_path=alignment_path,
                )

        if alg in {"kabsch", "simple"}:
            if alg == "simple":
                rotation, translation, rmsd = alignment_utils.simple_align_structures(
                    reference_df,
                    mobile_df,
                    atom_selection=atom_selection,
                    chain_id=chain_id,
                )
            else:  # 'kabsch'
                ref_coords = self._select_alignment_coordinates(
                    reference_df, atom_selection=atom_selection, chain_id=chain_id
                )
                mob_coords = self._select_alignment_coordinates(
                    mobile_df, atom_selection=atom_selection, chain_id=chain_id
                )
                if len(ref_coords) != len(mob_coords):
                    raise ValueError(
                        "Kabsch alignment requires equal numbers of selected coordinates"
                    )
                rotation, translation, rmsd = alignment_utils.kabsch_alignment(
                    ref_coords.to_numpy(), mob_coords.to_numpy()
                )

            rotation = np.asarray(rotation, dtype=float).T
            translation = np.asarray(translation, dtype=float).reshape(-1)
            if translation.size != 3:
                raise ValueError("Translation vector must have length 3")

        if rotation is None or translation is None:
            raise RuntimeError("Alignment did not yield transformation parameters")

        aligned_id = structure_id
        if apply_transform:
            aligned_id = self._store_transformed_frame(
                structure_id,
                rotation=rotation,
                translation=translation,
                save_as=target_id,
            )

        return AlignmentResult(
            structure_id=structure_id,
            aligned_id=aligned_id,
            rotation=rotation,
            translation=translation,
            rmsd=rmsd if rmsd is not None else float("nan"),
            algorithm=alg,
            alignment_path=alignment_path,
        )

    def _select_alignment_coordinates(
        self,
        df: pd.DataFrame,
        *,
        atom_selection: str,
        chain_id: Optional[str],
    ) -> pd.DataFrame:
        df_reset = df.reset_index() if isinstance(df.index, pd.MultiIndex) else df.copy()

        if chain_id is not None:
            df_reset = df_reset[df_reset['auth_chain_id'] == chain_id]

        if atom_selection == 'CA':
            mask = df_reset['atom_name'] == 'CA'
        elif atom_selection == 'backbone':
            mask = df_reset['atom_name'].isin(['N', 'CA', 'C', 'O'])
        elif atom_selection == 'all':
            if 'group' in df_reset.columns:
                mask = df_reset['group'] == 'ATOM'
            else:
                mask = pd.Series(True, index=df_reset.index)
        else:
            raise ValueError(f"Unknown atom selection: {atom_selection}")

        coords = df_reset.loc[mask, ['x', 'y', 'z']].dropna()
        if coords.empty:
            raise ValueError(
                "No coordinates found for alignment using selection "
                f"'{atom_selection}' and chain '{chain_id or 'all'}'"
            )
        return coords.reset_index(drop=True)

    def _store_transformed_frame(
        self,
        structure_id: str,
        rotation: Optional[np.ndarray],
        translation: Optional[np.ndarray],
        *,
        save_as: Optional[str],
        source_df: Optional[pd.DataFrame] = None,
    ) -> str:
        if source_df is None:
            source_df = self._load_structure(structure_id)
        if source_df is None:
            raise ValueError(f"Structure {structure_id} not found")

        transformed = source_df.copy()
        coords = transformed[['x', 'y', 'z']].to_numpy(dtype=float, copy=True)

        if rotation is not None:
            if rotation.shape != (3, 3):
                raise ValueError("Rotation matrix must be 3x3")
            coords = coords @ rotation

        if translation is not None:
            translation = np.asarray(translation, dtype=float)
            if translation.shape != (3,):
                raise ValueError("Translation vector must have length 3")
            coords = coords + translation

        transformed.loc[:, 'x'] = coords[:, 0]
        transformed.loc[:, 'y'] = coords[:, 1]
        transformed.loc[:, 'z'] = coords[:, 2]

        target_id = save_as or structure_id
        transformed_reset = transformed.reset_index()
        transformed_reset['structure_id'] = target_id
        canonical = self._processor._ensure_canonical(transformed_reset, target_id)  # noqa: SLF001
        self._processor._set_frame(target_id, canonical)  # noqa: SLF001
        return target_id


__all__ = ["StructureAlignmentEngine", "AlignmentResult"]
