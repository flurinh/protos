"""Processor for small-molecule descriptors (SMILES, InChI, metadata)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional, Union

from protos.io.core.base_processor import BaseProcessor
from protos.io.paths import sanitize_storage_name
from protos.io.ingest.utils.ligand_utils import (
    calculate_molecular_properties,
    is_drug_like,
    validate_smiles,
)


class MoleculeProcessor(BaseProcessor):
    """Manage ligand descriptors (SMILES, InChI) and related metadata."""

    processor_type = "molecule"

    def __init__(self, name: str = "molecule_processor"):
        super().__init__(name=name)
        self.records_dir = Path(self.get_subdirectory_path("records_dir"))
        self.datasets_dir = Path(self.get_subdirectory_path("datasets_dir"))
        self.records_dir.mkdir(parents=True, exist_ok=True)
        self.datasets_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # BaseProcessor API
    # ------------------------------------------------------------------
    def load_entity(self, name: str) -> Optional[Dict[str, Any]]:
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info is None:
            return None

        record: Dict[str, Any] = {}
        if entity_info.file_path:
            record_path = Path(self.paths.data_root) / entity_info.file_path
            if record_path.exists():
                with open(record_path, "r", encoding="utf-8") as handle:
                    record = json.load(handle)

        record.setdefault("name", entity_info.original_id)
        record.setdefault("metadata", entity_info.metadata or {})
        return record

    def save_entity(
        self,
        name: str,
        data: Union[str, Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        if isinstance(data, str):
            record: Dict[str, Any] = {"smiles": data, "kind": "smiles"}
        else:
            record = dict(data)

        safe_name = sanitize_storage_name(name)
        record.setdefault("name", safe_name)
        record.setdefault("kind", "smiles")

        record_path = self.records_dir / f"{safe_name}.json"
        with open(record_path, "w", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2, ensure_ascii=False)

        rel_path = record_path.relative_to(self.paths.data_root)
        entity_metadata = {"kind": record.get("kind", "smiles")}
        if metadata:
            entity_metadata.update(metadata)

        self.entity_registry.register_entity(
            name=safe_name,
            format_type=self.processor_type,
            file_path=str(rel_path),
            metadata=entity_metadata,
        )
        return safe_name

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------
    def register_smiles_map(
        self,
        smiles_map: Dict[str, str],
        *,
        dataset_name: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> str:
        entity_names = []
        for label, smiles in smiles_map.items():
            entity_name = self.save_entity(label, {"smiles": smiles, "kind": "smiles"}, metadata)
            entity_names.append(entity_name)

        dataset_id = dataset_name or self._sanitize_filename("smiles_dataset")
        dataset_meta = {
            "kind": "smiles",
            "entity_count": len(entity_names),
        }
        if metadata:
            dataset_meta.update(metadata)

        if self.dataset_manager.dataset_exists(dataset_id):
            self.dataset_manager.delete_dataset(dataset_id)

        self.create_dataset(dataset_id, entity_names, dataset_meta)
        return dataset_id

    def list_ligands(self) -> list[str]:
        return self.entity_registry.list_entities(self.processor_type)

    # ------------------------------------------------------------------
    # Analysis helpers (compatibility surface for MCP tools)
    # ------------------------------------------------------------------
    def calculate_properties(self, smiles: str) -> Optional[Dict[str, Any]]:
        return calculate_molecular_properties(smiles)

    def filter_drug_like(self, entity_names: list[str], strict: bool = False) -> list[str]:
        drug_like_entities: list[str] = []
        for name in entity_names:
            smiles = self._resolve_smiles(name)
            if smiles and is_drug_like(smiles, strict=strict):
                drug_like_entities.append(name)
        return drug_like_entities

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _resolve_smiles(self, name_or_smiles: str) -> Optional[str]:
        record = self.load_entity(name_or_smiles)
        if record and isinstance(record, dict):
            payload = record.get('metadata', {})
            if 'smiles' in record:
                return record['smiles']
            if 'smiles' in payload:
                return payload['smiles']
        # Fallback: treat incoming value as SMILES if no entity exists
        valid, canonical = validate_smiles(name_or_smiles)
        if valid:
            return canonical or name_or_smiles
        return None
