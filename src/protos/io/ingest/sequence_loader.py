"""Loader for sequence data (FASTA) focusing on local ingestion."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any, Dict, Optional, List

import io
import requests

from protos.io.core.base_loader import BaseLoader
from protos.io.formats.fasta_utils import read_fasta, write_fasta
from protos.processing.sequence import SequenceProcessor


class SequenceLoader(BaseLoader):
    loader_type = "sequence"

    def __init__(self, name: str = "sequence_loader", *, processor: Optional[SequenceProcessor] = None) -> None:
        super().__init__(name=name)
        self._processor: Optional[SequenceProcessor] = processor

    def parse_identifier(self, identifier: str) -> Dict[str, Any]:
        path = Path(identifier)
        if path.exists():
            return {
                "source": "local",
                "path": path,
                "name": path.stem,
            }

        token = identifier.strip()
        if token.lower().startswith("uniprot:"):
            ids = [part.strip() for part in token.split(":", 1)[1].split(",") if part.strip()]
            if not ids:
                raise ValueError("No UniProt identifiers provided")
            return {"source": "uniprot", "ids": ids, "name": ids[0]}

        # Bare UniProt accession (e.g., P00533)
        if self._looks_like_uniprot_id(token):
            return {"source": "uniprot", "ids": [token], "name": token}

        raise ValueError(f"Unsupported sequence identifier: {identifier}")

    def fetch_entity(
        self,
        identifier: str,
        source: Optional[str] = None,
        **kwargs,
    ) -> Optional[Path]:
        info = self.parse_identifier(identifier)
        source = info["source"]
        if source == "local":
            return info["path"]
        if source == "uniprot":
            return self._fetch_from_uniprot(info["ids"], info.get("name"))
        raise ValueError(f"Unsupported source for sequence loader: {source}")

    def register_sequence_records(
        self,
        records: List[Dict[str, Any]],
        *,
        dataset_name: Optional[str] = None,
        dataset_metadata: Optional[Dict[str, Any]] = None,
        overwrite: bool = False,
    ) -> Dict[str, Any]:
        """Register sequences provided in-memory, returning registry results."""

        if not records:
            return {"entities": [], "dataset": None}

        processor = self._get_processor()
        registered_entities: List[str] = []

        for record in records:
            name = record.get("name")
            sequence = record.get("sequence")
            metadata = record.get("metadata") or {}

            if not name or sequence is None:
                raise ValueError("Sequence record requires 'name' and 'sequence'")

            entity_exists = processor.entity_exists(name)

            if entity_exists and not overwrite:
                if metadata:
                    processor.entity_registry.update_metadata(
                        name,
                        processor.processor_type,
                        metadata,
                    )
            else:
                processor.save_entity(name, sequence, metadata=metadata)

            registered_entities.append(name)

        dataset_created = None
        if dataset_name:
            ds_metadata = dataset_metadata.copy() if dataset_metadata else {}
            ds_metadata.setdefault("size", len(registered_entities))
            processor.create_dataset(dataset_name, registered_entities, ds_metadata)
            dataset_created = dataset_name

        return {"entities": registered_entities, "dataset": dataset_created}

    def download_and_register(
        self,
        identifier: str,
        name: Optional[str] = None,
        source: str = "local",
        materialize_entities: bool = True,
        metadata: Optional[Dict[str, Any]] = None,
        **kwargs,
    ) -> Optional[str]:
        info = self.parse_identifier(identifier)
        source_path = self.fetch_entity(identifier, source=info["source"])
        if source_path is None:
            return None

        sequences = read_fasta(str(source_path))
        if not sequences:
            return None

        processor = self._get_processor()

        if len(sequences) == 1 and materialize_entities:
            seq_id, sequence = next(iter(sequences.items()))
            entity_name = name or seq_id
            processor.save_entity(entity_name, sequence, metadata=metadata)
            saved_name = entity_name
        else:
            dataset_name = name or source_path.stem
            processor.save_sequences(
                sequences,
                output_file=dataset_name,
                dataset_name=dataset_name,
                metadata=metadata,
                materialize_entities=materialize_entities,
            )
            saved_name = dataset_name

        if info["source"] == "local" and source_path.exists():
            try:
                source_path.unlink()
            except OSError:
                pass

        return saved_name

    def _get_processor(self) -> SequenceProcessor:
        if self._processor is None:
            self._processor = SequenceProcessor(name=f"{self.name}_processor")
        return self._processor

    @staticmethod
    def _looks_like_uniprot_id(identifier: str) -> bool:
        return bool(identifier) and len(identifier) >= 6 and identifier[0].isalpha()

    def _fetch_from_uniprot(self, ids: List[str], target_name: Optional[str]) -> Path:
        processor = self._get_processor()
        safe_name = processor._sanitize_filename(target_name or ids[0])
        output_path = Path(processor.path_fasta_dir) / f"{safe_name}.fasta"

        sequences: Dict[str, str] = {}
        for uniprot_id in ids:
            seq_dict = self._download_uniprot_fasta(uniprot_id)
            sequences.update(seq_dict)

        write_fasta(sequences, str(output_path))
        return output_path

    @staticmethod
    def _download_uniprot_fasta(uniprot_id: str) -> Dict[str, str]:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
        except requests.RequestException as exc:
            raise RuntimeError(f"Failed to fetch UniProt entry '{uniprot_id}'") from exc
        return SequenceLoader._parse_fasta_text(response.text)

    @staticmethod
    def _parse_fasta_text(text: str) -> Dict[str, str]:
        sequences: Dict[str, str] = {}
        current_id: Optional[str] = None
        buffer: List[str] = []

        for line in io.StringIO(text):
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(buffer)
                header = stripped[1:].split()[0]
                current_id = header
                buffer = []
            else:
                buffer.append(stripped)

        if current_id is not None:
            sequences[current_id] = ''.join(buffer)

        if not sequences:
            raise ValueError("No sequences found in UniProt response")
        return sequences
