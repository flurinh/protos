"""Sequence exporter implementation."""

from __future__ import annotations

from typing import Dict, List, Optional, Union

from protos.io.core.base_exporter import BaseExporter
from protos.io.formats.fasta_utils import write_fasta
from protos.io.paths import (
    get_sequence_entity_path,
    get_sequence_dataset_path,
    to_data_relative_path,
)


class SequenceExporter(BaseExporter):
    """Serialize sequences managed by SequenceProcessor into FASTA files."""

    def export_entity(
        self,
        name: str,
        *,
        export_name: Optional[str] = None,
        format: Optional[str] = None,
        overwrite: bool = False,
        sequence_ids: Optional[List[str]] = None,
    ) -> Dict[str, Union[str, int, bool, List[str]]]:
        export_format = (format or "fasta").lower()
        target_name = export_name or name
        if export_name and "." in export_name:
            stem, ext = export_name.rsplit(".", 1)
            if ext.lower() in {"fasta", "fa"}:
                target_name = stem
                export_format = ext.lower()

        if export_format not in {"fasta", "fa"}:
            raise ValueError("SequenceExporter supports only 'fasta' (or 'fa') outputs")
        target_path = get_sequence_entity_path(
            self.processor.paths,
            target_name,
            extension=export_format,
        )

        if target_path.exists() and not overwrite:
            raise FileExistsError(
                f"Entity export already exists at {target_path}. "
                "Pass overwrite=True to regenerate it."
            )

        data = self.processor.load_entity(name)
        if data is None:
            raise ValueError(f"Sequence '{name}' not found")

        if isinstance(data, str):
            sequences: Dict[str, str] = {name: data}
        else:
            sequences = dict(data)

        if sequence_ids is not None:
            filtered = {sid: seq for sid, seq in sequences.items() if sid in sequence_ids}
            if not filtered:
                raise ValueError("No sequences matched the requested sequence_ids")
            sequences = filtered

        target_path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(sequences, str(target_path))

        self.logger.info(
            "Exported %d sequence(s) for '%s' to %s",
            len(sequences),
            name,
            target_path,
        )

        return {
            "entity_name": name,
            "export_name": target_name,
            "format": export_format,
            "file_path": str(target_path),
            "artifact_path": to_data_relative_path(self.processor.paths, target_path),
            "sequence_count": len(sequences),
            "sequence_ids": list(sequences.keys()),
        }

    def export_dataset(
        self,
        dataset_name: str,
        *,
        format: Optional[str] = None,
        overwrite: bool = False,
        sequence_ids: Optional[List[str]] = None,
        export_name: Optional[str] = None,
        materialize_entities: bool = False,
    ) -> Dict[str, Union[str, int, bool, List[str]]]:
        export_format = (format or "fasta").lower()
        target_name = export_name or dataset_name
        if export_name and "." in export_name:
            stem, ext = export_name.rsplit(".", 1)
            if ext.lower() in {"fasta", "fa"}:
                target_name = stem
                export_format = ext.lower()

        if export_format not in {"fasta", "fa"}:
            raise ValueError("SequenceExporter supports only 'fasta' (or 'fa') outputs")

        sequences = self.processor.load_dataset(dataset_name)
        if sequence_ids is not None:
            sequences = {sid: sequences[sid] for sid in sequence_ids if sid in sequences}
            if not sequences:
                raise ValueError("No sequences matched the requested sequence_ids")

        target_path = get_sequence_dataset_path(
            self.processor.paths,
            target_name,
            extension=export_format,
        )

        if target_path.exists() and not overwrite:
            raise FileExistsError(
                f"Dataset export already exists at {target_path}. "
                "Use overwrite=True to regenerate it."
            )

        target_path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(sequences, str(target_path))

        if materialize_entities:
            for seq_id, sequence in sequences.items():
                self.processor.save_entity(seq_id, sequence)

        self.logger.info(
            "Exported dataset '%s' (%d sequences) to %s",
            dataset_name,
            len(sequences),
            target_path,
        )

        return {
            "dataset_name": dataset_name,
            "export_name": target_name,
            "format": export_format,
            "file_path": str(target_path),
            "artifact_path": to_data_relative_path(self.processor.paths, target_path),
            "sequence_count": len(sequences),
            "sequence_ids": list(sequences.keys()),
            "materialized_entities": materialize_entities,
        }
