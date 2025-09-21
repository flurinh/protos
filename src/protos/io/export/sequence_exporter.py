"""Sequence exporter implementation."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional

import shutil

from protos.io.core.base_exporter import BaseExporter
from protos.io.formats.fasta_utils import read_fasta, write_fasta


class SequenceExporter(BaseExporter):
    """Serialize sequences managed by SequenceProcessor into FASTA files."""

    def export_entity(
        self,
        name: str,
        out_path: Path,
        format: Optional[str] = None,
        overwrite: bool = False,
        sequence_ids: Optional[list[str]] = None,
    ) -> Path:
        export_format = (format or "fasta").lower()
        if export_format not in {"fasta", "fa"}:
            raise ValueError("SequenceExporter supports only 'fasta' (or 'fa') outputs")

        out_path = Path(out_path)
        if out_path.exists() and not overwrite:
            raise FileExistsError(
                f"File {out_path} already exists. Use overwrite=True to replace it."
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

        out_path.parent.mkdir(parents=True, exist_ok=True)
        write_fasta(sequences, str(out_path))
        self.logger.info("Exported %d sequence(s) for '%s' to %s", len(sequences), name, out_path)
        return out_path

    def export_dataset(
        self,
        dataset_name: str,
        output_dir: Path,
        format: Optional[str] = None,
        overwrite: bool = False,
        name_pattern: Optional[str] = None,
        sequence_ids: Optional[list[str]] = None,
    ) -> Dict[str, Path]:
        info = self.processor.get_dataset_info(dataset_name)
        metadata = info.get('metadata', {}) if info else {}
        artifact_rel = metadata.get('artifact_path')

        if artifact_rel:
            src_path = Path(self.processor.paths.data_root) / artifact_rel
            if not src_path.exists():
                raise FileNotFoundError(f"Dataset artifact missing at {src_path}")

            export_format = (format or 'fasta').lower()
            filename = (
                name_pattern.format(name=dataset_name)
                if name_pattern
                else f"{dataset_name}.{export_format}"
            )
            target = Path(output_dir) / filename
            if target.exists() and not overwrite:
                raise FileExistsError(
                    f"File {target} already exists. Use overwrite=True to replace it."
                )
            target.parent.mkdir(parents=True, exist_ok=True)
            sequences = read_fasta(str(src_path))
            if sequence_ids is not None:
                sequences = {sid: seq for sid, seq in sequences.items() if sid in sequence_ids}
                if not sequences:
                    raise ValueError("No sequences matched the requested sequence_ids")
            write_fasta(sequences, str(target))
            self.logger.info(
                "Exported dataset '%s' to %s", dataset_name, target
            )
            return {dataset_name: target}

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        entity_names = self.dataset_manager.get_dataset_entities(dataset_name)
        if sequence_ids is not None:
            entity_names = [name for name in entity_names if name in sequence_ids]
            if not entity_names:
                raise ValueError("No entities matched the requested sequence_ids")

        exported_paths: Dict[str, Path] = {}
        for entity_name in entity_names:
            filename: str
            export_format = (format or 'fasta').lower()
            if name_pattern:
                filename = name_pattern.format(name=entity_name)
            else:
                filename = f"{entity_name}.{export_format}"

            target_path = output_dir / filename
            exported_paths[entity_name] = self.export_entity(
                entity_name,
                target_path,
                format=export_format,
                overwrite=overwrite,
            )

        self.logger.info(
            "Exported %d sequences from dataset '%s'",
            len(exported_paths),
            dataset_name,
        )
        return exported_paths
