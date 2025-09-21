"""Base exporter class for the Protos framework.

Exporters encapsulate filesystem-aware write operations so processors can remain
focused on in-memory data manipulation. They coordinate with processors to
retrieve canonical data, but own all path validation, overwrite protection, and
batch export mechanics.
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Optional

from protos.io.paths import get_protos_paths


class BaseExporter(ABC):
    """Reusable export helper that wraps a processor instance."""

    def __init__(self, processor, name: str = "exporter") -> None:
        self.processor = processor
        self.name = name
        self.paths = get_protos_paths()
        self.dataset_manager = processor.dataset_manager
        self.processor_type = processor.processor_type
        self.logger = self._setup_logger()

    def _setup_logger(self) -> logging.Logger:
        logger = logging.getLogger(f"{self.__class__.__name__}.{self.name}")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger

    @abstractmethod
    def export_entity(
        self,
        name: str,
        out_path: Path,
        format: Optional[str] = None,
        overwrite: bool = False,
    ) -> Path:
        """Write a single entity to disk. Subclasses must implement."""

    def export_dataset(
        self,
        dataset_name: str,
        output_dir: Path,
        format: Optional[str] = None,
        overwrite: bool = False,
        name_pattern: Optional[str] = None,
    ) -> Dict[str, Path]:
        """Export every entity in a dataset using shared path logic."""

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        entity_names = self.dataset_manager.get_dataset_entities(dataset_name)
        exported_paths: Dict[str, Path] = {}

        for entity_name in entity_names:
            try:
                if name_pattern:
                    filename = name_pattern.format(name=entity_name)
                else:
                    if format:
                        filename = f"{entity_name}.{format}"
                    else:
                        filename = entity_name

                target_path = output_dir / filename
                exported_paths[entity_name] = self.export_entity(
                    entity_name,
                    target_path,
                    format=format,
                    overwrite=overwrite,
                )
            except Exception as exc:  # pragma: no cover - user-facing logging
                self.logger.error(
                    "Failed to export entity '%s' from dataset '%s': %s",
                    entity_name,
                    dataset_name,
                    exc,
                )

        self.logger.info(
            "Exported %d/%d entities from dataset '%s'",
            len(exported_paths),
            len(entity_names),
            dataset_name,
        )
        return exported_paths
