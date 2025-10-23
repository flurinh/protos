"""A tiny runtime entrypoint for testing ConfigurableRuntimeAdapter.

Writes a small property table annotating the first sequence in a dataset and
returns metadata pointing to the table for ingestion tests.
"""

from __future__ import annotations

from typing import Any, Dict, List

import pandas as pd

from protos.models.model_specs import ArtifactBundle, RuntimeResult
from protos.processing.property import PropertyProcessor


def run(
    *,
    card,  # ModelCard
    request,  # ModelRequest
    inputs: List[ArtifactBundle],
    paths,  # ProtosPaths
    work_dir: str | None = None,
    outputs_dir: str | None = None,
    inputs_dir: str | None = None,
) -> Dict[str, Any]:
    # Find sequence dataset bundle and pick first sequence id
    seq_bundle = next((b for b in inputs if b.spec.name == "sequence_dataset"), None)
    seqs = (seq_bundle.metadata.get("sequences") if seq_bundle else {}) or {}
    first_id = next(iter(seqs.keys()), "TESTSEQ")

    rows = [
        {
            "entity_name": first_id,
            "scope": [{"format": "sequence", "name": first_id}],
            "score": 0.42,
            "model": card.name,
        }
    ]

    table_name = "dummy_runtime_properties"
    prop = PropertyProcessor()
    df = pd.DataFrame(rows)
    prop.create_property_table(table_name, df, metadata={"source": card.name}, allow_create=True)

    # Return metadata so ModelManager.ingest_outputs can confirm registration
    return {
        "outputs": {"property_table": table_name},
        "artifacts": [],
        "metadata": {
            "property_table": table_name,
            "work_dir": work_dir,
            "outputs_dir": outputs_dir,
            "inputs_dir": inputs_dir,
        },
    }

