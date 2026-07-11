# ModelManager

Import: `from protos.models.model_manager import ModelManager`

`ModelManager` resolves processor-managed artifacts and either prepares an
external command or invokes a registered Python runtime. It is orchestration,
not a scheduler: `prepare()` does not submit an external job and does not run
the returned command.

## Core objects

The dataclasses in `protos.models.model_specs` are:

| Type | Purpose |
| --- | --- |
| `ArtifactSpec` | required/optional input or output description |
| `ExecutionSpec` | execution mode, entrypoint, and environment metadata |
| `ModelCard` | model name/version plus execution and artifact specs |
| `ArtifactBundle` | a resolved path and metadata for one artifact |
| `PreparedJob` | external command, working directory, and artifacts |
| `RuntimeResult` | outputs from an in-process adapter |
| `ModelInvocation` | common wrapper containing either a job or runtime result |
| `ModelBatch` | normalized batch metadata; it is not an executed batch |

`ModelInvocation.is_external()` and `.is_runtime()` report which payload is
present.

## Registered cards

The current default card names are held in `manager.cards`. There is no
`ModelManager.list_models()` method.

```python
from protos.models.model_manager import ModelManager

manager = ModelManager()
assert "boltz2" in manager.cards
assert manager.cards["boltz2"].execution.mode == "external_config"

batch = manager.prepare_batch(
    "boltz2",
    [
        {
            "entity": "protein_1",
            "dataset_name": "prediction_inputs",
            "config": {"use_msa_server": True},
        }
    ],
    batch_name="demo_batch",
)
assert batch.to_dict()["input_count"] == 1
```

`prepare_batch()` only validates and normalizes these dictionaries. It does not
check that the entities/datasets exist and does not call `prepare()` for them.

## Preparing an external Boltz job

```python
from protos.models.model_manager import ModelManager
from protos.processing.sequence import SequenceProcessor

sequences = SequenceProcessor()
sequences.save_entity("demo", "ACDEFGHIK")
sequences.create_dataset("demo_sequences", ["demo"])

manager = ModelManager()
invocation = manager.prepare(
    "boltz2",
    inputs={"sequence_dataset": "demo_sequences", "entity": "demo"},
)

assert invocation.is_external()
assert invocation.job is not None
assert invocation.job.command[:2] == ["boltz", "predict"]
assert invocation.job.working_dir.is_dir()
```

This writes the Boltz YAML, FASTA, and run metadata beneath
`<data_root>/models/boltz2/`. It does not prove that the `boltz` executable,
weights, GPU, or container runtime is available. Execution and output ingestion
are separate deployment-specific steps.

## Default model-card status

The manager registers cards/adapters for `boltz2`, `boltzgen`, three embedding
cards, `lambda`, `ligand_mpnn`, `unidock`, `equibind`, `pocketdta`,
`graphscoredta`, and `pocket2mol`. Registration means the adapter can assemble
or attempt the declared invocation; it is not a claim that every external model
is installed or validated on the current machine.

Runtime embedding cards load model weights when prepared. Lambda depends on the
pinned Lambda implementation and its model data. External cards require their
own executables, repositories, containers, checkpoints, and input formats.

Custom cards can be installed with `register_model(card, adapter)`. Passing
`None` selects the configurable adapter based on `card.execution.mode`.
