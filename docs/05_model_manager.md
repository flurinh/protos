# Model Manager

The Model Manager orchestrates model invocations in Protos. It provides a declarative system for specifying model requirements (via Model Cards), assembling input data from processors, and preparing execution packages for both in-process and external model runs.

## Architecture Overview

```
ModelCard (what a model needs)
    ↓
ModelManager.prepare()
    ↓
assemble_inputs() → ArtifactBundles from processors
    ↓
Adapter._prepare() → ModelInvocation
    ↓
PreparedJob (external) or RuntimeResult (in-process)
```

---

## Model Cards

A **ModelCard** is a declarative specification of what a model needs and how it runs.

**Location:** `protos.models.model_specs`

### ModelCard Structure

```python
@dataclass
class ModelCard:
    name: str                           # Unique identifier (e.g., "boltz2")
    version: str                        # Version string
    description: str                    # Human-readable description
    execution: ExecutionSpec            # How to execute the model
    input_spec: List[ArtifactSpec]      # Required/optional inputs
    output_spec: List[ArtifactSpec]     # Expected outputs
    metadata: Dict[str, Any]            # Additional configuration
```

### Example Model Card

```python
from protos.models.model_specs import ModelCard, ExecutionSpec, ArtifactSpec

boltz_card = ModelCard(
    name="boltz2",
    version="2",
    description="Boltz structure prediction",
    execution=ExecutionSpec(
        mode="external_config",
        entrypoint="boltz predict"
    ),
    input_spec=[
        ArtifactSpec(
            name="sequence_dataset",
            kind="sequence",
            provider="sequence_dataset",
            format="fasta"
        )
    ],
    output_spec=[
        ArtifactSpec(
            name="predicted_structure",
            kind="structure",
            provider="boltz",
            format="cif"
        )
    ]
)
```

---

## Artifact Specifications

**ArtifactSpec** describes an input or output artifact:

```python
@dataclass
class ArtifactSpec:
    name: str                    # Identifier for this artifact
    kind: str                    # Type: "sequence", "structure", "embedding", etc.
    provider: str                # How to obtain: "sequence_dataset", "grn_table", etc.
    format: Optional[str]        # File format: "fasta", "cif", "csv", "npz"
    optional: bool = False       # Whether the artifact is required
    params: Dict[str, Any]       # Additional parameters
```

### Built-in Artifact Providers

| Provider | Description |
|----------|-------------|
| `sequence_dataset` | Load sequences from SequenceProcessor |
| `structure_entity` | Load structure from StructureProcessor |
| `grn_table` | Load GRN table from GRNProcessor |
| `embedding_dataset` | Load embeddings from EmbeddingProcessor |
| `graph_entity` | Load graph from GraphProcessor |
| `ligand_file` | Load ligand file |
| `file_path` | Generic file path resolution |
| `json_payload` | JSON configuration object |

---

## Execution Specifications

**ExecutionSpec** describes how a model runs:

```python
@dataclass
class ExecutionSpec:
    mode: str                           # "external_config" or "runtime"
    entrypoint: Optional[str]           # Command or callable path
    environment: Dict[str, Any]         # Environment requirements
    expected_config: Optional[str]      # Config format (yaml/json)
```

### Execution Modes

| Mode | Description |
|------|-------------|
| `external_config` | Runs as external process (e.g., `boltz predict`) |
| `runtime` | Runs in-process via Python callable |

---

## Data Packaging

The ModelManager assembles artifacts from processors based on the ModelCard's input specifications.

### Assembly Flow

```python
def assemble_inputs(card: ModelCard, request: ModelRequest) -> List[ArtifactBundle]:
    for spec in card.input_spec:
        provider = artifact_providers[spec.provider]
        bundle = provider(spec, request)
        bundles.append(bundle)
    return bundles
```

### ArtifactBundle

Resolved artifact ready for use:

```python
@dataclass
class ArtifactBundle:
    spec: ArtifactSpec          # Original specification
    path: Path                  # Path to artifact file
    metadata: Dict[str, Any]    # Additional metadata
```

---

## Model Manager Usage

### Basic Usage

```python
from protos.models.model_manager import ModelManager

manager = ModelManager()

# List available models
models = manager.list_models()
for name in models:
    card = manager.cards[name]
    print(f"{name}: {card.description}")

# Prepare a model invocation
invocation = manager.prepare(
    "boltz2",
    inputs={"sequence_dataset": "my_proteins", "entity": "protein_1"},
    config={"recycling_steps": 10},
    metadata={"experiment": "test_run"}
)
```

### ModelRequest

Input to the prepare method:

```python
class ModelRequest:
    inputs: Dict[str, Any]    # Named inputs matching ArtifactSpec
    config: Dict[str, Any]    # Model configuration
    metadata: Dict[str, Any]  # Additional metadata
```

### ModelInvocation

Output from prepare:

```python
@dataclass
class ModelInvocation:
    model: str                              # Model name
    card: ModelCard                         # Full model card
    job: Optional[PreparedJob]              # For external execution
    runtime: Optional[RuntimeResult]        # For in-process execution
    inputs: List[ArtifactBundle]            # Resolved inputs
    outputs: List[ArtifactBundle]           # Resolved outputs
    metadata: Dict[str, Any]                # Combined metadata

    def is_external(self) -> bool:
        return self.job is not None

    def is_runtime(self) -> bool:
        return self.runtime is not None
```

---

## Job Submission

### PreparedJob (External Execution)

For models that run as external processes:

```python
@dataclass
class PreparedJob:
    command: List[str]              # Command to execute
    working_dir: Path               # Working directory
    artifacts: List[ArtifactBundle] # Input artifacts
    metadata: Dict[str, Any]        # Job metadata
```

### Example: External Job

```python
invocation = manager.prepare("boltz2", inputs={
    "sequence_dataset": "my_proteins",
    "entity": "protein_1"
})

if invocation.is_external():
    job = invocation.job
    print(f"Command: {' '.join(job.command)}")
    print(f"Working dir: {job.working_dir}")

    # Execute externally
    import subprocess
    subprocess.run(job.command, cwd=job.working_dir)
```

### RuntimeResult (In-Process Execution)

For models that run within the Python process:

```python
@dataclass
class RuntimeResult:
    outputs: Dict[str, Any]         # Output data
    artifacts: List[ArtifactBundle] # Output artifacts
    metadata: Dict[str, Any]        # Execution metadata
```

### Example: Runtime Execution

```python
invocation = manager.prepare("embedding_esm2_t12_35m", inputs={
    "sequence_dataset": "my_proteins"
})

if invocation.is_runtime():
    result = invocation.runtime
    embeddings = result.outputs.get("embeddings")
    print(f"Generated {len(embeddings)} embeddings")
```

---

## Built-in Models

The ModelManager registers several default models:

### Boltz Structure Prediction

```python
invocation = manager.prepare(
    "boltz2",
    inputs={
        "sequence_dataset": "my_proteins",
        "entity": "protein_1"
    },
    config={
        "mutations": [{"position": 48, "mutant": "R", "original": "K"}],
        "ligand_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "use_msa_server": True,
        "recycling_steps": 10,
        "diffusion_samples": 5
    }
)

# Execute the prepared job
job = invocation.job
print(f"Run: {' '.join(job.command)}")
```

### BoltzGen Protein Design

```python
invocation = manager.prepare(
    "boltzgen",
    config={
        "job_name": "design_binder",
        "designed_proteins": [{
            "id": "B",
            "min_length": 80,
            "max_length": 140
        }],
        "target_structure": "target.cif",
        "ligands": [{"id": "LIG", "ccd": "ATP"}]
    }
)
```

### Embedding Generation

```python
invocation = manager.prepare(
    "embedding_esm2_t12_35m",
    inputs={"sequence_dataset": "my_proteins"}
)

# Outputs saved to data/models/embedding_esm2_t12_35m/run_*/outputs/
output_dir = invocation.runtime.metadata["outputs_dir"]
```

### Lambda Property Prediction

```python
invocation = manager.prepare(
    "lambda",
    inputs={
        "sequence_dataset": "gpcrdb_sequences",
        "grn_table": "class_a_gpcr",
        "embedding_dataset": "gpcrdb_embeddings"  # optional
    },
    config={
        "run_id": "my_run",
        "binding_config": "binding_domain.json"
    }
)
```

---

## Registering Custom Models

### Define Model Card

```python
my_card = ModelCard(
    name="my_model",
    version="1.0",
    description="My custom model",
    execution=ExecutionSpec(
        mode="runtime",
        entrypoint="my_package.models:run_prediction"
    ),
    input_spec=[
        ArtifactSpec(
            name="structure",
            kind="structure",
            provider="structure_entity",
            format="pkl"
        )
    ],
    output_spec=[
        ArtifactSpec(
            name="predictions",
            kind="property",
            provider="property_processor",
            format="csv"
        )
    ]
)
```

### Implement Runtime Function

```python
# my_package/models.py
from protos.models.model_specs import RuntimeResult, ArtifactBundle

def run_prediction(card, request, inputs, outputs_dir):
    """Runtime entrypoint for my_model."""

    # Get input artifact
    structure_bundle = next(b for b in inputs if b.spec.name == "structure")
    structure_path = structure_bundle.path

    # Run prediction
    predictions = my_prediction_function(structure_path)

    # Save outputs
    output_path = Path(outputs_dir) / "predictions.csv"
    predictions.to_csv(output_path)

    return RuntimeResult(
        outputs={"predictions": predictions},
        artifacts=[],
        metadata={"output_path": str(output_path)}
    )
```

### Register Model

```python
manager = ModelManager()
manager.register_model(my_card)  # Uses ConfigurableRuntimeAdapter

# Or with custom adapter
class MyAdapter(RuntimeAdapter):
    def run_runtime(self, card, request, inputs):
        # Custom execution logic
        ...

manager.register_model(my_card, MyAdapter(manager))
```

---

## Model Adapters

Adapters translate ModelCards into executable invocations.

### Adapter Hierarchy

```
ModelAdapterBase (abstract)
├── ExternalJobAdapter (external processes)
│   ├── BoltzAdapter
│   ├── BoltzGenAdapter
│   └── ConfigurableExternalAdapter
├── RuntimeAdapter (in-process)
│   ├── LambdaAdapter
│   └── ConfigurableRuntimeAdapter
└── PlaceholderAdapter (stub for unimplemented)
```

### Custom External Adapter

```python
class MyExternalAdapter(ExternalJobAdapter):
    def build_job(self, card, request, inputs) -> PreparedJob:
        ctx = ModelRunContext(self.paths, card)
        ctx.create()

        # Package inputs
        packaged = ctx.package_inputs(inputs)

        # Build command
        command = ["my_model", "--input", packaged["structure"]]

        return PreparedJob(
            command=command,
            working_dir=ctx.work_dir,
            artifacts=inputs,
            metadata=ctx.as_metadata()
        )
```

---

## ModelRunContext

Helper for managing job directories:

```python
class ModelRunContext:
    work_dir: Path      # data/models/<model>/<run_id>/
    inputs_dir: Path    # .../inputs/
    outputs_dir: Path   # .../outputs/
    config_path: Path   # Optional config file

    def create(self):
        """Create directory structure."""

    def package_inputs(self, bundles) -> Dict[str, str]:
        """Copy inputs to inputs_dir, return name->path mapping."""

    def as_metadata(self) -> Dict[str, str]:
        """Return paths as metadata dict."""
```

---

## Complete Example

```python
from protos import SequenceLoader, SequenceProcessor
from protos.processing.grn import GRNProcessor
from protos.models.model_manager import ModelManager

# 1. Load and prepare data
loader = SequenceLoader()
loader.download_batch(
    ["P00533", "P04637"],
    dataset_name="test_proteins"
)

# 2. Annotate with GRN
grn_proc = GRNProcessor()
seq_proc = SequenceProcessor()
sequences = seq_proc.load_dataset("test_proteins")

annotations, _ = grn_proc.annotate_sequences(
    sequences,
    reference_table="class_a_gpcr",
    protein_family="gpcr_a"
)
grn_proc.record_table("test_grn", annotations)

# 3. Prepare model
manager = ModelManager()
invocation = manager.prepare(
    "boltz2",
    inputs={
        "sequence_dataset": "test_proteins",
        "entity": "P00533"
    },
    config={
        "mutations": [{"position": 790, "mutant": "M", "original": "T"}],
        "recycling_steps": 10
    }
)

# 4. Execute
if invocation.is_external():
    job = invocation.job
    print(f"Working directory: {job.working_dir}")
    print(f"Command: {' '.join(job.command)}")
    print(f"Config file: {job.metadata.get('context', {}).get('config_path')}")

    # Run externally
    # subprocess.run(job.command, cwd=job.working_dir)
```

---

## Directory Structure

Model runs create this structure under the data root:

```
~/protos_data/
└── models/
    ├── model_registry.json          # Registered models
    ├── boltz2/
    │   └── protein_1_K790M/         # Job directory
    │       ├── inputs/              # Packaged inputs
    │       ├── outputs/             # Model outputs
    │       ├── config.yaml          # Generated config
    │       └── metadata.json        # Job metadata
    └── embedding_esm2_t12_35m/
        └── run_20240115_103045/
            ├── inputs/
            └── outputs/
                └── embeddings.npz
```
