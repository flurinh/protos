# Protos Models: AI Integration for Structural Biology

The Protos Models system provides a unified interface for integrating state-of-the-art AI models into your structural biology workflows. Models run in isolated Docker containers to avoid dependency conflicts, with automatic backend selection and seamless format conversion.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Available Models](#available-models)
4. [Quick Start](#quick-start)
5. [Docker Backend](#docker-backend)
6. [Using Models](#using-models)
7. [Model Management](#model-management)
8. [Format System](#format-system)
9. [Custom Models](#custom-models)
10. [API Reference](#api-reference)

## Overview

The Protos Models system follows the framework's core principles:

- **Zero Configuration**: Models work out of the box with sensible defaults
- **Dependency Isolation**: Each model runs in its own Docker container
- **Automatic Backend Selection**: Seamlessly switches between local and Docker execution
- **Human-Readable Names**: All models use descriptive names, not cryptic IDs
- **Entity Integration**: Predictions are tracked in the entity registry
- **Unified Interface**: All models share the same API regardless of backend

### Key Features

- **Docker-Based Isolation**: No dependency conflicts between models
- **Pre-configured Models**: ESM-2, Ankh, Boltz, Lambda, and more
- **Automatic Downloads**: Fetches model weights from official sources
- **Format Validation**: Ensures data compatibility before prediction
- **Batch Processing**: Efficient processing of multiple entities
- **Service Management**: Start/stop models as needed
- **Model Registry**: Discover and manage available models

## Architecture

The models system consists of several key components:

```
protos.models/
├── base_model.py           # Abstract base class for all models
├── standard_model.py       # Generic implementation for standard models
├── model_definitions.py    # Configuration for all standard models
├── model_registry.py       # Model discovery and management
├── model_downloader.py     # Handles weight downloading
├── model_installer.py      # Manages dependencies
└── examples/              # Usage examples
```

### Data Flow

```
Entity → Processors → Format Adapters → Model → Predictions → Entity Registry
         (Structure,    (Convert to      (AI      (Store &      (Track
          Sequence,      model input)    Model)    Return)       results)
          GRN, etc.)
```

### Directory Structure

```
data/models/
├── model_registry.json      # Central model registry
├── .cache/                  # Download cache
├── .install_registry.json   # Installation tracking
├── esm2/                    # ESM-2 models
│   ├── config.json         # Model configuration
│   ├── weights/            # Model weights
│   ├── cache/              # Cached computations
│   └── predictions/        # Saved predictions
├── ankh/                   # Ankh models
├── boltz1/                 # Boltz models
└── lambda/                 # Lambda models
```

## Available Models

### Protein Language Models

#### ESM-2 (Meta)
- **Description**: State-of-the-art protein language model
- **Variants**: 8M, 35M, 150M, 650M, 3B parameters
- **Input**: Protein sequences
- **Output**: Embeddings, contacts, attention maps
- **Use Cases**: Feature extraction, variant effect prediction, homology detection

#### Ankh (ElnaggarLab)
- **Description**: Efficient protein language model
- **Variants**: Base, Large
- **Input**: Protein sequences
- **Output**: Embeddings
- **Use Cases**: Fast embedding generation, long sequences

#### ProstT5 (Rostlab)
- **Description**: Structure-informed protein language model
- **Input**: Protein sequences
- **Output**: Embeddings with structural awareness
- **Use Cases**: Structure-aware sequence analysis

#### RITA (LightOn)
- **Description**: Family-specific protein language model
- **Variants**: Small, Extra Large
- **Input**: Protein sequences
- **Output**: Embeddings
- **Use Cases**: Protein family analysis

### Structure Prediction Models

#### Boltz-1
- **Description**: Fast structure prediction
- **Input**: Sequences (MSA optional)
- **Output**: 3D structures
- **Requirements**: GPU required

#### ESMFold (Meta)
- **Description**: End-to-end structure prediction
- **Input**: Sequences
- **Output**: 3D structures
- **Requirements**: GPU required, high memory

### Property Prediction Models

#### Lambda
- **Description**: Graph-based property prediction
- **Input**: Structures, sequences, GRN, embeddings
- **Output**: Property predictions
- **Use Cases**: GPCR/opsin analysis, custom properties

### Generative Models

#### ProtGPT2
- **Description**: Protein sequence generation
- **Input**: Seed sequences or prompts
- **Output**: Generated sequences
- **Use Cases**: De novo protein design

## Quick Start

### Basic Usage (Docker Backend - Recommended)

```python
from protos.models import predict, ModelBackend

# Quick prediction using Docker (no dependency conflicts!)
result = predict(
    "esm2",
    inputs={"sequence": "MKTAYIAKQRQISFVK"},
    backend=ModelBackend.DOCKER
)

embeddings = result["embeddings"]
```

### Using the Unified Client

```python
from protos.models import UnifiedModelClient, ModelBackend

# Create client with Docker backend
client = UnifiedModelClient(backend=ModelBackend.DOCKER)

# Make predictions
result = client.predict(
    model_name="esm2",
    entity_name="BACR_HALSA",  # Or use inputs={"sequence": "..."}
    model_variant="esm2_t33_650M"
)

# Check service status
services = client.service_manager.list_services()
print(services)  # {'esm2': ServiceStatus.READY, ...}

# Cleanup when done
client.shutdown()
```

### Automatic Backend Selection

```python
# Let Protos choose the best backend
from protos.models import UnifiedModelClient, ModelBackend

client = UnifiedModelClient(backend=ModelBackend.AUTO)

# Will use Docker if available, otherwise local
result = client.predict("esm2", entity_name="BACR_HALSA")
```

## Docker Backend

The Docker backend runs each model in an isolated container with its specific dependencies. This eliminates version conflicts between models (e.g., different PyTorch versions).

### Setup

1. **Install Docker**: Follow instructions at https://docs.docker.com/get-docker/
2. **Build model images**:
   ```bash
   cd /path/to/protos
   docker-compose -f docker/docker-compose.yml build
   ```

### Starting Models

```bash
# Start a specific model
docker-compose -f docker/docker-compose.yml up -d esm2

# Start all models
docker-compose -f docker/docker-compose.yml up -d

# Check status
docker ps --filter 'label=protos.type=model-service'
```

### Model Services

Each model runs as a REST API service:

| Model | Port | Image | GPU Memory |
|-------|------|-------|------------|
| ESM-2 | 8001 | protos/esm2:latest | 8 GB |
| Ankh | 8002 | protos/ankh:latest | 16 GB |
| Lambda | 8003 | protos/lambda:v2 | 12 GB |
| Boltz-1 | 8004 | protos/boltz:1.0 | 24 GB |
| ESMFold | 8005 | protos/esmfold:latest | 16 GB |

### Benefits

- **No Dependency Conflicts**: Each model has its own environment
- **Easy Updates**: Update one model without affecting others
- **Reproducibility**: Consistent environment across machines
- **GPU Support**: Works with nvidia-docker
- **Scalability**: Run models on different machines

### Local Backend (Fallback)

If Docker is not available, models can run locally:

```python
from protos.models import StandardModel

# Traditional local execution (requires all dependencies)
model = StandardModel("esm2", "esm2_t33_650M")
model.load_model()
prediction = model.predict("BACR_HALSA")
```

Note: Local execution requires installing all model dependencies, which may conflict.

## Model Definitions

Models are defined using configuration dictionaries in `model_definitions.py`. Each definition includes:

### Core Information
- `name`: Short identifier
- `full_name`: Display name
- `version`: Model version
- `description`: What the model does

### Technical Specifications
- `framework`: PyTorch, TensorFlow, or JAX
- `input_formats`: Required input types
- `output_format`: What the model produces
- `max_sequence_length`: Input limits

### Sources
- Download URLs for weights
- File sizes
- Checksums for verification
- Authentication requirements

### Requirements
- GPU memory needs
- CPU support
- Batch size recommendations

### Example Definition

```python
ModelDefinition(
    name="esm2",
    full_name="ESM-2",
    version="2.0",
    description="Evolutionary Scale Modeling protein language model",
    framework=ModelFramework.PYTORCH,
    input_formats=[InputFormat.SEQUENCE],
    output_format=OutputFormat.EMBEDDING,
    sources={
        "esm2_t33_650M": ModelSource(
            url="https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt",
            size_mb=2480
        )
    },
    pip_dependencies=["fair-esm"],
    requirements=ModelRequirements(
        min_gpu_memory_gb=8,
        supports_cpu=True
    )
)
```

## Using Models

### Installation and Setup

```python
from protos.models import ModelInstaller, ModelDownloader

# Check and install dependencies
installer = ModelInstaller()
installer.install_model("esm2")

# Download model weights
downloader = ModelDownloader()
downloader.download_model("esm2", "esm2_t33_650M", target_dir)
```

### Making Predictions

```python
# Single entity
model = StandardModel("esm2", "esm2_t33_650M")
model.load_model()
prediction = model.predict("P62988")  # UniProt ID

# Multiple entities
predictions = model.predict(["P62988", "EGFR_HUMAN", "1UBQ"])

# Access results
for pred in predictions:
    print(f"{pred.entity_name}: {pred.prediction.keys()}")
```

### Working with Different Input Formats

The system automatically handles format conversions:

```python
# Model needs sequence, but entity only has structure?
# Protos extracts sequence automatically
model = StandardModel("esm2")
model.load_model()

# This works even if "1UBQ" is only registered as a structure
prediction = model.predict("1UBQ")
```

### Saving and Loading Predictions

Predictions are automatically saved and tracked:

```python
# Predictions saved to: data/models/esm2/predictions/
# Format: entity_name_timestamp.json

# Load previous prediction
old_pred = model.load_prediction("BACR_HALSA")

# List all predictions for an entity
predictions = model.list_predictions("BACR_HALSA")
```

## Model Management

### Model Registry

The registry tracks all available models:

```python
from protos.models import ModelRegistry

registry = ModelRegistry()

# Discover models
models = registry.discover_models()

# Load a model
model = registry.load_model("standard/esm2_small")

# Find compatible models for an entity
compatible = registry.get_models_for_entity(
    "BACR_HALSA", 
    output_format="embedding"
)
```

### Model Information

```python
# Get model details
info = model.get_model_info()
print(f"Model: {info['full_name']}")
print(f"Version: {info['version']}")
print(f"Framework: {info['framework']}")
print(f"Device: {info['device']}")

# Get definition
from protos.models import get_model_definition
definition = get_model_definition("esm2")
print(f"Citation: {definition.citation}")
print(f"Paper: {definition.paper_url}")
```

## Format Adapters

Format adapters convert between Protos data formats and model inputs:

### Built-in Adapters

```python
from protos.models.model_definitions import StandardAdapters

# Sequence tokenization
tokenizer = StandardAdapters.sequence_to_tokens("esm")

# Structure to graph
graph_builder = StandardAdapters.structure_to_graph(
    method="knn", 
    k=10, 
    distance_threshold=8.0
)

# Embedding aggregation
aggregator = StandardAdapters.embedding_to_features("mean")
```

### Custom Adapters

```python
def my_sequence_adapter(sequence_data):
    """Custom sequence preprocessing."""
    # Your logic here
    return processed_data

model.register_input_adapter("sequence", my_sequence_adapter)
```

## Custom Models

### Adding a New Standard Model

1. Add definition to `model_definitions.py`:

```python
STANDARD_MODELS["my_model"] = ModelDefinition(
    name="my_model",
    full_name="My Custom Model",
    version="1.0",
    description="Description",
    framework=ModelFramework.PYTORCH,
    input_formats=[InputFormat.SEQUENCE],
    output_format=OutputFormat.PROPERTY,
    sources={
        "default": ModelSource(
            url="https://example.com/model.pth",
            size_mb=100
        )
    },
    pip_dependencies=["torch"],
    requirements=ModelRequirements(supports_cpu=True)
)
```

2. Use it like any standard model:

```python
model = StandardModel("my_model")
model.load_model()
predictions = model.predict("entity_name")
```

### Creating a Custom Model Class

For complex models requiring custom logic:

```python
from protos.models import BaseModel, ModelConfig

class MyCustomModel(BaseModel):
    def load_model(self, checkpoint_path=None):
        # Load your model
        self.model = load_my_model(checkpoint_path)
        self.is_loaded = True
    
    def _preprocess_input(self, input_data):
        # Convert Protos data to model format
        return processed_data
    
    def _predict_single(self, input_data):
        # Make prediction
        return self.model(input_data)

# Register the model type
from protos.models import register_model_type
register_model_type("my_custom", MyCustomModel)
```

## API Reference

### Core Classes

#### StandardModel
```python
StandardModel(
    model_name: str,
    model_variant: Optional[str] = None,
    paths: Optional[ProtosPaths] = None,
    device: Optional[str] = None
)
```

Main class for using standard models.

**Methods:**
- `load_model(checkpoint_path=None)`: Load model weights
- `predict(entity_names, save=True, batch_size=None)`: Make predictions
- `predict_dataset(dataset_name, save=True, batch_size=None)`: Predict on dataset
- `create_dataset(name, entities, metadata=None)`: Create entity dataset
- `get_model_info()`: Get model information

#### ModelRegistry
```python
ModelRegistry(paths: Optional[ProtosPaths] = None)
```

Manages model discovery and loading.

**Methods:**
- `discover_models()`: Find all available models
- `list_models(model_type=None)`: List models
- `load_model(identifier, device=None)`: Load a model
- `create_model(type, name, config)`: Create new model
- `get_models_for_entity(entity, output_format=None)`: Find compatible models

#### ModelDownloader
```python
ModelDownloader(paths: Optional[ProtosPaths] = None)
```

Handles model weight downloading.

**Methods:**
- `download_model(name, variant, target_dir)`: Download weights
- `list_downloaded_models()`: Show downloaded models
- `estimate_download_size(name, variant=None)`: Get download size
- `check_disk_space(required_mb)`: Verify available space

#### ModelInstaller
```python
ModelInstaller(paths: Optional[ProtosPaths] = None)
```

Manages model dependencies.

**Methods:**
- `check_requirements(model_name)`: Check if requirements met
- `install_model(model_name, force=False)`: Install dependencies
- `create_environment(model_name, env_name=None)`: Create conda env
- `list_installed_models()`: Show installed models

### Utility Functions

```python
# List all available models
models = list_available_models()

# Get model definition
definition = get_model_definition("esm2")

# Find models by capability
seq_models = get_models_by_input(InputFormat.SEQUENCE)
emb_models = get_models_by_output(OutputFormat.EMBEDDING)
```

### Enumerations

```python
class ModelFramework(Enum):
    PYTORCH = "pytorch"
    TENSORFLOW = "tensorflow"
    JAX = "jax"

class InputFormat(Enum):
    SEQUENCE = "sequence"
    STRUCTURE = "structure"
    EMBEDDING = "embedding"
    MSA = "msa"
    GRAPH = "graph"
    GRN = "grn"
    PROPERTY = "property"

class OutputFormat(Enum):
    EMBEDDING = "embedding"
    STRUCTURE = "structure"
    PROPERTY = "property"
    LOGITS = "logits"
    ATTENTION = "attention"
    GRAPH = "graph"
```

## Best Practices

1. **Choose the Right Model Size**: Start with smaller variants for testing
2. **Use Batch Processing**: More efficient for multiple entities
3. **Cache Predictions**: Reuse saved predictions when possible
4. **Monitor Memory**: Check GPU memory for large models
5. **Create Datasets**: Group related entities for organization

## Troubleshooting

### Common Issues

1. **Out of Memory**: Use smaller model variant or reduce batch size
2. **Missing Dependencies**: Run `ModelInstaller.install_model()`
3. **Download Failed**: Check internet connection, retry with `force=True`
4. **GPU Not Found**: Models fall back to CPU if supported

### Debug Mode

```python
import logging
logging.basicConfig(level=logging.DEBUG)

# Now model operations will show detailed logs
```

## Examples

See the `examples/` directory for complete examples:
- `model_integration_example.py`: Basic model usage
- `standard_models_example.py`: Using standard models
- Jupyter notebooks with interactive demos

## Contributing

To add support for a new model:

1. Add definition to `model_definitions.py`
2. Test with `StandardModel`
3. Add examples and documentation
4. Submit pull request

The models system is designed to be extensible while maintaining consistency with Protos principles.