"""Standard models definitions for Protos.

This module defines configurations for standard protein AI models,
including their requirements, inputs, outputs, and download sources.
"""

from typing import Dict, List, Optional, Any, Callable
from dataclasses import dataclass, field
from enum import Enum


class ModelFramework(Enum):
    """Supported models frameworks."""
    PYTORCH = "pytorch"
    TENSORFLOW = "tensorflow"
    JAX = "jax"
    CUSTOM = "custom"


class InputFormat(Enum):
    """Standard input formats."""
    SEQUENCE = "sequence"
    STRUCTURE = "structure"
    EMBEDDING = "embedding"
    MSA = "msa"  # Multiple sequence alignment
    GRAPH = "graph"
    GRN = "grn"
    PROPERTY = "property"


class OutputFormat(Enum):
    """Standard output formats."""
    EMBEDDING = "embedding"
    STRUCTURE = "structure"
    PROPERTY = "property"
    LOGITS = "logits"
    ATTENTION = "attention"
    GRAPH = "graph"
    SEQUENCE = "sequence"  # For generative models


@dataclass
class ModelSource:
    """Model source information."""
    url: str  # GitHub URL or direct download link
    format: str = "weights"  # weights, checkpoint, onnx, etc.
    size_mb: Optional[float] = None
    checksum: Optional[str] = None
    requires_auth: bool = False
    auth_instructions: Optional[str] = None


@dataclass
class ModelRequirements:
    """Model computational requirements."""
    min_gpu_memory_gb: Optional[float] = None
    recommended_gpu_memory_gb: Optional[float] = None
    supports_cpu: bool = True
    supports_mps: bool = False  # Apple Silicon
    batch_size_recommendations: Dict[str, int] = field(default_factory=dict)


@dataclass
class ModelDefinition:
    """Complete definition of a standard models."""
    # Basic info
    name: str
    full_name: str
    version: str
    description: str
    
    # Technical specs
    framework: ModelFramework
    input_formats: List[InputFormat]
    output_format: OutputFormat
    max_sequence_length: Optional[int] = None
    
    # Source and installation
    sources: Dict[str, ModelSource] = field(default_factory=dict)
    pip_dependencies: List[str] = field(default_factory=list)
    conda_dependencies: List[str] = field(default_factory=list)
    setup_commands: List[str] = field(default_factory=list)
    
    # Requirements
    requirements: ModelRequirements = field(default_factory=ModelRequirements)
    
    # Processing specs
    preprocessing_config: Dict[str, Any] = field(default_factory=dict)
    postprocessing_config: Dict[str, Any] = field(default_factory=dict)
    
    # Model-specific parameters
    model_config: Dict[str, Any] = field(default_factory=dict)
    default_params: Dict[str, Any] = field(default_factory=dict)
    
    # Citations
    paper_url: Optional[str] = None
    citation: Optional[str] = None


# Standard models definitions
STANDARD_MODELS = {
    "esm2": ModelDefinition(
        name="esm2",
        full_name="ESM-2",
        version="2.0",
        description="Evolutionary Scale Modeling protein language models from Meta",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE],
        output_format=OutputFormat.EMBEDDING,
        max_sequence_length=1024,
        sources={
            "esm2_t6_8M": ModelSource(
                url="https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t6_8M_UR50D.pt",
                size_mb=30
            ),
            "esm2_t12_35M": ModelSource(
                url="https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t12_35M_UR50D.pt",
                size_mb=134
            ),
            "esm2_t30_150M": ModelSource(
                url="https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t30_150M_UR50D.pt",
                size_mb=572
            ),
            "esm2_t33_650M": ModelSource(
                url="https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt",
                size_mb=2480
            ),
            "esm2_t36_3B": ModelSource(
                url="https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt",
                size_mb=11400
            )
        },
        pip_dependencies=["fair-esm"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=2,
            recommended_gpu_memory_gb=8,
            supports_cpu=True,
            batch_size_recommendations={
                "8M": 32,
                "35M": 16,
                "150M": 8,
                "650M": 4,
                "3B": 1
            }
        ),
        preprocessing_config={
            "tokenizer": "esm",
            "add_special_tokens": True,
            "truncation": True,
            "max_length": 1024
        },
        model_config={
            "repr_layers": [-1],  # Last layer
            "include": ["mean", "per_tok", "contacts"]
        },
        paper_url="https://www.biorxiv.org/content/10.1101/2022.07.20.500902v1",
        citation="Lin et al. 2022"
    ),
    
    "ankh": ModelDefinition(
        name="ankh",
        full_name="Ankh",
        version="1.0",
        description="Large protein language models optimized for efficiency",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE],
        output_format=OutputFormat.EMBEDDING,
        max_sequence_length=1024,
        sources={
            "base": ModelSource(
                url="https://huggingface.co/ElnaggarLab/ankh-base",
                format="huggingface",
                size_mb=600
            ),
            "large": ModelSource(
                url="https://huggingface.co/ElnaggarLab/ankh-large", 
                format="huggingface",
                size_mb=2400
            )
        },
        pip_dependencies=["transformers", "tokenizers"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=4,
            recommended_gpu_memory_gb=16,
            supports_cpu=True,
            supports_mps=True
        ),
        preprocessing_config={
            "tokenizer": "ankh",
            "max_length": 1024
        },
        paper_url="https://arxiv.org/abs/2301.06568",
        citation="Elnaggar et al. 2023"
    ),
    
    "boltz1": ModelDefinition(
        name="boltz1",
        full_name="Boltz-1",
        version="1.0",
        description="Structure prediction models",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE, InputFormat.MSA],
        output_format=OutputFormat.STRUCTURE,
        sources={
            "models": ModelSource(
                url="https://github.com/jwohlwend/boltz/releases/download/v1.0.0/boltz1.ckpt",
                format="checkpoint",
                size_mb=800
            )
        },
        pip_dependencies=["boltz"],
        conda_dependencies=["pytorch", "pytorch-cuda"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=8,
            recommended_gpu_memory_gb=24,
            supports_cpu=False
        ),
        preprocessing_config={
            "msa_required": False,
            "confidence_threshold": 0.7
        },
        paper_url="https://github.com/jwohlwend/boltz",
        citation="Wohlwend et al. 2024"
    ),
    
    "lambda": ModelDefinition(
        name="lambda",
        full_name="Lambda",
        version="1.0",
        description="Graph-based property prediction for proteins",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.STRUCTURE, InputFormat.SEQUENCE, InputFormat.GRAPH, InputFormat.GRN],
        output_format=OutputFormat.PROPERTY,
        sources={
            "gpcr_model": ModelSource(
                url="https://github.com/user/lambda/releases/download/v1.0/gpcr_model.pth",
                size_mb=150
            ),
            "opsin_model": ModelSource(
                url="https://github.com/user/lambda/releases/download/v1.0/opsin_model.pth",
                size_mb=150
            )
        },
        pip_dependencies=["torch-geometric", "torch-scatter", "torch-sparse"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=4,
            recommended_gpu_memory_gb=8,
            supports_cpu=True
        ),
        preprocessing_config={
            "graph_construction": {
                "method": "knn",
                "k": 10,
                "distance_threshold": 8.0
            },
            "node_features": ["residue_type", "coordinates", "embeddings"],
            "edge_features": ["distance", "orientation"]
        },
        model_config={
            "architecture": "GraphTransformer",
            "hidden_dim": 256,
            "num_layers": 4,
            "num_heads": 8
        }
    ),
    
    "esmfold": ModelDefinition(
        name="esmfold",
        full_name="ESMFold",
        version="1.0",
        description="End-to-end structure prediction using ESM-2",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE],
        output_format=OutputFormat.STRUCTURE,
        max_sequence_length=512,
        sources={
            "models": ModelSource(
                url="https://dl.fbaipublicfiles.com/fair-esm/models/esmfold_3B_v1.pt",
                size_mb=11000,
                requires_auth=False
            )
        },
        pip_dependencies=["fair-esm[esmfold]", "openmm"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=16,
            recommended_gpu_memory_gb=40,
            supports_cpu=False
        ),
        preprocessing_config={
            "chunk_size": 256,
            "overlap": 128
        },
        paper_url="https://www.science.org/doi/10.1126/science.add2187",
        citation="Lin et al. 2023"
    ),
    
    "protgpt2": ModelDefinition(
        name="protgpt2",
        full_name="ProtGPT2",
        version="1.0",
        description="Generative protein language models",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE],
        output_format=OutputFormat.SEQUENCE,
        sources={
            "models": ModelSource(
                url="https://huggingface.co/nferruz/ProtGPT2",
                format="huggingface",
                size_mb=1500
            )
        },
        pip_dependencies=["transformers"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=4,
            recommended_gpu_memory_gb=8,
            supports_cpu=True
        ),
        preprocessing_config={
            "generation_params": {
                "max_length": 512,
                "temperature": 1.0,
                "top_k": 950,
                "repetition_penalty": 1.2
            }
        },
        paper_url="https://www.nature.com/articles/s41467-022-32007-7",
        citation="Ferruz et al. 2022"
    ),
    
    "prostt5": ModelDefinition(
        name="prostt5",
        full_name="ProstT5",
        version="1.0",
        description="Structure-informed protein language models",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE],
        output_format=OutputFormat.EMBEDDING,
        max_sequence_length=1024,
        sources={
            "models": ModelSource(
                url="https://huggingface.co/Rostlab/ProstT5",
                format="huggingface",
                size_mb=3000
            )
        },
        pip_dependencies=["transformers", "sentencepiece"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=8,
            recommended_gpu_memory_gb=16,
            supports_cpu=True
        ),
        preprocessing_config={
            "add_special_tokens": True,
            "tokenizer": "t5"
        },
        paper_url="https://www.biorxiv.org/content/10.1101/2023.07.23.550085v1",
        citation="Heinzinger et al. 2023"
    ),
    
    "rita": ModelDefinition(
        name="rita",
        full_name="RITA",
        version="1.0",
        description="Protein family-specific language models",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.SEQUENCE],
        output_format=OutputFormat.EMBEDDING,
        sources={
            "base": ModelSource(
                url="https://huggingface.co/lightonai/RITA_s",
                format="huggingface",
                size_mb=400
            ),
            "large": ModelSource(
                url="https://huggingface.co/lightonai/RITA_xl",
                format="huggingface", 
                size_mb=2400
            )
        },
        pip_dependencies=["transformers"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=4,
            recommended_gpu_memory_gb=8,
            supports_cpu=True
        ),
        preprocessing_config={
            "family_specific": True,
            "max_length": 1024
        },
        paper_url="https://arxiv.org/abs/2205.05789",
        citation="Hesslow et al. 2022"
    ),

    "rfdiffusion2": ModelDefinition(
        name="rfdiffusion2",
        full_name="RFdiffusion2",
        version="2.0",
        description="All-atom protein structure generation via diffusion",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.STRUCTURE],
        output_format=OutputFormat.STRUCTURE,
        max_sequence_length=None,  # Variable length generation
        sources={
            "sif_container": ModelSource(
                url="https://files.ipd.uw.edu/pub/rfdiffusion2/sifs/rfdiffusion.sif",
                format="singularity",
                size_mb=1750
            ),
            "weights_173": ModelSource(
                url="https://files.ipd.uw.edu/pub/rfdiffusion2/model_weights/RFD_173.pt",
                format="weights",
                size_mb=800
            ),
            "weights_140": ModelSource(
                url="https://files.ipd.uw.edu/pub/rfdiffusion2/model_weights/RFD_140.pt",
                format="weights",
                size_mb=800
            ),
            "ligandmpnn": ModelSource(
                url="https://files.ipd.uw.edu/pub/rfdiffusion2/third_party_model_weights/ligand_mpnn/",
                format="weights",
                size_mb=200
            )
        },
        pip_dependencies=[],
        conda_dependencies=[],
        setup_commands=["python setup.py"],
        requirements=ModelRequirements(
            min_gpu_memory_gb=16,
            recommended_gpu_memory_gb=40,
            supports_cpu=False,
            batch_size_recommendations={
                "small": 4,
                "medium": 2,
                "large": 1
            }
        ),
        preprocessing_config={
            "input_format": "pdb",
            "contigs": "required",  # e.g., "A1-100/0 B1-50"
            "hotspot_residues": "optional",
            "partial_diffusion": True,
        },
        postprocessing_config={
            "output_format": "pdb",
            "num_designs": 10,
            "ligandmpnn_integration": True,
        },
        model_config={
            "diffusion_steps": 50,
            "noise_scale": 1.0,
            "partial_T": 20,  # For partial diffusion
        },
        default_params={
            "inference.num_designs": 10,
            "inference.output_prefix": "design",
            "diffuser.partial_T": 20,
        },
        paper_url="https://www.biorxiv.org/content/10.1101/2024.11.20.624467v1",
        citation="Baker Lab 2024"
    ),

    "boltzgen": ModelDefinition(
        name="boltzgen",
        full_name="BoltzGen",
        version="1.0",
        description="All-atom generative model for protein binder design",
        framework=ModelFramework.PYTORCH,
        input_formats=[InputFormat.STRUCTURE],
        output_format=OutputFormat.STRUCTURE,
        sources={
            "docker": ModelSource(
                url="docker://protos/boltzgen:latest",
                format="docker",
                size_mb=5000
            )
        },
        pip_dependencies=[],
        requirements=ModelRequirements(
            min_gpu_memory_gb=8,
            recommended_gpu_memory_gb=24,
            supports_cpu=False
        ),
        preprocessing_config={
            "config_format": "yaml",
            "design_regions": "required",
            "constraints": "optional",
        },
        model_config={
            "num_designs": 24,
            "pipeline_steps": ["design", "inverse_folding", "folding", "design_folding", "analysis", "filtering"],
        },
        paper_url="https://www.biorxiv.org/content/10.1101/2025.11.20.689494v1",
        citation="BoltzGen Team 2025"
    )
}


# Format adapter function signatures
FormatAdapter = Callable[[Any], Any]


# Standard format adapters
class StandardAdapters:
    """Collection of standard format adapters."""
    
    @staticmethod
    def sequence_to_tokens(tokenizer_type: str = "esm"):
        """Create adapter for sequence to tokens."""
        def adapter(sequence: str) -> Dict[str, Any]:
            if tokenizer_type == "esm":
                # ESM tokenization logic
                return {"input_ids": sequence}  # Simplified
            elif tokenizer_type == "ankh":
                # Ankh tokenization
                return {"input_ids": sequence}
            else:
                raise ValueError(f"Unknown tokenizer: {tokenizer_type}")
        return adapter
    
    @staticmethod
    def structure_to_graph(method: str = "knn", **kwargs):
        """Create adapter for structure to graph."""
        def adapter(structure_df) -> Dict[str, Any]:
            # Convert structure DataFrame to graph
            # This would contain actual implementation
            return {
                "node_features": None,
                "edge_index": None,
                "edge_features": None
            }
        return adapter
    
    @staticmethod
    def embedding_to_features(aggregation: str = "mean"):
        """Create adapter for embeddings to features."""
        def adapter(embeddings) -> Any:
            if aggregation == "mean":
                return embeddings.mean(axis=0)
            elif aggregation == "max":
                return embeddings.max(axis=0)
            else:
                return embeddings
        return adapter
    
    @staticmethod
    def grn_to_features(encoding: str = "onehot"):
        """Create adapter for GRN to features."""
        def adapter(grn_data) -> Any:
            # Convert GRN annotations to features
            return grn_data  # Simplified
        return adapter


def get_model_definition(model_name: str) -> ModelDefinition:
    """Get a models definition by name."""
    if model_name not in STANDARD_MODELS:
        raise ValueError(f"Unknown models: {model_name}. Available: {list(STANDARD_MODELS.keys())}")
    return STANDARD_MODELS[model_name]


def list_available_models() -> List[str]:
    """List all available standard models."""
    return list(STANDARD_MODELS.keys())


def get_models_by_input(input_format: InputFormat) -> List[str]:
    """Get models that accept a specific input format."""
    models = []
    for name, definition in STANDARD_MODELS.items():
        if input_format in definition.input_formats:
            models.append(name)
    return models


def get_models_by_output(output_format: OutputFormat) -> List[str]:
    """Get models that produce a specific output format."""
    models = []
    for name, definition in STANDARD_MODELS.items():
        if definition.output_format == output_format:
            models.append(name)
    return models