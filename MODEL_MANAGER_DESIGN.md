# Model Manager Design for Protos

## Overview

The Model Manager acts as an intermediary between Protos processors and external structure prediction/design models (Boltz-2, RFdiffusion, AlphaFold, etc.). It handles:
1. Input data formatting from processor outputs
2. Model-specific configuration generation
3. Output parsing and registration back into Protos

## Architecture

```
Processors → ModelManager → Model Adapters → External Models
                    ↓                              ↓
             Entity Registry ← Output Parser ← Model Output
```

## Core Components

### 1. ModelManager (Main Interface)

```python
class ModelManager:
    """Central interface for preparing data for external models."""
    
    def __init__(self, data_root: Optional[Path] = None):
        self.paths = ProtosPaths(user_data_root=data_root)
        self.adapters = {}
        self._register_default_adapters()
    
    def prepare_input(
        self,
        model_name: str,
        entity_name: str,
        entity_format: str = "sequence",
        dataset_name: Optional[str] = None,
        config: Optional[Dict] = None
    ) -> ModelInput:
        """Prepare input data for a specific model."""
        adapter = self.get_adapter(model_name)
        return adapter.prepare_input(entity_name, entity_format, dataset_name, config)
    
    def parse_output(
        self,
        model_name: str,
        output_path: Path,
        register_as: Optional[str] = None,
        dataset_name: Optional[str] = None
    ) -> ProcessedOutput:
        """Parse model output and optionally register in Protos."""
        adapter = self.get_adapter(model_name)
        return adapter.parse_output(output_path, register_as, dataset_name)
```

### 2. Model Adapters (Model-Specific Logic)

```python
class ModelAdapter(ABC):
    """Base class for model-specific adapters."""
    
    @abstractmethod
    def prepare_input(
        self,
        entity_name: str,
        entity_format: str,
        dataset_name: Optional[str],
        config: Optional[Dict]
    ) -> ModelInput:
        """Prepare input files/configs for the model."""
        pass
    
    @abstractmethod
    def parse_output(
        self,
        output_path: Path,
        register_as: Optional[str],
        dataset_name: Optional[str]
    ) -> ProcessedOutput:
        """Parse model output and convert to Protos format."""
        pass


class Boltz2Adapter(ModelAdapter):
    """Adapter for Boltz-2 structure prediction."""
    
    def prepare_input(self, entity_name: str, entity_format: str, 
                     dataset_name: Optional[str], config: Optional[Dict]) -> ModelInput:
        # Load entity data
        if entity_format == "sequence":
            seq_proc = SequenceProcessor()
            sequences = seq_proc.load_entity(entity_name)
            
        # Generate YAML configuration
        yaml_config = self._generate_yaml(sequences, config)
        
        # Save to model input directory
        input_dir = self.paths.get_model_input_path("boltz2", entity_name)
        yaml_path = input_dir / "config.yaml"
        fasta_path = input_dir / "sequences.fasta"
        
        # Return ModelInput object
        return ModelInput(
            model="boltz2",
            config_path=yaml_path,
            data_paths={"sequences": fasta_path},
            metadata={"entity": entity_name, "format": entity_format}
        )
    
    def _generate_yaml(self, sequences: Dict[str, str], config: Dict) -> Dict:
        """Generate Boltz-2 YAML configuration."""
        yaml_data = {
            "sequences": [],
            "output": config.get("output_dir", "./output"),
            "num_recycling_iters": config.get("recycling", 3),
            "use_msa": config.get("use_msa", True),
            "device": config.get("device", "cuda")
        }
        
        # Handle mutations if present
        if "mutations" in config:
            for mutation in config["mutations"]:
                mutated_seq = apply_mutation(sequences, mutation)
                yaml_data["sequences"].append({
                    "name": f"{entity_name}_{mutation}",
                    "sequence": mutated_seq,
                    "chain": mutation.get("chain", "A")
                })
        else:
            # Single sequence
            for name, seq in sequences.items():
                yaml_data["sequences"].append({
                    "name": name,
                    "sequence": seq,
                    "chain": "A"
                })
        
        return yaml_data


class RFdiffusionAdapter(ModelAdapter):
    """Adapter for RFdiffusion protein design."""
    
    def prepare_input(self, entity_name: str, entity_format: str,
                     dataset_name: Optional[str], config: Optional[Dict]) -> ModelInput:
        # Load structure if needed
        if entity_format == "structure":
            struct_proc = StructureProcessor()
            structure = struct_proc.load_entity(entity_name)
            
        # Prepare input configuration
        config_dict = {
            "contigs": config.get("contigs", ["A1-100"]),
            "ppi.hotspot_res": config.get("hotspots", []),
            "inference.num_designs": config.get("num_designs", 10),
            "inference.ckpt_path": config.get("checkpoint", None)
        }
        
        # Save configuration
        input_dir = self.paths.get_model_input_path("rfdiffusion", entity_name)
        config_path = input_dir / "config.json"
        
        return ModelInput(
            model="rfdiffusion",
            config_path=config_path,
            data_paths={"template": structure_path} if entity_format == "structure" else {},
            metadata={"entity": entity_name, "format": entity_format}
        )
```

### 3. Data Classes

```python
@dataclass
class ModelInput:
    """Standardized model input container."""
    model: str
    config_path: Path
    data_paths: Dict[str, Path]
    metadata: Dict[str, Any]
    
    def get_command(self) -> List[str]:
        """Generate command line for model execution."""
        if self.model == "boltz2":
            return ["boltz", "predict", str(self.config_path)]
        elif self.model == "rfdiffusion":
            return ["python", "run_inference.py", f"--config={self.config_path}"]


@dataclass
class ProcessedOutput:
    """Container for parsed model output."""
    entity_name: str
    entity_format: str
    output_data: Union[pd.DataFrame, Dict]
    metadata: Dict[str, Any]
    registered: bool = False
```

### 4. Integration with Processors

```python
class ModelManagerIntegration:
    """Helper methods for processor integration."""
    
    @staticmethod
    def prepare_mutation_screen(
        seq_proc: SequenceProcessor,
        dataset_name: str,
        mutations: List[Dict],
        model: str = "boltz2"
    ) -> List[ModelInput]:
        """Prepare multiple mutations for structure prediction."""
        manager = ModelManager()
        inputs = []
        
        # Load sequences from dataset
        sequences = seq_proc.load_dataset(dataset_name)
        
        for seq_name in sequences:
            for mutation in mutations:
                config = {
                    "mutations": [mutation],
                    "output_dir": f"predictions/{seq_name}/{mutation['name']}"
                }
                
                model_input = manager.prepare_input(
                    model_name=model,
                    entity_name=seq_name,
                    entity_format="sequence",
                    config=config
                )
                inputs.append(model_input)
        
        return inputs
    
    @staticmethod
    def register_predictions(
        manager: ModelManager,
        output_dir: Path,
        model: str,
        dataset_prefix: str = "predicted"
    ) -> List[str]:
        """Register all predictions from a model run."""
        registered = []
        
        for output_file in output_dir.glob("*.pdb"):
            result = manager.parse_output(
                model_name=model,
                output_path=output_file,
                register_as=output_file.stem,
                dataset_name=f"{dataset_prefix}_{model}"
            )
            if result.registered:
                registered.append(result.entity_name)
        
        return registered
```

## Usage Examples

### Example 1: Preparing Mutations for Boltz-2

```python
# Initialize processors
seq_proc = SequenceProcessor()
manager = ModelManager()

# Define mutations
mutations = [
    {"position": 42, "original": "A", "mutant": "V", "name": "A42V"},
    {"position": 58, "original": "L", "mutant": "I", "name": "L58I"}
]

# Prepare inputs for each mutation
for mutation in mutations:
    config = {
        "mutations": [mutation],
        "recycling": 5,
        "use_msa": True,
        "device": "cuda"
    }
    
    model_input = manager.prepare_input(
        model_name="boltz2",
        entity_name="GPCR_sequence",
        entity_format="sequence",
        config=config
    )
    
    # Generated files:
    # data/models/boltz2/GPCR_sequence_A42V/config.yaml
    # data/models/boltz2/GPCR_sequence_A42V/sequences.fasta
    print(f"Prepared: {model_input.config_path}")
```

### Example 2: RFdiffusion Design

```python
# Design a binder for a GPCR structure
manager = ModelManager()

config = {
    "contigs": ["A1-280/0 B1-50"],  # Keep receptor, design 50-residue binder
    "hotspots": [125, 130, 200],     # Key interface residues
    "num_designs": 100,
    "checkpoint": "models/rf_diffusion_weights.pt"
}

model_input = manager.prepare_input(
    model_name="rfdiffusion",
    entity_name="5d5a",
    entity_format="structure",
    config=config
)

# Execute model (external)
# subprocess.run(model_input.get_command())

# Parse and register outputs
output_dir = Path("rfdiffusion_outputs/5d5a_binders")
results = manager.parse_output(
    model_name="rfdiffusion",
    output_path=output_dir,
    register_as="5d5a_binders",
    dataset_name="designed_binders"
)
```

### Example 3: Complete Workflow

```python
def run_structure_prediction_pipeline(
    sequence_dataset: str,
    grn_positions: List[str],
    mutations_per_position: List[str] = ["A", "V", "L", "I", "F"]
):
    """Run structure predictions for mutations at GRN positions."""
    
    # Initialize
    seq_proc = SequenceProcessor()
    grn_proc = GRNProcessor()
    manager = ModelManager()
    
    # Load GRN table
    grn_table = grn_proc.load_table(f"{sequence_dataset}_grn")
    
    # Generate mutations at specified GRN positions
    mutation_configs = []
    for grn_pos in grn_positions:
        # Find sequence position for GRN
        seq_positions = grn_proc.get_sequence_positions_for_grn(grn_pos)
        
        for seq_name, seq_pos in seq_positions.items():
            original = grn_table.loc[seq_name, grn_pos][0]  # Original amino acid
            
            for mutant in mutations_per_position:
                if mutant != original:
                    mutation_configs.append({
                        "sequence": seq_name,
                        "mutation": {
                            "position": seq_pos,
                            "original": original,
                            "mutant": mutant,
                            "name": f"{original}{seq_pos}{mutant}_GRN{grn_pos}"
                        }
                    })
    
    # Prepare all inputs
    inputs = []
    for config in mutation_configs:
        model_input = manager.prepare_input(
            model_name="boltz2",
            entity_name=config["sequence"],
            entity_format="sequence",
            config={"mutations": [config["mutation"]]}
        )
        inputs.append(model_input)
    
    # Save batch configuration
    batch_file = manager.save_batch_config(inputs, "grn_mutation_screen")
    print(f"Prepared {len(inputs)} predictions")
    print(f"Batch config: {batch_file}")
    
    return inputs
```

## Directory Structure

```
data/
├── models/
│   ├── boltz2/
│   │   ├── <entity_name>_<config>/
│   │   │   ├── config.yaml
│   │   │   ├── sequences.fasta
│   │   │   └── metadata.json
│   ├── rfdiffusion/
│   │   ├── <entity_name>_<config>/
│   │   │   ├── config.json
│   │   │   ├── template.pdb
│   │   │   └── metadata.json
│   └── alphafold/
│       └── ...
└── predictions/
    ├── boltz2/
    │   └── <entity_name>/
    │       ├── predicted_structure.pdb
    │       ├── confidence.json
    │       └── metadata.json
    └── rfdiffusion/
        └── <entity_name>/
            ├── design_*.pdb
            └── trajectories/
```

## Key Features

1. **Zero Configuration**: Follows Protos patterns - works out of the box
2. **Entity Registry Integration**: Automatically registers predictions
3. **Batch Processing**: Support for mutation screens and design campaigns
4. **Format Agnostic**: Works with sequences, structures, or any processor output
5. **Extensible**: Easy to add new model adapters
6. **Traceable**: Full provenance tracking through metadata

## Implementation Notes

1. Model adapters handle all model-specific formatting
2. ModelManager provides unified interface regardless of model
3. Outputs automatically registered in entity registry
4. Support for both single predictions and batch campaigns
5. Integration with existing processor workflows
6. Model execution is external - manager only prepares inputs/parses outputs