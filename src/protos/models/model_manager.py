#!/usr/bin/env python3
"""Model Manager for preparing inputs and parsing outputs for external structure prediction/design models."""

from __future__ import annotations

import json
import shutil
import yaml
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import pandas as pd

from protos.io.paths import ProtosPaths
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor


@dataclass
class ModelInput:
    """Standardized model input container."""
    
    model: str
    config_path: Path
    data_paths: Dict[str, Path]
    metadata: Dict[str, Any]
    
    def get_command(self) -> List[str]:
        """Generate command line for model execution."""
        commands = {
            "boltz2": ["boltz", "predict", str(self.config_path)],
            "rfdiffusion": ["python", "run_inference.py", f"--config={self.config_path}"],
            "alphafold": ["python", "run_alphafold.py", f"--fasta_paths={self.data_paths.get('sequences')}"],
            "esmfold": ["python", "-m", "esm.fold", f"--input={self.data_paths.get('sequences')}"],
        }
        return commands.get(self.model, [])
    
    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return {
            "model": self.model,
            "config_path": str(self.config_path),
            "data_paths": {k: str(v) for k, v in self.data_paths.items()},
            "metadata": self.metadata
        }


@dataclass
class ProcessedOutput:
    """Container for parsed model output."""
    
    entity_name: str
    entity_format: str
    output_data: Union[pd.DataFrame, Dict]
    metadata: Dict[str, Any]
    registered: bool = False
    registration_id: Optional[str] = None


@dataclass
class BatchConfig:
    """Configuration for batch model runs."""
    
    name: str
    model: str
    inputs: List[ModelInput] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def save(self, path: Path) -> None:
        """Save batch configuration to file."""
        config_data = {
            "name": self.name,
            "model": self.model,
            "metadata": self.metadata,
            "inputs": [inp.to_dict() for inp in self.inputs]
        }
        with open(path, 'w') as f:
            json.dump(config_data, f, indent=2)


class ModelAdapter(ABC):
    """Base class for model-specific adapters."""
    
    def __init__(self, paths: ProtosPaths):
        self.paths = paths
        
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
    
    def get_model_dir(self, model_name: str) -> Path:
        """Get model-specific directory."""
        model_dir = Path(self.paths.data_root) / "models" / model_name
        model_dir.mkdir(parents=True, exist_ok=True)
        return model_dir
    
    def get_input_dir(self, model_name: str, entity_name: str, config_id: Optional[str] = None) -> Path:
        """Get input directory for specific entity/config."""
        if config_id:
            dir_name = f"{entity_name}_{config_id}"
        else:
            dir_name = entity_name
        input_dir = self.get_model_dir(model_name) / dir_name
        input_dir.mkdir(parents=True, exist_ok=True)
        return input_dir


class Boltz2Adapter(ModelAdapter):
    """Adapter for Boltz-2 structure prediction."""
    
    def prepare_input(
        self,
        entity_name: str,
        entity_format: str,
        dataset_name: Optional[str],
        config: Optional[Dict]
    ) -> ModelInput:
        """Prepare Boltz-2 input files."""
        config = config or {}
        
        # Load entity data
        sequences = {}
        if entity_format == "sequence":
            seq_proc = SequenceProcessor()
            if dataset_name:
                dataset_sequences = seq_proc.load_dataset(dataset_name)
                if entity_name in dataset_sequences:
                    sequences[entity_name] = dataset_sequences[entity_name]
            else:
                entity_data = seq_proc.load_entity(entity_name)
                if entity_data:
                    sequences[entity_name] = entity_data
                    
        elif entity_format == "structure":
            struct_proc = StructureProcessor()
            struct_df = struct_proc.load_entity(entity_name)
            if struct_df is not None:
                # Extract sequences from structure
                chain_seqs = struct_proc.collect_chain_sequences([entity_name])
                for struct_chains in chain_seqs.values():
                    for chain_id, chain_data in struct_chains.items():
                        seq_name = f"{entity_name}_chain_{chain_id}"
                        sequences[seq_name] = chain_data["sequence"]
        
        if not sequences:
            raise ValueError(f"No sequences found for entity {entity_name}")
        
        # Apply mutations if specified
        if "mutations" in config:
            sequences = self._apply_mutations(sequences, config["mutations"])
            config_id = "_".join([m["name"] for m in config["mutations"]])
        else:
            config_id = "wild_type"
        
        # Generate YAML configuration
        yaml_config = self._generate_yaml(sequences, config)
        
        # Save files
        input_dir = self.get_input_dir("boltz2", entity_name, config_id)
        yaml_path = input_dir / "config.yaml"
        fasta_path = input_dir / "sequences.fasta"
        
        # Write YAML
        with open(yaml_path, 'w') as f:
            yaml.dump(yaml_config, f, default_flow_style=False)
        
        # Write FASTA
        with open(fasta_path, 'w') as f:
            for seq_name, seq in sequences.items():
                f.write(f">{seq_name}\n{seq}\n")
        
        # Save metadata
        metadata_path = input_dir / "metadata.json"
        metadata = {
            "entity": entity_name,
            "format": entity_format,
            "dataset": dataset_name,
            "config_id": config_id,
            "mutations": config.get("mutations", [])
        }
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        return ModelInput(
            model="boltz2",
            config_path=yaml_path,
            data_paths={"sequences": fasta_path},
            metadata=metadata
        )
    
    def _apply_mutations(self, sequences: Dict[str, str], mutations: List[Dict]) -> Dict[str, str]:
        """Apply mutations to sequences."""
        mutated_sequences = {}
        
        for seq_name, seq in sequences.items():
            mutated_seq = list(seq)
            mutation_names = []
            
            for mutation in mutations:
                pos = mutation["position"] - 1  # Convert to 0-based
                if pos >= len(seq):
                    continue
                    
                original = mutation.get("original", seq[pos])
                mutant = mutation["mutant"]
                
                if seq[pos] == original:
                    mutated_seq[pos] = mutant
                    mutation_names.append(mutation["name"])
            
            if mutation_names:
                new_name = f"{seq_name}_{'_'.join(mutation_names)}"
                mutated_sequences[new_name] = "".join(mutated_seq)
            else:
                mutated_sequences[seq_name] = seq
                
        return mutated_sequences
    
    def _generate_yaml(self, sequences: Dict[str, str], config: Dict) -> Dict:
        """Generate Boltz-2 YAML configuration."""
        yaml_data = {
            "sequences": [],
            "output": {
                "name": config.get("output_name", "prediction"),
                "dir": config.get("output_dir", "./output")
            },
            "sampling": {
                "seed": config.get("seed", 42),
                "num_steps": config.get("num_steps", 200),
                "diffusion_samples": config.get("num_samples", 1)
            },
            "model": {
                "checkpoint": config.get("checkpoint", "boltz2.ckpt")
            },
            "data": {
                "crop": config.get("crop_size", 256),
                "msa": {
                    "num_msa": config.get("num_msa", 64),
                    "max_msa": config.get("max_msa", 512)
                }
            },
            "recycling_steps": config.get("recycling", 3),
            "device": config.get("device", "cuda")
        }
        
        # Add sequences
        chain_id = ord('A')
        for name, seq in sequences.items():
            yaml_data["sequences"].append({
                "protein": {
                    "id": name,
                    "sequence": seq,
                    "chain": chr(chain_id)
                }
            })
            chain_id += 1
        
        return yaml_data
    
    def parse_output(
        self,
        output_path: Path,
        register_as: Optional[str],
        dataset_name: Optional[str]
    ) -> ProcessedOutput:
        """Parse Boltz-2 output and optionally register."""
        # Boltz-2 outputs PDB files
        if output_path.suffix not in ['.pdb', '.cif']:
            raise ValueError(f"Expected PDB/CIF output, got {output_path.suffix}")
        
        # Load structure
        struct_proc = StructureProcessor()
        
        # Register if requested
        if register_as:
            # Copy to structure directory
            target_path = Path(self.paths.data_root) / "structure" / "predicted" / output_path.name
            target_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(output_path, target_path)
            
            # Register entity
            struct_proc.register_entity(
                register_as,
                "structure",
                {"path": str(target_path), "source": "boltz2_prediction"}
            )
            
            # Add to dataset if specified
            if dataset_name:
                struct_proc.add_to_dataset(dataset_name, [register_as], create=True)
            
            registered = True
            registration_id = register_as
        else:
            registered = False
            registration_id = None
        
        # Read confidence if available
        confidence_path = output_path.parent / "confidence_metrics.json"
        metadata = {}
        if confidence_path.exists():
            with open(confidence_path) as f:
                metadata["confidence"] = json.load(f)
        
        return ProcessedOutput(
            entity_name=register_as or output_path.stem,
            entity_format="structure",
            output_data={"path": str(output_path)},
            metadata=metadata,
            registered=registered,
            registration_id=registration_id
        )


class RFdiffusionAdapter(ModelAdapter):
    """Adapter for RFdiffusion protein design."""
    
    def prepare_input(
        self,
        entity_name: str,
        entity_format: str,
        dataset_name: Optional[str],
        config: Optional[Dict]
    ) -> ModelInput:
        """Prepare RFdiffusion input files."""
        config = config or {}
        config_id = config.get("run_name", "default")
        
        # Get input directory
        input_dir = self.get_input_dir("rfdiffusion", entity_name, config_id)
        
        # Handle structure input (template-based design)
        data_paths = {}
        if entity_format == "structure":
            struct_proc = StructureProcessor()
            struct_df = struct_proc.load_entity(entity_name)
            if struct_df is not None:
                # Export as CIF for RFdiffusion (convert to PDB if needed)
                template_path = input_dir / f"{entity_name}_template.cif"
                struct_proc.export_entity(entity_name, template_path, format="cif")
                data_paths["template"] = template_path
        
        # Prepare configuration
        config_dict = {
            "inference": {
                "input_pdb": str(data_paths.get("template", "")),
                "output_prefix": str(input_dir / "design"),
                "num_designs": config.get("num_designs", 10),
                "design_startnum": config.get("start_num", 0),
                "ckpt_override_path": config.get("checkpoint", None),
                "seed": config.get("seed", None),
            },
            "contigmap": {
                "contigs": config.get("contigs", ["A1-100/0 B1-50"]),
            },
            "ppi": {
                "hotspot_res": config.get("hotspots", []),
            },
            "diffuser": {
                "T": config.get("timesteps", 50),
            },
            "potentials": {
                "guiding_potentials": config.get("potentials", []),
                "guide_scale": config.get("guide_scale", 1.0),
            }
        }
        
        # Save configuration
        config_path = input_dir / "config.json"
        with open(config_path, 'w') as f:
            json.dump(config_dict, f, indent=2)
        
        # Save metadata
        metadata = {
            "entity": entity_name,
            "format": entity_format,
            "dataset": dataset_name,
            "config_id": config_id,
            "design_params": config
        }
        metadata_path = input_dir / "metadata.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        return ModelInput(
            model="rfdiffusion",
            config_path=config_path,
            data_paths=data_paths,
            metadata=metadata
        )
    
    def parse_output(
        self,
        output_path: Path,
        register_as: Optional[str],
        dataset_name: Optional[str]
    ) -> ProcessedOutput:
        """Parse RFdiffusion output."""
        # RFdiffusion outputs multiple designs
        if output_path.is_dir():
            design_files = list(output_path.glob("design_*.pdb"))
        else:
            design_files = [output_path]
        
        struct_proc = StructureProcessor()
        registered_designs = []
        
        for design_file in design_files:
            if register_as:
                # Register each design
                design_name = f"{register_as}_{design_file.stem}"
                
                # Copy to structure directory
                target_path = Path(self.paths.data_root) / "structure" / "designed" / design_file.name
                target_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(design_file, target_path)
                
                # Register entity
                struct_proc.register_entity(
                    design_name,
                    "structure",
                    {"path": str(target_path), "source": "rfdiffusion"}
                )
                
                # Add to dataset
                if dataset_name:
                    struct_proc.add_to_dataset(dataset_name, [design_name], create=True)
                
                registered_designs.append(design_name)
        
        # Load trajectory info if available
        traj_file = output_path / "inference_outputs.pkl"
        metadata = {"num_designs": len(design_files)}
        if traj_file.exists():
            metadata["trajectory_file"] = str(traj_file)
        
        return ProcessedOutput(
            entity_name=register_as or "rfdiffusion_output",
            entity_format="structure",
            output_data={"designs": [str(f) for f in design_files]},
            metadata=metadata,
            registered=len(registered_designs) > 0,
            registration_id=registered_designs
        )


class ModelManager:
    """Central interface for preparing data for external models."""
    
    def __init__(self, data_root: Optional[Path] = None):
        """Initialize model manager."""
        self.paths = ProtosPaths(data_root=str(data_root) if data_root else None)
        self.adapters: Dict[str, ModelAdapter] = {}
        self._register_default_adapters()
    
    def _register_default_adapters(self) -> None:
        """Register built-in model adapters."""
        self.register_adapter("boltz2", Boltz2Adapter(self.paths))
        self.register_adapter("rfdiffusion", RFdiffusionAdapter(self.paths))
    
    def register_adapter(self, model_name: str, adapter: ModelAdapter) -> None:
        """Register a model adapter."""
        self.adapters[model_name] = adapter
    
    def get_adapter(self, model_name: str) -> ModelAdapter:
        """Get adapter for a specific model."""
        if model_name not in self.adapters:
            raise ValueError(f"No adapter registered for model: {model_name}")
        return self.adapters[model_name]
    
    def list_models(self) -> List[str]:
        """List available models."""
        return list(self.adapters.keys())
    
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
        output_path: Union[str, Path],
        register_as: Optional[str] = None,
        dataset_name: Optional[str] = None
    ) -> ProcessedOutput:
        """Parse model output and optionally register in Protos."""
        adapter = self.get_adapter(model_name)
        output_path = Path(output_path)
        return adapter.parse_output(output_path, register_as, dataset_name)
    
    def prepare_batch(
        self,
        model_name: str,
        entity_configs: List[Dict[str, Any]],
        batch_name: str
    ) -> BatchConfig:
        """Prepare a batch of model inputs."""
        batch = BatchConfig(name=batch_name, model=model_name)
        
        for config in entity_configs:
            model_input = self.prepare_input(
                model_name=model_name,
                entity_name=config["entity"],
                entity_format=config.get("format", "sequence"),
                dataset_name=config.get("dataset"),
                config=config.get("config", {})
            )
            batch.inputs.append(model_input)
        
        # Save batch configuration
        batch_dir = Path(self.paths.data_root) / "models" / "batches"
        batch_dir.mkdir(parents=True, exist_ok=True)
        batch_file = batch_dir / f"{batch_name}.json"
        batch.save(batch_file)
        
        return batch
    
    def get_model_inputs(self, model_name: str) -> List[Path]:
        """List all prepared inputs for a model."""
        model_dir = Path(self.paths.data_root) / "models" / model_name
        if not model_dir.exists():
            return []
        
        inputs = []
        for entity_dir in model_dir.iterdir():
            if entity_dir.is_dir():
                config_files = list(entity_dir.glob("config.*"))
                inputs.extend(config_files)
        
        return inputs
    
    def get_model_outputs(self, model_name: str) -> List[Path]:
        """List all outputs for a model."""
        predictions_dir = Path(self.paths.data_root) / "predictions" / model_name
        if not predictions_dir.exists():
            return []
        
        outputs = []
        for output_file in predictions_dir.rglob("*"):
            if output_file.is_file() and output_file.suffix in [".pdb", ".cif", ".json"]:
                outputs.append(output_file)
        
        return outputs


def prepare_mutation_screen(
    seq_proc: SequenceProcessor,
    dataset_name: str,
    grn_positions: List[str],
    mutations: List[str],
    grn_table_name: Optional[str] = None
) -> List[ModelInput]:
    """Helper function to prepare mutation screen at GRN positions."""
    from protos.processing.grn import GRNProcessor
    
    manager = ModelManager()
    grn_proc = GRNProcessor()
    
    # Load GRN table
    if grn_table_name:
        grn_table = grn_proc.load_table(grn_table_name)
    else:
        grn_table = grn_proc.load_table(f"{dataset_name}_grn")
    
    # Generate mutation configs
    inputs = []
    sequences = seq_proc.load_dataset(dataset_name)
    
    for seq_name in sequences.keys():
        if seq_name not in grn_table.index:
            continue
            
        for grn_pos in grn_positions:
            if grn_pos not in grn_table.columns:
                continue
                
            # Get original residue
            grn_value = grn_table.loc[seq_name, grn_pos]
            if grn_value == '-' or not isinstance(grn_value, str) or len(grn_value) < 2:
                continue
                
            original = grn_value[0]
            seq_pos = int(grn_value[1:])
            
            # Create mutations
            for mutant in mutations:
                if mutant == original:
                    continue
                    
                config = {
                    "mutations": [{
                        "position": seq_pos,
                        "original": original,
                        "mutant": mutant,
                        "name": f"{original}{seq_pos}{mutant}_GRN{grn_pos}"
                    }],
                    "output_name": f"{seq_name}_{original}{seq_pos}{mutant}"
                }
                
                model_input = manager.prepare_input(
                    model_name="boltz2",
                    entity_name=seq_name,
                    entity_format="sequence",
                    dataset_name=dataset_name,
                    config=config
                )
                inputs.append(model_input)
    
    return inputs