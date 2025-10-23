#!/usr/bin/env python3
"""Example usage of ModelManager for preparing models inputs."""

from pathlib import Path
import sys

# Add src to path
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos
from protos.models.model_manager import ModelManager, prepare_mutation_screen
from protos.processing.sequence import SequenceProcessor
from protos.processing.structure import StructureProcessor
from protos.processing.grn import GRNProcessor


def ensure_data_root() -> Path:
    """Set up data root directory."""
    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def example_boltz2_single_sequence():
    """Example: Prepare single sequence for Boltz-2 prediction."""
    print("\n=== Example 1: Boltz-2 Single Sequence ===")
    
    manager = ModelManager()
    
    # Prepare input for a single sequence
    config = {
        "recycling": 5,
        "num_samples": 3,
        "device": "cuda",
        "crop_size": 384
    }
    
    # Assuming we have a sequence entity registered
    try:
        model_input = manager.prepare_input(
            model_name="boltz2",
            entity_name="GPCR_A2A",
            entity_format="sequence",
            config=config
        )
        
        print(f"✓ Prepared Boltz-2 input:")
        print(f"  Config: {model_input.config_path}")
        print(f"  FASTA: {model_input.data_paths['sequences']}")
        print(f"  Command: {' '.join(model_input.get_command())}")
    except Exception as e:
        print(f"✗ Could not prepare input: {e}")
        print("  (This is expected if GPCR_A2A sequence is not registered)")


def example_boltz2_mutations():
    """Example: Prepare mutations for Boltz-2 structure prediction."""
    print("\n=== Example 2: Boltz-2 Mutations ===")
    
    manager = ModelManager()
    
    # Define mutations at specific positions
    mutations = [
        {"position": 91, "original": "V", "mutant": "A", "name": "V91A"},
        {"position": 91, "original": "V", "mutant": "F", "name": "V91F"},
        {"position": 219, "original": "T", "mutant": "A", "name": "T219A"}
    ]
    
    # Create mock sequence for demo
    test_sequence = "MVLQGNAVRSLLLLLLLLGPGCGRGPPPPGPPPGPPPGPPPGLPGPPPGPPPGVPPPGPPPG" * 6
    test_sequence = test_sequence[:300]  # Reasonable length
    
    # Prepare each mutation
    inputs = []
    for mutation in mutations:
        config = {
            "mutations": [mutation],
            "recycling": 3,
            "num_samples": 1,
            "output_name": f"prediction_{mutation['name']}"
        }
        
        # For demo, we'll use a mock sequence
        print(f"\n  Preparing mutation {mutation['name']}...")
        
        # In real usage, this would use an actual registered sequence
        # model_input = manager.prepare_input(
        #     model_name="boltz2",
        #     entity_name="GPCR_A2A",
        #     entity_format="sequence",
        #     config=config
        # )
        
        print(f"  ✓ Would prepare: {mutation['name']} at position {mutation['position']}")
        print(f"    Original: {mutation['original']} → Mutant: {mutation['mutant']}")


def example_rfdiffusion_design():
    """Example: Prepare RFdiffusion binder design."""
    print("\n=== Example 3: RFdiffusion Binder Design ===")
    
    manager = ModelManager()
    
    # Design configuration
    config = {
        "contigs": ["A1-380/0 B1-60"],  # Keep chain A (receptor), design 60-residue binder
        "hotspots": [125, 130, 219, 225],  # Key interface residues
        "num_designs": 50,
        "timesteps": 50,
        "seed": 42,
        "run_name": "binder_design_v1"
    }
    
    print("  Design parameters:")
    print(f"    Receptor: chain A (380 residues)")
    print(f"    Binder: chain B (60 residues to design)")
    print(f"    Hotspots: {config['hotspots']}")
    print(f"    Number of designs: {config['num_designs']}")
    
    # In real usage with a registered structure:
    # model_input = manager.prepare_input(
    #     model_name="rfdiffusion",
    #     entity_name="5d5a",
    #     entity_format="structure",
    #     config=config
    # )


def example_grn_mutation_screen():
    """Example: Prepare mutation screen at GRN positions."""
    print("\n=== Example 4: GRN-based Mutation Screen ===")
    
    # This example shows how to prepare mutations at conserved GRN positions
    grn_positions = ["3.50", "7.50"]  # DRY motif and NPxxY motif
    amino_acids = ["A", "V", "L", "I", "F", "W", "Y"]  # Hydrophobic mutations
    
    print(f"  Target GRN positions: {grn_positions}")
    print(f"  Mutation scan: {amino_acids}")
    
    # Calculate number of mutations
    # For each sequence, at each GRN position, test each amino acid
    num_sequences = 3  # Example
    mutations_per_position = len(amino_acids) - 1  # Excluding original
    total_mutations = num_sequences * len(grn_positions) * mutations_per_position
    
    print(f"\n  Total structure predictions to prepare: ~{total_mutations}")
    print(f"  (Exact number depends on original amino acids at each position)")
    
    # In real usage:
    # inputs = prepare_mutation_screen(
    #     seq_proc=SequenceProcessor(),
    #     dataset_name="gpcr_sequences",
    #     grn_positions=grn_positions,
    #     mutations=amino_acids,
    #     grn_table_name="gpcr_sequences_grn"
    # )


def example_batch_preparation():
    """Example: Prepare batch of models runs."""
    print("\n=== Example 5: Batch Preparation ===")
    
    manager = ModelManager()
    
    # Define batch of predictions
    entity_configs = [
        {
            "entity": "GPCR_A2A",
            "format": "sequence",
            "config": {
                "mutations": [{"position": 91, "original": "V", "mutant": "A", "name": "V91A"}],
                "recycling": 3
            }
        },
        {
            "entity": "GPCR_B2AR",
            "format": "sequence",
            "config": {
                "mutations": [{"position": 113, "original": "Y", "mutant": "F", "name": "Y113F"}],
                "recycling": 3
            }
        },
        {
            "entity": "GPCR_D2R",
            "format": "sequence",
            "config": {
                "recycling": 5,  # Wild-type, higher quality
                "num_samples": 5
            }
        }
    ]
    
    # This would prepare all inputs and save batch configuration
    # batch = manager.prepare_batch(
    #     model_name="boltz2",
    #     entity_configs=entity_configs,
    #     batch_name="gpcr_mutation_study_2024"
    # )
    
    print("  Batch configuration:")
    print(f"    Model: boltz2")
    print(f"    Predictions: {len(entity_configs)}")
    print(f"    - 2 mutations")
    print(f"    - 1 wild-type (high quality)")


def example_output_parsing():
    """Example: Parse models outputs and register in Protos."""
    print("\n=== Example 6: Output Parsing and Registration ===")
    
    manager = ModelManager()
    
    # Example Boltz-2 output parsing
    print("  Boltz-2 output:")
    print("    Input: predictions/boltz2/GPCR_A2A_V91A/prediction.pdb")
    print("    Parse: Extract structure and confidence")
    print("    Register: As 'GPCR_A2A_V91A_predicted' in 'boltz2_predictions' dataset")
    
    # In real usage:
    # result = manager.parse_output(
    #     model_name="boltz2",
    #     output_path="predictions/boltz2/GPCR_A2A_V91A/prediction.pdb",
    #     register_as="GPCR_A2A_V91A_predicted",
    #     dataset_name="boltz2_predictions"
    # )
    
    print("\n  RFdiffusion output:")
    print("    Input: predictions/rfdiffusion/5d5a_binder/")
    print("    Parse: All design_*.pdb files")
    print("    Register: Each as '5d5a_binder_design_X' in 'designed_binders' dataset")


def example_model_listing():
    """Example: List available models and prepared inputs."""
    print("\n=== Example 7: Model and Input Management ===")
    
    manager = ModelManager()
    
    # List available models
    print("  Available models:")
    for model in manager.list_models():
        print(f"    - {model}")
    
    # Check for prepared inputs (would show actual files if they exist)
    print("\n  Checking for prepared inputs:")
    for model in manager.list_models():
        inputs = manager.get_model_inputs(model)
        if inputs:
            print(f"    {model}: {len(inputs)} prepared inputs")
        else:
            print(f"    {model}: No inputs prepared yet")


def main():
    """Run all examples."""
    print("ModelManager Examples")
    print("====================")
    print("This demonstrates how to prepare inputs for structure prediction/design models.")
    print("Note: These are examples - actual usage requires registered entities in Protos.")
    
    # Set up data directory
    data_root = ensure_data_root()
    print(f"\nData root: {data_root}")
    
    # Run examples
    example_boltz2_single_sequence()
    example_boltz2_mutations()
    example_rfdiffusion_design()
    example_grn_mutation_screen()
    example_batch_preparation()
    example_output_parsing()
    example_model_listing()
    
    print("\n" + "="*50)
    print("ModelManager is ready for use!")
    print("\nKey features:")
    print("- Prepare inputs for Boltz-2, RFdiffusion, and other models")
    print("- Handle mutations, batch predictions, and design tasks")
    print("- Parse outputs and register back into Protos")
    print("- Full integration with processor ecosystem")


if __name__ == "__main__":
    main()