"""Demo of the format system capabilities without full model loading."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

from protos.models import (
    list_available_models, get_model_definition,
    SequenceFormat, StructureFormat, EmbeddingFormat,
    GRNFormat, PropertyFormat,
    FormatValidator, FormatAdapter, FormatConverter,
    validate_model_compatibility, suggest_format_conversions
)


def demo_complete_workflow():
    """Demonstrate a complete workflow from data to model input."""
    print("=" * 70)
    print("COMPLETE WORKFLOW: Entity → Format → Model")
    print("=" * 70)
    
    # Simulate having data for a protein
    entity_name = "BACR_HALSA"
    print(f"\nEntity: {entity_name} (Bacteriorhodopsin)")
    
    # Step 1: Create format objects for different data types
    print("\n1. Creating format objects from raw data:")
    
    # Structure data (simplified)
    structure_df = pd.DataFrame({
        'atom_name': ['CA'] * 7,
        'auth_comp_id': ['MET', 'LYS', 'THR', 'ALA', 'TYR', 'ILE', 'ALA'],
        'auth_seq_id': [1, 2, 3, 4, 5, 6, 7],
        'auth_chain_id': ['A'] * 7,
        'x': [0.0, 3.8, 7.6, 11.4, 15.2, 19.0, 22.8],
        'y': [0.0, 1.0, -0.5, 1.5, -1.0, 0.5, 0.0],
        'z': [0.0, 0.5, 1.0, 0.5, 0.0, -0.5, 0.0],
        'pdb_id': ['BACR'] * 7
    })
    
    structure = StructureFormat(
        coordinates=structure_df,
        pdb_id="BACR",
        chains=['A'],
        metadata={"organism": "Halobacterium salinarum"}
    )
    print(f"  ✓ Structure: {len(structure.coordinates)} atoms")
    
    # Sequence
    sequence = SequenceFormat(
        sequence="MKTAYIA",
        sequence_id=entity_name,
        metadata={"length": 7}
    )
    print(f"  ✓ Sequence: {sequence.sequence}")
    
    # GRN annotations
    grn = GRNFormat(
        grn_positions=pd.Series({
            '1.50': 'M', '2.50': 'K', '3.50': 'T',
            '4.50': 'A', '5.50': 'Y', '6.50': 'I', '7.50': 'A'
        }),
        sequence="MKTAYIA",
        grn_system="ballesteros_weinstein",
        metadata={"family": "microbial_rhodopsin"}
    )
    print(f"  ✓ GRN: {len(grn.grn_positions)} positions")
    
    # Embeddings (mock)
    embeddings = EmbeddingFormat(
        embeddings=np.random.randn(7, 1280),
        embedding_type="per_residue",
        model_name="esm2",
        layer=-1
    )
    print(f"  ✓ Embeddings: shape {embeddings.embeddings.shape}")
    
    # Properties
    properties = PropertyFormat(
        properties={
            "lambda_max": 568,
            "photocycle": "fast",
            "pump_type": "proton"
        },
        property_types={
            "lambda_max": "regression",
            "photocycle": "classification",
            "pump_type": "classification"
        }
    )
    print(f"  ✓ Properties: {list(properties.properties.keys())}")
    
    # Step 2: Validate all formats
    print("\n2. Validating formats:")
    validator = FormatValidator()
    
    formats = {
        "structure": structure,
        "sequence": sequence,
        "grn": grn,
        "embedding": embeddings,
        "property": properties
    }
    
    all_valid = True
    for fmt_name, fmt_obj in formats.items():
        # Direct object validation
        is_valid = fmt_obj.validate()
        print(f"  {'✓' if is_valid else '✗'} {fmt_name}: {is_valid}")
        all_valid &= is_valid
    
    # Step 3: Check model compatibility
    print("\n3. Checking model compatibility:")
    available_formats = list(formats.keys())
    
    models_to_check = [
        ("ESM-2", ["sequence"]),
        ("Lambda", ["structure", "sequence", "grn"]),
        ("Property Predictor", ["embedding"]),
        ("ESMFold", ["sequence"]),
        ("Multi-modal", ["structure", "sequence", "embedding", "grn"])
    ]
    
    for model_name, required in models_to_check:
        compatible = validate_model_compatibility(available_formats, required)
        if compatible:
            print(f"  ✓ {model_name}: Compatible")
        else:
            missing = set(required) - set(available_formats)
            print(f"  ✗ {model_name}: Missing {missing}")
    
    # Step 4: Convert between formats
    print("\n4. Format conversions:")
    converter = FormatConverter()
    
    # Structure to sequence
    extracted_seq = converter.structure_to_sequence(structure, chain='A')
    print(f"  Structure → Sequence: {extracted_seq.sequence}")
    
    # Structure to graph
    graph = converter.structure_to_graph(structure, method="knn", k=3)
    print(f"  Structure → Graph: {graph.node_features.shape[0]} nodes, {graph.edge_index.shape[1]} edges")
    
    # Step 5: Adapt for specific models
    print("\n5. Model-specific adaptations:")
    adapter = FormatAdapter()
    
    # For ESM-2
    esm_data = {"sequence": sequence}
    esm_input = adapter.adapt_for_model(esm_data, "esm2", {})
    print(f"  ESM-2 input: {esm_input}")
    
    # For Lambda
    lambda_data = {
        "structure": structure,
        "sequence": sequence,
        "grn": grn
    }
    lambda_config = {
        "preprocessing_params": {
            "method": "knn",
            "k": 5,
            "threshold": 8.0
        }
    }
    try:
        lambda_input = adapter.adapt_for_model(lambda_data, "lambda", lambda_config)
        print(f"  Lambda input: Graph format ready")
    except Exception as e:
        print(f"  Lambda input: Would create graph (PyTorch Geometric not installed)")
    
    return formats


def demo_model_capabilities():
    """Show capabilities of different models."""
    print("\n\n" + "=" * 70)
    print("MODEL CAPABILITIES")
    print("=" * 70)
    
    # Show detailed info for key models
    models_to_show = ["esm2", "lambda", "boltz1"]
    
    for model_name in models_to_show:
        definition = get_model_definition(model_name)
        print(f"\n{definition.full_name}:")
        print(f"  Purpose: {definition.description}")
        print(f"  Inputs: {[f.value for f in definition.input_formats]}")
        print(f"  Output: {definition.output_format.value}")
        
        if definition.requirements.min_gpu_memory_gb:
            print(f"  GPU Memory: {definition.requirements.min_gpu_memory_gb} GB")
        
        if definition.max_sequence_length:
            print(f"  Max length: {definition.max_sequence_length}")
        
        # Show variants
        if len(definition.sources) > 1:
            print(f"  Variants: {list(definition.sources.keys())}")


def demo_format_suggestions():
    """Show format conversion suggestions."""
    print("\n\n" + "=" * 70)
    print("FORMAT CONVERSION SUGGESTIONS")
    print("=" * 70)
    
    scenarios = [
        {
            "name": "Have structure only",
            "available": ["structure"],
            "want_model": "ESM-2",
            "required": ["sequence"]
        },
        {
            "name": "Have sequence, want structure prediction",
            "available": ["sequence"],
            "want_model": "Boltz-1",
            "required": ["sequence", "msa"]
        },
        {
            "name": "Have structure + sequence, want Lambda",
            "available": ["structure", "sequence"],
            "want_model": "Lambda",
            "required": ["structure", "sequence", "grn"]
        }
    ]
    
    for scenario in scenarios:
        print(f"\nScenario: {scenario['name']}")
        print(f"  Available: {scenario['available']}")
        print(f"  Want to use: {scenario['want_model']} (needs {scenario['required']})")
        
        missing = set(scenario['required']) - set(scenario['available'])
        if missing:
            print(f"  Missing: {missing}")
            
            # Get suggestions
            suggestions = suggest_format_conversions(scenario['available'], scenario['required'])
            if suggestions:
                print(f"  Suggested conversions:")
                for from_fmt, to_fmt in suggestions:
                    print(f"    - {from_fmt} → {to_fmt}")
            else:
                print(f"  Manual preparation needed for: {missing}")
        else:
            print(f"  ✓ All formats available!")


def main():
    """Run the demo."""
    # Run demos
    formats = demo_complete_workflow()
    demo_model_capabilities()
    demo_format_suggestions()
    
    print("\n\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("\nThe Protos format system provides:")
    print("✓ Standardized format objects for all biological data types")
    print("✓ Automatic validation of data integrity")
    print("✓ Format conversions (structure↔sequence, structure→graph)")
    print("✓ Model compatibility checking")
    print("✓ Model-specific data adaptation")
    print("✓ Clear suggestions for missing formats")
    print("\nThis ensures reliable data flow from Protos processors to AI models!")


if __name__ == "__main__":
    main()