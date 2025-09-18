"""Examples of using the format system for model inputs and outputs.

This demonstrates:
1. Format validation
2. Format conversion
3. Model compatibility checking
4. Adapting data for specific models
"""

import numpy as np
import pandas as pd
from protos.models import (
    SequenceFormat, StructureFormat, EmbeddingFormat,
    GraphFormat, GRNFormat, PropertyFormat,
    FormatValidator, FormatAdapter, FormatConverter,
    validate_model_compatibility, suggest_format_conversions
)


def example_format_objects():
    """Show how to create and use format objects."""
    print("=== Format Objects ===\n")
    
    # Sequence format
    seq = SequenceFormat(
        sequence="MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQD",
        sequence_id="BACR_HALSA",
        metadata={"organism": "Halobacterium salinarum", "length": 45}
    )
    print(f"Sequence format valid: {seq.validate()}")
    print(f"Sequence ID: {seq.sequence_id}")
    print(f"Length: {len(seq.sequence)}")
    
    # Structure format (mock data)
    structure_data = pd.DataFrame({
        'atom_name': ['CA', 'CB', 'CA', 'CB'],
        'auth_comp_id': ['MET', 'MET', 'LYS', 'LYS'],
        'auth_seq_id': [1, 1, 2, 2],
        'auth_chain_id': ['A', 'A', 'A', 'A'],
        'x': [10.0, 11.0, 12.0, 13.0],
        'y': [20.0, 21.0, 22.0, 23.0],
        'z': [30.0, 31.0, 32.0, 33.0]
    })
    
    struct = StructureFormat(
        coordinates=structure_data,
        pdb_id="1ABC",
        chains=['A'],
        metadata={"resolution": 2.0, "method": "X-RAY"}
    )
    print(f"\nStructure format valid: {struct.validate()}")
    print(f"PDB ID: {struct.pdb_id}")
    print(f"Number of atoms: {len(struct.coordinates)}")
    
    # Embedding format
    emb = EmbeddingFormat(
        embeddings=np.random.randn(45, 1280),  # 45 residues, 1280-dim
        embedding_type="per_residue",
        model_name="esm2",
        layer=-1
    )
    print(f"\nEmbedding shape: {emb.embeddings.shape}")
    print(f"Pooled shape: {emb.get_pooled('mean').shape}")
    
    # Property format
    prop = PropertyFormat(
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
    print(f"\nProperties: {list(prop.properties.keys())}")
    print()


def example_format_validation():
    """Show format validation."""
    print("=== Format Validation ===\n")
    
    validator = FormatValidator()
    
    # Valid sequence
    valid, error = validator.validate_input("MKTAYIAKQRQ", "sequence")
    print(f"Valid sequence: {valid}")
    
    # Invalid sequence
    valid, error = validator.validate_input("MKTAYIAKQRQ123", "sequence")
    print(f"Invalid sequence: {valid}, Error: {error}")
    
    # Valid structure
    structure_df = pd.DataFrame({
        'atom_name': ['CA'], 'auth_comp_id': ['ALA'],
        'auth_seq_id': [1], 'auth_chain_id': ['A'],
        'x': [1.0], 'y': [2.0], 'z': [3.0]
    })
    valid, error = validator.validate_input(structure_df, "structure")
    print(f"Valid structure: {valid}")
    
    # Invalid structure (missing columns)
    invalid_df = pd.DataFrame({'atom_name': ['CA'], 'x': [1.0]})
    valid, error = validator.validate_input(invalid_df, "structure")
    print(f"Invalid structure: {valid}, Error: {error}")
    
    # Valid embedding
    embedding = np.random.randn(100, 768)
    valid, error = validator.validate_input(embedding, "embedding")
    print(f"Valid embedding: {valid}")
    print()


def example_format_conversion():
    """Show format conversions."""
    print("=== Format Conversion ===\n")
    
    converter = FormatConverter()
    
    # Structure to sequence
    structure_data = pd.DataFrame({
        'atom_name': ['CA', 'CA', 'CA'],
        'auth_comp_id': ['MET', 'LYS', 'THR'],
        'auth_seq_id': [1, 2, 3],
        'auth_chain_id': ['A', 'A', 'A'],
        'x': [1.0, 2.0, 3.0],
        'y': [1.0, 2.0, 3.0],
        'z': [1.0, 2.0, 3.0]
    })
    
    struct = StructureFormat(
        coordinates=structure_data,
        pdb_id="TEST",
        chains=['A']
    )
    
    # Convert to sequence
    seq = converter.structure_to_sequence(struct, chain='A')
    print(f"Extracted sequence: {seq.sequence}")
    print(f"Sequence ID: {seq.sequence_id}")
    
    # Structure to graph
    graph = converter.structure_to_graph(struct, method="knn", k=2)
    print(f"\nGraph nodes: {graph.node_features.shape[0]}")
    print(f"Graph edges: {graph.edge_index.shape[1]}")
    print()


def example_model_compatibility():
    """Check model compatibility with entity formats."""
    print("=== Model Compatibility ===\n")
    
    # Entity has these formats
    entity_formats = ["structure", "sequence", "grn"]
    
    # Check compatibility with different models
    models = {
        "ESM-2": ["sequence"],
        "Lambda": ["structure", "sequence", "grn"],
        "Boltz": ["sequence", "msa"],
        "ESMFold": ["sequence"]
    }
    
    for model_name, required_formats in models.items():
        compatible = validate_model_compatibility(entity_formats, required_formats)
        print(f"{model_name}: {'✓ Compatible' if compatible else '✗ Incompatible'}")
        
        if not compatible:
            # Suggest conversions
            conversions = suggest_format_conversions(entity_formats, required_formats)
            if conversions:
                print(f"  Possible conversions: {conversions}")
            else:
                missing = set(required_formats) - set(entity_formats)
                print(f"  Missing formats: {missing}")
    print()


def example_model_adaptation():
    """Adapt data for specific models."""
    print("=== Model Adaptation ===\n")
    
    adapter = FormatAdapter()
    
    # Prepare data for ESM-2
    seq_data = {"sequence": SequenceFormat(sequence="MKTAYIAK", sequence_id="test")}
    esm_input = adapter.adapt_for_model(seq_data, "esm2", {})
    print(f"ESM-2 input: {esm_input}")
    
    # Prepare data for Lambda (with structure)
    struct = StructureFormat(
        coordinates=pd.DataFrame({
            'atom_name': ['CA', 'CA'],
            'auth_comp_id': ['MET', 'LYS'],
            'auth_seq_id': [1, 2],
            'auth_chain_id': ['A', 'A'],
            'x': [0.0, 3.8], 'y': [0.0, 0.0], 'z': [0.0, 0.0]
        }),
        pdb_id="TEST",
        chains=['A']
    )
    
    struct_data = {"structure": struct}
    lambda_config = {
        "preprocessing_params": {
            "method": "knn",
            "k": 5,
            "threshold": 8.0
        }
    }
    
    try:
        # This would create a PyTorch Geometric Data object
        lambda_input = adapter.adapt_for_model(struct_data, "lambda", lambda_config)
        print(f"\nLambda input type: Graph Data")
    except ImportError:
        print(f"\nLambda input: Would create Graph Data (requires torch_geometric)")
    
    print()


def example_output_adaptation():
    """Adapt model outputs back to Protos format."""
    print("=== Output Adaptation ===\n")
    
    adapter = FormatAdapter()
    
    # ESM-2 output (mock)
    esm_output = {
        "embeddings": np.random.randn(10, 1280),
        "contacts": np.random.randn(10, 10),
        "attentions": np.random.randn(12, 20, 10, 10)  # layers, heads, seq, seq
    }
    
    adapted = adapter.adapt_from_model(esm_output, "embedding", "esm2")
    if isinstance(adapted, EmbeddingFormat):
        print(f"ESM-2 output adapted to EmbeddingFormat")
        print(f"Shape: {adapted.embeddings.shape}")
        print(f"Type: {adapted.embedding_type}")
    
    # Property prediction output
    property_output = {
        "lambda_max": 470.5,
        "photocycle": "slow",
        "confidence": {"lambda_max": 0.95, "photocycle": 0.87}
    }
    
    adapted_prop = adapter.adapt_from_model(property_output, "property", "lambda")
    print(f"\nProperty output: {adapted_prop}")
    print()


def main():
    """Run all examples."""
    print("=== Protos Model Format System Examples ===\n")
    
    example_format_objects()
    example_format_validation()
    example_format_conversion()
    example_model_compatibility()
    example_model_adaptation()
    example_output_adaptation()
    
    print("=== Summary ===\n")
    print("The format system provides:")
    print("- Standardized format objects for all data types")
    print("- Validation to ensure data compatibility")
    print("- Conversion between different formats")
    print("- Model-specific adaptations")
    print("- Automatic compatibility checking")
    print("\nThis ensures reliable data flow between Protos and AI models!")


if __name__ == "__main__":
    main()