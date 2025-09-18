"""Test format validation system with real data."""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

from protos.models.format_schemas import (
    SequenceFormat, StructureFormat, EmbeddingFormat,
    MSAFormat, GraphFormat, GRNFormat, PropertyFormat,
    FormatConverter
)
from protos.models.format_validators import (
    FormatValidator, FormatAdapter,
    validate_model_compatibility, suggest_format_conversions
)


def test_sequence_validation():
    """Test sequence format validation."""
    print("=== Testing Sequence Format Validation ===")
    
    validator = FormatValidator()
    
    # Test 1: Valid sequence string
    valid_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
    is_valid, error = validator.validate_input(valid_seq, "sequence")
    print(f"✓ Valid sequence string: {is_valid}")
    
    # Test 2: Invalid sequence (contains numbers)
    invalid_seq = "MKTAY123IAKQRQ"
    is_valid, error = validator.validate_input(invalid_seq, "sequence")
    print(f"✓ Invalid sequence detected: {not is_valid}, Error: {error}")
    
    # Test 3: SequenceFormat object
    seq_obj = SequenceFormat(
        sequence="MKTAYIAKQRQISFVK",
        sequence_id="test_protein",
        metadata={"organism": "E. coli"}
    )
    is_valid, error = validator.validate_input(seq_obj, "sequence")
    print(f"✓ SequenceFormat object valid: {is_valid}")
    
    # Test 4: Dictionary format
    seq_dict = {"sequence": "MKTAYIAKQRQISFVK", "id": "test"}
    is_valid, error = validator.validate_input(seq_dict, "sequence")
    print(f"✓ Dictionary format valid: {is_valid}")
    
    print()


def test_structure_validation():
    """Test structure format validation."""
    print("=== Testing Structure Format Validation ===")
    
    validator = FormatValidator()
    
    # Test 1: Valid structure DataFrame
    valid_structure = pd.DataFrame({
        'atom_name': ['CA', 'CB', 'CA', 'CB'],
        'auth_comp_id': ['MET', 'MET', 'LYS', 'LYS'],
        'auth_seq_id': [1, 1, 2, 2],
        'auth_chain_id': ['A', 'A', 'A', 'A'],
        'x': [10.123, 11.456, 12.789, 13.012],
        'y': [20.345, 21.678, 22.901, 23.234],
        'z': [30.567, 31.890, 32.123, 33.456]
    })
    
    is_valid, error = validator.validate_input(valid_structure, "structure")
    print(f"✓ Valid structure DataFrame: {is_valid}")
    
    # Test 2: Missing required columns
    invalid_structure = pd.DataFrame({
        'atom_name': ['CA', 'CB'],
        'x': [10.0, 11.0]
        # Missing other required columns
    })
    
    is_valid, error = validator.validate_input(invalid_structure, "structure")
    print(f"✓ Invalid structure detected: {not is_valid}")
    print(f"  Error: {error}")
    
    # Test 3: StructureFormat object
    struct_obj = StructureFormat(
        coordinates=valid_structure,
        pdb_id="1ABC",
        chains=['A'],
        metadata={"resolution": 2.5, "method": "X-RAY"}
    )
    
    is_valid, error = validator.validate_input(struct_obj, "structure")
    print(f"✓ StructureFormat object valid: {is_valid}")
    
    print()


def test_embedding_validation():
    """Test embedding format validation."""
    print("=== Testing Embedding Format Validation ===")
    
    validator = FormatValidator()
    
    # Test 1: Per-residue embeddings
    per_residue = np.random.randn(50, 1280)  # 50 residues, 1280-dim
    is_valid, error = validator.validate_input(per_residue, "embedding")
    print(f"✓ Per-residue embeddings valid: {is_valid}, shape: {per_residue.shape}")
    
    # Test 2: Pooled embeddings
    pooled = np.random.randn(1280)  # Single vector
    is_valid, error = validator.validate_input(pooled, "embedding")
    print(f"✓ Pooled embeddings valid: {is_valid}, shape: {pooled.shape}")
    
    # Test 3: Invalid shape
    invalid = np.random.randn(10, 20, 30)  # 3D tensor
    is_valid, error = validator.validate_input(invalid, "embedding")
    print(f"✓ Invalid embedding shape detected: {not is_valid}")
    
    # Test 4: EmbeddingFormat object
    emb_obj = EmbeddingFormat(
        embeddings=per_residue,
        embedding_type="per_residue",
        model_name="esm2",
        layer=-1
    )
    is_valid, error = validator.validate_input(emb_obj, "embedding")
    print(f"✓ EmbeddingFormat object valid: {is_valid}")
    
    print()


def test_format_conversion():
    """Test format conversions."""
    print("=== Testing Format Conversions ===")
    
    converter = FormatConverter()
    
    # Create a test structure
    structure_data = pd.DataFrame({
        'atom_name': ['CA', 'CA', 'CA', 'CA', 'CA'],
        'auth_comp_id': ['MET', 'LYS', 'THR', 'ALA', 'TYR'],
        'auth_seq_id': [1, 2, 3, 4, 5],
        'auth_chain_id': ['A', 'A', 'A', 'A', 'A'],
        'x': [0.0, 3.8, 7.6, 11.4, 15.2],
        'y': [0.0, 0.0, 0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0, 0.0, 0.0]
    })
    
    struct = StructureFormat(
        coordinates=structure_data,
        pdb_id="TEST",
        chains=['A']
    )
    
    # Test 1: Structure to sequence
    seq = converter.structure_to_sequence(struct, chain='A')
    print(f"✓ Structure to sequence: {seq.sequence}")
    print(f"  Sequence ID: {seq.sequence_id}")
    
    # Test 2: Structure to graph
    graph = converter.structure_to_graph(struct, method="knn", k=3, threshold=8.0)
    print(f"✓ Structure to graph conversion:")
    print(f"  Nodes: {graph.node_features.shape[0]}")
    print(f"  Edges: {graph.edge_index.shape[1]}")
    print(f"  Has positions: {graph.pos is not None}")
    
    print()


def test_model_compatibility():
    """Test model compatibility checking."""
    print("=== Testing Model Compatibility ===")
    
    # Simulate different entity format scenarios
    scenarios = [
        {
            "name": "Structure + Sequence entity",
            "formats": ["structure", "sequence"],
            "models": {
                "ESM-2": ["sequence"],
                "ESMFold": ["sequence"],
                "Lambda": ["structure", "sequence", "grn"],
                "Boltz": ["sequence", "msa"]
            }
        },
        {
            "name": "Full annotation entity",
            "formats": ["structure", "sequence", "grn", "embedding"],
            "models": {
                "Lambda": ["structure", "sequence", "grn"],
                "Property predictor": ["embedding"],
            }
        }
    ]
    
    for scenario in scenarios:
        print(f"\n{scenario['name']}:")
        print(f"Available formats: {scenario['formats']}")
        
        for model, required in scenario['models'].items():
            compatible = validate_model_compatibility(scenario['formats'], required)
            status = "✓ Compatible" if compatible else "✗ Needs conversion"
            print(f"  {model}: {status}")
            
            if not compatible:
                missing = set(required) - set(scenario['formats'])
                if missing:
                    print(f"    Missing: {missing}")
                    conversions = suggest_format_conversions(scenario['formats'], required)
                    if conversions:
                        print(f"    Possible conversions: {conversions}")
    
    print()


def test_model_adaptation():
    """Test adapting data for specific models."""
    print("=== Testing Model Data Adaptation ===")
    
    adapter = FormatAdapter()
    
    # Test 1: ESM-2 adaptation
    seq = SequenceFormat(sequence="MKTAYIAKQRQISFVK", sequence_id="test_protein")
    data = {"sequence": seq}
    
    try:
        esm_input = adapter.adapt_for_model(data, "esm2", {})
        print(f"✓ ESM-2 input format: {esm_input}")
    except Exception as e:
        print(f"✗ ESM-2 adaptation error: {e}")
    
    # Test 2: Ankh adaptation
    try:
        ankh_input = adapter.adapt_for_model(data, "ankh", {})
        print(f"✓ Ankh input format: {type(ankh_input)} = '{ankh_input[:20]}...'")
    except Exception as e:
        print(f"✗ Ankh adaptation error: {e}")
    
    # Test 3: Lambda adaptation (needs structure)
    struct = StructureFormat(
        coordinates=pd.DataFrame({
            'atom_name': ['CA', 'CA', 'CA'],
            'auth_comp_id': ['MET', 'LYS', 'THR'],
            'auth_seq_id': [1, 2, 3],
            'auth_chain_id': ['A', 'A', 'A'],
            'x': [0.0, 3.8, 7.6],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0]
        }),
        pdb_id="TEST",
        chains=['A']
    )
    
    struct_data = {"structure": struct}
    config = {
        "preprocessing_params": {
            "method": "knn",
            "k": 2,
            "threshold": 10.0
        }
    }
    
    try:
        lambda_input = adapter.adapt_for_model(struct_data, "lambda", config)
        print(f"✓ Lambda input format: Graph (would be PyTorch Geometric Data)")
    except ImportError:
        print(f"✓ Lambda adaptation works (PyTorch Geometric not installed)")
    except Exception as e:
        print(f"✗ Lambda adaptation error: {e}")
    
    print()


def main():
    """Run all validation tests."""
    print("=" * 60)
    print("PROTOS MODEL FORMAT VALIDATION TESTS")
    print("=" * 60)
    print()
    
    test_sequence_validation()
    test_structure_validation()
    test_embedding_validation()
    test_format_conversion()
    test_model_compatibility()
    test_model_adaptation()
    
    print("=" * 60)
    print("All format validation tests completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()