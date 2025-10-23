"""Test integration between Protos processors and models formats."""

import os
import sys
from pathlib import Path
import numpy as np

# Set up paths
os.environ["PROTOS_DATA_ROOT"] = str(Path.cwd() / "data")

from protos.processing.structure.structure_processor import StructureProcessor
from protos.processing.sequence.sequence_processor import SequenceProcessor
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.property.property_processor import PropertyProcessor
from protos.processing.embedding.embedding_processor import EmbeddingProcessor

from protos.models.format_schemas import (
    SequenceFormat, StructureFormat, EmbeddingFormat,
    GRNFormat, PropertyFormat, FormatConverter
)
from protos.models.format_validators import (
    FormatValidator, FormatAdapter
)


def test_structure_processor_to_model():
    """Test converting CifBaseProcessor data to models formats."""
    print("=== Testing Structure Processor → Model Format ===")
    
    # Initialize processor
    struct_proc = StructureProcessor(name="test_structures")
    
    # Check if we have any test structures
    entities = struct_proc.list_entities()
    print(f"Available structures: {len(entities)}")
    
    if entities:
        # Load first available structure
        test_pdb = entities[0]
        print(f"\nUsing structure: {test_pdb}")
        
        try:
            # Load structure data
            structure_df = struct_proc.load_entity(test_pdb)
            
            if structure_df is not None and not structure_df.empty:
                print(f"  Loaded {len(structure_df)} atoms")
                
                # Get unique chains
                chains = structure_df['auth_chain_id'].unique()
                print(f"  Chains: {chains}")
                
                # Create StructureFormat
                struct_format = StructureFormat(
                    coordinates=structure_df,
                    pdb_id=test_pdb,
                    chains=list(chains)
                )
                
                # Validate
                validator = FormatValidator()
                is_valid, error = validator.validate_input(struct_format, "structure")
                print(f"  StructureFormat valid: {is_valid}")
                
                # Convert to sequence
                if 'A' in chains:
                    converter = FormatConverter()
                    seq_format = converter.structure_to_sequence(struct_format, chain='A')
                    print(f"  Extracted sequence: {seq_format.sequence[:50]}...")
                    print(f"  Sequence length: {len(seq_format.sequence)}")
                    
                    # Convert to graph
                    graph_format = converter.structure_to_graph(
                        struct_format, 
                        method="knn", 
                        k=10,
                        threshold=8.0
                    )
                    print(f"  Graph conversion: {graph_format.node_features.shape[0]} nodes, "
                          f"{graph_format.edge_index.shape[1]} edges")
                
        except Exception as e:
            print(f"  Error loading structure: {e}")
    else:
        print("No structures available. Creating mock data...")
        
        # Create mock structure data
        import pandas as pd
        mock_structure = pd.DataFrame({
            'atom_name': ['CA'] * 10,
            'auth_comp_id': ['ALA', 'GLY', 'SER', 'THR', 'VAL', 
                             'ILE', 'LEU', 'MET', 'PHE', 'TRP'],
            'auth_seq_id': list(range(1, 11)),
            'auth_chain_id': ['A'] * 10,
            'x': np.linspace(0, 30, 10),
            'y': np.linspace(0, 10, 10),
            'z': np.linspace(0, 5, 10),
            'pdb_id': ['MOCK'] * 10
        })
        
        struct_format = StructureFormat(
            coordinates=mock_structure,
            pdb_id="MOCK",
            chains=['A']
        )
        
        validator = FormatValidator()
        is_valid, error = validator.validate_input(struct_format, "structure")
        print(f"  Mock StructureFormat valid: {is_valid}")
    
    print()


def test_sequence_processor_to_model():
    """Test converting SeqProcessor data to models formats."""
    print("=== Testing Sequence Processor → Model Format ===")
    
    # Initialize processor
    seq_proc = SequenceProcessor(name="test_sequences")
    
    # Check for available sequences
    entities = seq_proc.list_entities()
    print(f"Available sequences: {len(entities)}")
    
    if entities:
        test_seq_id = entities[0]
        print(f"\nUsing sequence: {test_seq_id}")
        
        try:
            # Load sequence
            sequence = seq_proc.load_entity(test_seq_id)
            
            if sequence:
                print(f"  Sequence length: {len(sequence)}")
                print(f"  First 50 AA: {sequence[:50]}...")
                
                # Create SequenceFormat
                seq_format = SequenceFormat(
                    sequence=sequence,
                    sequence_id=test_seq_id
                )
                
                # Validate
                validator = FormatValidator()
                is_valid, error = validator.validate_input(seq_format, "sequence")
                print(f"  SequenceFormat valid: {is_valid}")
                
                # Test models adaptation
                adapter = FormatAdapter()
                
                # For ESM-2
                try:
                    esm_input = adapter.adapt_for_model(
                        {"sequence": seq_format}, 
                        "esm2", 
                        {}
                    )
                    print(f"  ESM-2 format: {esm_input}")
                except Exception as e:
                    print(f"  ESM-2 adaptation error: {e}")
                
        except Exception as e:
            print(f"  Error loading sequence: {e}")
    else:
        print("No sequences available. Creating mock data...")
        
        # Create mock sequence
        mock_seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEK"
        seq_format = SequenceFormat(
            sequence=mock_seq,
            sequence_id="mock_protein"
        )
        
        validator = FormatValidator()
        is_valid, error = validator.validate_input(seq_format, "sequence")
        print(f"  Mock SequenceFormat valid: {is_valid}")
    
    print()


def test_grn_processor_to_model():
    """Test converting GRNBaseProcessor data to models formats."""
    print("=== Testing GRN Processor → Model Format ===")
    
    # Initialize processor
    grn_proc = GRNProcessor(name="test_grn")
    
    # Check for GRN tables
    tables = grn_proc.list_entities()
    print(f"Available GRN entities: {len(tables)}")
    
    # Create mock GRN data
    import pandas as pd
    mock_grn = pd.Series({
        '1.50': 'A',
        '2.50': 'L', 
        '3.50': 'V',
        '4.50': 'I',
        '5.50': 'F',
        '6.50': 'W',
        '7.50': 'Y'
    }, name='mock_protein')
    
    grn_format = GRNFormat(
        grn_positions=mock_grn,
        sequence="ALVIFW",
        grn_system="ballesteros_weinstein"
    )
    
    validator = FormatValidator()
    is_valid, error = validator.validate_input(grn_format, "grn")
    print(f"  GRNFormat valid: {is_valid}")
    print(f"  GRN positions: {list(grn_format.grn_positions.index)}")
    print(f"  Conserved positions: {grn_format.get_conserved_positions(0.8)}")
    
    print()


def test_property_processor_to_model():
    """Test converting PropertyProcessor data to models formats."""
    print("=== Testing Property Processor → Model Format ===")
    
    # Initialize processor
    prop_proc = PropertyProcessor(name="test_properties")
    
    # Create mock property data
    mock_properties = {
        "lambda_max": 470,
        "photocycle": "slow",
        "pump_type": "channel",
        "conductance": 85.5
    }
    
    prop_format = PropertyFormat(
        properties=mock_properties,
        property_types={
            "lambda_max": "regression",
            "photocycle": "classification",
            "pump_type": "classification",
            "conductance": "regression"
        },
        metadata={"source": "literature", "confidence": "high"}
    )
    
    validator = FormatValidator()
    is_valid, error = validator.validate_input(prop_format, "property")
    print(f"  PropertyFormat valid: {is_valid}")
    print(f"  Properties: {list(prop_format.properties.keys())}")
    print(f"  DataFrame representation:")
    print(prop_format.to_dataframe())
    
    print()


def test_full_pipeline():
    """Test complete pipeline from processor to models input."""
    print("=== Testing Full Pipeline: Processor → Model Input ===")
    
    # Simulate having all data types for one entity
    entity_name = "test_protein"
    
    # Mock data for each format
    import pandas as pd
    
    # Structure
    structure_df = pd.DataFrame({
        'atom_name': ['CA'] * 5,
        'auth_comp_id': ['ALA', 'GLY', 'SER', 'THR', 'VAL'],
        'auth_seq_id': [1, 2, 3, 4, 5],
        'auth_chain_id': ['A'] * 5,
        'x': [0.0, 3.8, 7.6, 11.4, 15.2],
        'y': [0.0, 0.0, 0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0, 0.0, 0.0],
        'pdb_id': [entity_name] * 5
    })
    
    # Create all formats
    formats = {
        "structure": StructureFormat(
            coordinates=structure_df,
            pdb_id=entity_name,
            chains=['A']
        ),
        "sequence": SequenceFormat(
            sequence="AGSTV",
            sequence_id=entity_name
        ),
        "grn": GRNFormat(
            grn_positions=pd.Series({'1.50': 'A', '2.50': 'G'}, name=entity_name),
            sequence="AGSTV"
        ),
        "embedding": EmbeddingFormat(
            embeddings=np.random.randn(5, 320),  # 5 residues, 320-dim
            embedding_type="per_residue",
            model_name="mock_model"
        )
    }
    
    # Validate all formats
    validator = FormatValidator()
    adapter = FormatAdapter()
    
    print(f"\nEntity: {entity_name}")
    print("Available formats:", list(formats.keys()))
    
    # Check models compatibility
    models = {
        "ESM-2": ["sequence"],
        "Lambda": ["structure", "sequence", "grn"],
        "Property Predictor": ["embedding"]
    }
    
    for model_name, required in models.items():
        print(f"\n{model_name} (requires: {required}):")
        
        # Check if we have all required formats
        missing = set(required) - set(formats.keys())
        if missing:
            print(f"  ✗ Missing formats: {missing}")
        else:
            print(f"  ✓ All formats available")
            
            # Prepare input
            model_input = {fmt: formats[fmt] for fmt in required}
            
            # Validate
            all_valid = True
            for fmt, data in model_input.items():
                is_valid, error = validator.validate_input(data, fmt)
                if not is_valid:
                    print(f"  ✗ {fmt} validation failed: {error}")
                    all_valid = False
            
            if all_valid:
                print(f"  ✓ All inputs validated")
                
                # Try adaptation
                try:
                    if model_name == "Lambda":
                        # Lambda needs special config
                        config = {
                            "preprocessing_params": {
                                "method": "knn",
                                "k": 3,
                                "threshold": 10.0
                            }
                        }
                        adapted = adapter.adapt_for_model(model_input, "lambda", config)
                        print(f"  ✓ Adapted to models format")
                except Exception as e:
                    print(f"  ! Adaptation note: {e}")
    
    print()


def main():
    """Run all integration tests."""
    print("=" * 60)
    print("PROTOS PROCESSOR → MODEL FORMAT INTEGRATION TESTS")
    print("=" * 60)
    print()
    
    test_structure_processor_to_model()
    test_sequence_processor_to_model()
    test_grn_processor_to_model()
    test_property_processor_to_model()
    test_full_pipeline()
    
    print("=" * 60)
    print("Integration tests completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()