"""
Test PropertyProcessor with ProtosPaths integration.

Verifies that PropertyProcessor:
1. Can be instantiated without path parameters
2. Uses ProtosPaths for all path operations
3. Implements abstract methods correctly
"""

import pytest
import pandas as pd
from pathlib import Path

from protos.processing.property import PropertyProcessor


class TestPropertyProcessorProtosPaths:
    """Test PropertyProcessor with ProtosPaths."""
    
    def test_zero_configuration_init(self):
        """Test that PropertyProcessor can be initialized without any path parameters."""
        # Should work with zero configuration
        processor = PropertyProcessor()
        
        assert processor is not None
        assert processor.processor_type == "property"
        assert processor.data_path.exists()
    
    def test_no_path_parameters_allowed(self):
        """Test that path parameters are not accepted."""
        # The old constructor had processor_data_dir parameter
        # This should no longer work
        with pytest.raises(TypeError):
            processor = PropertyProcessor(
                name="pp"
            )
    
    def test_path_properties_use_protospath(self):
        """Test that path properties use ProtosPaths."""
        processor = PropertyProcessor()
        
        # Check that subdirectories are created by ProtosPaths
        datasets_dir = processor.get_subdirectory_path('datasets_dir')
        assert datasets_dir.exists()
        assert datasets_dir.name == "datasets"
        
        tables_dir = processor.get_subdirectory_path('tables_dir')
        assert tables_dir.exists()
        assert tables_dir.name == "tables"
        
        # Both should be under the processor's data path
        assert datasets_dir.parent == processor.data_path
        assert tables_dir.parent == processor.data_path
    
    def test_load_entity_implementation(self):
        """Test that load_entity is properly implemented."""
        processor = PropertyProcessor()
        
        # PropertyProcessor returns None for load_entity (entities are rows in tables)
        # This is expected behavior as documented in the abstract method
        props = processor.load_entity('unknown_entity')
        assert props is None
    
    def test_save_entity_implementation(self):
        """Test that save_entity is properly implemented."""
        processor = PropertyProcessor()
        
        # PropertyProcessor manages entities as rows in tables
        # save_entity just logs a warning - this is expected behavior
        test_data = {
            'lambda_max': 500,
            'expression_level': 'high',
            'temperature': 37.5
        }
        
        # This should work but log a warning
        processor.save_entity('opsin_123', test_data, metadata={'dataset': 'test_opsins'})
    
    def test_property_assignment(self):
        """Test property assignment functionality."""
        processor = PropertyProcessor()
        
        # Assign a single property
        processor.assign_property(
            entity_identifier='protein_abc',
            property_name='molecular_weight',
            property_value=25000.5,
            dataset_name='protein_properties'
        )
        
        # Get the property back using table API
        df = processor.get_property_table('protein_properties')
        assert df.loc['protein_abc', 'molecular_weight'] == 25000.5
        
        # Get all properties for entity
        all_props = processor.get_entity_properties('protein_abc')
        assert all_props is not None
        assert all_props['molecular_weight'] == 25000.5
    
    def test_dataset_operations(self):
        """Test property dataset operations."""
        processor = PropertyProcessor()
        
        # Create test data using table-based API
        test_props = {
            'entity1': {'prop_a': 1, 'prop_b': 'x'},
            'entity2': {'prop_a': 2, 'prop_b': 'y'},
            'entity3': {'prop_a': 3, 'prop_b': 'z'}
        }
        
        # Create property table
        df = processor.create_property_table('test_dataset', test_props)
        assert len(df) == 3
        
        # Create new processor and load dataset
        processor2 = PropertyProcessor()
        loaded_df = processor2.load_property_dataset('test_dataset')
        assert isinstance(loaded_df, pd.DataFrame)
        assert len(loaded_df) == 3
        
        # List datasets
        datasets = processor2.list_datasets()
        assert 'test_dataset' in datasets