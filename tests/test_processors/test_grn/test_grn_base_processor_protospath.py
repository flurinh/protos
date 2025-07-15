"""
Test GRNProcessor with ProtosPaths integration.

Verifies that GRNProcessor:
1. Can be instantiated without path parameters
2. Uses ProtosPaths for all path operations
3. Implements abstract methods correctly
"""

import pytest
import pandas as pd
from pathlib import Path

from protos.processing.grn import GRNProcessor


class TestGRNBaseProcessorProtosPaths:
    """Test GRNProcessor with ProtosPaths."""
    
    def test_zero_configuration_init(self):
        """Test that GRNProcessor can be initialized without any path parameters."""
        # Should work with zero configuration
        processor = GRNProcessor()
        
        assert processor is not None
        assert processor.processor_type == "grn"
        assert processor.data_path.exists()
    
    def test_no_path_parameters_allowed(self):
        """Test that path parameters are not accepted."""
        # The old constructor had path and processor_data_dir parameters
        # These should no longer work
        with pytest.raises(TypeError):
            processor = GRNProcessor(
                path="/some/path"  # This parameter should not exist
            )
    
    def test_path_properties_use_protospath(self):
        """Test that path properties use ProtosPaths."""
        processor = GRNProcessor()
        
        # Check path properties
        assert processor.path_grn_dir.exists()
        assert processor.path_grn_dir.name == "tables"
        
        assert processor.path_ref_dir.exists()
        assert processor.path_ref_dir.name == "ref"
        
        assert processor.path_config_dir.exists()
        assert processor.path_config_dir.name == "configs"
        
        # All paths should be under the GRN data directory
        grn_root = processor.data_path
        assert processor.path_grn_dir.parent == grn_root
        assert processor.path_ref_dir.parent == grn_root
        assert processor.path_config_dir.parent == grn_root
    
    def test_load_entity_implementation(self):
        """Test that load_entity is properly implemented."""
        processor = GRNProcessor()
        
        # Create test data
        test_data = pd.DataFrame({
            '3.50': ['A', 'V', 'L'],
            '7.53': ['F', 'Y', 'W']
        }, index=['protein1', 'protein2', 'protein3'])
        
        processor.data = test_data
        
        # Test loading an entity
        entity = processor.load_entity('protein2')
        assert entity is not None
        assert isinstance(entity, pd.Series)
        assert entity['3.50'] == 'V'
        assert entity['7.53'] == 'Y'
        
        # Test loading non-existent entity
        entity = processor.load_entity('protein99')
        assert entity is None
    
    def test_save_entity_implementation(self):
        """Test that save_entity is properly implemented."""
        processor = GRNProcessor()
        
        # Create test entity
        entity_data = pd.Series({
            '3.50': 'A',
            '7.53': 'F',
            'species': 'Human'
        })
        
        # Save entity
        processor.save_entity('test_protein', entity_data)
        
        # Verify it was saved
        assert processor.data is not None
        assert 'test_protein' in processor.data.index
        assert processor.data.loc['test_protein', '3.50'] == 'A'
        assert processor.data.loc['test_protein', '7.53'] == 'F'
        
        # Verify entity was registered
        assert processor.entity_exists('test_protein')
    
    def test_dataset_operations(self):
        """Test dataset save/load operations."""
        processor = GRNProcessor()
        
        # Create test dataset
        test_data = pd.DataFrame({
            '3.50': ['A', 'V', 'L'],
            '7.53': ['F', 'Y', 'W'],
            'species': ['Human', 'Mouse', 'Rat']
        }, index=['prot1', 'prot2', 'prot3'])
        
        processor.data = test_data
        processor.ids = list(test_data.index)
        processor.grns = ['3.50', '7.53']
        processor.dataset = 'test_grn_dataset'
        
        # Save dataset
        saved_path = processor.save_grn_table()
        assert saved_path is not None
        assert Path(saved_path).exists()
        
        # Create new processor and load dataset
        processor2 = GRNProcessor()
        processor2.load_grn_table('test_grn_dataset')
        
        # Verify data was loaded correctly
        assert len(processor2.data) == 3
        assert set(processor2.ids) == {'prot1', 'prot2', 'prot3'}
        assert '3.50' in processor2.grns
        assert '7.53' in processor2.grns
    
    def test_load_grn_table_return_only(self):
        """Test load_grn_table with return_only parameter."""
        processor = GRNProcessor()
        
        # Create and save test data
        test_data = pd.DataFrame({
            '3.50': ['A', 'V'],
            '7.53': ['F', 'Y']
        }, index=['p1', 'p2'])
        
        processor.data = test_data
        processor.dataset = 'test_return_only'
        processor.save_grn_table()
        
        # Test return_only=True
        processor2 = GRNProcessor()
        loaded_data = processor2.load_grn_table('test_return_only', return_only=True)
        
        # Should return data without setting instance state
        assert loaded_data is not None
        assert len(loaded_data) == 2
        assert processor2.data is None or processor2.data.empty
        assert processor2.dataset is None
        
        # Test return_only=False (default)
        result = processor2.load_grn_table('test_return_only', return_only=False)
        
        # Should set instance state and return None
        assert result is None
        assert processor2.data is not None
        assert len(processor2.data) == 2
        assert processor2.dataset == 'test_return_only'