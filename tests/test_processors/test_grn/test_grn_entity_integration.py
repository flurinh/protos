"""
Tests for GRNBaseProcessor entity integration using real data.
"""

import pytest
import os
import shutil
import pandas as pd
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.io.data_access import GlobalRegistry, generate_entity_id


@pytest.fixture
def setup_test_environment(request):
    """Use the test-data directory directly for GRN tests."""
    # Set the global data root to our test-data directory
    test_data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/tests/test-data")
    ProtosPaths.set_data_root(str(test_data_dir))
    
    def teardown():
        # Clear the data root after test
        ProtosPaths.set_data_root(None)
    
    request.addfinalizer(teardown)
    return test_data_dir


class TestGRNEntityIntegration:
    """Test entity integration in GRNBaseProcessor."""
    
    def test_save_grn_table_with_entity_ids(self, setup_test_environment):
        """Test that saving a GRN table includes entity IDs."""
        # Create processor
        processor = GRNBaseProcessor(name="test_grn")
        
        # Load a real GRN table from ref directory
        processor.load_grn_table("ref/mo_ref")
        assert processor.data is not None
        assert len(processor.ids) > 0
        
        # Save with entity IDs
        saved_path = processor.save_grn_table("test_with_entities", include_entity_ids=True)
        assert os.path.exists(saved_path)
        
        # Load the saved table and verify entity_id column
        saved_df = pd.read_csv(saved_path)
        assert 'entity_id' in saved_df.columns
        assert saved_df.columns[0] == 'entity_id'  # Should be first column
        
        # Verify each entity was registered
        global_registry = GlobalRegistry()
        for idx, row in saved_df.iterrows():
            if 'protein_id' in saved_df.columns:
                seq_id = row['protein_id']
            else:
                # Handle case where first non-entity_id column is the ID
                seq_id = saved_df.iloc[idx, 1]  # Second column
            
            entity_id = row['entity_id']
            
            # Check entity in registry
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            assert entity_info is not None
            assert entity_info['original_id'] == str(seq_id)
            assert 'grn' in entity_info['formats']
    
    def test_load_grn_table_registers_entities(self, setup_test_environment):
        """Test that loading a GRN table registers all entities."""
        # Create processor
        processor = GRNBaseProcessor(name="test_load")
        
        # Get initial entity count
        global_registry = GlobalRegistry()
        initial_count = len(global_registry.entity_registry.list_entities(format_type="grn"))
        
        # Load table with entity registration
        processor.load_grn_table("mo_ref", register_entities=True)
        
        # Check entities were registered
        final_count = len(global_registry.entity_registry.list_entities(format_type="grn"))
        assert final_count > initial_count
        assert final_count - initial_count == len(processor.ids)
        
        # Verify each sequence has an entity
        for seq_id in processor.ids:
            entity_id = generate_entity_id(str(seq_id))
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            assert entity_info is not None
            assert entity_info['original_id'] == str(seq_id)
    
    def test_load_grn_entity(self, setup_test_environment):
        """Test loading a single GRN entity."""
        # Create processor and load table
        processor = GRNBaseProcessor(name="test_entity")
        processor.load_grn_table("ref/mo_ref")
        
        # Get a test sequence ID
        test_seq_id = processor.ids[0]
        
        # Load entity by sequence ID
        entity_data = processor.load_grn_entity(test_seq_id)
        assert entity_data is not None
        assert isinstance(entity_data, pd.Series)
        assert entity_data.name == test_seq_id
        
        # Load by entity hash ID
        entity_id = generate_entity_id(str(test_seq_id))
        entity_data2 = processor.load_grn_entity(entity_id)
        assert entity_data2 is not None
        assert entity_data.equals(entity_data2)
    
    def test_list_grn_entities(self, setup_test_environment):
        """Test listing GRN entities returns sequence IDs not hashes."""
        # Create processor and load table
        processor = GRNBaseProcessor(name="test_list")
        processor.load_grn_table("ref/mo_ref")
        
        # List entities
        entities = processor.list_grn_entities()
        
        # Should return sequence IDs, not hashes
        assert len(entities) > 0
        for entity in entities:
            # Hash IDs are 10 alphanumeric characters
            assert not (len(entity) == 10 and entity.isalnum())
            # Should match loaded IDs
            assert entity in processor.ids
    
    def test_get_entity_grn_positions(self, setup_test_environment):
        """Test getting GRN positions for an entity."""
        # Create processor and load table
        processor = GRNBaseProcessor(name="test_positions")
        processor.load_grn_table("ref/mo_ref")
        
        # Get a test sequence with some GRN positions
        test_seq_id = processor.ids[0]
        
        # Get GRN positions
        positions = processor.get_entity_grn_positions(test_seq_id)
        assert isinstance(positions, list)
        assert len(positions) > 0
        
        # Verify these are actual GRN positions (not metadata columns)
        for pos in positions:
            assert pos in processor.grns
            # Should not be metadata columns
            assert pos not in ['entity_id', 'family', 'species', 'name', 'grn_system', 'id']
            # Should have a residue (not '-')
            assert processor.data.loc[test_seq_id, pos] != '-'
    
    def test_grn_table_multiple_entities(self, setup_test_environment):
        """Test that GRN tables correctly handle multiple entities."""
        # Create processor
        processor = GRNBaseProcessor(name="test_multi")
        
        # Create a test GRN table with multiple sequences
        test_data = pd.DataFrame({
            '1.50': ['A', 'V', 'L'],
            '2.50': ['G', 'G', 'P'],
            '3.50': ['I', 'V', 'F']
        }, index=['SEQ1', 'SEQ2', 'SEQ3'])
        
        processor.data = test_data
        processor.ids = test_data.index.tolist()
        processor.grns = test_data.columns.tolist()
        
        # Save with entity IDs
        saved_path = processor.save_grn_table("test_multi_entities")
        
        # Verify each row got its own entity ID
        saved_df = pd.read_csv(saved_path)
        assert 'entity_id' in saved_df.columns
        
        # Check all entity IDs are unique
        entity_ids = saved_df['entity_id'].tolist()
        assert len(entity_ids) == len(set(entity_ids))
        
        # Verify registry has all entities
        global_registry = GlobalRegistry()
        for seq_id in ['SEQ1', 'SEQ2', 'SEQ3']:
            entity_id = generate_entity_id(seq_id)
            assert global_registry.entity_registry.entity_exists(entity_id)