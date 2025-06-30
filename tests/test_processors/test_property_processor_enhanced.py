#!/usr/bin/env python3
"""
Comprehensive tests for the Enhanced PropertyProcessor

This test suite rigorously tests all aspects of the PropertyProcessor:
1. Entity mapping and association across all processor types
2. Property assignment and retrieval 
3. Dataset management and organization
4. Secondary selection features (GRN positions, atom selections)
5. Cross-format entity integration
6. File I/O and data persistence
7. Error handling and edge cases

Test Design Principles:
- Test entity resolution from all processor types
- Verify property assignment accuracy and consistency
- Test dataset operations (create, load, save, filter)
- Validate advanced selection criteria
- Ensure proper integration with entity registry
- Test batch operations and performance
"""

import pytest
import pandas as pd
import numpy as np
import json
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock

from protos.processing.property.property_processor_enhanced import PropertyProcessor
from protos.io.data_access import generate_entity_id


@pytest.fixture
def temp_data_dir():
    """Create temporary directory for test data."""
    temp_dir = tempfile.mkdtemp()
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)


@pytest.fixture
def property_processor(temp_data_dir):
    """Create PropertyProcessor instance with temporary directory."""
    # Mock ProtosPaths to use temp directory
    with patch('protos.core.base_processor.ProtosPaths') as mock_paths:
        mock_paths_instance = MagicMock()
        mock_paths_instance.data_root = str(temp_data_dir)
        mock_paths_instance.get_processor_path.return_value = str(temp_data_dir / 'property')
        mock_paths.return_value = mock_paths_instance
        
        processor = PropertyProcessor(name="test_property_processor")
        processor.data_path = temp_data_dir / 'property'
        processor.data_path.mkdir(parents=True, exist_ok=True)
        
        # Re-initialize data directories with correct paths
        processor.data_dirs = {
            'datasets': processor.data_path / 'datasets',
            'metadata': processor.data_path / 'metadata', 
            'assignments': processor.data_path / 'assignments',
            'cache': processor.data_path / 'cache'
        }
        
        for dir_path in processor.data_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        yield processor


@pytest.fixture
def sample_entities():
    """Sample entity data for testing."""
    return {
        'structure_entities': {
            '36c2c0da93': {'original_id': '1ubq', 'type': 'structure'},
            'a1b2c3d4e5': {'original_id': '2abc', 'type': 'structure'},
        },
        'sequence_entities': {
            '7e77394211': {'original_id': 'TEST_PROTEIN', 'type': 'sequence'},
            'f1e2d3c4b5': {'original_id': 'OPSIN_001', 'type': 'sequence'},
        },
        'grn_entities': {
            'b3c4d5e6f7': {'original_id': 'GRN_ENTRY_1', 'type': 'grn'},
            'c4d5e6f7g8': {'original_id': 'GRN_ENTRY_2', 'type': 'grn'},
        }
    }


class TestPropertyProcessorBasics:
    """Test basic PropertyProcessor functionality."""
    
    def test_initialization(self, property_processor):
        """Test PropertyProcessor initialization."""
        assert property_processor.name == "test_property_processor"
        assert property_processor.data_path.exists()
        
        # Check all required directories are created
        for dir_name, dir_path in property_processor.data_dirs.items():
            assert dir_path.exists(), f"Directory {dir_name} not created"
        
        # Check initial state
        assert len(property_processor.property_datasets) == 0
        assert len(property_processor.entity_properties) == 0
        assert len(property_processor.property_metadata) == 0
    
    def test_entity_id_resolution(self, property_processor):
        """Test entity identifier resolution."""
        # Test direct entity ID (10 characters)
        entity_id = "36c2c0da93"
        resolved = property_processor._resolve_entity_identifier(entity_id)
        assert resolved == entity_id
        
        # Test identifier that needs to be hashed
        identifier = "1ubq"
        resolved = property_processor._resolve_entity_identifier(identifier)
        expected = generate_entity_id(identifier)
        assert resolved == expected
        assert len(resolved) == 10


class TestPropertyAssignment:
    """Test property assignment functionality."""
    
    def test_basic_property_assignment(self, property_processor):
        """Test basic property assignment to entities."""
        entity_id = "36c2c0da93"
        property_name = "lambda_max"
        property_value = 500.0
        dataset_name = "test_properties"
        
        # Assign property
        result_entity_id = property_processor.assign_property(
            entity_identifier=entity_id,
            property_name=property_name,
            property_value=property_value,
            dataset_name=dataset_name
        )
        
        assert result_entity_id == entity_id
        
        # Verify property was assigned
        retrieved_value = property_processor.get_entity_property(
            entity_id, property_name, dataset_name
        )
        assert retrieved_value == property_value
        
        # Check dataset was created
        assert dataset_name in property_processor.property_datasets
        
        # Check entity properties structure
        assert entity_id in property_processor.entity_properties
        assert dataset_name in property_processor.entity_properties[entity_id]
        assert property_name in property_processor.entity_properties[entity_id][dataset_name]
    
    def test_property_assignment_with_metadata(self, property_processor):
        """Test property assignment with metadata."""
        entity_id = "36c2c0da93"
        property_name = "resolution"
        property_value = 1.8
        dataset_name = "structural_properties"
        metadata = {"method": "X-ray crystallography", "source": "PDB"}
        
        property_processor.assign_property(
            entity_identifier=entity_id,
            property_name=property_name,
            property_value=property_value,
            dataset_name=dataset_name,
            metadata=metadata
        )
        
        # Check metadata is stored
        prop_entry = property_processor.entity_properties[entity_id][dataset_name][property_name]
        assert prop_entry['metadata'] == metadata
        assert 'assigned_at' in prop_entry
        assert prop_entry['value'] == property_value
    
    def test_property_overwrite_behavior(self, property_processor):
        """Test property overwrite behavior."""
        entity_id = "36c2c0da93"
        property_name = "lambda_max"
        dataset_name = "test_properties"
        
        # Initial assignment
        property_processor.assign_property(
            entity_id, property_name, 500.0, dataset_name
        )
        
        # Try to assign again without overwrite
        result_entity_id = property_processor.assign_property(
            entity_id, property_name, 520.0, dataset_name, overwrite=False
        )
        
        # Should keep original value
        value = property_processor.get_entity_property(entity_id, property_name, dataset_name)
        assert value == 500.0
        
        # Now overwrite
        property_processor.assign_property(
            entity_id, property_name, 520.0, dataset_name, overwrite=True
        )
        
        value = property_processor.get_entity_property(entity_id, property_name, dataset_name)
        assert value == 520.0
    
    def test_multiple_properties_same_entity(self, property_processor):
        """Test assigning multiple properties to the same entity."""
        entity_id = "36c2c0da93"
        dataset_name = "opsin_properties"
        
        properties = {
            "lambda_max": 500.0,
            "extinction_coefficient": 40000,
            "protein_type": "rhodopsin",
            "organism": "bovine"
        }
        
        # Assign all properties
        for prop_name, prop_value in properties.items():
            property_processor.assign_property(
                entity_id, prop_name, prop_value, dataset_name
            )
        
        # Retrieve all properties
        all_props = property_processor.get_entity_properties(entity_id, dataset_name)
        
        assert len(all_props) == len(properties)
        for prop_name, expected_value in properties.items():
            assert all_props[prop_name] == expected_value
    
    def test_batch_property_assignment(self, property_processor):
        """Test batch property assignment."""
        dataset_name = "batch_test"
        
        assignments = [
            {
                'entity_identifier': '36c2c0da93',
                'property_name': 'lambda_max',
                'property_value': 500.0,
                'metadata': {'method': 'spectroscopy'}
            },
            {
                'entity_identifier': 'a1b2c3d4e5',
                'property_name': 'lambda_max', 
                'property_value': 480.0,
                'metadata': {'method': 'spectroscopy'}
            },
            {
                'entity_identifier': '36c2c0da93',
                'property_name': 'expression_level',
                'property_value': 'high'
            }
        ]
        
        entity_ids = property_processor.assign_properties_batch(
            assignments, dataset_name
        )
        
        assert len(entity_ids) == 3
        assert '36c2c0da93' in entity_ids
        assert 'a1b2c3d4e5' in entity_ids
        
        # Verify assignments
        assert property_processor.get_entity_property('36c2c0da93', 'lambda_max', dataset_name) == 500.0
        assert property_processor.get_entity_property('a1b2c3d4e5', 'lambda_max', dataset_name) == 480.0
        assert property_processor.get_entity_property('36c2c0da93', 'expression_level', dataset_name) == 'high'


class TestPropertyRetrieval:
    """Test property retrieval functionality."""
    
    def test_get_entity_property_specific_dataset(self, property_processor):
        """Test retrieving specific property from specific dataset."""
        entity_id = "36c2c0da93"
        
        # Set up properties in multiple datasets
        property_processor.assign_property(entity_id, "lambda_max", 500.0, "dataset1")
        property_processor.assign_property(entity_id, "lambda_max", 520.0, "dataset2")
        property_processor.assign_property(entity_id, "temperature", 37.0, "dataset1")
        
        # Test specific dataset retrieval
        assert property_processor.get_entity_property(entity_id, "lambda_max", "dataset1") == 500.0
        assert property_processor.get_entity_property(entity_id, "lambda_max", "dataset2") == 520.0
        assert property_processor.get_entity_property(entity_id, "temperature", "dataset1") == 37.0
        assert property_processor.get_entity_property(entity_id, "temperature", "dataset2") is None
    
    def test_get_entity_property_across_datasets(self, property_processor):
        """Test retrieving property across all datasets."""
        entity_id = "36c2c0da93"
        
        property_processor.assign_property(entity_id, "unique_prop", "value1", "dataset1")
        property_processor.assign_property(entity_id, "another_prop", "value2", "dataset2")
        
        # Search across all datasets
        assert property_processor.get_entity_property(entity_id, "unique_prop") == "value1"
        assert property_processor.get_entity_property(entity_id, "another_prop") == "value2"
        assert property_processor.get_entity_property(entity_id, "nonexistent") is None
    
    def test_get_all_entity_properties(self, property_processor):
        """Test retrieving all properties for an entity."""
        entity_id = "36c2c0da93"
        
        # Properties in different datasets
        property_processor.assign_property(entity_id, "prop1", "value1", "dataset1")
        property_processor.assign_property(entity_id, "prop2", "value2", "dataset1")
        property_processor.assign_property(entity_id, "prop3", "value3", "dataset2")
        
        # Get properties from specific dataset
        dataset1_props = property_processor.get_entity_properties(entity_id, "dataset1")
        assert len(dataset1_props) == 2
        assert dataset1_props["prop1"] == "value1"
        assert dataset1_props["prop2"] == "value2"
        
        # Get properties from all datasets
        all_props = property_processor.get_entity_properties(entity_id)
        assert len(all_props) == 3
        assert "prop1" in all_props
        assert "prop2" in all_props
        assert "prop3" in all_props
    
    def test_property_not_found(self, property_processor):
        """Test behavior when property or entity not found."""
        # Nonexistent entity
        assert property_processor.get_entity_property("nonexistent", "prop") is None
        assert property_processor.get_entity_properties("nonexistent") == {}
        
        # Nonexistent property on existing entity
        entity_id = "36c2c0da93"
        property_processor.assign_property(entity_id, "existing_prop", "value", "dataset")
        
        assert property_processor.get_entity_property(entity_id, "nonexistent_prop") is None


class TestDatasetManagement:
    """Test dataset management functionality."""
    
    def test_dataset_creation_and_statistics(self, property_processor):
        """Test dataset creation and statistics."""
        dataset_name = "test_dataset"
        
        # Assign properties to create dataset
        entities = ["36c2c0da93", "a1b2c3d4e5", "b3c4d5e6f7"]
        for i, entity_id in enumerate(entities):
            property_processor.assign_property(
                entity_id, "lambda_max", 500 + i * 10, dataset_name
            )
            property_processor.assign_property(
                entity_id, "organism", f"organism_{i}", dataset_name
            )
        
        # Check dataset statistics
        stats = property_processor.get_dataset_statistics(dataset_name)
        
        assert stats['dataset_name'] == dataset_name
        assert stats['entity_count'] == 3
        assert stats['property_count'] == 2
        assert 'created_at' in stats
        assert 'modified_at' in stats
        
        # Check property statistics
        assert 'properties' in stats
        assert 'lambda_max' in stats['properties']
        assert 'organism' in stats['properties']
        
        lambda_stats = stats['properties']['lambda_max']
        assert lambda_stats['type'] == 'float64'
        assert lambda_stats['non_null_count'] == 3
        assert lambda_stats['mean'] == 510.0  # (500 + 510 + 520) / 3
    
    def test_get_dataset_as_dataframe(self, property_processor):
        """Test getting dataset as pandas DataFrame."""
        dataset_name = "dataframe_test"
        
        # Create test data
        test_data = {
            "36c2c0da93": {"lambda_max": 500.0, "expression": "high"},
            "a1b2c3d4e5": {"lambda_max": 480.0, "expression": "medium"},
            "b3c4d5e6f7": {"lambda_max": 520.0, "expression": "low"}
        }
        
        for entity_id, properties in test_data.items():
            for prop_name, prop_value in properties.items():
                property_processor.assign_property(
                    entity_id, prop_name, prop_value, dataset_name
                )
        
        # Get as DataFrame
        df = property_processor.get_dataset_properties(dataset_name)
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert len(df.columns) == 2
        assert 'lambda_max' in df.columns
        assert 'expression' in df.columns
        
        # Check index is entity IDs
        assert all(entity_id in df.index for entity_id in test_data.keys())
        
        # Check values
        assert df.loc["36c2c0da93", "lambda_max"] == 500.0
        assert df.loc["a1b2c3d4e5", "expression"] == "medium"
    
    def test_filter_entities_by_property(self, property_processor):
        """Test filtering entities by property values."""
        dataset_name = "filter_test"
        
        # Create test data with various property values
        test_data = [
            ("entity1", {"lambda_max": 500.0, "type": "rhodopsin", "count": 1}),
            ("entity2", {"lambda_max": 480.0, "type": "opsin", "count": 2}),
            ("entity3", {"lambda_max": 520.0, "type": "rhodopsin", "count": 3}),
            ("entity4", {"lambda_max": 460.0, "type": "opsin", "count": 4}),
        ]
        
        for entity_id, properties in test_data:
            for prop_name, prop_value in properties.items():
                property_processor.assign_property(
                    entity_id, prop_name, prop_value, dataset_name
                )
        
        # Test various filters
        
        # Equal filter
        rhodopsins = property_processor.filter_entities_by_property(
            dataset_name, {"type": "rhodopsin"}
        )
        assert len(rhodopsins) == 2
        assert "entity1" in rhodopsins
        assert "entity3" in rhodopsins
        
        # Greater than filter
        high_lambda = property_processor.filter_entities_by_property(
            dataset_name, {"lambda_max": {"gt": 490}}
        )
        assert len(high_lambda) == 2
        assert "entity1" in high_lambda
        assert "entity3" in high_lambda
        
        # Less than filter
        low_lambda = property_processor.filter_entities_by_property(
            dataset_name, {"lambda_max": {"lt": 490}}
        )
        assert len(low_lambda) == 2
        assert "entity2" in low_lambda
        assert "entity4" in low_lambda
        
        # In filter
        specific_types = property_processor.filter_entities_by_property(
            dataset_name, {"type": {"in": ["rhodopsin"]}}
        )
        assert len(specific_types) == 2
        
        # Combined filters
        specific_rhodopsins = property_processor.filter_entities_by_property(
            dataset_name, {"type": "rhodopsin", "lambda_max": {"gt": 510}}
        )
        assert len(specific_rhodopsins) == 1
        assert "entity3" in specific_rhodopsins
    
    def test_list_datasets(self, property_processor):
        """Test listing available datasets."""
        # Initially no datasets
        datasets = property_processor.list_datasets()
        assert len(datasets) == 0
        
        # Create some datasets
        dataset_names = ["dataset1", "dataset2", "dataset3"]
        for ds_name in dataset_names:
            property_processor.assign_property("test_entity", "test_prop", "value", ds_name)
        
        # List datasets
        datasets = property_processor.list_datasets()
        assert len(datasets) == 3
        
        # Check dataset structure
        for dataset_info in datasets:
            assert 'id' in dataset_info
            assert 'entity_count' in dataset_info
            assert 'property_count' in dataset_info
            assert 'type' in dataset_info
            assert dataset_info['type'] == 'property_dataset'
            assert dataset_info['id'] in dataset_names
    
    def test_list_entities(self, property_processor):
        """Test listing entities with properties."""
        # Initially no entities
        entities = property_processor.list_entities()
        assert len(entities) == 0
        
        # Add entities to different datasets
        dataset1_entities = ["entity1", "entity2", "entity3"]
        dataset2_entities = ["entity2", "entity3", "entity4"]  # Some overlap
        
        for entity in dataset1_entities:
            property_processor.assign_property(entity, "prop1", "value1", "dataset1")
        
        for entity in dataset2_entities:
            property_processor.assign_property(entity, "prop2", "value2", "dataset2")
        
        # List entities from specific dataset
        ds1_entities = property_processor.list_entities("dataset1")
        assert len(ds1_entities) == 3
        assert all(entity in ds1_entities for entity in dataset1_entities)
        
        ds2_entities = property_processor.list_entities("dataset2")
        assert len(ds2_entities) == 3
        assert all(entity in ds2_entities for entity in dataset2_entities)
        
        # List all entities
        all_entities = property_processor.list_entities()
        unique_entities = set(dataset1_entities + dataset2_entities)
        assert len(all_entities) == len(unique_entities)
        assert all(entity in all_entities for entity in unique_entities)


class TestFileIO:
    """Test file input/output operations."""
    
    def test_save_and_load_dataset_json(self, property_processor, temp_data_dir):
        """Test saving and loading dataset in JSON format."""
        dataset_name = "save_load_test"
        
        # Create test dataset
        test_data = {
            "entity1": {"lambda_max": 500.0, "type": "rhodopsin"},
            "entity2": {"lambda_max": 480.0, "type": "opsin"},
        }
        
        for entity_id, properties in test_data.items():
            for prop_name, prop_value in properties.items():
                property_processor.assign_property(
                    entity_id, prop_name, prop_value, dataset_name
                )
        
        # Save dataset
        property_processor.save_property_dataset(dataset_name, 'json')
        
        # Check file was created
        json_file = property_processor.data_dirs['datasets'] / dataset_name / f"{dataset_name}.json"
        assert json_file.exists()
        
        # Create new processor instance and load dataset
        new_processor = PropertyProcessor(name="new_processor")
        new_processor.data_dirs = property_processor.data_dirs
        
        success = new_processor.load_property_dataset(dataset_name)
        assert success
        
        # Verify loaded data
        assert dataset_name in new_processor.property_datasets
        for entity_id in test_data.keys():
            for prop_name, expected_value in test_data[entity_id].items():
                loaded_value = new_processor.get_entity_property(entity_id, prop_name, dataset_name)
                assert loaded_value == expected_value
    
    def test_save_dataset_csv(self, property_processor):
        """Test saving dataset in CSV format."""
        dataset_name = "csv_test"
        
        # Create test dataset
        for i in range(3):
            entity_id = f"entity_{i}"
            property_processor.assign_property(entity_id, "lambda_max", 500 + i * 10, dataset_name)
            property_processor.assign_property(entity_id, "organism", f"organism_{i}", dataset_name)
        
        # Save as CSV
        property_processor.save_property_dataset(dataset_name, 'csv')
        
        # Check file was created
        csv_file = property_processor.data_dirs['datasets'] / dataset_name / f"{dataset_name}.csv"
        assert csv_file.exists()
        
        # Load and verify CSV content
        df = pd.read_csv(csv_file, index_col='entity_id')
        assert len(df) == 3
        assert 'lambda_max' in df.columns
        assert 'organism' in df.columns
        assert df.loc['entity_0', 'lambda_max'] == 500.0
    
    def test_create_dataset_from_csv_file(self, property_processor, temp_data_dir):
        """Test creating dataset from CSV file."""
        # Create test CSV file
        test_data = pd.DataFrame({
            'entity_id': ['36c2c0da93', 'a1b2c3d4e5', 'b3c4d5e6f7'],
            'lambda_max': [500.0, 480.0, 520.0],
            'organism': ['bovine', 'human', 'mouse'],
            'type': ['rhodopsin', 'opsin', 'rhodopsin']
        })
        
        csv_file = temp_data_dir / 'test_properties.csv'
        test_data.to_csv(csv_file, index=False)
        
        # Create dataset from file
        dataset_name = "from_csv_test"
        count = property_processor.create_property_dataset_from_file(
            str(csv_file), dataset_name, 'entity_id'
        )
        
        assert count == 3
        assert dataset_name in property_processor.property_datasets
        
        # Verify data was loaded correctly
        for _, row in test_data.iterrows():
            entity_id = row['entity_id']
            for col in ['lambda_max', 'organism', 'type']:
                value = property_processor.get_entity_property(entity_id, col, dataset_name)
                assert value == row[col]
    
    def test_create_dataset_from_json_file(self, property_processor, temp_data_dir):
        """Test creating dataset from JSON file."""
        # Create test JSON file
        test_data = [
            {'entity_id': '36c2c0da93', 'lambda_max': 500.0, 'type': 'rhodopsin'},
            {'entity_id': 'a1b2c3d4e5', 'lambda_max': 480.0, 'type': 'opsin'},
        ]
        
        json_file = temp_data_dir / 'test_properties.json'
        with open(json_file, 'w') as f:
            json.dump(test_data, f)
        
        # Create dataset from file
        dataset_name = "from_json_test"
        count = property_processor.create_property_dataset_from_file(
            str(json_file), dataset_name, 'entity_id'
        )
        
        assert count == 2
        
        # Verify data
        for item in test_data:
            entity_id = item['entity_id']
            for key, value in item.items():
                if key != 'entity_id':
                    loaded_value = property_processor.get_entity_property(entity_id, key, dataset_name)
                    assert loaded_value == value


class TestSecondarySelection:
    """Test advanced secondary selection features."""
    
    def test_grn_position_assignment(self, property_processor):
        """Test assigning properties to specific GRN positions."""
        entity_id = "36c2c0da93"
        dataset_name = "grn_properties"
        
        # Assign properties to specific GRN positions
        grn_properties = {
            "3.50": {"polarity": "hydrophobic", "importance": 0.9},
            "7.45": {"polarity": "polar", "importance": 0.8},
        }
        
        for grn_pos, properties in grn_properties.items():
            for prop_name, prop_value in properties.items():
                property_processor.assign_grn_property(
                    entity_id, grn_pos, prop_name, prop_value, dataset_name
                )
        
        # Verify assignments
        polarity_3_50 = property_processor.get_entity_property(
            entity_id, "3.50.polarity", dataset_name
        )
        assert polarity_3_50 == "hydrophobic"
        
        importance_7_45 = property_processor.get_entity_property(
            entity_id, "7.45.importance", dataset_name
        )
        assert importance_7_45 == 0.8
        
        # Check secondary selection was stored
        assert entity_id in property_processor.secondary_selections
        selections = property_processor.secondary_selections[entity_id]
        assert any(sel['type'] == 'grn_position' and sel['grn_position'] == '3.50' for sel in selections)
    
    def test_atom_selection_assignment(self, property_processor):
        """Test assigning properties to specific atoms."""
        entity_id = "36c2c0da93"
        dataset_name = "atom_properties"
        
        # Define atom selections
        atom_selectors = [
            {"chain": "A", "residue": 100},
            {"chain": "B", "atom_name": "CA"},
        ]
        
        for i, selector in enumerate(atom_selectors):
            property_processor.assign_atom_property(
                entity_id, selector, "b_factor", 20.0 + i * 5, dataset_name
            )
        
        # Check properties were assigned with correct names
        all_props = property_processor.get_entity_properties(entity_id, dataset_name)
        
        # Should have properties named like "chainAresidue100.b_factor"
        prop_names = list(all_props.keys())
        assert any("chainA" in name and "b_factor" in name for name in prop_names)
        assert any("chainB" in name and "b_factor" in name for name in prop_names)
        
        # Check secondary selections
        assert entity_id in property_processor.secondary_selections
        selections = property_processor.secondary_selections[entity_id]
        assert any(sel['type'] == 'atom_selection' for sel in selections)
    
    def test_secondary_selection_metadata(self, property_processor):
        """Test that secondary selection metadata is properly stored."""
        entity_id = "test_entity"
        dataset_name = "test_dataset"
        
        # Assign property with secondary selection
        secondary_selection = {
            'type': 'custom_selection',
            'criteria': {'custom_field': 'custom_value'}
        }
        
        property_processor.assign_property(
            entity_id, "test_prop", "test_value", dataset_name,
            secondary_selection=secondary_selection
        )
        
        # Check that secondary selection is stored in property entry
        prop_entry = property_processor.entity_properties[entity_id][dataset_name]["test_prop"]
        assert 'secondary_selection' in prop_entry
        assert prop_entry['secondary_selection'] == secondary_selection
        
        # Check it's also stored in secondary_selections
        assert entity_id in property_processor.secondary_selections
        stored_selections = property_processor.secondary_selections[entity_id]
        assert secondary_selection in stored_selections


class TestEntityIntegration:
    """Test integration with entity registry and cross-format support."""
    
    @patch('protos.processing.property.property_processor_enhanced.GlobalRegistry')
    def test_entity_registry_resolution(self, mock_global_registry, property_processor):
        """Test entity resolution through entity registry."""
        # Mock entity registry
        mock_registry = MagicMock()
        mock_registry.resolve_identifier.return_value = "36c2c0da93"
        mock_global_registry.return_value.entity_registry = mock_registry
        
        # Test resolution
        resolved_id = property_processor._resolve_entity_identifier("1ubq", "structure")
        assert resolved_id == "36c2c0da93"
        
        # Verify registry was called
        mock_registry.resolve_identifier.assert_called_with("1ubq", format_type="structure")
    
    @patch('protos.processing.property.property_processor_enhanced.GlobalRegistry')
    def test_entity_format_discovery(self, mock_global_registry, property_processor):
        """Test discovering available formats for an entity."""
        # Mock entity data
        mock_entity_data = {
            'formats': {
                'structure': {'file_path': '/path/to/structure'},
                'sequence': {'file_path': '/path/to/sequence'},
                'grn': {'metadata': {'positions': 5}}
            }
        }
        
        mock_registry = MagicMock()
        mock_registry.get_entity.return_value = mock_entity_data
        mock_global_registry.return_value.entity_registry = mock_registry
        
        # Test format discovery
        formats = property_processor.get_entity_formats("36c2c0da93")
        
        assert len(formats) == 3
        assert 'structure' in formats
        assert 'sequence' in formats
        assert 'grn' in formats
    
    def test_cross_format_property_assignment(self, property_processor, sample_entities):
        """Test assigning properties to entities across different formats."""
        dataset_name = "cross_format_test"
        
        # Assign properties to entities of different types
        entity_properties = [
            ("36c2c0da93", "structure", "resolution", 1.8),
            ("7e77394211", "sequence", "length", 280),
            ("b3c4d5e6f7", "grn", "positions", 7),
        ]
        
        for entity_id, entity_type, prop_name, prop_value in entity_properties:
            property_processor.assign_property(
                entity_id, prop_name, prop_value, dataset_name,
                metadata={'entity_type': entity_type}
            )
        
        # Verify all properties were assigned
        df = property_processor.get_dataset_properties(dataset_name)
        assert len(df) == 3
        
        # Check that entities from different formats can coexist
        all_entities = property_processor.list_entities(dataset_name)
        assert "36c2c0da93" in all_entities  # structure
        assert "7e77394211" in all_entities  # sequence  
        assert "b3c4d5e6f7" in all_entities  # grn


class TestErrorHandling:
    """Test error handling and edge cases."""
    
    def test_invalid_dataset_operations(self, property_processor):
        """Test operations on invalid or nonexistent datasets."""
        # Get statistics for nonexistent dataset
        stats = property_processor.get_dataset_statistics("nonexistent_dataset")
        assert stats == {}
        
        # Get DataFrame for nonexistent dataset
        df = property_processor.get_dataset_properties("nonexistent_dataset")
        assert df.empty
        
        # Filter nonexistent dataset
        entities = property_processor.filter_entities_by_property(
            "nonexistent_dataset", {"prop": "value"}
        )
        assert entities == []
        
        # Load nonexistent dataset
        success = property_processor.load_property_dataset("nonexistent_dataset")
        assert not success
    
    def test_invalid_file_operations(self, property_processor, temp_data_dir):
        """Test file operations with invalid files."""
        # Try to create dataset from nonexistent file
        with pytest.raises(FileNotFoundError):
            property_processor.create_property_dataset_from_file(
                "/nonexistent/file.csv", "test_dataset", "entity_id"
            )
        
        # Try to create dataset from file with missing entity column
        test_csv = temp_data_dir / "no_entity_column.csv"
        pd.DataFrame({"prop1": [1, 2, 3], "prop2": ["a", "b", "c"]}).to_csv(test_csv, index=False)
        
        with pytest.raises(ValueError, match="Entity column.*not found"):
            property_processor.create_property_dataset_from_file(
                str(test_csv), "test_dataset", "entity_id"
            )
        
        # Try to create dataset from unsupported file format
        unsupported_file = temp_data_dir / "data.xml"
        unsupported_file.write_text("<data></data>")
        
        with pytest.raises(ValueError, match="Unsupported file format"):
            property_processor.create_property_dataset_from_file(
                str(unsupported_file), "test_dataset", "entity_id"
            )
    
    def test_property_filtering_edge_cases(self, property_processor):
        """Test property filtering with edge cases."""
        dataset_name = "filter_edge_cases"
        
        # Create dataset with NaN values and mixed types
        property_processor.assign_property("entity1", "numeric_prop", 100, dataset_name)
        property_processor.assign_property("entity2", "numeric_prop", np.nan, dataset_name)
        property_processor.assign_property("entity1", "string_prop", "test", dataset_name)
        property_processor.assign_property("entity2", "string_prop", None, dataset_name)
        
        # Filter by nonexistent property (should return empty list)
        result = property_processor.filter_entities_by_property(
            dataset_name, {"nonexistent_prop": "value"}
        )
        assert result == []
        
        # Filter with valid property should work
        result = property_processor.filter_entities_by_property(
            dataset_name, {"numeric_prop": 100}
        )
        assert len(result) == 1
        assert "entity1" in result


class TestPerformanceAndScaling:
    """Test performance with larger datasets."""
    
    def test_large_dataset_operations(self, property_processor):
        """Test operations with larger number of entities and properties."""
        dataset_name = "large_dataset"
        num_entities = 1000
        num_properties = 10
        
        # Create large dataset
        assignments = []
        for i in range(num_entities):
            entity_id = f"entity_{i:04d}"
            for j in range(num_properties):
                assignments.append({
                    'entity_identifier': entity_id,
                    'property_name': f'prop_{j}',
                    'property_value': np.random.rand(),
                })
        
        # Batch assign properties
        entity_ids = property_processor.assign_properties_batch(assignments, dataset_name)
        
        assert len(entity_ids) == num_entities * num_properties
        
        # Test retrieval performance
        df = property_processor.get_dataset_properties(dataset_name)
        assert len(df) == num_entities
        assert len(df.columns) == num_properties
        
        # Test filtering performance
        entities = property_processor.filter_entities_by_property(
            dataset_name, {"prop_0": {"gt": 0.5}}
        )
        assert len(entities) > 0  # Should find some entities
        
        # Test statistics calculation
        stats = property_processor.get_dataset_statistics(dataset_name)
        assert stats['entity_count'] == num_entities
        assert stats['property_count'] == num_properties
    
    def test_multiple_datasets_performance(self, property_processor):
        """Test performance with multiple datasets."""
        num_datasets = 10
        entities_per_dataset = 50
        
        # Create multiple datasets
        for i in range(num_datasets):
            dataset_name = f"dataset_{i}"
            for j in range(entities_per_dataset):
                entity_id = f"entity_{i}_{j}"
                property_processor.assign_property(
                    entity_id, "value", i * entities_per_dataset + j, dataset_name
                )
        
        # Test cross-dataset operations
        all_datasets = property_processor.list_datasets()
        assert len(all_datasets) == num_datasets
        
        all_entities = property_processor.list_entities()
        assert len(all_entities) == num_datasets * entities_per_dataset
        
        # Test entity with properties in multiple datasets
        shared_entity = "shared_entity"
        for i in range(num_datasets):
            dataset_name = f"dataset_{i}"
            property_processor.assign_property(
                shared_entity, f"prop_{i}", f"value_{i}", dataset_name
            )
        
        all_props = property_processor.get_entity_properties(shared_entity)
        assert len(all_props) == num_datasets


if __name__ == "__main__":
    pytest.main([__file__, "-v"])