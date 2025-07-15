#!/usr/bin/env python3
"""
Tests for the Enhanced PropertyProcessor

This test suite tests the PropertyProcessor functionality:
1. Table-based property management
2. Dataset operations
3. Property assignment and retrieval
4. Cross-format entity support
5. Data persistence

Test Design Principles:
- NO tmp_path or temporary directories
- NO mocking of ProtosPaths
- Use real test-data directory via conftest.py
- Test with real biological data
- Follow ProtosPaths principles throughout
"""

import pytest
import pandas as pd
import numpy as np
import json
from pathlib import Path

from protos.processing.property import PropertyProcessor


class TestPropertyProcessorBasics:
    """Test basic PropertyProcessor functionality."""
    
    def test_initialization(self):
        """Test PropertyProcessor initialization with ProtosPaths."""
        # Create processor - ProtosPaths handles everything
        processor = PropertyProcessor(name="test_properties")
        
        assert processor.name == "test_properties"
        assert processor.processor_type == "property"
        assert hasattr(processor, 'property_tables')
        assert isinstance(processor.property_tables, dict)
        
        # Verify ProtosPaths created the directories
        tables_dir = processor.get_subdirectory_path('tables_dir')
        datasets_dir = processor.get_subdirectory_path('datasets_dir')
        
        assert tables_dir.exists()
        assert datasets_dir.exists()
    
    def test_create_property_table(self):
        """Test creating a property table."""
        processor = PropertyProcessor(name="test_create_table")
        
        # Create test data
        data = {
            'protein_1': {'lambda_max': 500, 'expression': 'high', 'organism': 'E.coli'},
            'protein_2': {'lambda_max': 520, 'expression': 'medium', 'organism': 'Human'},
            'protein_3': {'lambda_max': 480, 'expression': 'low', 'organism': 'Mouse'}
        }
        
        # Create table
        df = processor.create_property_table('opsin_properties', data)
        
        # Verify table structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert len(df.columns) == 3
        assert df.index.name == 'entity_id'
        assert list(df.columns) == ['lambda_max', 'expression', 'organism']
        
        # Verify data
        assert df.loc['protein_1', 'lambda_max'] == 500
        assert df.loc['protein_2', 'organism'] == 'Human'
        
        # Verify table was saved
        csv_path = processor.get_subdirectory_path('tables_dir') / 'opsin_properties.csv'
        assert csv_path.exists()
    
    def test_load_property_table(self):
        """Test loading a property table."""
        processor = PropertyProcessor(name="test_load_table")
        
        # First create a table
        data = {
            'entity_A': {'prop1': 10, 'prop2': 'value_A'},
            'entity_B': {'prop1': 20, 'prop2': 'value_B'}
        }
        processor.create_property_table('test_dataset', data)
        
        # Clear memory cache to test loading from disk
        processor.property_tables.clear()
        
        # Load table
        loaded_df = processor.load_property_table('test_dataset')
        
        assert isinstance(loaded_df, pd.DataFrame)
        assert len(loaded_df) == 2
        assert loaded_df.loc['entity_A', 'prop1'] == 10
        assert loaded_df.loc['entity_B', 'prop2'] == 'value_B'


class TestPropertyOperations:
    """Test property manipulation operations."""
    
    def test_add_property_column(self):
        """Test adding a new property column."""
        processor = PropertyProcessor(name="test_add_column")
        
        # Create initial table
        data = {
            'protein_1': {'lambda_max': 500},
            'protein_2': {'lambda_max': 520}
        }
        processor.create_property_table('test_props', data)
        
        # Add new property column with scalar value
        processor.add_property_column('test_props', 'validated', True)
        
        df = processor.get_property_table('test_props')
        assert 'validated' in df.columns
        assert df['validated'].all()  # All True
        
        # Add column with dict values
        expression_levels = {'protein_1': 'high', 'protein_2': 'low'}
        processor.add_property_column('test_props', 'expression', expression_levels)
        
        df = processor.get_property_table('test_props')
        assert df.loc['protein_1', 'expression'] == 'high'
        assert df.loc['protein_2', 'expression'] == 'low'
    
    def test_filter_by_property(self):
        """Test filtering entities by property value."""
        processor = PropertyProcessor(name="test_filter")
        
        # Create test data
        data = {
            f'protein_{i}': {
                'lambda_max': 450 + i * 10,
                'expression': ['low', 'medium', 'high'][i % 3]
            }
            for i in range(10)
        }
        processor.create_property_table('filter_test', data)
        
        # Filter by lambda_max > 500
        filtered = processor.filter_by_property(
            'filter_test', 
            'lambda_max', 
            lambda x: x > 500
        )
        
        assert len(filtered) == 4  # proteins 6,7,8,9
        assert all(filtered['lambda_max'] > 500)
        
        # Filter by expression level
        high_expr = processor.filter_by_property(
            'filter_test',
            'expression',
            lambda x: x == 'high'
        )
        
        assert len(high_expr) == 3  # proteins 2,5,8
        assert all(high_expr['expression'] == 'high')
    
    def test_merge_property_tables(self):
        """Test merging multiple property tables."""
        processor = PropertyProcessor(name="test_merge")
        
        # Create first table
        data1 = {
            'protein_A': {'lambda_max': 500, 'method': 'spectroscopy'},
            'protein_B': {'lambda_max': 520, 'method': 'spectroscopy'}
        }
        processor.create_property_table('spectral_data', data1)
        
        # Create second table with some overlap
        data2 = {
            'protein_A': {'expression': 'high', 'organism': 'E.coli'},
            'protein_C': {'expression': 'low', 'organism': 'Human'}
        }
        processor.create_property_table('expression_data', data2)
        
        # Merge tables
        merged = processor.merge_property_tables(
            ['spectral_data', 'expression_data'],
            'combined_data',
            how='outer'
        )
        
        # Verify merge
        assert len(merged) == 3  # A, B, C
        assert len(merged.columns) == 4  # lambda_max, method, expression, organism
        
        # Check specific values
        assert merged.loc['protein_A', 'lambda_max'] == 500
        assert merged.loc['protein_A', 'expression'] == 'high'
        assert pd.isna(merged.loc['protein_B', 'expression'])
        assert merged.loc['protein_C', 'organism'] == 'Human'


class TestEntityProperties:
    """Test entity-centric property operations."""
    
    def test_get_entity_properties(self):
        """Test retrieving all properties for an entity."""
        processor = PropertyProcessor(name="test_entity_props")
        
        # Create multiple datasets with same entity
        data1 = {'protein_X': {'lambda_max': 480, 'quantum_yield': 0.65}}
        processor.create_property_table('optical_props', data1)
        
        data2 = {'protein_X': {'expression': 'high', 'stability': 'stable'}}
        processor.create_property_table('biochem_props', data2)
        
        # Get all properties for entity
        all_props = processor.get_entity_properties('protein_X')
        
        assert len(all_props) == 4
        assert all_props['lambda_max'] == 480
        assert all_props['quantum_yield'] == 0.65
        assert all_props['expression'] == 'high'
        assert all_props['stability'] == 'stable'
        
        # Get properties from specific dataset
        optical_only = processor.get_entity_properties('protein_X', 'optical_props')
        assert len(optical_only) == 2
        assert 'expression' not in optical_only
    
    def test_list_datasets(self):
        """Test listing available datasets."""
        processor = PropertyProcessor(name="test_list_datasets")
        
        # Create several datasets
        processor.create_property_table('dataset1', {'e1': {'p1': 1}})
        processor.create_property_table('dataset2', {'e2': {'p2': 2}})
        processor.create_property_table('dataset3', {'e3': {'p3': 3}})
        
        datasets = processor.list_datasets()
        
        assert len(datasets) >= 3
        assert 'dataset1' in datasets
        assert 'dataset2' in datasets
        assert 'dataset3' in datasets


class TestCompatibilityAPI:
    """Test compatibility with old PropertyProcessor API."""
    
    def test_assign_property_compatibility(self):
        """Test single property assignment (compatibility method)."""
        processor = PropertyProcessor(name="test_compat")
        
        # Use old API to assign properties
        processor.assign_property('protein_1', 'lambda_max', 500, 'compat_dataset')
        processor.assign_property('protein_1', 'expression', 'high', 'compat_dataset')
        processor.assign_property('protein_2', 'lambda_max', 520, 'compat_dataset')
        
        # Verify using new API
        df = processor.get_property_table('compat_dataset')
        
        assert len(df) == 2
        assert df.loc['protein_1', 'lambda_max'] == 500
        assert df.loc['protein_1', 'expression'] == 'high'
        assert df.loc['protein_2', 'lambda_max'] == 520
    
    def test_load_property_dataset_alias(self):
        """Test dataset loading alias."""
        processor = PropertyProcessor(name="test_alias")
        
        # Create dataset
        data = {'entity1': {'prop1': 'value1'}}
        processor.create_property_table('alias_test', data)
        
        # Use alias method
        df = processor.load_property_dataset('alias_test')
        
        assert isinstance(df, pd.DataFrame)
        assert df.loc['entity1', 'prop1'] == 'value1'


class TestDataPersistence:
    """Test data persistence and file operations."""
    
    def test_save_with_metadata(self):
        """Test saving property table with metadata."""
        processor = PropertyProcessor(name="test_metadata")
        
        # Create table with metadata
        data = {
            'opsin_1': {'lambda_max': 568, 'type': 'proton_pump'},
            'opsin_2': {'lambda_max': 470, 'type': 'channel'}
        }
        
        metadata = {
            'source': 'literature_review',
            'date': '2024-01-15',
            'curator': 'test_user',
            'notes': 'Microbial opsin properties'
        }
        
        processor.create_property_table('opsin_metadata', data, metadata)
        
        # Check metadata file was created
        metadata_path = processor.get_subdirectory_path('datasets_dir') / 'opsin_metadata.json'
        assert metadata_path.exists()
        
        # Load and verify metadata
        with open(metadata_path, 'r') as f:
            saved_data = json.load(f)
        
        assert saved_data['name'] == 'opsin_metadata'
        assert saved_data['metadata']['source'] == 'literature_review'
        assert saved_data['metadata']['curator'] == 'test_user'
        assert saved_data['shape'] == [2, 2]
        assert saved_data['properties'] == ['lambda_max', 'type']
        assert saved_data['entities'] == ['opsin_1', 'opsin_2']
    
    def test_persistence_across_instances(self):
        """Test that data persists across processor instances."""
        # Create first processor and save data
        proc1 = PropertyProcessor(name="test_persist_1")
        data = {'entity_A': {'value': 100}}
        proc1.create_property_table('persistent_data', data)
        
        # Create new processor instance
        proc2 = PropertyProcessor(name="test_persist_2")
        
        # Should be able to load the data
        df = proc2.load_property_table('persistent_data')
        assert df.loc['entity_A', 'value'] == 100
        
        # List datasets should show it
        datasets = proc2.list_datasets()
        assert 'persistent_data' in datasets


class TestRealWorldScenarios:
    """Test real-world usage scenarios."""
    
    def test_opsin_property_workflow(self):
        """Test a complete opsin property management workflow."""
        processor = PropertyProcessor(name="opsin_workflow")
        
        # Create opsin spectral properties
        spectral_data = {
            'bacteriorhodopsin': {
                'lambda_max': 568,
                'quantum_yield': 0.64,
                'photocycle': 'fast'
            },
            'channelrhodopsin2': {
                'lambda_max': 470,
                'quantum_yield': 0.30,
                'photocycle': 'slow'
            },
            'halorhodopsin': {
                'lambda_max': 578,
                'quantum_yield': 0.55,
                'photocycle': 'fast'
            }
        }
        
        processor.create_property_table(
            'opsin_spectral', 
            spectral_data,
            metadata={'experiment': 'spectroscopy', 'temperature': '20C'}
        )
        
        # Add functional properties
        processor.add_property_column(
            'opsin_spectral',
            'pump_type',
            {
                'bacteriorhodopsin': 'H+',
                'channelrhodopsin2': 'cation_channel',
                'halorhodopsin': 'Cl-'
            }
        )
        
        # Filter fast photocycle opsins
        fast_opsins = processor.filter_by_property(
            'opsin_spectral',
            'photocycle',
            lambda x: x == 'fast'
        )
        
        assert len(fast_opsins) == 2
        assert 'bacteriorhodopsin' in fast_opsins.index
        assert 'halorhodopsin' in fast_opsins.index
        
        # Create expression data in separate table
        expression_data = {
            'bacteriorhodopsin': {'expression_system': 'E.coli', 'yield_mg_L': 50},
            'channelrhodopsin2': {'expression_system': 'HEK293', 'yield_mg_L': 10}
        }
        
        processor.create_property_table('opsin_expression', expression_data)
        
        # Merge spectral and expression data
        complete_data = processor.merge_property_tables(
            ['opsin_spectral', 'opsin_expression'],
            'opsin_complete',
            how='outer'
        )
        
        # Verify complete dataset
        assert len(complete_data) == 3
        assert 'lambda_max' in complete_data.columns
        assert 'expression_system' in complete_data.columns
        
        # Check specific opsin has both types of data
        br_props = processor.get_entity_properties('bacteriorhodopsin', 'opsin_complete')
        assert br_props['lambda_max'] == 568
        assert br_props['pump_type'] == 'H+'
        assert br_props['expression_system'] == 'E.coli'
        assert br_props['yield_mg_L'] == 50