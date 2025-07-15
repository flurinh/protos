"""
Test PropertyProcessor with real property data.

This test suite verifies that PropertyProcessor correctly implements
the Protos data management principles:
1. Uses ProtosPaths exclusively (no path parameters)
2. Works with real property data files
3. Saves files with human-readable names
4. Properly implements entity and dataset operations
"""

import pytest
import pandas as pd
import json
from pathlib import Path

from protos.processing.property.property_processor import PropertyProcessor


class TestPropertyProcessorRealData:
    """Test PropertyProcessor with real property data."""
    
    @pytest.fixture
    def sample_property_data(self):
        """Create sample property data for testing."""
        return pd.DataFrame({
            'entity_id': ['BACR_HALSA', 'ChR2', 'NpHR', 'C1C2', 'ReaChR'],
            'lambda_max': [568, 470, 590, 460, 590],
            'photocycle_time_ms': [15.0, 20.0, 4000.0, 10.0, 50.0],
            'pump_type': ['proton', 'channel', 'chloride', 'channel', 'channel'],
            'expression_system': ['E.coli', 'HEK293', 'E.coli', 'HEK293', 'HEK293']
        })
    
    @pytest.fixture
    def experimental_properties(self):
        """Create experimental measurement data."""
        return {
            'BACR_HALSA': {
                'lambda_max': 568,
                'quantum_yield': 0.65,
                'extinction_coefficient': 52000,
                'photocycle_type': 'BR-like',
                'measured_date': '2024-01-15'
            },
            'ChR2': {
                'lambda_max': 470,
                'conductance_pS': 50,
                'channel_type': 'non-selective cation',
                'inactivation_tau_ms': 10,
                'measured_date': '2024-01-20'
            }
        }
    
    def test_zero_configuration_with_real_data(self):
        """Test that PropertyProcessor works with zero configuration."""
        # Initialize with no parameters
        processor = PropertyProcessor()
        
        # Should work immediately
        assert processor is not None
        assert processor.processor_type == "property"
        
        # Check that paths are created
        assert processor.data_path.exists()
        assert processor.data_dirs['datasets'].exists()
        assert processor.data_dirs['tables'].exists()
    
    def test_save_and_load_individual_properties(self):
        """Test saving and loading properties for individual entities."""
        processor = PropertyProcessor()
        
        # Save properties for BACR_HALSA
        bacr_props = {
            'lambda_max': 568,
            'quantum_yield': 0.65,
            'extinction_coefficient': 52000,
            'organism': 'Halobacterium salinarum'
        }
        
        processor.save_entity('BACR_HALSA', bacr_props, 
                            metadata={'dataset': 'microbial_opsins'})
        
        # Load back and verify
        loaded_props = processor.load_entity('BACR_HALSA')
        assert loaded_props is not None
        assert 'microbial_opsins' in loaded_props
        
        # Properties are stored with metadata, so check the 'value' field
        assert loaded_props['microbial_opsins']['lambda_max']['value'] == 568
        assert loaded_props['microbial_opsins']['quantum_yield']['value'] == 0.65
        
        # Verify entity is registered with human-readable name
        assert 'BACR_HALSA' in processor.entity_registry.list_entities('property')
    
    def test_batch_property_assignment(self, experimental_properties):
        """Test batch assignment of properties."""
        processor = PropertyProcessor()
        
        # Prepare batch assignments
        assignments = []
        for entity_name, props in experimental_properties.items():
            for prop_name, prop_value in props.items():
                assignments.append({
                    'entity_identifier': entity_name,
                    'property_name': prop_name,
                    'property_value': prop_value,
                    'metadata': {'source': 'experimental_measurements'}
                })
        
        # Batch assign
        entity_ids = processor.assign_properties_batch(
            assignments, 
            'experimental_data'
        )
        
        assert len(entity_ids) > 0
        
        # Verify assignments
        bacr_lambda = processor.get_entity_property('BACR_HALSA', 'lambda_max', 'experimental_data')
        assert bacr_lambda == 568
        
        chr2_conductance = processor.get_entity_property('ChR2', 'conductance_pS', 'experimental_data')
        assert chr2_conductance == 50
    
    def test_property_dataset_creation_from_dataframe(self, sample_property_data):
        """Test creating a property dataset from a pandas DataFrame."""
        processor = PropertyProcessor()
        
        # Save the DataFrame as a CSV file first
        csv_path = processor.data_dirs['tables'] / 'opsin_properties.csv'
        sample_property_data.to_csv(csv_path, index=False)
        
        # Create dataset from CSV
        count = processor.create_property_dataset_from_file(
            str(csv_path),
            'opsin_characterization',
            entity_column='entity_id'
        )
        
        assert count == 5  # All 5 entities
        
        # Load dataset and verify
        df = processor.get_dataset_properties('opsin_characterization')
        assert len(df) == 5
        assert 'lambda_max' in df.columns
        assert 'pump_type' in df.columns
        
        # Check specific values
        assert df.loc['BACR_HALSA', 'lambda_max'] == 568
        assert df.loc['ChR2', 'pump_type'] == 'channel'
    
    def test_property_filtering(self, sample_property_data):
        """Test filtering entities by property values."""
        processor = PropertyProcessor()
        
        # Create dataset first
        csv_path = processor.data_dirs['tables'] / 'opsin_filter_test.csv'
        sample_property_data.to_csv(csv_path, index=False)
        processor.create_property_dataset_from_file(
            str(csv_path),
            'filter_test',
            entity_column='entity_id'
        )
        
        # Filter by exact value
        channels = processor.filter_entities_by_property(
            'filter_test',
            {'pump_type': 'channel'}
        )
        assert len(channels) == 3
        assert 'ChR2' in channels
        assert 'C1C2' in channels
        assert 'ReaChR' in channels  # Also a channel
        
        # Filter by numeric comparison
        red_shifted = processor.filter_entities_by_property(
            'filter_test',
            {'lambda_max': {'gt': 500}}
        )
        assert len(red_shifted) == 3  # BACR_HALSA (568), NpHR (590), ReaChR (590)
        assert 'BACR_HALSA' in red_shifted
        assert 'NpHR' in red_shifted
        assert 'ReaChR' in red_shifted
    
    def test_dataset_persistence(self):
        """Test that datasets persist across processor instances."""
        # First processor - create dataset
        processor1 = PropertyProcessor()
        
        # Assign properties
        processor1.assign_property(
            'test_protein_1',
            'stability_tm',
            65.5,
            'thermal_stability',
            metadata={'method': 'DSF', 'buffer': 'PBS'}
        )
        
        processor1.assign_property(
            'test_protein_2',
            'stability_tm',
            72.3,
            'thermal_stability',
            metadata={'method': 'DSF', 'buffer': 'PBS'}
        )
        
        # Save dataset
        processor1.save_property_dataset('thermal_stability', 'both')
        
        # New processor instance
        processor2 = PropertyProcessor()
        
        # Load dataset
        success = processor2.load_property_dataset('thermal_stability')
        assert success
        
        # Verify data
        tm1 = processor2.get_entity_property('test_protein_1', 'stability_tm', 'thermal_stability')
        assert tm1 == 65.5
        
        tm2 = processor2.get_entity_property('test_protein_2', 'stability_tm', 'thermal_stability')
        assert tm2 == 72.3
        
        # Verify file organization
        # JSON metadata should be in datasets/
        json_file = processor2.data_dirs['datasets'] / 'thermal_stability.json'
        assert json_file.exists()
        
        # CSV data should be in tables/
        csv_file = processor2.data_dirs['tables'] / 'thermal_stability.csv'
        assert csv_file.exists()
    
    def test_dataset_statistics(self, sample_property_data):
        """Test dataset statistics calculation."""
        processor = PropertyProcessor()
        
        # Create dataset
        csv_path = processor.data_dirs['tables'] / 'stats_test.csv'
        sample_property_data.to_csv(csv_path, index=False)
        processor.create_property_dataset_from_file(
            str(csv_path),
            'stats_dataset',
            entity_column='entity_id'
        )
        
        # Get statistics
        stats = processor.get_dataset_statistics('stats_dataset')
        
        assert stats['entity_count'] == 5
        assert stats['property_count'] == 4  # Excluding entity_id
        
        # Check property statistics
        assert 'lambda_max' in stats['properties']
        lambda_stats = stats['properties']['lambda_max']
        assert lambda_stats['type'] == 'int64'
        assert lambda_stats['mean'] == pytest.approx((568 + 470 + 590 + 460 + 590) / 5)
        assert lambda_stats['min'] == 460
        assert lambda_stats['max'] == 590
    
    def test_secondary_selection_grn_property(self):
        """Test assigning properties to specific GRN positions."""
        processor = PropertyProcessor()
        
        # Assign property to specific GRN position
        entity_id = processor.assign_grn_property(
            'BACR_HALSA',
            '3.50',  # GRN position
            'residue_accessibility',
            0.15,
            'grn_structural_properties'
        )
        
        # Retrieve property
        value = processor.get_entity_property(
            'BACR_HALSA',
            '3.50.residue_accessibility',
            'grn_structural_properties'
        )
        assert value == 0.15
    
    def test_integration_with_entity_registry(self):
        """Test that PropertyProcessor properly integrates with EntityRegistry."""
        processor = PropertyProcessor()
        
        # Save entity with properties
        processor.save_entity(
            'novel_opsin_2024',
            {
                'lambda_max': 515,
                'discovered_date': '2024-03-15',
                'source_organism': 'uncultured marine bacterium'
            },
            metadata={'dataset': 'novel_opsins', 'publication': 'doi:10.1234/test'}
        )
        
        # Check entity is registered
        entity_info = processor.entity_registry.find_entity('novel_opsin_2024')
        assert entity_info is not None
        # EntityInfo has format_type attribute
        assert hasattr(entity_info, 'format_type')
        assert entity_info.format_type == 'property'
        
        # List all property entities
        property_entities = processor.entity_registry.list_entities('property')
        assert 'novel_opsin_2024' in property_entities
    
    def test_human_readable_file_organization(self):
        """Test that all files use human-readable names."""
        processor = PropertyProcessor()
        
        # Create some test data
        processor.assign_property(
            'human_GPCR_A2A',
            'binding_affinity_nM',
            25.3,
            'gpcr_ligand_binding'
        )
        
        processor.save_property_dataset('gpcr_ligand_binding', 'both')
        
        # Check file organization
        # Dataset JSON should be in datasets/ (not in subdirectory)
        dataset_json = processor.data_dirs['datasets'] / 'gpcr_ligand_binding.json'
        assert dataset_json.exists()
        
        # CSV should be in tables/
        csv_file = processor.data_dirs['tables'] / 'gpcr_ligand_binding.csv'
        assert csv_file.exists()
        
        # Load and check content uses human names
        with open(dataset_json, 'r') as f:
            dataset_data = json.load(f)
        
        # Should have human-readable entity names in the data
        assert 'entities' in dataset_data
        assert 'human_GPCR_A2A' in str(dataset_data)
        
        # Check CSV also uses human names
        import pandas as pd
        df = pd.read_csv(csv_file, index_col=0)
        assert 'human_GPCR_A2A' in df.index
        
        # Entity properties cache should exist
        cache_file = processor.data_dirs['cache'] / 'entity_properties.json'
        assert cache_file.exists()
    
    def test_list_operations(self):
        """Test listing entities and datasets."""
        processor = PropertyProcessor()
        
        # Create some test data
        test_entities = ['protein_A', 'protein_B', 'protein_C']
        for entity in test_entities:
            processor.assign_property(
                entity,
                'test_value',
                42,
                'list_test_dataset'
            )
        
        # List entities in dataset
        entities = processor.list_entities('list_test_dataset')
        assert len(entities) >= 3
        for entity in test_entities:
            assert entity in entities
        
        # List all datasets
        datasets = processor.list_property_datasets()
        assert 'list_test_dataset' in datasets
        
        # Get dataset info as list
        dataset_list = processor.list_datasets()
        list_test = next((d for d in dataset_list if d['id'] == 'list_test_dataset'), None)
        assert list_test is not None
        assert list_test['entity_count'] >= 3
        assert list_test['type'] == 'property_dataset'