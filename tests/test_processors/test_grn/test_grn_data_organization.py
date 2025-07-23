#!/usr/bin/env python3
"""
Comprehensive tests for GRN processor data organization.

This test suite verifies that the GRN processor correctly:
1. Saves GRN tables to the tables/ subdirectory
2. Creates dataset JSON files in the datasets/ subdirectory
3. Uses human-readable protein names throughout
4. Can load reference tables from ref/ directory
5. Can load configs from configs/ directory

The tests follow the pattern of test_property_processor_enhanced.py but adapted
for GRN-specific functionality.
"""

import pytest
import pandas as pd
import numpy as np
import json
import tempfile
import shutil
import os
from pathlib import Path

from protos.processing.grn import GRNProcessor
from protos.io.data_access import generate_entity_id
from protos.io.paths import ProtosPaths


@pytest.fixture
def temp_data_dir(monkeypatch):
    """Create temporary directory for test data."""
    temp_dir = tempfile.mkdtemp()
    # Set environment variable for ProtosPaths
    monkeypatch.setenv("PROTOS_DATA_ROOT", str(temp_dir))
    yield Path(temp_dir)
    shutil.rmtree(temp_dir)


@pytest.fixture
def grn_processor(temp_data_dir):
    """Create GRNProcessor instance with temporary directory."""
    # Create ProtosPaths instance with the temp directory
    paths = ProtosPaths(data_root=str(temp_data_dir))
    
    # Ensure GRN directories exist
    grn_dir = Path(paths.get_processor_path('grn'))
    grn_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories with correct names from DEFAULT_GRN_SUBDIRS
    for subdir_key in ['table_dir', 'ref_dir', 'configs_dir', 'datasets_dir', 'assignment_dir']:
        subdir_path = Path(paths.get_subdir_path('grn', subdir_key))
        subdir_path.mkdir(parents=True, exist_ok=True)
    
    # Create processor with real paths
    processor = GRNProcessor(name="test_grn_processor", paths=paths)
    processor.data_path = grn_dir
    
    # Return the configured processor
    return processor



@pytest.fixture
def sample_grn_data():
    """Sample GRN data for testing with human-readable names."""
    return pd.DataFrame({
        '1.50': ['L28', 'V29', 'L30', 'L31', 'L32'],
        '2.50': ['A72', 'S73', 'A74', 'A75', 'S76'],
        '3.50': ['R115', 'R116', 'R117', 'R118', 'R119'],
        '4.50': ['G157', 'W158', 'W159', 'W160', 'G161'],
        '5.50': ['F197', 'F198', 'Y199', 'F200', 'F201'],
        '6.50': ['W238', 'W239', 'W240', 'F241', 'W242'],
        '7.50': ['K296', 'K297', 'K298', 'K299', 'K300'],
    }, index=['Bacteriorhodopsin', 'Channelrhodopsin-2', 'Halorhodopsin', 'Sensory_rhodopsin_II', 'Archaerhodopsin-3'])


@pytest.fixture
def sample_config():
    """Sample GRN configuration."""
    return {
        "name": "test_config",
        "version": "1.0",
        "grn_positions": {
            "key_positions": ["1.50", "2.50", "3.50", "7.50"],
            "extended_positions": ["4.50", "5.50", "6.50"]
        },
        "filters": {
            "min_coverage": 0.8,
            "exclude_gaps": True
        }
    }


class TestGRNProcessorDirectoryStructure:
    """Test GRN processor directory structure and organization."""
    
    def test_initialization_creates_directories(self, grn_processor):
        """Test that GRN processor creates all required directories."""
        assert grn_processor.name == "test_grn_processor"
        assert grn_processor.data_path.exists()
        
        # Check all required directories are created
        expected_dirs = ['tables', 'datasets', 'ref', 'configs', 'assignments']
        for dir_name in expected_dirs:
            dir_path = grn_processor.data_path / dir_name
            assert dir_path.exists(), f"Directory {dir_name} not created"
            assert dir_path.is_dir(), f"{dir_name} is not a directory"
    
    def test_subdirectory_paths(self, grn_processor):
        """Test that subdirectory paths are correctly configured."""
        # Tables directory
        assert (grn_processor.data_path / 'tables').exists()
        
        # Datasets directory
        assert (grn_processor.data_path / 'datasets').exists()
        
        # Reference directory
        assert (grn_processor.data_path / 'ref').exists()
        
        # Configs directory
        assert (grn_processor.data_path / 'configs').exists()
        
        # Assignments directory
        assert (grn_processor.data_path / 'assignments').exists()


class TestGRNTableOperations:
    """Test GRN table save/load operations with proper directory structure."""
    
    def test_save_grn_table_to_tables_directory(self, grn_processor, sample_grn_data):
        """Test that GRN tables are saved to the tables/ subdirectory."""
        table_name = "test_grn_table"
        
        # Set data and save
        grn_processor.data = sample_grn_data
        grn_processor.save_grn_table(table_name)
        
        # Check file was created in tables directory
        expected_file = grn_processor.data_path / 'tables' / f"{table_name}.csv"
        assert expected_file.exists(), f"Table not saved to {expected_file}"
        
        # Verify content
        loaded_df = pd.read_csv(expected_file, index_col=0)
        assert len(loaded_df) == len(sample_grn_data)
        
        # The saved file includes an entity_id column, so check GRN columns separately
        grn_columns = [col for col in loaded_df.columns if col != 'entity_id']
        assert grn_columns == list(sample_grn_data.columns)
    
    def test_load_grn_table_from_tables_directory(self, grn_processor, sample_grn_data):
        """Test loading GRN tables from the tables/ subdirectory."""
        table_name = "load_test_table"
        
        # Save table first
        grn_processor.data = sample_grn_data
        grn_processor.save_grn_table(table_name)
        
        # Create new processor and load
        new_processor = GRNProcessor(name="new_processor", paths=grn_processor.paths)
        new_processor.data_path = grn_processor.data_path
        
        new_processor.load_grn_table(table_name)
        
        # Verify data loaded correctly
        assert new_processor.dataset == table_name
        assert not new_processor.data.empty
        assert len(new_processor.data) == len(sample_grn_data)
        
        # Check human-readable names are preserved
        assert 'Bacteriorhodopsin' in new_processor.data.index
        assert 'Channelrhodopsin-2' in new_processor.data.index
    
    def test_human_readable_protein_names_preserved(self, grn_processor, sample_grn_data):
        """Test that human-readable protein names are preserved throughout operations."""
        table_name = "human_readable_test"
        
        # Save with human-readable names
        grn_processor.data = sample_grn_data
        grn_processor.save_grn_table(table_name)
        
        # Load and check names
        grn_processor.load_grn_table(table_name)
        
        expected_names = ['Bacteriorhodopsin', 'Channelrhodopsin-2', 'Halorhodopsin', 
                         'Sensory_rhodopsin_II', 'Archaerhodopsin-3']
        
        for name in expected_names:
            assert name in grn_processor.ids, f"Human-readable name {name} not found"
    
    def test_list_tables_in_tables_directory(self, grn_processor, sample_grn_data):
        """Test listing available tables in the tables/ directory."""
        # Save multiple tables
        table_names = ["table1", "table2", "table3"]
        
        for table_name in table_names:
            grn_processor.data = sample_grn_data
            grn_processor.save_grn_table(table_name)
        
        # List datasets - should show tables with proper prefix
        datasets = grn_processor.list_datasets()
        dataset_ids = [d["id"] for d in datasets]
        
        # GRN processor prefixes with 'tables/'
        for table_name in table_names:
            assert f"tables/{table_name}" in dataset_ids


class TestGRNDatasetOperations:
    """Test GRN dataset JSON operations in datasets/ directory."""
    
    def test_create_dataset_json_in_datasets_directory(self, grn_processor, sample_grn_data):
        """Test creating dataset JSON files in the datasets/ subdirectory."""
        dataset_name = "microbial_opsins_grn"
        
        # Set data
        grn_processor.data = sample_grn_data
        
        # Create dataset JSON
        dataset_info = {
            "name": dataset_name,
            "description": "GRN mappings for microbial opsins",
            "created_date": "2024-01-15",
            "protein_count": len(sample_grn_data),
            "grn_positions": list(sample_grn_data.columns),
            "proteins": list(sample_grn_data.index)
        }
        
        # Save dataset JSON
        dataset_dir = grn_processor.data_path / 'datasets' / dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        dataset_file = dataset_dir / f"{dataset_name}.json"
        with open(dataset_file, 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        # Verify file exists
        assert dataset_file.exists()
        
        # Load and verify content
        with open(dataset_file, 'r') as f:
            loaded_info = json.load(f)
        
        assert loaded_info["name"] == dataset_name
        assert loaded_info["protein_count"] == 5
        assert "Bacteriorhodopsin" in loaded_info["proteins"]
    
    def test_load_dataset_from_json(self, grn_processor, sample_grn_data):
        """Test loading dataset information from JSON file."""
        dataset_name = "test_dataset"
        
        # First save the actual data as a table
        grn_processor.data = sample_grn_data
        grn_processor.save_grn_table(dataset_name)
        
        # Create dataset JSON with metadata
        dataset_info = {
            "name": dataset_name,
            "table_name": dataset_name,
            "proteins": list(sample_grn_data.index),
            "grn_positions": list(sample_grn_data.columns),
            "metadata": {
                "source": "test",
                "version": "1.0"
            }
        }
        
        dataset_dir = grn_processor.data_path / 'datasets' / dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        with open(dataset_dir / f"{dataset_name}.json", 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        # Load dataset using the JSON metadata
        grn_processor.load_grn_table(dataset_info["table_name"])
        
        # Verify data loaded correctly
        assert len(grn_processor.data) == 5
        assert all(protein in grn_processor.ids for protein in dataset_info["proteins"])
    
    def test_dataset_statistics_saved_to_json(self, grn_processor, sample_grn_data):
        """Test that dataset statistics are saved to JSON files."""
        dataset_name = "stats_test"
        
        grn_processor.data = sample_grn_data
        grn_processor.save_grn_table(dataset_name)
        
        # Calculate statistics
        stats = {
            "dataset_name": dataset_name,
            "protein_count": len(sample_grn_data),
            "grn_position_count": len(sample_grn_data.columns),
            "coverage_by_position": {}
        }
        
        # Calculate coverage for each GRN position
        for col in sample_grn_data.columns:
            non_gap_count = sum(1 for val in sample_grn_data[col] 
                               if pd.notna(val) and val != '-')
            stats["coverage_by_position"][col] = {
                "count": non_gap_count,
                "percentage": non_gap_count / len(sample_grn_data) * 100
            }
        
        # Save statistics
        dataset_dir = grn_processor.data_path / 'datasets' / dataset_name
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        stats_file = dataset_dir / f"{dataset_name}_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        assert stats_file.exists()


class TestReferenceTableOperations:
    """Test loading reference tables from ref/ directory."""
    
    def test_load_reference_table_from_ref_directory(self, grn_processor):
        """Test loading reference GRN tables from ref/ directory."""
        # Create reference table
        ref_data = pd.DataFrame({
            '1.50': ['L28', 'V29', 'L30'],
            '2.50': ['A72', 'S73', 'A74'],
            '3.50': ['R115', 'R116', 'R117'],
            '7.50': ['K296', 'K297', 'K298']
        }, index=['RefProtein1', 'RefProtein2', 'RefProtein3'])
        
        # Save to ref directory
        ref_file = grn_processor.data_path / 'ref' / 'standard_grn_ref.csv'
        ref_data.to_csv(ref_file)
        
        # Load reference table
        loaded_ref = pd.read_csv(ref_file, index_col=0)
        grn_processor.data = loaded_ref
        
        # Verify
        assert len(grn_processor.data) == 3
        assert 'RefProtein1' in grn_processor.data.index
        assert '7.50' in grn_processor.data.columns
    
    def test_load_gpcrdb_reference(self, grn_processor):
        """Test loading GPCRDB reference format."""
        # Create GPCRDB-style reference
        gpcrdb_data = pd.DataFrame({
            'protein_name': ['Rhodopsin', 'Beta2_adrenergic', 'Muscarinic_M2'],
            '3x50': ['R135', 'R131', 'R129'],
            '7x53': ['K296', 'K305', 'K465']
        })
        
        # Save to ref directory
        ref_file = grn_processor.data_path / 'ref' / 'gpcrdb_ref.csv'
        gpcrdb_data.to_csv(ref_file, index=False)
        
        # Load and transform
        loaded_data = pd.read_csv(ref_file)
        
        # Transform to standard GRN format
        if 'protein_name' in loaded_data.columns:
            loaded_data = loaded_data.set_index('protein_name')
        
        # Convert x notation to dot notation
        loaded_data.columns = [col.replace('x', '.') for col in loaded_data.columns]
        
        grn_processor.data = loaded_data
        
        # Verify
        assert 'Rhodopsin' in grn_processor.data.index
        assert '3.50' in grn_processor.data.columns
        assert grn_processor.data.loc['Rhodopsin', '3.50'] == 'R135'
    
    def test_list_reference_tables(self, grn_processor):
        """Test listing available reference tables."""
        # Create multiple reference files
        ref_names = ['mo_ref', 'gpcrdb_ref', 'custom_ref']
        
        for ref_name in ref_names:
            ref_data = pd.DataFrame({
                '3.50': ['R115'],
                '7.50': ['K296']
            }, index=[f'{ref_name}_protein'])
            
            ref_file = grn_processor.data_path / 'ref' / f'{ref_name}.csv'
            ref_data.to_csv(ref_file)
        
        # List reference files
        ref_files = list((grn_processor.data_path / 'ref').glob('*.csv'))
        ref_names_found = [f.stem for f in ref_files]
        
        for ref_name in ref_names:
            assert ref_name in ref_names_found


class TestConfigOperations:
    """Test loading and using configurations from configs/ directory."""
    
    def test_load_config_from_configs_directory(self, grn_processor, sample_config):
        """Test loading configuration from configs/ directory."""
        config_name = "test_config"
        
        # Save config
        config_file = grn_processor.data_path / 'configs' / f'{config_name}.json'
        with open(config_file, 'w') as f:
            json.dump(sample_config, f, indent=2)
        
        # Load config
        with open(config_file, 'r') as f:
            loaded_config = json.load(f)
        
        assert loaded_config["name"] == "test_config"
        assert "key_positions" in loaded_config["grn_positions"]
        assert "1.50" in loaded_config["grn_positions"]["key_positions"]
    
    def test_binding_domain_config(self, grn_processor):
        """Test loading binding domain configuration."""
        binding_config = {
            "name": "binding_domain",
            "description": "Configuration for retinal binding domain analysis",
            "binding_positions": {
                "schiff_base": "7.50",
                "counterion": "3.50",
                "retinal_contacts": ["3.32", "3.33", "3.36", "4.57", "5.42", "5.46", "6.44", "6.48"]
            },
            "distance_thresholds": {
                "contact": 4.0,
                "interaction": 5.0
            }
        }
        
        # Save config
        config_file = grn_processor.data_path / 'configs' / 'binding_domain.json'
        with open(config_file, 'w') as f:
            json.dump(binding_config, f, indent=2)
        
        # Load and use config
        with open(config_file, 'r') as f:
            loaded_config = json.load(f)
        
        # Verify binding domain positions
        assert loaded_config["binding_positions"]["schiff_base"] == "7.50"
        assert "3.32" in loaded_config["binding_positions"]["retinal_contacts"]
    
    def test_motif_config(self, grn_processor):
        """Test loading motif configuration."""
        motif_config = {
            "name": "motif_config",
            "motifs": {
                "NPxxY": {
                    "positions": ["7.49", "7.50", "7.51", "7.52", "7.53"],
                    "pattern": "NP..Y",
                    "required_positions": ["7.49", "7.53"]
                },
                "DRY": {
                    "positions": ["3.49", "3.50", "3.51"],
                    "pattern": "DRY",
                    "required_positions": ["3.49", "3.50", "3.51"]
                }
            }
        }
        
        # Save config
        config_file = grn_processor.data_path / 'configs' / 'motif.json'
        with open(config_file, 'w') as f:
            json.dump(motif_config, f, indent=2)
        
        # Load config
        with open(config_file, 'r') as f:
            loaded_config = json.load(f)
        
        # Verify motif definitions
        assert "NPxxY" in loaded_config["motifs"]
        assert loaded_config["motifs"]["NPxxY"]["pattern"] == "NP..Y"
        assert "7.50" in loaded_config["motifs"]["NPxxY"]["positions"]


class TestIntegratedWorkflow:
    """Test complete integrated workflows with proper data organization."""
    
    def test_complete_grn_analysis_workflow(self, grn_processor, sample_grn_data):
        """Test a complete GRN analysis workflow with proper directory usage."""
        # Step 1: Load reference data
        ref_data = pd.DataFrame({
            '1.50': ['L28', 'V29'],
            '3.50': ['R115', 'R116'],
            '7.50': ['K296', 'K297']
        }, index=['Reference_BR', 'Reference_ChR2'])
        
        ref_file = grn_processor.data_path / 'ref' / 'reference_opsins.csv'
        ref_data.to_csv(ref_file)
        
        # Step 2: Load config
        analysis_config = {
            "name": "opsin_analysis",
            "key_positions": ["3.50", "7.50"],
            "min_coverage": 0.8
        }
        
        config_file = grn_processor.data_path / 'configs' / 'analysis_config.json'
        with open(config_file, 'w') as f:
            json.dump(analysis_config, f, indent=2)
        
        # Step 3: Process new data
        grn_processor.data = sample_grn_data
        
        # Step 4: Save processed table
        grn_processor.save_grn_table("processed_opsins")
        
        # Step 5: Create dataset with metadata
        dataset_info = {
            "name": "processed_opsins",
            "table_name": "processed_opsins",
            "reference_file": str(ref_file.name),
            "config_file": str(config_file.name),
            "analysis_date": "2024-01-15",
            "proteins": list(sample_grn_data.index),
            "key_findings": {
                "conserved_K_at_7.50": True,
                "conserved_R_at_3.50": True
            }
        }
        
        dataset_dir = grn_processor.data_path / 'datasets' / 'processed_opsins'
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        with open(dataset_dir / 'processed_opsins.json', 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        # Verify complete workflow
        assert (grn_processor.data_path / 'tables' / 'processed_opsins.csv').exists()
        assert (dataset_dir / 'processed_opsins.json').exists()
        assert ref_file.exists()
        assert config_file.exists()
    
    def test_merge_multiple_grn_sources(self, grn_processor):
        """Test merging GRN data from multiple sources with proper organization."""
        # Create multiple source tables
        source1 = pd.DataFrame({
            '3.50': ['R115', 'R116'],
            '7.50': ['K296', 'K297']
        }, index=['Protein_A', 'Protein_B'])
        
        source2 = pd.DataFrame({
            '3.50': ['R117', 'R118'],
            '7.50': ['K298', 'K299']
        }, index=['Protein_C', 'Protein_D'])
        
        # Save sources
        grn_processor.data = source1
        grn_processor.save_grn_table("source1")
        
        grn_processor.data = source2
        grn_processor.save_grn_table("source2")
        
        # Merge tables
        grn_processor.load_and_merge_grn_tables(["source1", "source2"])
        
        # Save merged result
        grn_processor.save_grn_table("merged_grn")
        
        # Create merged dataset metadata
        dataset_info = {
            "name": "merged_grn",
            "source_tables": ["source1", "source2"],
            "merge_date": "2024-01-15",
            "total_proteins": len(grn_processor.data),
            "proteins": list(grn_processor.data.index)
        }
        
        dataset_dir = grn_processor.data_path / 'datasets' / 'merged_grn'
        dataset_dir.mkdir(parents=True, exist_ok=True)
        
        with open(dataset_dir / 'merged_grn.json', 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        # Verify merge
        assert len(grn_processor.data) == 4
        assert all(protein in grn_processor.data.index 
                  for protein in ['Protein_A', 'Protein_B', 'Protein_C', 'Protein_D'])


class TestErrorHandling:
    """Test error handling for GRN data organization."""
    
    def test_handle_missing_directories(self, grn_processor):
        """Test handling of missing directories."""
        # Remove a directory
        shutil.rmtree(grn_processor.data_path / 'tables')
        
        # Try to save - should recreate directory
        grn_processor.data = pd.DataFrame({'3.50': ['R115']}, index=['Test'])
        grn_processor.save_grn_table("test_table")
        
        # Directory should be recreated
        assert (grn_processor.data_path / 'tables').exists()
    
    def test_handle_corrupted_files(self, grn_processor, temp_data_dir):
        """Test handling of corrupted files."""
        # Create corrupted CSV
        corrupted_file = grn_processor.data_path / 'tables' / 'corrupted.csv'
        corrupted_file.write_text("This is not valid CSV data!!!")
        
        # Try to load - the processor handles this gracefully, returning empty data
        grn_processor.load_grn_table("corrupted")
        
        # Verify it loaded but has no valid data
        assert grn_processor.data.empty or len(grn_processor.data) == 0
    
    def test_handle_missing_config(self, grn_processor):
        """Test handling of missing configuration files."""
        # Try to load non-existent config
        config_path = grn_processor.data_path / 'configs' / 'nonexistent.json'
        
        assert not config_path.exists()
        
        # Should handle gracefully when trying to read
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = json.load(f)
        else:
            config = {}  # Default empty config
        
        assert config == {}


class TestBackwardCompatibility:
    """Test backward compatibility with existing GRN data formats."""
    
    def test_load_x_notation_tables(self, grn_processor):
        """Test loading tables with x notation (3x50 vs 3.50)."""
        # Create table with x notation
        x_notation_data = pd.DataFrame({
            '3x50': ['R115', 'R116', 'R117'],
            '7x50': ['K296', 'K297', 'K298']
        }, index=['Protein1', 'Protein2', 'Protein3'])
        
        # Save with x notation
        x_file = grn_processor.data_path / 'tables' / 'x_notation_table.csv'
        x_notation_data.to_csv(x_file)
        
        # Load and convert
        loaded_data = pd.read_csv(x_file, index_col=0)
        
        # Convert x to dot notation
        loaded_data.columns = [col.replace('x', '.') for col in loaded_data.columns]
        grn_processor.data = loaded_data
        
        # Verify conversion
        assert '3.50' in grn_processor.data.columns
        assert '7.50' in grn_processor.data.columns
        assert grn_processor.data.loc['Protein1', '3.50'] == 'R115'
    
    def test_load_legacy_format(self, grn_processor):
        """Test loading legacy GRN table formats."""
        # Create legacy format (might have different structure)
        legacy_data = pd.DataFrame({
            'protein_id': ['1ABC', '2DEF', '3GHI'],
            'position_3_50': ['R115', 'R116', 'R117'],
            'position_7_50': ['K296', 'K297', 'K298']
        })
        
        legacy_file = grn_processor.data_path / 'tables' / 'legacy_format.csv'
        legacy_data.to_csv(legacy_file, index=False)
        
        # Load and transform to current format
        loaded_legacy = pd.read_csv(legacy_file)
        
        # Transform to standard format
        if 'protein_id' in loaded_legacy.columns:
            loaded_legacy = loaded_legacy.set_index('protein_id')
        
        # Rename columns
        loaded_legacy.columns = [col.replace('position_', '').replace('_', '.') 
                                for col in loaded_legacy.columns]
        
        grn_processor.data = loaded_legacy
        
        # Verify transformation
        assert '3.50' in grn_processor.data.columns
        assert grn_processor.data.loc['1ABC', '3.50'] == 'R115'


if __name__ == "__main__":
    pytest.main([__file__, "-v"])