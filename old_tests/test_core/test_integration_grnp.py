"""
Tests for integrating BaseProcessor with GRNProcessor.

This old_tests the combined functionality of the BaseProcessor class with
the existing GRNProcessor in a real-world scenario.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil
from datetime import datetime

from protos.core.base_processor import BaseProcessor
from protos.processing.grn.grn_processor import GRNProcessor


# Fixtures for testing
@pytest.fixture
def test_data_root():
    """Create a test data directory structure."""
    # Create a unique test directory in the protos/old_tests/data path
    base_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "temp")
    os.makedirs(base_dir, exist_ok=True)
    
    # Create a unique test directory
    dir_name = f"test_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
    data_dir = os.path.join(base_dir, dir_name)
    
    # Create main directories
    for subdir in ["grn/reference", "grn/opsin", "structures/mmcif", 
                  "embeddings", "properties"]:
        os.makedirs(os.path.join(data_dir, subdir), exist_ok=True)
    
    yield data_dir
    
    # Clean up
    shutil.rmtree(data_dir, ignore_errors=True)


@pytest.fixture
def sample_grn_table():
    """Sample GRN table for testing."""
    data = {
        "protein_id": ["opsin1", "opsin2", "opsin3"],
        "3.50": ["M125", "M123", "M127"],
        "6.48": ["Y261", "Y259", "Y263"],
        "7.49": ["P296", "P294", "P298"]
    }
    return pd.DataFrame(data)


class GRNBaseProcessor(BaseProcessor):
    """GRN processor that utilizes the BaseProcessor functionality."""
    
    def __init__(self, name="grn_processor", data_root=None, processor_data_dir="grn"):
        super().__init__(name=name, data_root=data_root, processor_data_dir=processor_data_dir)
        self.ids = []
        self.grns = []
        self.grn_dict = {}
    
    def load_grn_table(self, dataset_id, **kwargs):
        """Load a GRN table using BaseProcessor."""
        data = self.load_data(dataset_id, file_format="csv", **kwargs)
        
        # Handle index settings
        if "protein_id" in data.columns and data.index.name != "protein_id":
            data = data.set_index("protein_id")
        
        # Set internal data
        self.data = data
        self.ids = list(data.index)
        self.grns = [col for col in data.columns if col != "protein_id"]
        
        return self.data
    
    def save_grn_table(self, dataset_id, **kwargs):
        """Save a GRN table using BaseProcessor."""
        # Reset index if protein_id is the index name
        if self.data is not None and self.data.index.name == "protein_id":
            data_to_save = self.data.reset_index()
        else:
            data_to_save = self.data
            
        return self.save_data(dataset_id, data_to_save, file_format="csv", **kwargs)
    
    def get_grn_dict(self):
        """Get GRN dictionary mapping proteins to GRN positions."""
        if self.data is None:
            return {}
        
        self.grn_dict = {}
        for protein_id in self.ids:
            self.grn_dict[protein_id] = {}
            for grn in self.grns:
                value = self.data.loc[protein_id, grn]
                if pd.notna(value) and value != '':
                    self.grn_dict[protein_id][grn] = value
        
        return self.grn_dict
    
    def filter_by_ids(self, ids_to_keep):
        """Filter data to include only the specified IDs."""
        if self.data is None:
            return
        
        # Filter data - use a more explicit approach to avoid pandas dtype issues
        valid_ids = [id for id in self.ids if id in ids_to_keep]
        if valid_ids:
            self.data = self.data.loc[valid_ids]
            self.ids = valid_ids
        else:
            self.logger.warning(f"No valid IDs found in {ids_to_keep}")
    
    def apply_interval(self, grn_interval):
        """Limit data to specific GRN positions."""
        if self.data is None:
            return
        
        # Keep only the specified GRNs
        valid_grns = [grn for grn in grn_interval if grn in self.data.columns]
        if not valid_grns:
            self.logger.warning("No valid GRNs in interval")
            return
        
        self.data = self.data[valid_grns]
        self.grns = valid_grns


class TestBaseProcessorIntegration:
    """Tests for integrating BaseProcessor with GRNProcessor."""
    
    def test_basic_integration(self, test_data_root, sample_grn_table):
        """Test basic integration of BaseProcessor with GRNProcessor."""
        # Save sample data to test directory
        grn_path = os.path.join(test_data_root, "grn", "reference", "test_grn.csv")
        sample_grn_table.to_csv(grn_path, index=False)
        
        # Create processor with test data root
        processor = GRNBaseProcessor(name="test_integration", data_root=test_data_root, processor_data_dir="grn")
        
        # Load the test data
        processor.load_grn_table("reference/test_grn")
        
        # Verify data loaded correctly
        assert processor.data is not None
        assert len(processor.ids) == 3
        assert len(processor.grns) == 3
        assert "3.50" in processor.grns
        
        # Test GRN dictionary
        grn_dict = processor.get_grn_dict()
        assert len(grn_dict) == 3
        assert "3.50" in grn_dict["opsin1"]
        assert grn_dict["opsin1"]["3.50"] == "M125"
        
        # Test filtering
        processor.filter_by_ids(["opsin1", "opsin3"])
        assert len(processor.ids) == 2
        assert "opsin1" in processor.ids
        assert "opsin3" in processor.ids
        assert "opsin2" not in processor.ids
        
        # Test GRN interval
        processor.apply_interval(["3.50", "7.49"])
        assert len(processor.grns) == 2
        assert "3.50" in processor.grns
        assert "7.49" in processor.grns
        assert "6.48" not in processor.grns
        
        # Test saving
        new_path = processor.save_grn_table("reference/new_test_grn")
        assert os.path.exists(new_path)
        
        # Test registry
        datasets = processor.list_datasets()
        assert len(datasets) >= 1
        
        # Get dataset info
        info = processor.get_dataset_info("reference/new_test_grn")
        assert info is not None
        assert info["format"] == "csv"
        assert info.get("directory") == "grn"
    
    @pytest.mark.skipif(not os.path.exists("/mnt/c/Users/hidbe/PycharmProjects/phd/protos/data/grn/opsin/bistable_ao.csv"), 
                       reason="Real data file not available")
    def test_real_data_integration(self, test_data_root):
        """Test integration with real opsin data."""
        # Copy real data file to test directory
        source_file = "/mnt/c/Users/hidbe/PycharmProjects/phd/protos/data/grn/opsin/bistable_ao.csv"
        target_file = os.path.join(test_data_root, "grn", "opsin", "bistable_ao.csv")
        
        # Only run if the real data file exists
        if os.path.exists(source_file):
            shutil.copy(source_file, target_file)
            
            # Create processor with test data root
            processor = GRNBaseProcessor(data_root=test_data_root, processor_data_dir="grn")
            
            # Load the real data
            processor.load_grn_table("opsin/bistable_ao")
            
            # Verify data loaded correctly
            assert processor.data is not None
            assert len(processor.ids) > 0
            assert len(processor.grns) > 0
            
            # Test filtering to key positions
            key_positions = ["1.50", "3.50", "6.48", "7.49"]
            valid_positions = [pos for pos in key_positions if pos in processor.grns]
            
            if valid_positions:
                processor.apply_interval(valid_positions)
                assert all(pos in processor.grns for pos in valid_positions)
                
                # Test GRN dictionary with real data
                grn_dict = processor.get_grn_dict()
                assert isinstance(grn_dict, dict)
                assert len(grn_dict) > 0