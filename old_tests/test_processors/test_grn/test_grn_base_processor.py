"""
Tests for the GRNBaseProcessor class.

These old_tests verify that the BaseProcessor-integrated GRN processor
works correctly with real world data.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil
from datetime import datetime

from protos.processing.grn.grn_utils import sort_grns_str
from protos.processing.grn.grn_base_processor import GRNBaseProcessor


# Fixtures for testing
@pytest.fixture
def test_data_dir():
    """Create a temporary test data directory."""
    data_dir = tempfile.mkdtemp(prefix="protos_test_")
    
    # Set up the directory structure
    for subdir in ["grn/reference", "grn/opsin", "structures", "embeddings"]:
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


@pytest.fixture
def grn_processor(test_data_dir, sample_grn_table):
    """Create a GRNBaseProcessor with test data."""
    # Save the sample data
    save_path = os.path.join(test_data_dir, "grn", "reference", "test_grn.csv")
    sample_grn_table.to_csv(save_path, index=False)
    
    # Create the processor
    processor = GRNBaseProcessor(
        name="test_processor",
        data_root=test_data_dir,
        processor_data_dir="grn",
        dataset=None,  # Don't preload any data
        preload=False
    )
    
    return processor


class TestGRNBaseProcessor:
    """Tests for the GRNBaseProcessor."""
    
    def test_initialization(self, grn_processor):
        """Test that the processor initializes correctly."""
        assert grn_processor.name == "test_processor"
        assert grn_processor.processor_data_dir == "grn"
        assert isinstance(grn_processor.dataset_registry, dict)
    
    def test_load_grn_table(self, grn_processor):
        """Test loading a GRN table."""
        # Load the test data
        data = grn_processor.load_grn_table("reference/test_grn")
        
        # Check that the data is loaded correctly
        assert isinstance(data, pd.DataFrame)
        assert len(data) == 3
        assert len(grn_processor.ids) == 3
        assert len(grn_processor.grns) == 3
        assert "3.50" in grn_processor.grns
        
        # Check dataset registry
        datasets = grn_processor.list_datasets()
        assert len(datasets) >= 1
        dataset_info = grn_processor.get_dataset_info("reference/test_grn")
        assert dataset_info is not None
        assert dataset_info["type"] == "grn_table"
        assert dataset_info["protein_count"] == 3
    
    def test_save_grn_table(self, grn_processor):
        """Test saving a GRN table."""
        # Load the test data
        grn_processor.load_grn_table("reference/test_grn")
        
        # Save it to a new file
        file_path = grn_processor.save_grn_table("reference/new_test_grn")
        
        # Check that the file exists
        assert os.path.exists(file_path)
        
        # Check that it's registered
        dataset_info = grn_processor.get_dataset_info("reference/new_test_grn")
        assert dataset_info is not None
        assert dataset_info["format"] == "csv"
        
        # Load the new file and check contents
        new_data = grn_processor.load_grn_table("reference/new_test_grn")
        assert len(new_data) == 3
        assert "3.50" in new_data.columns
    
    def test_filter_by_ids(self, grn_processor):
        """Test filtering by IDs."""
        # Load data
        grn_processor.load_grn_table("reference/test_grn")
        original_count = len(grn_processor.ids)
        
        # Filter by IDs
        grn_processor.filter_by_ids(["opsin1", "opsin3"])
        
        # Check results
        assert len(grn_processor.ids) == 2
        assert "opsin1" in grn_processor.ids
        assert "opsin3" in grn_processor.ids
        assert "opsin2" not in grn_processor.ids
    
    def test_apply_interval(self, grn_processor):
        """Test applying a GRN interval."""
        # Load data
        grn_processor.load_grn_table("reference/test_grn")
        original_grn_count = len(grn_processor.grns)
        
        # Apply interval
        grn_processor.apply_interval(["3.50", "7.49"])
        
        # Check results
        assert len(grn_processor.grns) == 2
        assert "3.50" in grn_processor.grns
        assert "7.49" in grn_processor.grns
        assert "6.48" not in grn_processor.grns
    
    def test_get_seq_dict(self, grn_processor):
        """Test getting sequences from a GRN table."""
        # Load data
        grn_processor.load_grn_table("reference/test_grn")
        
        # Get sequences
        seq_dict = grn_processor.get_seq_dict()
        
        # Check results
        assert len(seq_dict) == 3
        assert "opsin1" in seq_dict
        assert "opsin2" in seq_dict
        assert "opsin3" in seq_dict
        assert all(isinstance(seq, str) for seq in seq_dict.values())
        assert all(len(seq) == 3 for seq in seq_dict.values())  # Each sequence has 3 amino acids
        
        # Check specific sequences
        assert seq_dict["opsin1"] == "MYP"
        assert seq_dict["opsin2"] == "MYP"
        assert seq_dict["opsin3"] == "MYP"
    
    def test_get_grn_dict(self, grn_processor):
        """Test getting a GRN dictionary."""
        # Load data
        grn_processor.load_grn_table("reference/test_grn")
        
        # Get GRN dictionary
        grn_dict = grn_processor.get_grn_dict()
        
        # Check results
        assert len(grn_dict) == 3
        assert "opsin1" in grn_dict
        assert "opsin2" in grn_dict
        assert "opsin3" in grn_dict
        assert all(isinstance(grns, list) for grns in grn_dict.values())
        assert all(len(grns) == 3 for grns in grn_dict.values())  # Each protein has 3 GRNs
        
        # Check specific positions
        for protein_id in grn_dict:
            assert "3.50" in grn_dict[protein_id]
            assert "6.48" in grn_dict[protein_id]
            assert "7.49" in grn_dict[protein_id]
    
    def test_filter_occurances(self, grn_processor, sample_grn_table):
        """Test filtering by occurrences with a modified dataset."""
        # Create a dataset with missing values
        mod_table = sample_grn_table.copy()
        mod_table.loc[0, "6.48"] = "-"  # Make 6.48 occur in only 2/3 proteins
        mod_path = os.path.join(grn_processor.data_path, "reference", "mod_test_grn.csv")
        mod_table.to_csv(mod_path, index=False)
        
        # Load the modified data
        grn_processor.load_grn_table("reference/mod_test_grn")
        assert len(grn_processor.grns) == 3
        
        # Filter by occurrences with threshold=3 (only keep columns with 3 proteins)
        grn_processor.filter_data_by_occurances(3)
        assert len(grn_processor.grns) == 2
        assert "3.50" in grn_processor.grns
        assert "7.49" in grn_processor.grns
        assert "6.48" not in grn_processor.grns  # Should be filtered out (only 2 occurrences)
    
    def test_multiple_datasets(self, grn_processor, sample_grn_table):
        """Test working with multiple datasets."""
        # Create a second dataset
        mod_table = sample_grn_table.copy()
        mod_table["protein_id"] = ["opsin4", "opsin5", "opsin6"]  # Different proteins
        mod_path = os.path.join(grn_processor.data_path, "reference", "second_test_grn.csv")
        mod_table.to_csv(mod_path, index=False)
        
        # Load both datasets
        merged_data = grn_processor.load_and_merge_grn_tables(["reference/test_grn", "reference/second_test_grn"])
        
        # Check results
        assert len(grn_processor.ids) == 6
        assert "opsin1" in grn_processor.ids
        assert "opsin4" in grn_processor.ids
        assert grn_processor.dataset.startswith("merged_")
        
        # Check registry
        dataset_info = grn_processor.get_dataset_info(grn_processor.dataset)
        assert dataset_info is not None
        assert "merged_from" in dataset_info
        assert len(dataset_info["merged_from"]) == 2
    
    def test_reset_data(self, grn_processor):
        """Test resetting data after modifications."""
        # Load data
        grn_processor.load_grn_table("reference/test_grn")
        original_id_count = len(grn_processor.ids)
        original_grn_count = len(grn_processor.grns)
        
        # Modify the data
        grn_processor.filter_by_ids(["opsin1"])
        assert len(grn_processor.ids) == 1
        
        # Reset
        grn_processor.reset_data()
        
        # Check that data is restored
        assert len(grn_processor.ids) == original_id_count
        assert len(grn_processor.grns) == original_grn_count
    
    def test_backwards_compatibility(self, grn_processor, sample_grn_table):
        """Test backwards compatibility with legacy path parameter."""
        # Create a processor with legacy path parameter
        legacy_processor = GRNBaseProcessor(
            name="legacy_processor",
            path=os.path.join(grn_processor.data_root, "grn"),
            preload=False
        )
        
        # Check that it works
        assert legacy_processor.processor_data_dir == "grn"
        
        # Load data
        legacy_processor.load_grn_table("reference/test_grn")
        assert len(legacy_processor.ids) == 3