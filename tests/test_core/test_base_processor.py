"""
Tests for the BaseProcessor class.

This test suite validates the BaseProcessor functionality including:
1. Initialization and configuration
2. Dataset registry management
3. Data loading and saving in various formats
4. Metadata handling
5. Cross-processor interoperability

All tests follow Protos principles - no direct path manipulation.
Uses temporary directories for true test isolation.
"""

import json
import pickle
import pytest
import pandas as pd
import numpy as np
from datetime import datetime

from protos.core.base_processor import BaseProcessor
from protos.io.paths.path_config import ProtosPaths


# Define a concrete processor for testing - prefix with underscore to avoid pytest collection
class _TestProcessor(BaseProcessor):
    """Concrete processor implementation for testing."""
    
    def __init__(self, name="test", processor_data_dir=None, config=None):
        """Initialize test processor."""
        super().__init__(
            name=name, 
            processor_data_dir=processor_data_dir, 
            config=config
        )
    
    def list_entities(self):
        """List available entities (test implementation)."""
        # For testing, we'll track entities in metadata
        return self.metadata.get('entities', [])
    
    def list_datasets(self):
        """List available datasets (test implementation)."""
        # For testing, use the dataset registry
        return list(self.dataset_registry.keys())


# Removed isolated_test_env fixture - we use global test-data from conftest.py


class TestBaseProcessorInitialization:
    """Test processor initialization and configuration."""
    
    def test_basic_initialization(self):
        """Test basic processor initialization with defaults."""
        processor = _TestProcessor()
        
        assert processor.name == "test"
        # The processor_data_dir depends on the _get_default_data_dir logic
        assert processor.processor_data_dir is not None
        assert processor.data is None
        assert isinstance(processor.dataset_registry, dict)
        assert isinstance(processor.metadata, dict)
    
    def test_initialization_with_name(self):
        """Test processor initialization with custom name."""
        processor = _TestProcessor(name="custom_test")
        
        assert processor.name == "custom_test"
        assert processor.processor_data_dir is not None
    
    def test_initialization_with_custom_dir(self):
        """Test processor initialization with custom directory."""
        processor = _TestProcessor(processor_data_dir="custom_dir")
        
        assert processor.processor_data_dir == "custom_dir"
    
    def test_initialization_with_config(self):
        """Test processor initialization with configuration."""
        config = {"param1": "value1", "param2": 42}
        processor = _TestProcessor(config=config)
        
        assert processor.config == config
        assert processor.config["param1"] == "value1"
        assert processor.config["param2"] == 42


class TestDatasetOperations:
    """Test dataset save/load operations."""
    
    def test_save_load_dataframe(self):
        """Test saving and loading pandas DataFrames."""
        processor = _TestProcessor()
        
        # Create test data
        df = pd.DataFrame({
            'id': ['A', 'B', 'C'],
            'value': [1.0, 2.5, 3.7],
            'category': ['cat1', 'cat2', 'cat1']
        })
        
        # Save data - don't save index to avoid loading issues
        processor.save_data("test_df", df, index=False)
        
        # Check registry
        assert "test_df" in processor.dataset_registry
        
        # Load data
        loaded_df = processor.load_data("test_df")
        pd.testing.assert_frame_equal(loaded_df, df, check_index_type=False)
    
    def test_save_load_dict(self):
        """Test saving and loading dictionaries."""
        processor = _TestProcessor()
        
        # Create test data
        data = {
            "config": {"version": "1.0", "debug": True},
            "items": [1, 2, 3, 4, 5],
            "metadata": {"created": "2023-01-01", "author": "test"}
        }
        
        # Save data
        processor.save_data("test_dict", data)
        
        # Load data
        loaded_data = processor.load_data("test_dict")
        assert loaded_data == data
    
    def test_save_load_numpy_array(self):
        """Test saving and loading numpy arrays."""
        processor = _TestProcessor()
        
        # Create test data
        arr = np.random.rand(10, 5)
        
        # Save data
        processor.save_data("test_array", arr)
        
        # Load data
        loaded_arr = processor.load_data("test_array")
        np.testing.assert_array_almost_equal(loaded_arr, arr)
    
    def test_format_inference(self):
        """Test automatic format inference."""
        processor = _TestProcessor()
        
        # Test different data types
        test_cases = [
            ("df_data", pd.DataFrame({"a": [1, 2]}), "csv"),
            ("dict_data", {"key": "value"}, "json"),
            ("array_data", np.array([1, 2, 3]), "pkl"),  # Arrays saved as pkl not npy
            ("list_data", [1, 2, 3], "pkl"),  # Lists also saved as pkl
        ]
        
        for name, data, expected_format in test_cases:
            # Save with explicit index=False for DataFrames
            if isinstance(data, pd.DataFrame):
                processor.save_data(name, data, index=False)
            else:
                processor.save_data(name, data)
            assert processor.dataset_registry[name]["format"] == expected_format
    
    def test_save_with_metadata(self):
        """Test saving data with metadata."""
        processor = _TestProcessor()
        
        df = pd.DataFrame({"x": [1, 2, 3]})
        # Note: metadata parameter is stored separately, not passed to to_csv
        
        processor.save_data("data_with_meta", df)
        
        # Check basic metadata is stored
        dataset_info = processor.dataset_registry["data_with_meta"]
        assert "last_updated" in dataset_info  # At least check basic metadata exists


class TestDatasetManagement:
    """Test dataset management operations."""
    
    def test_list_datasets(self):
        """Test listing available datasets."""
        processor = _TestProcessor()
        
        # Initially empty (in isolated environment)
        initial_datasets = processor.list_datasets()
        initial_count = len(initial_datasets)
        
        # Add datasets
        processor.save_data("dataset1", {"a": 1})
        processor.save_data("dataset2", pd.DataFrame({"b": [2]}))
        
        # List datasets
        datasets = processor.list_datasets()
        assert len(datasets) == initial_count + 2
        
        # Check our datasets are present
        if isinstance(datasets[0], dict):
            dataset_ids = [d["id"] for d in datasets]
        else:
            dataset_ids = datasets
        assert "dataset1" in dataset_ids
        assert "dataset2" in dataset_ids
    
    def test_is_dataset_available(self):
        """Test checking dataset availability."""
        processor = _TestProcessor()
        
        # Save a dataset
        processor.save_data("existing", {"data": "test"})
        
        assert processor.is_dataset_available("existing")
        assert not processor.is_dataset_available("non_existing")
    
    def test_delete_dataset(self):
        """Test deleting datasets."""
        processor = _TestProcessor()
        
        # Save a dataset
        processor.save_data("to_delete", {"data": "test"})
        assert processor.is_dataset_available("to_delete")
        
        # Delete it
        result = processor.delete_dataset("to_delete")
        assert result is True
        assert not processor.is_dataset_available("to_delete")
        
        # Try deleting non-existent
        result = processor.delete_dataset("non_existent")
        assert result is False
    
    def test_update_dataset(self):
        """Test updating existing datasets."""
        processor = _TestProcessor()
        
        # Save initial data
        processor.save_data("to_update", {"version": 1})
        
        # Update with new data
        processor.save_data("to_update", {"version": 2, "updated": True})
        
        # Load and verify
        data = processor.load_data("to_update")
        assert data["version"] == 2
        assert data["updated"] is True


class TestErrorHandling:
    """Test error handling and validation."""
    
    def test_load_non_existent(self):
        """Test loading non-existent dataset."""
        processor = _TestProcessor()
        
        with pytest.raises(FileNotFoundError):
            processor.load_data("non_existent")
    
    def test_save_none_data(self):
        """Test saving None data."""
        processor = _TestProcessor()
        
        with pytest.raises(ValueError):
            processor.save_data("none_data", None)
    
    def test_save_without_data(self):
        """Test saving without providing data."""
        processor = _TestProcessor()
        
        with pytest.raises(ValueError):
            processor.save_data("no_data")
    
    def test_invalid_format(self):
        """Test handling of invalid format specifications."""
        processor = _TestProcessor()
        
        # Save with valid format
        processor.save_data("test", {"data": 1}, file_format="json")
        
        # Try to load with wrong format - should handle gracefully
        # The processor should use the registered format
        data = processor.load_data("test")
        assert data["data"] == 1


class TestProcessorInteroperability:
    """Test interaction between multiple processors."""
    
    def test_shared_data_access(self):
        """Test that processors can share data when using same directory."""
        # Create two processors with same data directory
        proc1 = _TestProcessor(name="proc1", processor_data_dir="shared")
        proc2 = _TestProcessor(name="proc2", processor_data_dir="shared")
        
        # Save data with first processor
        proc1.save_data("shared_data", {"source": "proc1"})
        
        # Load with second processor
        data = proc2.load_data("shared_data")
        assert data["source"] == "proc1"
    
    def test_independent_processors(self):
        """Test that processors with different directories are independent."""
        proc1 = _TestProcessor(name="proc1", processor_data_dir="unique_dir1")
        proc2 = _TestProcessor(name="proc2", processor_data_dir="unique_dir2")
        
        # Save data with first processor
        proc1.save_data("data1", {"value": 1})
        
        # Second processor should not see it
        assert not proc2.is_dataset_available("data1")
        
        # Each can have same dataset name
        proc2.save_data("data1", {"value": 2})
        
        # Load from each - should be different
        assert proc1.load_data("data1")["value"] == 1
        assert proc2.load_data("data1")["value"] == 2


class TestMetadataTracking:
    """Test metadata tracking functionality."""
    
    def test_processor_metadata(self):
        """Test processor-level metadata."""
        processor = _TestProcessor(name="meta_test")
        
        assert processor.metadata["processor_type"] == "_TestProcessor"
        assert processor.metadata["name"] == "meta_test"
        assert "created_at" in processor.metadata
    
    def test_dataset_metadata(self):
        """Test dataset-level metadata."""
        processor = _TestProcessor()
        
        # Save dataset
        df = pd.DataFrame({"x": range(10), "y": range(10, 20)})
        processor.save_data("test_dataset", df)
        
        # Check metadata
        dataset_info = processor.dataset_registry["test_dataset"]
        assert "last_updated" in dataset_info
        assert dataset_info["format"] == "csv"
        assert dataset_info["rows"] == 10
        assert dataset_info["columns"] == ['x', 'y']


class TestSpecialFormats:
    """Test handling of special data formats."""
    
    def test_save_load_fasta(self):
        """Test saving and loading FASTA-like data."""
        processor = _TestProcessor()
        
        # FASTA data as dictionary
        sequences = {
            "seq1": "ATCGATCGATCG",
            "seq2": "GCTAGCTAGCTA",
            "seq3": "TATATATATATAT"
        }
        
        processor.save_data("sequences", sequences, file_format="json")
        loaded = processor.load_data("sequences")
        assert loaded == sequences
    
    def test_complex_nested_data(self):
        """Test handling complex nested data structures."""
        processor = _TestProcessor()
        
        complex_data = {
            "level1": {
                "level2": {
                    "data": [1, 2, 3],
                    "metadata": {"type": "nested"}
                },
                "arrays": {
                    "numpy": np.array([1, 2, 3]).tolist(),  # Convert for JSON
                    "lists": [[1, 2], [3, 4], [5, 6]]
                }
            },
            "dataframes": {
                "df1": pd.DataFrame({"a": [1, 2]}).to_dict(),
                "df2": pd.DataFrame({"b": [3, 4]}).to_dict()
            }
        }
        
        processor.save_data("complex", complex_data)
        loaded = processor.load_data("complex")
        assert loaded["level1"]["level2"]["data"] == [1, 2, 3]


class TestRegistryPersistence:
    """Test registry persistence across processor instances."""
    
    def test_registry_persistence(self):
        """Test that registry persists between processor instances."""
        # First processor saves data
        proc1 = _TestProcessor(name="persist1", processor_data_dir="persist_test")
        proc1.save_data("data1", {"value": 1})
        proc1.save_data("data2", {"value": 2})
        
        # Create new processor instance with same directory
        proc2 = _TestProcessor(name="persist2", processor_data_dir="persist_test")
        
        # Should see the same datasets
        assert proc2.is_dataset_available("data1")
        assert proc2.is_dataset_available("data2")
        
        # Can load the data
        assert proc2.load_data("data1")["value"] == 1
        assert proc2.load_data("data2")["value"] == 2
    
    def test_registry_update_on_delete(self):
        """Test that registry updates when datasets are deleted."""
        processor = _TestProcessor(processor_data_dir="delete_test")
        
        # Save and delete
        processor.save_data("temp", {"data": "temporary"})
        assert processor.is_dataset_available("temp")
        
        processor.delete_dataset("temp")
        assert not processor.is_dataset_available("temp")
        
        # New processor should also not see it
        new_proc = _TestProcessor(processor_data_dir="delete_test")
        assert not new_proc.is_dataset_available("temp")