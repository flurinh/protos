"""
Tests for the BaseProcessor class.

This test suite validates the BaseProcessor functionality including:
1. Initialization and configuration
2. Path resolution with the new unified data system
3. Dataset registry management
4. Data loading and saving in various formats
5. Metadata handling
6. Cross-processor interoperability
"""

import os
import json
import pickle
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
from unittest.mock import patch, MagicMock

from protos.core.base_processor import BaseProcessor
from protos.io.paths.path_config import ProtosPaths, DataSource


# Define a concrete processor for testing - prefix with underscore to avoid pytest collection
class _TestProcessor(BaseProcessor):
    """Concrete processor implementation for testing."""
    
    def __init__(self, name="test", processor_data_dir=None, config=None):
        """Initialize test processor."""
        # Simply pass all parameters to parent
        super().__init__(
            name=name, 
            processor_data_dir=processor_data_dir, 
            config=config
        )


class TestBaseProcessorInitialization:
    """Test processor initialization and configuration."""
    
    def test_basic_initialization(self, test_data_root):
        """Test basic processor initialization with defaults."""
        processor = _TestProcessor()
        
        assert processor.name == "test"
        assert processor.data_root == str(test_data_root)
        assert processor.processor_data_dir == "__test"  # _TestProcessor -> __test via snake_case conversion
        assert processor.data_path == os.path.join(str(test_data_root), "test")  # But maps to "test" directory
        assert processor.data is None
        assert isinstance(processor.dataset_registry, dict)
        assert processor.metadata["processor_type"] == "_TestProcessor"
        assert processor.metadata["name"] == "test"
        
    def test_custom_initialization(self, test_data_root):
        """Test processor initialization with custom parameters."""
        custom_config = {"param1": "value1", "param2": "value2"}
        processor = _TestProcessor(
            name="custom",
            processor_data_dir="custom_dir",
            config=custom_config
        )
        
        assert processor.name == "custom"
        assert processor.data_root == str(test_data_root)
        assert processor.processor_data_dir == "custom_dir"
        assert processor.data_path == os.path.join(str(test_data_root), "custom_dir")
        assert processor.config == custom_config
        assert os.path.exists(processor.data_path)
        
    def test_processor_type_inference(self, test_data_root):
        """Test automatic processor type and directory name inference."""
        # Test simple processor
        class SimpleProcessor(BaseProcessor):
            pass
        
        processor = SimpleProcessor(name="test")
        assert processor.processor_data_dir == "simple"
        assert processor.metadata["processor_type"] == "SimpleProcessor"
        
        # Test complex processor name
        class ComplexProcessorWithLongName(BaseProcessor):
            pass
        
        processor = ComplexProcessorWithLongName(name="test")
        assert processor.processor_data_dir == "complex_processor_with_long_name"
        
    def test_path_initialization_modes(self, test_data_root):
        """Test different path initialization modes."""
        # Test that processor uses global configuration
        processor = _TestProcessor()
        assert os.path.isabs(processor.data_path)
        assert processor.data_root == str(test_data_root)


class TestDatasetRegistry:
    """Test dataset registry functionality."""
    
    def test_empty_registry(self, test_data_root):
        """Test that registry starts empty."""
        # Clean up any existing registry from previous tests
        registry_path = test_data_root / "test" / "registry.json"
        if registry_path.exists():
            registry_path.unlink()
        
        processor = _TestProcessor()
        assert processor.dataset_registry == {}
        
    def test_register_dataset(self, test_data_root):
        """Test registering a dataset in the registry."""
        processor = _TestProcessor()
        
        # Register a dataset
        metadata = {"type": "test", "columns": ["id", "name", "value"]}
        file_path = os.path.join(processor.data_path, "test_dataset.csv")
        processor._register_dataset("test_dataset", metadata, file_path)
        
        # Verify registration
        assert "test_dataset" in processor.dataset_registry
        registry_entry = processor.dataset_registry["test_dataset"]
        assert registry_entry["filename"] == "test_dataset.csv"
        assert registry_entry["type"] == "test"
        assert "last_updated" in registry_entry
        assert "columns" in registry_entry
        
        # Verify registry file was saved
        registry_path = os.path.join(processor.data_path, "registry.json")
        assert os.path.exists(registry_path)
        
        # Verify registry can be loaded
        with open(registry_path, 'r') as f:
            saved_registry = json.load(f)
        assert "test_dataset" in saved_registry
        
    def test_load_registry_from_file(self, test_data_root):
        """Test loading registry from file."""
        processor = _TestProcessor()
        
        # Create a registry file manually
        registry_data = {
            "existing_dataset": {
                "filename": "existing.csv",
                "type": "test",
                "last_updated": datetime.now().isoformat()
            }
        }
        registry_path = os.path.join(processor.data_path, "registry.json")
        with open(registry_path, 'w') as f:
            json.dump(registry_data, f)
        
        # Load registry
        loaded_registry = processor._load_dataset_registry()
        assert "existing_dataset" in loaded_registry
        assert loaded_registry["existing_dataset"]["type"] == "test"
        
    def test_registry_persistence(self, test_data_root):
        """Test that registry persists across processor instances."""
        # First processor saves dataset
        processor1 = _TestProcessor(processor_data_dir="persist_test")
        df = pd.DataFrame({"a": [1, 2, 3]})
        processor1.save_data("persistent_dataset", df, file_format="csv")
        
        # Second processor should see the dataset
        processor2 = _TestProcessor(processor_data_dir="persist_test")
        assert "persistent_dataset" in processor2.dataset_registry
        assert processor2.is_dataset_available("persistent_dataset")


class TestDataOperations:
    """Test data loading and saving operations."""
    
    def test_save_load_csv(self, test_data_root):
        """Test saving and loading CSV files."""
        processor = _TestProcessor()
        
        # Create test DataFrame
        df = pd.DataFrame({
            "id": [1, 2, 3],
            "name": ["a", "b", "c"],
            "value": [10.5, 20.5, 30.5]
        })
        
        # Save DataFrame
        file_path = processor.save_data("test_csv", df, file_format="csv", index=False)
        assert os.path.exists(file_path)
        assert file_path.endswith(".csv")
        
        # Load DataFrame
        loaded_df = processor.load_data("test_csv", file_format="csv")
        pd.testing.assert_frame_equal(loaded_df, df)
        
        # Verify registry entry
        assert "test_csv" in processor.dataset_registry
        assert processor.dataset_registry["test_csv"]["format"] == "csv"
        assert processor.dataset_registry["test_csv"]["rows"] == 3
        assert set(processor.dataset_registry["test_csv"]["columns"]) == set(df.columns)
        
    def test_save_load_json(self, test_data_root):
        """Test saving and loading JSON files."""
        processor = _TestProcessor()
        
        # Create test data
        data = {
            "metadata": {"version": "1.0", "type": "test"},
            "items": [{"id": 1, "value": "a"}, {"id": 2, "value": "b"}]
        }
        
        # Save JSON
        file_path = processor.save_data("test_json", data, file_format="json")
        assert os.path.exists(file_path)
        assert file_path.endswith(".json")
        
        # Load JSON
        loaded_data = processor.load_data("test_json", file_format="json")
        assert loaded_data == data
        
        # Verify registry entry
        assert "test_json" in processor.dataset_registry
        assert processor.dataset_registry["test_json"]["format"] == "json"
        
    def test_save_load_pickle(self, test_data_root):
        """Test saving and loading pickle files."""
        processor = _TestProcessor()
        
        # Create test data with numpy arrays
        data = {
            "array1": np.random.rand(10, 5),
            "array2": np.random.rand(15, 5),
            "metadata": {"type": "embeddings", "dimensions": 5}
        }
        
        # Save pickle
        file_path = processor.save_data("test_pickle", data, file_format="pkl")
        assert os.path.exists(file_path)
        assert file_path.endswith(".pkl")
        
        # Load pickle
        loaded_data = processor.load_data("test_pickle", file_format="pkl")
        assert list(loaded_data.keys()) == list(data.keys())
        np.testing.assert_array_equal(loaded_data["array1"], data["array1"])
        np.testing.assert_array_equal(loaded_data["array2"], data["array2"])
        assert loaded_data["metadata"] == data["metadata"]
        
    def test_format_inference(self, test_data_root):
        """Test automatic format inference from file extension."""
        processor = _TestProcessor()
        
        # Save files with different formats
        df = pd.DataFrame({"a": [1, 2, 3]})
        processor.save_data("test_infer_csv", df, file_format="csv")
        processor.save_data("test_infer_json", {"a": 1}, file_format="json")
        processor.save_data("test_infer_pkl", {"b": 2}, file_format="pkl")
        
        # Load without specifying format
        csv_data = processor.load_data("test_infer_csv")
        json_data = processor.load_data("test_infer_json")
        pkl_data = processor.load_data("test_infer_pkl")
        
        # Verify correct format was inferred
        assert isinstance(csv_data, pd.DataFrame)
        assert isinstance(json_data, dict)
        assert isinstance(pkl_data, dict)
        
    def test_save_without_data(self, test_data_root):
        """Test error when saving without data."""
        processor = _TestProcessor()
        
        with pytest.raises(ValueError):
            processor.save_data("test_no_data")
            
    def test_load_nonexistent(self, test_data_root):
        """Test error when loading nonexistent dataset."""
        processor = _TestProcessor()
        
        with pytest.raises(FileNotFoundError):
            processor.load_data("nonexistent_dataset")
            
    def test_unsupported_format(self, test_data_root):
        """Test error with unsupported file format."""
        processor = _TestProcessor()
        
        with pytest.raises(ValueError):
            processor.save_data("test_unsupported", {"a": 1}, file_format="xyz")


class TestDatasetManagement:
    """Test dataset management operations."""
    
    def test_list_datasets(self, test_data_root):
        """Test listing available datasets."""
        # Create a processor with a unique name to avoid registry conflicts
        processor = _TestProcessor(name="list_test", processor_data_dir="list_test")
        
        # Clear any existing registry entries
        processor.dataset_registry = {}
        processor._save_dataset_registry()
        
        # Initially empty
        datasets = processor.list_datasets()
        assert datasets == []
        
        # Add some datasets
        processor.save_data("dataset1", pd.DataFrame({"a": [1, 2]}))
        processor.save_data("dataset2", {"b": 2}, file_format="json")
        processor.save_data("dataset3", {"c": 3}, file_format="pkl")
        
        # List datasets
        datasets = processor.list_datasets()
        assert len(datasets) == 3
        
        # Verify dataset info
        dataset_ids = [d["id"] for d in datasets]
        assert set(dataset_ids) == {"dataset1", "dataset2", "dataset3"}
        
        # Check individual dataset info
        dataset1_info = next(d for d in datasets if d["id"] == "dataset1")
        assert dataset1_info["format"] == "csv"
        assert "saved_at" in dataset1_info
        
    def test_get_dataset_info(self, test_data_root):
        """Test getting information about a specific dataset."""
        processor = _TestProcessor()
        
        # Save a dataset
        df = pd.DataFrame({"id": [1, 2, 3], "value": [10, 20, 30]})
        processor.save_data("info_test", df)
        
        # Get dataset info
        info = processor.get_dataset_info("info_test")
        assert info["id"] == "info_test"
        assert info["format"] == "csv"
        assert info["rows"] == 3
        assert set(info["columns"]) == {"id", "value"}
        assert "saved_at" in info
        
        # Test nonexistent dataset
        assert processor.get_dataset_info("nonexistent") is None
        
    def test_delete_dataset(self, test_data_root):
        """Test deleting datasets."""
        processor = _TestProcessor()
        
        # Save a dataset
        processor.save_data("to_delete", pd.DataFrame({"a": [1, 2, 3]}))
        file_path = processor._get_dataset_path("to_delete", ".csv")
        
        # Verify it exists
        assert os.path.exists(file_path)
        assert processor.is_dataset_available("to_delete")
        
        # Delete it
        result = processor.delete_dataset("to_delete")
        assert result is True
        
        # Verify deletion
        assert not os.path.exists(file_path)
        assert not processor.is_dataset_available("to_delete")
        assert "to_delete" not in processor.dataset_registry
        
        # Try deleting nonexistent
        result = processor.delete_dataset("nonexistent")
        assert result is False
        
    def test_is_dataset_available(self, test_data_root):
        """Test checking dataset availability."""
        processor = _TestProcessor()
        
        # Save a dataset
        processor.save_data("available", pd.DataFrame({"a": [1]}))
        
        # Check availability
        assert processor.is_dataset_available("available") is True
        assert processor.is_dataset_available("not_available") is False
        
        # Test unregistered file detection
        unregistered_path = os.path.join(processor.data_path, "unregistered.csv")
        pd.DataFrame({"b": [2]}).to_csv(unregistered_path, index=False)
        
        # Should detect unregistered file
        assert processor.is_dataset_available("unregistered") is True


class _TestProcessorInteroperability:
    """Test interoperability between different processor types."""
    
    def test_cross_processor_data_sharing(self, test_data_root):
        """Test sharing data between different processor types."""
        # Create processors for different types
        grn_processor = _TestProcessor(
            name="grn_test", 
            processor_data_dir="grn"
        )
        struct_processor = _TestProcessor(
            name="struct_test", 
            processor_data_dir="structure"
        )
        
        # Save data in GRN processor
        grn_data = pd.DataFrame({
            "pdb_id": ["1ABC", "2DEF"],
            "grn": ["3.50", "3.51"],
            "res_id": [100, 101]
        })
        grn_processor.save_data("grn_mapping", grn_data)
        
        # Structure processor should be able to reference GRN data
        # by constructing the correct path
        grn_path = os.path.join(test_data_root, "grn", "grn_mapping.csv")
        assert os.path.exists(grn_path)
        
        # Load GRN data from structure processor using absolute path
        loaded_grn = pd.read_csv(grn_path)
        pd.testing.assert_frame_equal(loaded_grn, grn_data)
        
    def test_processor_isolation(self, test_data_root):
        """Test that processors maintain separate data directories."""
        # Create multiple processors
        processors = [
            _TestProcessor(name="proc1", processor_data_dir="proc1"),
            _TestProcessor(name="proc2", processor_data_dir="proc2"),
            _TestProcessor(name="proc3", processor_data_dir="proc3")
        ]
        
        # Save datasets with same name in each processor
        for i, proc in enumerate(processors):
            proc.save_data("common_name", pd.DataFrame({"value": [i]}))
        
        # Verify each processor has its own dataset
        for i, proc in enumerate(processors):
            data = proc.load_data("common_name")
            assert data["value"].iloc[0] == i
            
        # Verify separate directories
        assert all(os.path.exists(proc.data_path) for proc in processors)
        assert len(set(proc.data_path for proc in processors)) == 3


class TestMetadataHandling:
    """Test metadata handling in processors."""
    
    def test_processor_metadata(self, test_data_root):
        """Test processor metadata initialization and updates."""
        processor = _TestProcessor(
            name="metadata_test",
            config={"custom_param": "value"}
        )
        
        # Check initial metadata
        assert processor.metadata["processor_type"] == "_TestProcessor"
        assert processor.metadata["name"] == "metadata_test"
        assert "created_at" in processor.metadata
        
        # Update metadata
        processor.metadata["custom_field"] = "custom_value"
        processor.metadata["version"] = "1.0.0"
        
        assert processor.metadata["custom_field"] == "custom_value"
        assert processor.metadata["version"] == "1.0.0"
        
    def test_dataset_metadata_in_registry(self, test_data_root):
        """Test that dataset metadata is properly stored in registry."""
        processor = _TestProcessor()
        
        # Save dataset
        df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
        processor.save_data("metadata_dataset", df)
        
        # Check registry contains basic metadata
        registry_entry = processor.dataset_registry["metadata_dataset"]
        assert "filename" in registry_entry
        assert "format" in registry_entry
        assert "last_updated" in registry_entry
        assert registry_entry["format"] == "csv"
        assert registry_entry["rows"] == 3


class TestPathHandling:
    """Test path handling with the new unified data system."""
    
    def test_unified_data_root(self, test_data_root):
        """Test that processor uses unified data root correctly."""
        processor = _TestProcessor()
        
        # Verify paths use the single data root
        assert processor.data_root == str(test_data_root)
        assert processor.data_path.startswith(str(test_data_root))
        assert os.path.exists(processor.data_path)
        
    def test_backward_compatibility_warnings(self, test_data_root, caplog):
        """Test that backward compatibility is maintained with warnings."""
        # Create ProtosPaths with deprecated parameters
        paths = ProtosPaths(
            user_data_root=str(test_data_root),
            ref_data_root="/some/other/path"
        )
        
        # Should use user_data_root and warn about ref_data_root
        assert paths.data_root == str(test_data_root)
        assert "deprecated" in caplog.text.lower()
        
    def test_relative_to_absolute_path_conversion(self):
        """Test that paths are always absolute."""
        processor = _TestProcessor()
        
        # Should always be absolute paths from global configuration
        assert os.path.isabs(processor.data_root)
        assert os.path.isabs(processor.data_path)
        
    def test_environment_variable_expansion(self):
        """Test that global configuration can be set via environment variable."""
        # Save current global setting
        original = ProtosPaths.get_global_data_root()
        
        # Clear global setting so env var will be used
        ProtosPaths.set_data_root(None)
        
        # Set a test environment variable
        os.environ["PROTOS_DATA_ROOT"] = "/tmp/test_protos_env"
        
        try:
            processor = _TestProcessor()
            # Should use environment variable
            assert processor.data_root == "/tmp/test_protos_env"
        finally:
            # Clean up
            if "PROTOS_DATA_ROOT" in os.environ:
                del os.environ["PROTOS_DATA_ROOT"]
            # Restore original setting
            ProtosPaths.set_data_root(original)