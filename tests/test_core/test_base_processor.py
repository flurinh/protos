"""
Tests for the BaseProcessor class.

These tests cover:
1. Initialization and configuration
2. Path resolution
3. Dataset registry management
4. Data loading and saving
5. Helper methods
"""

import os
import pandas as pd
import numpy as np
import json
import pickle
import tempfile
import shutil
import pytest
from unittest.mock import patch, MagicMock
from datetime import datetime

from protos.core.base_processor import BaseProcessor


# Define a concrete processor for testing
class TestProcessor(BaseProcessor):
    """Concrete processor implementation for testing."""
    
    def __init__(self, name="test", data_root=None, processor_data_dir=None, config=None):
        """Initialize test processor."""
        super().__init__(name=name, data_root=data_root, processor_data_dir=processor_data_dir, config=config)
    
    # We will use the BaseProcessor's load_data and save_data methods directly, no need to override


# Fixtures for testing
@pytest.fixture
def temp_dir():
    """Temporary directory for file operations within the project structure."""
    # Create a temporary directory inside protos/tests/data
    base_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "temp")
    os.makedirs(base_dir, exist_ok=True)
    
    # Create a unique test directory
    dir_name = f"test_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{os.getpid()}"
    dir_path = os.path.join(base_dir, dir_name)
    os.makedirs(dir_path, exist_ok=True)
    
    yield dir_path
    
    # Clean up
    shutil.rmtree(dir_path, ignore_errors=True)


@pytest.fixture
def sample_df():
    """Sample DataFrame for testing."""
    return pd.DataFrame({
        "id": [1, 2, 3],
        "name": ["a", "b", "c"],
        "value": [10, 20, 30]
    })


@pytest.fixture
def test_data_root(temp_dir):
    """Create a realistic data directory structure for testing."""
    data_dir = os.path.join(temp_dir, "data")
    
    # Create main directories
    for subdir in ["grn/reference", "structures/mmcif", "structures/temp", 
                  "embeddings/esm2", "embeddings/ankh", 
                  "properties", "graphs/binding_domain"]:
        os.makedirs(os.path.join(data_dir, subdir), exist_ok=True)
    
    # Create temp dirs
    os.makedirs(os.path.join(temp_dir, "temp/tool_workspace"), exist_ok=True)
    os.makedirs(os.path.join(temp_dir, "temp/test_outputs"), exist_ok=True)
    
    return data_dir


@pytest.fixture
def processor_grn(test_data_root):
    """Create a GRN processor with appropriate data directory."""
    return TestProcessor(
        name="test_grn_processor", 
        data_root=test_data_root,
        processor_data_dir="grn"
    )

@pytest.fixture
def processor_struct(test_data_root):
    """Create a Structure processor with appropriate data directory."""
    return TestProcessor(
        name="test_struct_processor", 
        data_root=test_data_root,
        processor_data_dir="structures"
    )

@pytest.fixture
def processor_emb(test_data_root):
    """Create an Embedding processor with appropriate data directory."""
    return TestProcessor(
        name="test_emb_processor", 
        data_root=test_data_root,
        processor_data_dir="embeddings"
    )

@pytest.fixture
def processor_prop(test_data_root):
    """Create a Property processor with appropriate data directory."""
    return TestProcessor(
        name="test_prop_processor", 
        data_root=test_data_root,
        processor_data_dir="properties"
    )

@pytest.fixture
def processor(test_data_root):
    """Generic processor for basic tests."""
    return TestProcessor(name="test_processor", data_root=test_data_root)


# Tests
class TestBaseProcessor:
    """Tests for the BaseProcessor class."""
    
    def test_initialization(self, temp_dir):
        """Test processor initialization."""
        # Test with default parameters
        processor = TestProcessor()
        assert processor.name == "test"
        assert processor.data_root == "data"
        assert processor.processor_data_dir == "test"
        assert processor.data_path == os.path.join("data", "test")
        assert processor.data is None
        assert isinstance(processor.dataset_registry, dict)
        assert processor.metadata["processor_type"] == "TestProcessor"
        
        # Test with custom parameters
        custom_config = {"param1": "value1", "param2": "value2"}
        processor = TestProcessor(
            name="custom",
            data_root=temp_dir,
            processor_data_dir="custom_dir",
            config=custom_config
        )
        assert processor.name == "custom"
        assert processor.data_root == temp_dir
        assert processor.processor_data_dir == "custom_dir"
        assert processor.data_path == os.path.join(temp_dir, "custom_dir")
        assert processor.config == custom_config
        
        # Check that data directory was created
        assert os.path.exists(os.path.join(temp_dir, "custom_dir"))
    
    def test_default_data_dir(self):
        """Test the default data directory name generation."""
        # Test processing with different names
        class SimpleProcessor(BaseProcessor):
            pass
        
        class ComplexProcessorWithLongName(BaseProcessor):
            pass
        
        simple_proc = SimpleProcessor(name="test")
        assert simple_proc.processor_data_dir == "simple"
        
        complex_proc = ComplexProcessorWithLongName(name="test")
        assert complex_proc.processor_data_dir == "complex_processor_with_long_name"
    
    def test_dataset_registry(self, processor):
        """Test dataset registry functionality."""
        # Initially empty
        assert processor.dataset_registry == {}
        
        # Register a dataset
        metadata = {"type": "test", "columns": ["id", "name", "value"]}
        file_path = os.path.join(processor.data_path, "test_dataset.csv")
        processor._register_dataset("test_dataset", metadata, file_path)
        
        # Check registry
        assert "test_dataset" in processor.dataset_registry
        assert processor.dataset_registry["test_dataset"]["filename"] == "test_dataset.csv"
        assert "last_updated" in processor.dataset_registry["test_dataset"]
        assert processor.dataset_registry["test_dataset"]["type"] == "test"
        
        # Check that registry file was saved
        registry_path = os.path.join(processor.data_path, "registry.json")
        assert os.path.exists(registry_path)
        
        # Load registry from file
        processor.dataset_registry = {}  # Clear in-memory registry
        loaded_registry = processor._load_dataset_registry()
        assert "test_dataset" in loaded_registry
        assert loaded_registry["test_dataset"]["type"] == "test"
    
    def test_load_save_csv(self, processor_prop, sample_df):
        """Test loading and saving CSV files to the properties directory."""
        # Save DataFrame with index=False to avoid extra column
        file_path = processor_prop.save_data("test_properties", sample_df, file_format="csv", index=False)
        
        # Verify file exists in correct location
        assert os.path.exists(file_path)
        assert file_path.endswith(".csv")
        assert "properties" in file_path
        
        # Load DataFrame
        loaded_df = processor_prop.load_data("test_properties", file_format="csv")
        
        # Verify loaded data
        pd.testing.assert_frame_equal(loaded_df, sample_df)
        
        # Check that dataset is in registry
        assert "test_properties" in processor_prop.dataset_registry
        assert processor_prop.dataset_registry["test_properties"]["format"] == "csv"
        assert processor_prop.dataset_registry["test_properties"]["rows"] == 3
        assert set(processor_prop.dataset_registry["test_properties"]["columns"]) == set(["id", "name", "value"])
    
    def test_load_save_pickle(self, processor_emb):
        """Test loading and saving pickle files to the embeddings directory."""
        # Create test data similar to embeddings
        test_data = {
            "protein1": np.random.rand(10, 5),
            "protein2": np.random.rand(15, 5),
            "metadata": {"embedding_type": "test", "dimensions": 5}
        }
        
        # Save data
        file_path = processor_emb.save_data("test_embeddings", test_data, file_format="pkl")
        
        # Verify file exists in correct location
        assert os.path.exists(file_path)
        assert file_path.endswith(".pkl")
        assert "embeddings" in file_path
        
        # Load data
        loaded_data = processor_emb.load_data("test_embeddings", file_format="pkl")
        
        # Verify loaded data
        assert list(loaded_data.keys()) == list(test_data.keys())
        np.testing.assert_array_equal(loaded_data["protein1"], test_data["protein1"])
        np.testing.assert_array_equal(loaded_data["protein2"], test_data["protein2"])
        assert loaded_data["metadata"] == test_data["metadata"]
    
    def test_load_nonexistent_file(self, processor):
        """Test loading a nonexistent file."""
        with pytest.raises(FileNotFoundError):
            processor.load_data("nonexistent_dataset")
    
    def test_save_without_data(self, processor):
        """Test saving without providing data."""
        with pytest.raises(ValueError):
            processor.save_data("test_dataset")
    
    def test_infer_format(self, processor_grn, sample_df):
        """Test format inference with files in the GRN directory."""
        # Save with different formats
        processor_grn.save_data("test_grn_csv", sample_df, file_format="csv", index=False)
        processor_grn.save_data("test_grn_json", {"a": 1, "b": 2}, file_format="json")
        
        # Try loading without specifying format
        csv_data = processor_grn.load_data("test_grn_csv")
        json_data = processor_grn.load_data("test_grn_json")
        
        # Verify file paths
        assert "grn" in processor_grn.metadata["file_path"]
        
        # Verify correct format was inferred
        assert isinstance(csv_data, pd.DataFrame)
        assert isinstance(json_data, dict)
    
    def test_dataset_info(self, processor_struct, sample_df):
        """Test getting dataset information for files in the structures directory."""
        # Save a dataset
        processor_struct.save_data("test_structure_info", sample_df, file_format="csv", index=False)
        
        # Get dataset info
        info = processor_struct.get_dataset_info("test_structure_info")
        
        # Verify info contains expected fields
        assert info["id"] == "test_structure_info"
        assert info["format"] == "csv"
        assert info["rows"] == 3
        assert "saved_at" in info
        
        # Verify directory is correct
        assert info.get("directory", "") == "structures"
        
        # Test nonexistent dataset
        assert processor_struct.get_dataset_info("nonexistent") is None
    
    def test_list_datasets_across_directories(self, processor_grn, processor_struct, processor_emb, sample_df):
        """Test listing datasets in their respective directories."""
        # Save datasets to different processor directories
        processor_grn.save_data("grn_dataset", sample_df, file_format="csv", index=False)
        processor_struct.save_data("structure_dataset", {"a": 1}, file_format="json")
        processor_emb.save_data("embedding_dataset", {"test": np.random.rand(10, 5)}, file_format="pkl")
        
        # List datasets from each processor
        grn_datasets = processor_grn.list_datasets()
        struct_datasets = processor_struct.list_datasets()
        emb_datasets = processor_emb.list_datasets()
        
        # Verify each processor has its own dataset
        assert len(grn_datasets) == 1
        assert len(struct_datasets) == 1
        assert len(emb_datasets) == 1
        
        # Check dataset IDs
        assert grn_datasets[0]["id"] == "grn_dataset"
        assert struct_datasets[0]["id"] == "structure_dataset"
        assert emb_datasets[0]["id"] == "embedding_dataset"
        
        # Verify registry contains correct directories
        assert grn_datasets[0].get("directory", "") == "grn"
        assert struct_datasets[0].get("directory", "") == "structures"
        assert emb_datasets[0].get("directory", "") == "embeddings"
    
    def test_delete_dataset(self, processor_prop, sample_df):
        """Test deleting a dataset from the properties directory."""
        # Save a dataset in properties folder
        processor_prop.save_data("property_to_delete", sample_df, index=False)
        
        # Verify it exists
        file_path = processor_prop._get_dataset_path("property_to_delete", ".csv")
        assert os.path.exists(file_path)
        assert "property_to_delete" in processor_prop.dataset_registry
        assert "properties" in file_path
        
        # Delete it
        result = processor_prop.delete_dataset("property_to_delete")
        
        # Verify deletion
        assert result is True
        assert not os.path.exists(file_path)
        assert "property_to_delete" not in processor_prop.dataset_registry
        
        # Try deleting nonexistent dataset
        result = processor_prop.delete_dataset("nonexistent")
        assert result is False
    
    def test_is_dataset_available(self, processor_struct, sample_df):
        """Test checking dataset availability in the structures directory."""
        # Save a dataset
        processor_struct.save_data("available_structure", sample_df, index=False)
        
        # Check availability
        assert processor_struct.is_dataset_available("available_structure") is True
        assert processor_struct.is_dataset_available("nonexistent_structure") is False
        
        # Create file without registry entry
        unregistered_path = os.path.join(processor_struct.data_path, "unregistered.csv")
        with open(unregistered_path, 'w') as f:
            f.write("test,data\n1,2\n")
        
        # Verify file location is in structures directory
        assert "structures" in unregistered_path
        
        # Should detect unregistered file
        assert processor_struct.is_dataset_available("unregistered") is True