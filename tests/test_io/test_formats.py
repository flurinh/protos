"""
Tests for the protos.io.formats module.

These tests cover:
1. File format handlers
2. Format registry functionality
3. Temporary file management
4. External tool helpers
"""

import os
import pandas as pd
import numpy as np
import json
import pickle
import tempfile
import shutil
import pytest
from typing import Dict, List, Any

from protos.io.formats import (
    FormatHandler, CSVHandler, PickleHandler, JSONHandler, FASTAHandler,
    GRNTableHandler, PropertyTableHandler, EmbeddingHandler, GraphHandler,
    FormatRegistry, format_registry, TempFileManager, ExternalToolHelper
)


# Fixtures for commonly used test data
@pytest.fixture
def sample_sequences() -> Dict[str, str]:
    """Sample sequences for FASTA tests."""
    return {
        "seq1": "ACGTACGTACGT",
        "seq2": "TGCATGCATGCA",
        "seq3": "AAAAGGGGCCCCTTTT"
    }


@pytest.fixture
def sample_grn_table() -> pd.DataFrame:
    """Sample GRN table for tests."""
    data = {
        "3.50": ["M125", "M123", "M127"],
        "6.48": ["Y261", "Y259", "Y263"],
        "7.49": ["P296", "P294", "P298"]
    }
    index = ["opsin1", "opsin2", "opsin3"]
    df = pd.DataFrame(data, index=index)
    df.index.name = "protein_id"
    return df


@pytest.fixture
def sample_property_table() -> pd.DataFrame:
    """Sample property table for tests."""
    data = {
        "protein_id": ["opsin1", "opsin2", "opsin3"],
        "activation": ["light", "light", "chemical"],
        "spectral_sensitivity": [500, 460, None],
        "ligand_binding": ["all-trans-retinal", "11-cis-retinal", "none"]
    }
    df = pd.DataFrame(data)
    df.metadata = {
        "description": "Test dataset",
        "version": "1.0",
        "columns": {
            "activation": {"type": "categorical"},
            "spectral_sensitivity": {"type": "numeric", "unit": "nm"},
            "ligand_binding": {"type": "categorical"}
        }
    }
    return df


@pytest.fixture
def sample_embeddings() -> Dict[str, np.ndarray]:
    """Sample embeddings for tests."""
    return {
        "opsin1": np.random.rand(12, 10),
        "opsin2": np.random.rand(12, 10),
        "opsin3": np.random.rand(15, 10)
    }


@pytest.fixture
def sample_graph() -> Dict[str, Dict]:
    """Sample graph for tests."""
    return {
        "protein_id": {
            "nodes": [
                {"id": 0, "residue": "ALA", "position": 1, "grn": "1.50", "coords": [1.0, 2.0, 3.0]},
                {"id": 1, "residue": "GLY", "position": 2, "grn": None, "coords": [4.0, 5.0, 6.0]},
                {"id": 2, "residue": "SER", "position": 3, "grn": "1.52", "coords": [7.0, 8.0, 9.0]}
            ],
            "edges": [
                {"source": 0, "target": 1, "weight": 0.95, "type": "covalent"},
                {"source": 0, "target": 2, "weight": 0.82, "type": "spatial"},
                {"source": 1, "target": 2, "weight": 0.90, "type": "covalent"}
            ]
        }
    }


@pytest.fixture
def temp_dir():
    """Temporary directory for file operations."""
    dir_path = tempfile.mkdtemp(prefix="protos_test_")
    yield dir_path
    # Clean up
    shutil.rmtree(dir_path, ignore_errors=True)


# Format handler tests
class TestCSVHandler:
    """Tests for the CSVHandler."""

    def test_read_write(self, temp_dir):
        """Test reading and writing a CSV file."""
        # Create test data
        df = pd.DataFrame({
            "A": [1, 2, 3],
            "B": ["x", "y", "z"]
        })
        
        # Create handler
        handler = CSVHandler()
        
        # Write file
        file_path = os.path.join(temp_dir, "test.csv")
        handler.write(file_path, df, index=False)
        
        # Read file
        df_read = handler.read(file_path)
        print(df_read)
        
        # Check results
        pd.testing.assert_frame_equal(df_read, df)
    
    def test_validation(self):
        """Test data validation."""
        handler = CSVHandler()
        
        # Valid data should not raise errors
        handler._validate(pd.DataFrame({"A": [1, 2, 3]}))
        
        # Invalid data should raise errors
        with pytest.raises(ValueError):
            handler._validate("not a dataframe")
        with pytest.raises(ValueError):
            handler._validate([1, 2, 3])


class TestPickleHandler:
    """Tests for the PickleHandler."""

    def test_read_write_primitives(self, temp_dir):
        """Test reading and writing primitive types."""
        handler = PickleHandler()
        
        data_types = [
            123,
            "string",
            [1, 2, 3],
            {"a": 1, "b": 2},
            True,
            None
        ]
        
        for i, data in enumerate(data_types):
            file_path = os.path.join(temp_dir, f"test_{i}.pkl")
            handler.write(file_path, data)
            
            # Read and verify
            read_data = handler.read(file_path)
            assert read_data == data
    
    def test_read_write_complex(self, temp_dir, sample_embeddings):
        """Test reading and writing complex types like numpy arrays."""
        handler = PickleHandler()
        file_path = os.path.join(temp_dir, "embeddings.pkl")
        
        # Write
        handler.write(file_path, sample_embeddings)
        
        # Read
        read_data = handler.read(file_path)
        
        # Verify
        assert list(read_data.keys()) == list(sample_embeddings.keys())
        for key in sample_embeddings:
            np.testing.assert_array_equal(read_data[key], sample_embeddings[key])


class TestJSONHandler:
    """Tests for the JSONHandler."""

    def test_read_write(self, temp_dir):
        """Test reading and writing JSON files."""
        handler = JSONHandler()
        
        # Test data
        data = {
            "string": "value",
            "number": 123,
            "list": [1, 2, 3],
            "nested": {
                "a": True,
                "b": False
            }
        }
        
        # Write file
        file_path = os.path.join(temp_dir, "test.json")
        handler.write(file_path, data)
        
        # Read file
        read_data = handler.read(file_path)
        
        # Verify
        assert read_data == data
    
    def test_validation(self):
        """Test data validation."""
        handler = JSONHandler()
        
        # Valid data should not raise errors
        handler._validate({"a": 1, "b": 2})
        
        # Test with non-serializable data
        non_serializable = {"key": object()}
        with pytest.raises(ValueError):
            handler._validate(non_serializable)


class TestFASTAHandler:
    """Tests for the FASTAHandler."""

    def test_read_write(self, temp_dir, sample_sequences):
        """Test reading and writing FASTA files."""
        handler = FASTAHandler()
        
        # Write file
        file_path = os.path.join(temp_dir, "test.fasta")
        handler.write(file_path, sample_sequences)
        
        # Read file
        read_data = handler.read(file_path)
        
        # Verify
        assert read_data == sample_sequences
    
    def test_read_malformed(self, temp_dir):
        """Test reading a malformed FASTA file (without > prefix)."""
        handler = FASTAHandler()
        
        # Create malformed FASTA
        file_path = os.path.join(temp_dir, "malformed.fasta")
        with open(file_path, 'w') as f:
            f.write("ACGTACGTACGT\nACGTACGT")
        
        # Read file - should handle as a single unnamed sequence
        read_data = handler.read(file_path)
        
        # Verify
        assert "unnamed_sequence" in read_data
        # The implementation removes all whitespace, so the expected result is all characters together
        assert read_data["unnamed_sequence"] == "ACGTACGTACGTACGTACGT"
    
    def test_validation(self):
        """Test data validation."""
        handler = FASTAHandler()
        
        # Valid data should not raise errors
        handler._validate({"seq1": "ACGT", "seq2": "TGCA"})
        
        # Test with invalid data types
        with pytest.raises(ValueError):
            handler._validate("not a dict")
        
        with pytest.raises(ValueError):
            handler._validate({123: "ACGT"})
        
        with pytest.raises(ValueError):
            handler._validate({"seq1": 123})


class TestGRNTableHandler:
    """Tests for the GRNTableHandler."""

    def test_read_write(self, temp_dir, sample_grn_table):
        """Test reading and writing GRN tables."""
        handler = GRNTableHandler()
        
        # Write file
        file_path = os.path.join(temp_dir, "test_grn.csv")
        handler.write(file_path, sample_grn_table)
        
        # Debug: print the CSV content
        with open(file_path, 'r') as f:
            csv_content = f.read()
            print(f"\nCSV content:\n{csv_content}")
            
        # Read file
        read_data = handler.read(file_path)
        
        print("\nOriginal GRN table shape:", sample_grn_table.shape)
        print("Read GRN table shape:", read_data.shape)
        print("Original columns:", sample_grn_table.columns.tolist())
        print("Read columns:", read_data.columns.tolist())
        
        # Verify data matches
        pd.testing.assert_frame_equal(read_data, sample_grn_table)
    
    def test_validation(self, sample_grn_table):
        """Test data validation."""
        handler = GRNTableHandler()
        
        # Valid data should not raise errors
        handler._validate(sample_grn_table)
        
        # Test without protein_id
        invalid_df = pd.DataFrame({"col1": [1, 2, 3]})
        with pytest.raises(ValueError):
            handler._validate(invalid_df)
        
        # Test with invalid GRN format (needs amino acid + position)
        invalid_grn = sample_grn_table.copy()
        invalid_grn.at["opsin1", "3.50"] = "125"  # Missing amino acid
        with pytest.raises(ValueError):
            handler._validate(invalid_grn, skip_cell_validation=False)


class TestPropertyTableHandler:
    """Tests for the PropertyTableHandler."""

    def test_read_write(self, temp_dir, sample_property_table):
        """Test reading and writing property tables with metadata."""
        handler = PropertyTableHandler()
        
        # Write file
        file_path = os.path.join(temp_dir, "test_properties.csv")
        handler.write(file_path, sample_property_table, index=False)
        
        # Read file
        read_data = handler.read(file_path)
        
        # Check DataFrame content
        pd.testing.assert_frame_equal(read_data, sample_property_table.reset_index(drop=True))
        
        # Check if metadata was saved and loaded
        assert hasattr(read_data, "metadata")
        assert read_data.metadata == sample_property_table.metadata
        
        # Check if metadata file exists
        metadata_path = os.path.join(temp_dir, "test_properties_metadata.pkl")
        assert os.path.exists(metadata_path)
    
    def test_validation(self, sample_property_table):
        """Test data validation."""
        handler = PropertyTableHandler()
        
        # Valid data should not raise errors
        handler._validate(sample_property_table)
        
        # Test without protein_id
        invalid_df = pd.DataFrame({"col1": [1, 2, 3]})
        with pytest.raises(ValueError):
            handler._validate(invalid_df)


class TestEmbeddingHandler:
    """Tests for the EmbeddingHandler."""

    def test_read_write(self, temp_dir, sample_embeddings):
        """Test reading and writing embedding pickle files."""
        handler = EmbeddingHandler()
        
        # Write file
        file_path = os.path.join(temp_dir, "test_embeddings.pkl")
        handler.write(file_path, sample_embeddings)
        
        # Read file
        read_data = handler.read(file_path)
        
        # Verify
        assert list(read_data.keys()) == list(sample_embeddings.keys())
        for key in sample_embeddings:
            np.testing.assert_array_equal(read_data[key], sample_embeddings[key])
    
    def test_validation(self, sample_embeddings):
        """Test data validation."""
        handler = EmbeddingHandler()
        
        # Valid data should not raise errors
        handler._validate(sample_embeddings)
        
        # Test with non-dict data
        with pytest.raises(ValueError):
            handler._validate("not a dict")
        
        # Test with non-numpy array values
        invalid_emb = {"seq1": "not an array"}
        with pytest.raises(ValueError):
            handler._validate(invalid_emb)
        
        # Test with inconsistent dimensions
        invalid_dims = {
            "seq1": np.random.rand(10, 5),
            "seq2": np.random.rand(10, 10)  # Different dimension
        }
        with pytest.raises(ValueError):
            handler._validate(invalid_dims)


class TestGraphHandler:
    """Tests for the GraphHandler."""

    def test_read_write(self, temp_dir, sample_graph):
        """Test reading and writing graph JSON files."""
        handler = GraphHandler()
        
        # Write file
        file_path = os.path.join(temp_dir, "test_graph.json")
        handler.write(file_path, sample_graph)
        
        # Read file
        read_data = handler.read(file_path)
        
        # Verify
        assert read_data == sample_graph
    
    def test_validation(self, sample_graph):
        """Test data validation."""
        handler = GraphHandler()
        
        # Valid data should not raise errors
        handler._validate(sample_graph)
        
        # Test with non-dict data
        with pytest.raises(ValueError):
            handler._validate("not a dict")
        
        # Test with missing nodes or edges
        invalid_graph = {"protein_id": {"nodes": []}}  # Missing edges
        with pytest.raises(ValueError):
            handler._validate(invalid_graph)
        
        # Test with invalid nodes
        invalid_nodes = {
            "protein_id": {
                "nodes": [{"missing_id": "value"}],  # Missing id field
                "edges": [{"source": 0, "target": 1}]
            }
        }
        with pytest.raises(ValueError):
            handler._validate(invalid_nodes)
        
        # Test with invalid edges
        invalid_edges = {
            "protein_id": {
                "nodes": [{"id": 0}],
                "edges": [{"source": 0}]  # Missing target field
            }
        }
        with pytest.raises(ValueError):
            handler._validate(invalid_edges)


# Format registry tests
class TestFormatRegistry:
    """Tests for the FormatRegistry."""

    def test_get_handler(self):
        """Test getting handlers by format type."""
        registry = FormatRegistry()
        
        # Test all supported formats
        assert isinstance(registry.get_handler("csv"), CSVHandler)
        assert isinstance(registry.get_handler("pkl"), PickleHandler)
        assert isinstance(registry.get_handler("pickle"), PickleHandler)
        assert isinstance(registry.get_handler("json"), JSONHandler)
        assert isinstance(registry.get_handler("fasta"), FASTAHandler)
        assert isinstance(registry.get_handler("fas"), FASTAHandler)
        assert isinstance(registry.get_handler("grn"), GRNTableHandler)
        assert isinstance(registry.get_handler("emb"), EmbeddingHandler)
        assert isinstance(registry.get_handler("graph"), GraphHandler)
        
        # Test unsupported format
        with pytest.raises(ValueError):
            registry.get_handler("unsupported")
    
    def test_register_handler(self):
        """Test registering a custom handler."""
        registry = FormatRegistry()
        
        # Create a custom handler
        class CustomHandler(FormatHandler):
            def _read_impl(self, file_path, **kwargs):
                return "custom read"
            
            def _write_impl(self, file_path, data, **kwargs):
                pass
        
        # Register the handler
        registry.register_handler("custom", CustomHandler())
        
        # Get the handler
        handler = registry.get_handler("custom")
        
        # Verify
        assert isinstance(handler, CustomHandler)
    
    def test_infer_format(self, temp_dir):
        """Test inferring format from file path."""
        registry = FormatRegistry()
        
        # Test various extensions
        assert registry.infer_format("test.csv") == "csv"
        assert registry.infer_format("test.pkl") == "pkl"
        assert registry.infer_format("test.pickle") == "pickle"
        assert registry.infer_format("test.json") == "json"
        assert registry.infer_format("test.fasta") == "fasta"
        assert registry.infer_format("test.fas") == "fas"
        
        # Test GRN-specific inference
        assert registry.infer_format("grn_table.csv") == "grn"
        
        # Test unsupported extension
        with pytest.raises(ValueError):
            registry.infer_format("test.unsupported")
    
    def test_read_write(self, temp_dir, sample_sequences):
        """Test high-level read/write methods."""
        registry = FormatRegistry()
        
        # Write a file
        file_path = os.path.join(temp_dir, "test.fasta")
        registry.write(file_path, sample_sequences)
        
        # Read it back
        read_data = registry.read(file_path)
        
        # Verify
        assert read_data == sample_sequences


# Temporary file management tests
class TestTempFileManager:
    """Tests for the TempFileManager."""

    def test_temp_file_context(self):
        """Test creating a temporary file with context manager."""
        manager = TempFileManager()
        
        # Use context manager
        with manager.temp_file(suffix=".txt") as temp_path:
            # Check that file exists
            assert os.path.exists(temp_path)
            
            # Write to file
            with open(temp_path, 'w') as f:
                f.write("test content")
                
            # Read from file
            with open(temp_path, 'r') as f:
                content = f.read()
                assert content == "test content"
        
        # After context, file should be deleted
        assert not os.path.exists(temp_path)
    
    def test_temp_file_with_content(self):
        """Test creating a temporary file with initial content."""
        manager = TempFileManager()
        
        # Create file with content
        with manager.temp_file(suffix=".txt", content="initial content") as temp_path:
            # Check content
            with open(temp_path, 'r') as f:
                content = f.read()
                assert content == "initial content"
    
    def test_temp_dir_context(self):
        """Test creating a temporary directory with context manager."""
        manager = TempFileManager()
        
        # Use context manager
        with manager.temp_dir() as temp_dir:
            # Check that directory exists
            assert os.path.exists(temp_dir)
            assert os.path.isdir(temp_dir)
            
            # Create a file in the directory
            test_file = os.path.join(temp_dir, "test.txt")
            with open(test_file, 'w') as f:
                f.write("test content")
            
            # Check that file exists
            assert os.path.exists(test_file)
        
        # After context, directory should be deleted
        assert not os.path.exists(temp_dir)
    
    def test_explicit_temp_file(self):
        """Test creating a temporary file explicitly (without context)."""
        manager = TempFileManager()
        
        # Create file
        temp_path = manager.create_temp_file(suffix=".txt")
        
        try:
            # Check that file exists
            assert os.path.exists(temp_path)
            
            # Write to file
            with open(temp_path, 'w') as f:
                f.write("test content")
            
            # Check content
            with open(temp_path, 'r') as f:
                content = f.read()
                assert content == "test content"
        
        finally:
            # Clean up
            manager.cleanup(temp_path)
            assert not os.path.exists(temp_path)
    
    def test_explicit_temp_dir(self):
        """Test creating a temporary directory explicitly (without context)."""
        manager = TempFileManager()
        
        # Create directory
        temp_dir = manager.create_temp_dir()
        
        try:
            # Check that directory exists
            assert os.path.exists(temp_dir)
            assert os.path.isdir(temp_dir)
            
            # Create a file in the directory
            test_file = os.path.join(temp_dir, "test.txt")
            with open(test_file, 'w') as f:
                f.write("test content")
            
            # Check that file exists
            assert os.path.exists(test_file)
        
        finally:
            # Clean up
            manager.cleanup(temp_dir)
            assert not os.path.exists(temp_dir)
    
    def test_cleanup_all(self):
        """Test cleaning up all temporary files and directories."""
        manager = TempFileManager()
        
        # Create multiple files and directories
        temp_files = [manager.create_temp_file() for _ in range(3)]
        temp_dirs = [manager.create_temp_dir() for _ in range(2)]
        
        # Check tracking lists
        assert set(manager.get_temp_files()) == set(temp_files)
        assert set(manager.get_temp_dirs()) == set(temp_dirs)
        
        # Clean up all
        manager.cleanup()
        
        # Check that all are deleted
        for path in temp_files + temp_dirs:
            assert not os.path.exists(path)
        
        # Check tracking lists are empty
        assert manager.get_temp_files() == []
        assert manager.get_temp_dirs() == []


# External tool helper tests
class TestExternalToolHelper:
    """Tests for the ExternalToolHelper."""

    def test_tool_input(self, sample_sequences):
        """Test preparing input for an external tool."""
        helper = ExternalToolHelper("test_tool")
        
        # Use tool_input context manager
        with helper.tool_input(sample_sequences, format_type="fasta") as input_file:
            # Check that file exists and has the right extension
            assert os.path.exists(input_file)
            assert input_file.endswith(".fa")
            
            # Check file content
            read_data = format_registry.read(input_file, format_type="fasta")
            assert read_data == sample_sequences
    
    def test_tool_output(self):
        """Test preparing for output from an external tool."""
        helper = ExternalToolHelper("test_tool")
        
        # Use tool_output context manager
        with helper.tool_output(suffix=".json", expected_format="json") as (output_path, read_func):
            # Check that file exists
            assert os.path.exists(output_path)
            assert output_path.endswith(".json")
            
            # Write some data to the file
            data = {"result": "success"}
            with open(output_path, 'w') as f:
                json.dump(data, f)
            
            # Use the read function
            result = read_func()
            assert result == data
    
    def test_tool_workspace(self):
        """Test creating a workspace for an external tool."""
        helper = ExternalToolHelper("test_tool")
        
        # Use tool_workspace context manager
        with helper.tool_workspace() as workspace:
            # Check that directory exists
            assert os.path.exists(workspace)
            assert os.path.isdir(workspace)
            assert "_test_tool" in workspace
            
            # Create some files in the workspace
            file1 = os.path.join(workspace, "file1.txt")
            file2 = os.path.join(workspace, "file2.txt")
            
            with open(file1, 'w') as f:
                f.write("content1")
            
            with open(file2, 'w') as f:
                f.write("content2")
            
            # Check that files exist
            assert os.path.exists(file1)
            assert os.path.exists(file2)
        
        # After context, workspace should be deleted
        assert not os.path.exists(workspace)
    
    def test_create_temp_copy(self, temp_dir):
        """Test creating a temporary copy of a file."""
        helper = ExternalToolHelper("test_tool")
        
        # Create a source file
        source_path = os.path.join(temp_dir, "source.txt")
        with open(source_path, 'w') as f:
            f.write("source content")
        
        # Create a temporary copy
        temp_copy = helper.create_temp_copy(source_path)
        
        try:
            # Check that copy exists
            assert os.path.exists(temp_copy)
            
            # Check content
            with open(temp_copy, 'r') as f:
                content = f.read()
                assert content == "source content"
                
            # Modify copy (shouldn't affect source)
            with open(temp_copy, 'w') as f:
                f.write("modified content")
            
            # Check source is unchanged
            with open(source_path, 'r') as f:
                content = f.read()
                assert content == "source content"
        
        finally:
            # Clean up
            from protos.io.formats import temp_manager
            temp_manager.cleanup(temp_copy)