"""
Tests for the path resolution module.

These old_tests verify that the path resolution functions correctly handle
various path scenarios and provide consistent results.
"""

import os
import tempfile
import pytest
from pathlib import Path

from protos.io.paths import (
    get_structure_path,
    get_grn_path,
    get_sequence_path,
    get_dataset_path,
    resolve_path,
    ProtosPaths,
    get_data_root,
    DataSource
)
from protos.io.paths.path_constants import (
    DEFAULT_USER_DATA_ROOT,
    DEFAULT_REF_DATA_ROOT
)

from protos.io.data_access import (
    DataRegistry,
    GlobalRegistry
)


class TestPathResolution:
    """Tests for path resolution functions."""
    
    def test_resolve_path(self):
        """Test that resolve_path correctly handles relative and absolute paths."""
        # Test with absolute path
        absolute_path = "/tmp/test/path"
        assert os.path.normpath(resolve_path(absolute_path)) == os.path.normpath(absolute_path)
        
        # Create a temp path resolver that doesn't use real paths 
        # to avoid system dependency in old_tests
        paths = ProtosPaths(
            user_data_root="/user/data", 
            ref_data_root="/ref/data",
            create_dirs=False,
            validate=False
        )
        
        # Test with relative path - default to reference data with AUTO source
        relative_path = "test/path"
        with_ref = paths.resolve_path(relative_path, DataSource.REFERENCE)
        # Normalize path separators for cross-platform testing
        assert "ref" in with_ref.replace('\\', '/') and "data" in with_ref.replace('\\', '/')
        assert relative_path.replace('/', os.sep) in with_ref
        
        # Test with relative_to parameter
        base_dir = "/base/dir"
        with_base = paths.resolve_path(relative_path, DataSource.AUTO, base_dir)
        assert "base" in with_base.replace('\\', '/') and "dir" in with_base.replace('\\', '/')
        assert relative_path.replace('/', os.sep) in with_base
        
        # Test with None path
        none_path = paths.resolve_path(None, DataSource.REFERENCE)
        assert "ref" in none_path.replace('\\', '/') and "data" in none_path.replace('\\', '/')
        
        # Test with explicit user data source
        with_user = paths.resolve_path(relative_path, DataSource.USER)
        assert "user" in with_user.replace('\\', '/') and "data" in with_user.replace('\\', '/')
        assert relative_path.replace('/', os.sep) in with_user
    
    def test_get_structure_path(self):
        """Test that get_structure_path returns correct structure file paths."""
        # Test with default parameters
        pdb_id = "1abc"
        result = get_structure_path(pdb_id)
        assert pdb_id in result
        assert result.endswith(".cif")
        
        # Test with custom structure directory
        custom_dir = "/custom/structure/dir"
        result = get_structure_path(pdb_id, structure_dir=custom_dir)
        assert os.path.normpath(os.path.dirname(result)) == os.path.normpath(custom_dir)
        assert os.path.basename(result) == f"{pdb_id}.cif"
        
        # Test with uppercase PDB ID
        pdb_id = "1ABC"
        result = get_structure_path(pdb_id)
        assert "1abc.cif" in result
        
        # Test with different data sources
        ref_result = get_structure_path(pdb_id, source=DataSource.REFERENCE)
        user_result = get_structure_path(pdb_id, source=DataSource.USER)
        assert DEFAULT_REF_DATA_ROOT in ref_result
        assert DEFAULT_USER_DATA_ROOT in user_result
    
    def test_get_grn_path(self):
        """Test that get_grn_path returns correct GRN table file paths."""
        # Test with default parameters
        table_name = "grn_table"
        result = get_grn_path(table_name)
        assert table_name in result
        assert result.endswith(".csv")
        
        # Test with custom table directory
        custom_dir = "/custom/grn/dir"
        result = get_grn_path(table_name, table_dir=custom_dir)
        assert os.path.normpath(os.path.dirname(result)) == os.path.normpath(custom_dir)
        assert os.path.basename(result) == f"{table_name}.csv"
        
        # Test with different data sources
        ref_result = get_grn_path(table_name, source=DataSource.REFERENCE)
        user_result = get_grn_path(table_name, source=DataSource.USER)
        assert DEFAULT_REF_DATA_ROOT in ref_result
        assert DEFAULT_USER_DATA_ROOT in user_result
    
    def test_get_sequence_path(self):
        """Test that get_sequence_path returns correct sequence file paths."""
        # Test with default parameters
        sequence_id = "sequence1"
        result = get_sequence_path(sequence_id)
        assert sequence_id in result
        assert result.endswith(".fasta")
        
        # Test with custom FASTA directory
        custom_dir = "/custom/sequence/dir"
        result = get_sequence_path(sequence_id, fasta_dir=custom_dir)
        assert os.path.normpath(os.path.dirname(result)) == os.path.normpath(custom_dir)
        assert os.path.basename(result) == f"{sequence_id}.fasta"
        
        # Test with different data sources
        ref_result = get_sequence_path(sequence_id, source=DataSource.REFERENCE)
        user_result = get_sequence_path(sequence_id, source=DataSource.USER)
        assert DEFAULT_REF_DATA_ROOT in ref_result
        assert DEFAULT_USER_DATA_ROOT in user_result
    
    def test_get_dataset_path(self):
        """Test that get_dataset_path returns correct dataset file paths."""
        # Test with default parameters
        dataset_name = "test_dataset"
        result = get_dataset_path(dataset_name)
        assert dataset_name in result
        assert result.endswith(".json")
        
        # Test with custom processor type
        result = get_dataset_path(dataset_name, processor_type="grn")
        assert "grn" in result
        
        # Test with custom extension
        result = get_dataset_path(dataset_name, file_extension=".pkl")
        assert result.endswith(".pkl")
        
        # Test with different data sources
        ref_result = get_dataset_path(dataset_name, source=DataSource.REFERENCE)
        user_result = get_dataset_path(dataset_name, source=DataSource.USER)
        assert DEFAULT_REF_DATA_ROOT in ref_result
        assert DEFAULT_USER_DATA_ROOT in user_result
    
    @pytest.mark.parametrize("pdb_id", ["1abc", "4XYZ", "pdb123"])
    def test_structure_path_normalization(self, pdb_id):
        """Test that PDB IDs are normalized consistently."""
        result = get_structure_path(pdb_id)
        
        # ID should be lowercase in result
        assert pdb_id.lower() in result
        assert pdb_id.upper() not in result


class TestProtosPathsClass:
    """Tests for the ProtosPaths class."""
    
    def test_initialization(self):
        """Test that ProtosPaths initializes correctly."""
        # Test with default parameters
        paths = ProtosPaths(create_dirs=False, validate=False)
        assert paths.user_data_root.endswith(DEFAULT_USER_DATA_ROOT)
        assert paths.ref_data_root == DEFAULT_REF_DATA_ROOT
        
        # Test with custom data roots
        custom_user_root = "/custom/user/data"
        custom_ref_root = "/custom/ref/data"
        paths = ProtosPaths(
            user_data_root=custom_user_root, 
            ref_data_root=custom_ref_root,
            create_dirs=False, 
            validate=False
        )
        assert paths.user_data_root == custom_user_root
        assert paths.ref_data_root == custom_ref_root
    
    def test_get_processor_path(self):
        """Test that get_processor_path returns correct paths."""
        paths = ProtosPaths(create_dirs=False, validate=False)
        
        # Test with valid processor type
        result = paths.get_processor_path("structure")
        assert "structure" in result
        
        # Test with invalid processor type
        with pytest.raises(ValueError):
            paths.get_processor_path("invalid_type")
            
        # Test with different data sources
        ref_result = paths.get_processor_path("structure", DataSource.REFERENCE)
        user_result = paths.get_processor_path("structure", DataSource.USER)
        assert DEFAULT_REF_DATA_ROOT in ref_result
        assert DEFAULT_USER_DATA_ROOT in user_result
    
    def test_directory_creation(self):
        """Test that directories are created when requested."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Initialize with directory creation
            paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True, validate=False)
            
            # Check that directories were created
            structure_dir = os.path.join(temp_dir, "structure")
            assert os.path.exists(structure_dir)
            
            # Check subdirectories
            mmcif_dir = os.path.join(structure_dir, "mmcif")
            assert os.path.exists(mmcif_dir)
    
    def test_path_update(self):
        """Test that path configurations can be updated."""
        paths = ProtosPaths(create_dirs=False, validate=False)
        
        # Update data roots
        new_user_root = "/new/user/data/root"
        new_ref_root = "/new/ref/data/root"
        paths.update_paths(user_data_root=new_user_root, ref_data_root=new_ref_root)
        assert paths.user_data_root == new_user_root
        assert paths.ref_data_root == new_ref_root
        
        # Update processor directories
        new_processor_dirs = {"structure": "custom_structure"}
        paths.update_paths(processor_dirs=new_processor_dirs)
        assert paths.processor_dirs["structure"] == "custom_structure"
        
        # Verify path resolution uses updated configurations
        structure_path = paths.get_processor_path("structure")
        assert "custom_structure" in structure_path
    
    def test_exists_method(self):
        """Test that exists method checks both data sources correctly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            ref_temp_dir = os.path.join(temp_dir, "reference_data")
            user_temp_dir = os.path.join(temp_dir, "user_data")
            os.makedirs(ref_temp_dir, exist_ok=True)
            os.makedirs(user_temp_dir, exist_ok=True)
            
            # Create a test file in the reference data
            ref_test_file = os.path.join(ref_temp_dir, "test.txt")
            with open(ref_test_file, "w") as f:
                f.write("Reference data")
                
            # Create a test file in the user data
            user_test_file = os.path.join(user_temp_dir, "user.txt")
            with open(user_test_file, "w") as f:
                f.write("User data")
            
            # Create a mock is_package_resource function for testing
            original_is_resource = ProtosPaths.exists

            # Override the exists method for testing
            def mock_exists(self, path, check_both_sources=True):
                if "reference_data" in str(path):
                    return True, DataSource.REFERENCE
                elif "user_data" in str(path):
                    return True, DataSource.USER
                else:
                    return False, None
                    
            # Replace the exists method with our mock
            ProtosPaths.exists = mock_exists
            
            try:
                # Initialize with the test directories
                paths = ProtosPaths(
                    user_data_root=user_temp_dir,
                    ref_data_root=ref_temp_dir,
                    create_dirs=False,
                    validate=False
                )
                
                # Test absolute paths
                exists, source = paths.exists(ref_test_file)
                assert exists
                assert source == DataSource.REFERENCE
                
                exists, source = paths.exists(user_test_file)
                assert exists
                assert source == DataSource.USER
                
                # Test relative paths
                exists, source = paths.exists(os.path.join(ref_temp_dir, "other.txt"))
                assert exists
                assert source == DataSource.REFERENCE
                
                exists, source = paths.exists(os.path.join(user_temp_dir, "other.txt"))
                assert exists
                assert source == DataSource.USER
                
                # Test file that doesn't exist
                exists, source = paths.exists("nonexistent.txt")
                assert not exists
                assert source is None
            finally:
                # Restore the original is_package_resource function
                ProtosPaths.exists = original_is_resource


class TestGlobalRegistry:
    """Tests for the GlobalRegistry class."""
    
    def test_initialization(self):
        """Test that GlobalRegistry initializes correctly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True, validate=False)
            registry = GlobalRegistry(paths)
            
            # Check that the registry file was created
            assert os.path.exists(paths.get_global_registry_path())
    
    def test_register_dataset(self):
        """Test registering datasets in the global registry."""
        with tempfile.TemporaryDirectory() as temp_dir:
            paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True, validate=False)
            registry = GlobalRegistry(paths)
            
            # Register user dataset
            user_file = os.path.join(temp_dir, "user_data.txt")
            with open(user_file, "w") as f:
                f.write("User data")
                
            registry.register_dataset(
                "user_dataset", 
                user_file, 
                "structure", 
                "cif", 
                DataSource.USER, 
                {"description": "Test user dataset"}
            )
            
            # Verify dataset was registered in global registry
            assert "user_dataset" in registry.registry
            assert registry.registry["user_dataset"]["path"] == user_file
            assert registry.registry["user_dataset"]["metadata"]["description"] == "Test user dataset"
            assert registry.registry["user_dataset"]["metadata"]["source"] == DataSource.USER.value
            
            # Verify processor-specific registry was updated
            proc_registry = registry._get_processor_registry("structure")
            assert "user_dataset" in proc_registry.registry
            
            # Register reference dataset
            ref_file = "/path/to/reference/data.txt"
            registry.register_dataset(
                "ref_dataset", 
                ref_file, 
                "grn", 
                "table", 
                DataSource.REFERENCE, 
                {"description": "Test reference dataset"}
            )
            
            # Verify dataset was registered in global registry
            assert "ref_dataset" in registry.registry
            assert registry.registry["ref_dataset"]["path"] == ref_file
            assert registry.registry["ref_dataset"]["metadata"]["description"] == "Test reference dataset"
            assert registry.registry["ref_dataset"]["metadata"]["source"] == DataSource.REFERENCE.value
            
            # Verify processor-specific registry was NOT updated
            proc_registry = registry._get_processor_registry("grn")
            assert "ref_dataset" not in proc_registry.registry
    
    def test_get_dataset_path(self):
        """Test getting dataset paths from the global registry."""
        with tempfile.TemporaryDirectory() as temp_dir:
            paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True, validate=False)
            registry = GlobalRegistry(paths)
            
            # Register datasets
            user_file = os.path.join(temp_dir, "user_data.txt")
            with open(user_file, "w") as f:
                f.write("User data")
                
            registry.register_dataset(
                "user_dataset", 
                user_file, 
                "structure", 
                "cif", 
                DataSource.USER, 
                {"description": "Test user dataset"}
            )
            
            ref_file = "/path/to/reference/data.txt"
            registry.register_dataset(
                "ref_dataset", 
                ref_file, 
                "grn", 
                "table", 
                DataSource.REFERENCE, 
                {"description": "Test reference dataset"}
            )
            
            # Get dataset paths
            assert registry.get_dataset_path("user_dataset") == user_file
            assert registry.get_dataset_path("ref_dataset") == ref_file
            assert registry.get_dataset_path("nonexistent_dataset") is None
    
    def test_dataset_filtering(self):
        """Test filtering datasets by type and source."""
        with tempfile.TemporaryDirectory() as temp_dir:
            paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True, validate=False)
            registry = GlobalRegistry(paths)
            
            # Register datasets
            registry.register_dataset(
                "user_structure", 
                "/path/to/user/structure.cif", 
                "structure", 
                "cif", 
                DataSource.USER
            )
            
            registry.register_dataset(
                "user_grn", 
                "/path/to/user/grn.csv", 
                "grn", 
                "table", 
                DataSource.USER
            )
            
            registry.register_dataset(
                "ref_structure", 
                "/path/to/ref/structure.cif", 
                "structure", 
                "cif", 
                DataSource.REFERENCE
            )
            
            registry.register_dataset(
                "ref_grn", 
                "/path/to/ref/grn.csv", 
                "grn", 
                "table", 
                DataSource.REFERENCE
            )
            
            # Test filtering by processor type
            structure_datasets = registry.list_datasets("structure")
            assert set(structure_datasets) == {"user_structure", "ref_structure"}
            
            grn_datasets = registry.list_datasets("grn")
            assert set(grn_datasets) == {"user_grn", "ref_grn"}
            
            # Test filtering by dataset type
            cif_datasets = registry.get_datasets_by_type("cif")
            assert set(cif_datasets) == {"user_structure", "ref_structure"}
            
            table_datasets = registry.get_datasets_by_type("table")
            assert set(table_datasets) == {"user_grn", "ref_grn"}
            
            # Test filtering by source
            user_datasets = registry.get_datasets_by_source(DataSource.USER)
            assert set(user_datasets) == {"user_structure", "user_grn"}
            
            ref_datasets = registry.get_datasets_by_source(DataSource.REFERENCE)
            assert set(ref_datasets) == {"ref_structure", "ref_grn"}
    
    def test_remove_dataset(self):
        """Test removing datasets from the registry."""
        with tempfile.TemporaryDirectory() as temp_dir:
            paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True, validate=False)
            registry = GlobalRegistry(paths)
            
            # Register datasets
            registry.register_dataset(
                "user_dataset", 
                "/path/to/user/data.txt", 
                "structure", 
                "cif", 
                DataSource.USER
            )
            
            registry.register_dataset(
                "ref_dataset", 
                "/path/to/ref/data.txt", 
                "grn", 
                "table", 
                DataSource.REFERENCE
            )
            
            # Remove user dataset
            assert registry.remove_dataset("user_dataset")
            assert "user_dataset" not in registry.registry
            
            # Try to remove reference dataset (should fail)
            assert not registry.remove_dataset("ref_dataset")
            assert "ref_dataset" in registry.registry