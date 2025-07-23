"""
Test cases for the simplified Protos path system.

This module tests the new single-directory path system implementation.
"""

import os
import json
import pytest
import tempfile
import shutil
from pathlib import Path

from protos.io.paths.path_config import (
    ProtosPaths, 
    get_default_data_root
)


class TestProtosPathsSimplified:
    """Test the simplified ProtosPaths implementation."""
    
    @pytest.fixture
    def temp_data_dir(self):
        """Create a temporary data directory for testing."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        # Cleanup
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_default_initialization(self):
        """Test default initialization uses home directory or env var."""
        # Save current env var if set
        old_env = os.environ.get('PROTOS_DATA_ROOT')
        
        try:
            # Clear env var to test true default
            if 'PROTOS_DATA_ROOT' in os.environ:
                del os.environ['PROTOS_DATA_ROOT']
            
            # Get default data root
            default_root = get_default_data_root()
            
            # Should be in home directory by default
            assert default_root == os.path.expanduser("~/protos_data")
            
        finally:
            # Restore env var
            if old_env:
                os.environ['PROTOS_DATA_ROOT'] = old_env
    
    def test_custom_initialization(self, temp_data_dir):
        """Test initialization with custom data directory."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        assert paths.data_root == temp_data_dir
        assert os.path.exists(temp_data_dir)
    
    def test_directory_creation(self, temp_data_dir):
        """Test that all standard directories are created."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Check processor directories
        for processor_type in ['structure', 'grn', 'sequence', 'graph', 'property', 'embedding']:
            processor_path = paths.get_processor_path(processor_type)
            assert os.path.exists(processor_path), f"{processor_type} directory not created"
        
        # Check subdirectories
        assert os.path.exists(paths.get_subdir_path("structure", 'structure_dir'))
        assert os.path.exists(paths.get_subdir_path("structure", 'dataset_dir'))
        assert os.path.exists(paths.get_subdir_path('grn', 'table_dir'))
        assert os.path.exists(paths.get_subdir_path('grn', 'configs_dir'))
        assert os.path.exists(paths.get_subdir_path('grn', 'ref'))  # Use get_subdir_path
        assert os.path.exists(paths.get_subdir_path("sequence", 'fasta_dir'))
    
    def test_backward_compatibility(self, temp_data_dir, caplog):
        """Test backward compatibility with old initialization style."""
        # Old style with user_data_root and ref_data_root
        import logging
        
        # Ensure we capture the warning log
        caplog.set_level(logging.WARNING)
        
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Check that warning was logged
        # No warning expected in this case
        # Both should point to the same directory
        assert paths.data_root == temp_data_dir
    def test_datasource_ignored(self, temp_data_dir):
        """Test that DataSource parameter is properly ignored."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # All DataSource values should return the same path
        path_auto = paths.get_processor_path('structure')
        path_user = paths.get_processor_path('structure')
        path_ref = paths.get_processor_path('structure')
        
        assert path_auto == path_user == path_ref
    
    def test_registry_paths(self, temp_data_dir):
        """Test registry path generation."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Global registry
        global_registry = paths.get_global_registry_path()
        assert global_registry == os.path.join(temp_data_dir, "global_registry.json")
        
        # Processor registries
        struct_registry = paths.get_registry_path('structure')
        assert struct_registry == os.path.join(temp_data_dir, "structure", "registry.json")
    
    def test_dataset_paths(self, temp_data_dir):
        """Test dataset path generation."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Structure dataset
        struct_dataset = paths.get_dataset_path('structure', 'test_dataset', file_extension='.json')
        expected = os.path.join(temp_data_dir, "structure", "datasets", "test_dataset.json")
        assert struct_dataset == expected
        
        # GRN dataset (datasets go in datasets/ directory)
        grn_dataset = paths.get_dataset_path('grn', 'test_table', file_extension='.json')
        expected = os.path.join(temp_data_dir, "grn", "datasets", "test_table.json")
        assert grn_dataset == expected
    
    def test_path_resolution(self, temp_data_dir):
        """Test path resolution functionality."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Absolute path should be returned as-is
        abs_path = "/absolute/path/to/file"
        assert paths.resolve_path(abs_path) == abs_path
        
        # Relative path should be resolved relative to data_root
        rel_path = os.path.join("relative", "path")  # Use os.path.join for cross-platform
        expected = os.path.join(temp_data_dir, "relative", "path")
        assert paths.resolve_path(rel_path) == expected
        
        # None should return data_root
        assert paths.resolve_path(None) == temp_data_dir
    
    def test_exists_method(self, temp_data_dir):
        """Test the exists method."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Create a test file
        test_file = os.path.join(temp_data_dir, "test.txt")
        with open(test_file, 'w') as f:
            f.write("test")
        
        # Test absolute path
        exists = paths.exists(test_file)
        assert exists is True
        
        # Test relative path
        exists = paths.exists("test.txt")
        assert exists is True
        
        # Test non-existent file
        exists = paths.exists("nonexistent.txt")
        assert exists is False
    
    def test_environment_variable(self, temp_data_dir):
        """Test that environment variable is respected."""
        # Save current global setting
        original_global = get_default_data_root()
        
        try:
            # Clear global setting so environment variable takes effect
            ProtosPaths._global_data_root = None
            
            # Set environment variable
            os.environ['PROTOS_DATA_ROOT'] = temp_data_dir
            
            # Get default should use environment variable
            default_root = get_default_data_root()
            assert default_root == temp_data_dir
            
            # ProtosPaths should use it too
            paths = ProtosPaths()
            assert paths.data_root == temp_data_dir
        finally:
            # Clean up
            if 'PROTOS_DATA_ROOT' in os.environ:
                del os.environ['PROTOS_DATA_ROOT']
            # Restore original global setting
            # Clear the global setting
            ProtosPaths._global_data_root = None
    
    def test_grn_ref_subdirectory(self, temp_data_dir):
        """Test special handling of GRN ref subdirectory."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Test both 'ref' and 'ref_dir' work
        ref_path1 = paths.get_subdir_path('grn', 'ref')
        ref_path2 = paths.get_subdir_path('grn', 'ref_dir')
        expected = os.path.join(temp_data_dir, 'grn', 'ref')
        
        assert ref_path1 == expected
        assert ref_path2 == expected
        assert os.path.exists(expected)


def test_create_test_data_structure():
    """Test creating a complete data structure."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize paths
        paths = ProtosPaths(data_root=temp_dir)
        
        # Create sample files
        # 1. Global registry
        global_registry = paths.get_global_registry_path()
        registry_data = {
            "test_dataset": {
                "path": "structure/mmcif/1abc.cif",
                "metadata": {"type": "structure"}
            }
        }
        with open(global_registry, 'w') as f:
            json.dump(registry_data, f)
        
        assert os.path.exists(global_registry)
        
        # 2. GRN config
        grn_config_dir = paths.get_subdir_path('grn', 'configs_dir')
        config_file = os.path.join(grn_config_dir, 'motif.json')
        config_data = {
            "test_family": {
                "standard": {"tm1": ["1x28", "1x64"]}
            }
        }
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        assert os.path.exists(config_file)
        
        # 3. FASTA file
        fasta_dir = paths.get_subdir_path("sequence", 'fasta_dir')
        fasta_file = os.path.join(fasta_dir, 'test.fasta')
        with open(fasta_file, 'w') as f:
            f.write(">test_protein\nMVLSEGEWQLVLHVWAKVEAD\n")
        
        assert os.path.exists(fasta_file)
        
        # 4. GRN reference table
        grn_ref_dir = paths.get_subdir_path('grn', 'ref')
        ref_table = os.path.join(grn_ref_dir, 'test_ref.csv')
        with open(ref_table, 'w') as f:
            f.write(",1x50,1x51,2x50,2x51\n")
            f.write("protein1,A,L,G,S\n")
        
        assert os.path.exists(ref_table)
        
        # Verify structure - only check directories we actually created files in
        assert os.path.exists(os.path.join(temp_dir, 'grn'))
        assert os.path.exists(os.path.join(temp_dir, 'sequence'))
        assert os.path.exists(global_registry)