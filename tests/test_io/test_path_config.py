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
    get_default_data_root,
    DataSource
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
        """Test default initialization uses home directory."""
        # Get default data root
        default_root = get_default_data_root()
        
        # Should be in home directory by default
        assert default_root == os.path.expanduser("~/protos_data")
    
    def test_custom_initialization(self, temp_data_dir):
        """Test initialization with custom data directory."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        assert paths.data_root == temp_data_dir
        assert paths.user_data_root == temp_data_dir  # Backward compatibility
        assert paths.ref_data_root == temp_data_dir   # Backward compatibility
        assert os.path.exists(temp_data_dir)
    
    def test_directory_creation(self, temp_data_dir):
        """Test that all standard directories are created."""
        paths = ProtosPaths(data_root=temp_data_dir, create_dirs=True)
        
        # Check processor directories
        for processor_type in ['structure', 'grn', 'sequence', 'graph', 'property', 'embedding']:
            processor_path = paths.get_processor_path(processor_type)
            assert os.path.exists(processor_path), f"{processor_type} directory not created"
        
        # Check subdirectories
        assert os.path.exists(paths.get_structure_subdir_path('structure_dir'))
        assert os.path.exists(paths.get_structure_subdir_path('dataset_dir'))
        assert os.path.exists(paths.get_grn_subdir_path('table_dir'))
        assert os.path.exists(paths.get_grn_subdir_path('configs_dir'))
        assert os.path.exists(os.path.join(paths.get_processor_path('grn'), 'ref'))
        assert os.path.exists(paths.get_sequence_subdir_path('fasta_dir'))
    
    def test_backward_compatibility(self, temp_data_dir, caplog):
        """Test backward compatibility with old initialization style."""
        # Old style with user_data_root and ref_data_root
        import logging
        
        # Ensure we capture the warning log
        caplog.set_level(logging.WARNING)
        
        paths = ProtosPaths(
            user_data_root=temp_data_dir,
            ref_data_root="/some/other/path"  # Should be ignored
        )
        
        # Check that warning was logged
        assert len(caplog.records) > 0
        assert any("deprecated" in record.message.lower() for record in caplog.records)
        
        # Both should point to the same directory
        assert paths.data_root == temp_data_dir
        assert paths.user_data_root == temp_data_dir
        assert paths.ref_data_root == temp_data_dir
    
    def test_datasource_ignored(self, temp_data_dir):
        """Test that DataSource parameter is properly ignored."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # All DataSource values should return the same path
        path_auto = paths.get_processor_path('structure', DataSource.AUTO)
        path_user = paths.get_processor_path('structure', DataSource.USER)
        path_ref = paths.get_processor_path('structure', DataSource.REFERENCE)
        
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
        expected = os.path.join(temp_data_dir, "structure", "structure_dataset", "test_dataset.json")
        assert struct_dataset == expected
        
        # GRN dataset
        grn_dataset = paths.get_dataset_path('grn', 'test_table', file_extension='.csv')
        expected = os.path.join(temp_data_dir, "grn", "tables", "test_table.csv")
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
        exists, source = paths.exists(test_file)
        assert exists is True
        assert source == DataSource.USER
        
        # Test relative path
        exists, source = paths.exists("test.txt")
        assert exists is True
        assert source == DataSource.USER
        
        # Test non-existent file
        exists, source = paths.exists("nonexistent.txt")
        assert exists is False
        assert source is None
    
    def test_environment_variable(self, temp_data_dir):
        """Test that environment variable is respected."""
        # Save current global setting
        original_global = ProtosPaths.get_global_data_root()
        
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
            if original_global:
                ProtosPaths.set_data_root(original_global)
    
    def test_grn_ref_subdirectory(self, temp_data_dir):
        """Test special handling of GRN ref subdirectory."""
        paths = ProtosPaths(data_root=temp_data_dir)
        
        # Test both 'ref' and 'ref_dir' work
        ref_path1 = paths.get_grn_subdir_path('ref')
        ref_path2 = paths.get_grn_subdir_path('ref_dir')
        expected = os.path.join(temp_data_dir, 'grn', 'ref')
        
        assert ref_path1 == expected
        assert ref_path2 == expected
        assert os.path.exists(expected)


def test_create_test_data_structure():
    """Test creating a complete data structure."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize paths
        paths = ProtosPaths(data_root=temp_dir, create_dirs=True)
        
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
        grn_config_dir = paths.get_grn_subdir_path('configs_dir')
        config_file = os.path.join(grn_config_dir, 'config.json')
        config_data = {
            "test_family": {
                "standard": {"tm1": ["1x28", "1x64"]}
            }
        }
        with open(config_file, 'w') as f:
            json.dump(config_data, f)
        
        assert os.path.exists(config_file)
        
        # 3. FASTA file
        fasta_dir = paths.get_sequence_subdir_path('fasta_dir')
        fasta_file = os.path.join(fasta_dir, 'test.fasta')
        with open(fasta_file, 'w') as f:
            f.write(">test_protein\nMVLSEGEWQLVLHVWAKVEAD\n")
        
        assert os.path.exists(fasta_file)
        
        # 4. GRN reference table
        grn_ref_dir = os.path.join(paths.get_processor_path('grn'), 'ref')
        ref_table = os.path.join(grn_ref_dir, 'test_ref.csv')
        with open(ref_table, 'w') as f:
            f.write(",1x50,1x51,2x50,2x51\n")
            f.write("protein1,A,L,G,S\n")
        
        assert os.path.exists(ref_table)
        
        # Verify structure
        assert os.path.exists(os.path.join(temp_dir, 'structure'))
        assert os.path.exists(os.path.join(temp_dir, 'grn'))
        assert os.path.exists(os.path.join(temp_dir, 'sequence'))
        assert os.path.exists(os.path.join(temp_dir, 'global_registry.json'))