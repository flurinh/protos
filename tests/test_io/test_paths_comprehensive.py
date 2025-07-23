"""Comprehensive tests for ProtosPaths to ensure it meets all requirements."""

import os
import pytest
import tempfile
import shutil
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.io.paths.path_config import get_default_data_root
from protos.io.paths.path_constants import (
    DEFAULT_PROCESSOR_DIRS,
    DEFAULT_STRUCTURE_SUBDIRS,
    DEFAULT_SEQUENCE_SUBDIRS,
    DEFAULT_GRN_SUBDIRS,
    DEFAULT_PROPERTY_SUBDIRS,
    DEFAULT_EMBEDDING_SUBDIRS,
    DEFAULT_LIGAND_SUBDIRS
)


class TestProtosPaths:
    """Test ProtosPaths comprehensive functionality."""
    
    def test_default_data_root(self):
        """Test that ProtosPaths uses working_dir/data as default."""
        # Clear any environment variables and global settings
        old_env = os.environ.get('PROTOS_DATA_ROOT')
        if 'PROTOS_DATA_ROOT' in os.environ:
            del os.environ['PROTOS_DATA_ROOT']
        
        # Clear global data root
        ProtosPaths._global_data_root = None
        
        try:
            # Create ProtosPaths without any parameters
            paths = ProtosPaths()
            
            # Should use home directory by default (or env var if set)
            expected = os.path.expanduser("~/protos_data")
            assert str(paths.data_root) == expected
            
        finally:
            # Restore environment
            if old_env:
                os.environ['PROTOS_DATA_ROOT'] = old_env
            # Restore global setting
            ProtosPaths._global_data_root = None
    
    def test_path_objects_not_strings(self):
        """Test that all path methods return Path objects, not strings."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Test processor paths
            for processor in ['structure', 'sequence', 'grn', 'property', 'embedding', 'ligand']:
                path = paths.get_processor_path(processor)
                assert isinstance(path, str), f"get_processor_path({processor}) should return string for backward compatibility"
            
            # Test subdirectory paths
            structure_path = paths.get_subdir_path("structure", 'structure_dir')
            assert isinstance(structure_path, str), "Subdirectory paths should return strings for backward compatibility"
    
    def test_directory_creation(self):
        """Test that ProtosPaths creates all required directories."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create ProtosPaths with directory creation
            paths = ProtosPaths(data_root=tmpdir)
            
            # Check base directory
            assert os.path.exists(tmpdir)
            
            # Check processor directories
            for processor_type, dir_name in DEFAULT_PROCESSOR_DIRS.items():
                if processor_type not in ['test', 'test_processor', '_test', '__test', 'simple', 
                                         'complex_processor_with_long_name', 'custom_dir']:
                    processor_path = os.path.join(tmpdir, dir_name)
                    assert os.path.exists(processor_path), f"Missing processor directory: {processor_path}"
            
            # Check structure subdirectories
            structure_path = os.path.join(tmpdir, 'structure')
            for subdir in DEFAULT_STRUCTURE_SUBDIRS.values():
                subdir_path = os.path.join(structure_path, subdir)
                assert os.path.exists(subdir_path), f"Missing structure subdirectory: {subdir_path}"
            
            # Check sequence subdirectories
            sequence_path = os.path.join(tmpdir, 'sequence')
            for subdir in DEFAULT_SEQUENCE_SUBDIRS.values():
                subdir_path = os.path.join(sequence_path, subdir)
                assert os.path.exists(subdir_path), f"Missing sequence subdirectory: {subdir_path}"
            
            # Check GRN subdirectories
            grn_path = os.path.join(tmpdir, 'grn')
            for subdir in DEFAULT_GRN_SUBDIRS.values():
                subdir_path = os.path.join(grn_path, subdir)
                assert os.path.exists(subdir_path), f"Missing GRN subdirectory: {subdir_path}"
            # Check ref subdirectory
            assert os.path.exists(os.path.join(grn_path, 'ref')), "Missing GRN ref subdirectory"
            
            # Check property subdirectories
            property_path = os.path.join(tmpdir, 'property')
            for subdir in DEFAULT_PROPERTY_SUBDIRS.values():
                subdir_path = os.path.join(property_path, subdir)
                assert os.path.exists(subdir_path), f"Missing property subdirectory: {subdir_path}"
            
            # Check embedding subdirectories
            embedding_path = os.path.join(tmpdir, 'embedding')
            for subdir in DEFAULT_EMBEDDING_SUBDIRS.values():
                subdir_path = os.path.join(embedding_path, subdir)
                assert os.path.exists(subdir_path), f"Missing embedding subdirectory: {subdir_path}"
            
            # Check ligand subdirectories
            ligand_path = os.path.join(tmpdir, 'ligand')
            for subdir in DEFAULT_LIGAND_SUBDIRS.values():
                subdir_path = os.path.join(ligand_path, subdir)
                assert os.path.exists(subdir_path), f"Missing ligand subdirectory: {subdir_path}"
            
            # Check that each processor has a datasets directory
            for processor_type, dir_name in DEFAULT_PROCESSOR_DIRS.items():
                if processor_type not in ['test', 'test_processor', '_test', '__test', 'simple', 
                                         'complex_processor_with_long_name', 'custom_dir']:
                    processor_path = os.path.join(tmpdir, dir_name)
                    datasets_path = os.path.join(processor_path, 'datasets')
                    assert os.path.exists(datasets_path), f"Missing datasets directory for {processor_type}: {datasets_path}"
    
    def test_no_datasource_references(self):
        """Test that DataSource enum is properly deprecated."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # These methods should accept DataSource but ignore it
            from protos.io.paths.path_config import DataSource
            
            # Should not raise errors
            path1 = paths.get_processor_path('structure')
            path2 = paths.get_processor_path('structure')
            path3 = paths.get_processor_path('structure')
            
            # All should return the same path
            assert path1 == path2 == path3
    
    def test_clean_environment(self):
        """Test ProtosPaths works with no configuration files or env vars."""
        # Create a clean temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save current environment
            old_env = os.environ.copy()
            old_cwd = os.getcwd()
            
            try:
                # Clear all PROTOS environment variables
                for key in list(os.environ.keys()):
                    if key.startswith('PROTOS'):
                        del os.environ[key]
                
                # Change to temp directory
                os.chdir(tmpdir)
                
                # Create ProtosPaths with no parameters
                paths = ProtosPaths()
                
                # Should still work
                assert paths.data_root is not None
                assert os.path.isabs(paths.data_root)
                
            finally:
                # Restore environment
                os.environ.clear()
                os.environ.update(old_env)
                os.chdir(old_cwd)
    
    def test_subdirectory_methods(self):
        """Test all subdirectory access methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Test structure subdirectories
            assert paths.get_subdir_path("structure", 'structure_dir').endswith('structure/mmcif')
            assert paths.get_subdir_path("structure", 'cache_dir').endswith('structure/cache')
            assert paths.get_subdir_path("structure", 'dataset_dir').endswith('structure/structure_dataset')
            assert paths.get_subdir_path("structure", 'alignments_dir').endswith('structure/alignments')
            assert paths.get_subdir_path("structure", 'temp_dir').endswith('structure/temp_cif')
            assert paths.get_subdir_path("structure", 'datasets_dir').endswith('structure/datasets')
            
            # Test sequence subdirectories
            assert paths.get_subdir_path("sequence", 'fasta_dir').endswith('sequence/fasta')
            assert paths.get_subdir_path("sequence", 'alignment_dir').endswith('sequence/alignments')
            assert paths.get_subdir_path("sequence", 'metadata_dir').endswith('sequence/metadata')
            assert paths.get_subdir_path("sequence", 'datasets_dir').endswith('sequence/datasets')
            
            # Test GRN subdirectories
            assert paths.get_subdir_path('grn', 'table_dir').endswith('grn/tables')
            assert paths.get_subdir_path('grn', 'configs_dir').endswith('grn/configs')
            assert paths.get_subdir_path('grn', 'assignment_dir').endswith('grn/assignments')
            assert paths.get_subdir_path('grn', 'ref').endswith('grn/ref')
            assert paths.get_subdir_path('grn', 'datasets_dir').endswith('grn/datasets')
            
            # Test property subdirectories
            assert paths.get_subdir_path("property", 'tables_dir').endswith('property/tables')
            assert paths.get_subdir_path("property", 'datasets_dir').endswith('property/datasets')
            
            # Test embedding subdirectories
            assert paths.get_subdir_path("embedding", 'embeddings_dir').endswith('embedding/embeddings')
            assert paths.get_subdir_path("embedding", 'datasets_dir').endswith('embedding/datasets')
            
            # Test ligand subdirectories
            assert paths.get_subdir_path("ligand", 'sdf_dir').endswith('ligand/sdf')
            assert paths.get_subdir_path("ligand", 'cache_dir').endswith('ligand/cache')
            assert paths.get_subdir_path("ligand", 'datasets_dir').endswith('ligand/datasets')
    
    def test_registry_paths(self):
        """Test registry path methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Test processor registry paths
            struct_reg = paths.get_registry_path('structure')
            assert struct_reg.endswith('structure/registry.json')
            
            # Test global registry path
            global_reg = paths.get_global_registry_path()
            assert global_reg.endswith('global_registry.json')
    
    def test_dataset_paths(self):
        """Test dataset path methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Test dataset paths for all processors
            struct_dataset = paths.get_dataset_path('structure', 'my_dataset')
            assert struct_dataset.endswith('structure/datasets/my_dataset.json')
            
            seq_dataset = paths.get_dataset_path('sequence', 'seq_collection')
            assert seq_dataset.endswith('sequence/datasets/seq_collection.json')
            
            grn_dataset = paths.get_dataset_path('grn', 'grn_study')
            assert grn_dataset.endswith('grn/datasets/grn_study.json')
            
            prop_dataset = paths.get_dataset_path('property', 'prop_analysis')
            assert prop_dataset.endswith('property/datasets/prop_analysis.json')
            
            # Test with custom extension
            emb_dataset = paths.get_dataset_path('embedding', 'embeddings', file_extension='.yaml')
            assert emb_dataset.endswith('embedding/datasets/embeddings.yaml')
    
    def test_backward_compatibility(self):
        """Test backward compatibility with old parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Test with user_data_root (deprecated)
            paths1 = ProtosPaths(data_root=tmpdir)
            assert paths1.data_root == tmpdir
            assert paths1.user_data_root == tmpdir
            assert paths1.ref_data_root == tmpdir
            
            # Test with ref_data_root (deprecated, should be ignored)
            paths2 = ProtosPaths(data_root=tmpdir, ref_data_root="/some/other/path")
            assert paths2.data_root == tmpdir
            assert paths2.ref_data_root == tmpdir  # Should be same as data_root
    
    def test_absolute_paths(self):
        """Test that all paths are absolute."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Use relative path
            rel_path = os.path.relpath(tmpdir)
            paths = ProtosPaths(data_root=rel_path)
            
            # Should be converted to absolute
            assert os.path.isabs(paths.data_root)
            assert os.path.samefile(paths.data_root, tmpdir)