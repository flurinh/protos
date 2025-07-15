"""
Tests for new BaseProcessor implementation with ProtosPaths integration.
"""

import pytest
import tempfile
import json
from pathlib import Path
from typing import Any, Optional
import pandas as pd
import numpy as np

from protos.core.base_processor import BaseProcessor
from protos.io.paths import ProtosPaths


class TestProcessor(BaseProcessor):
    """Test processor for testing BaseProcessor functionality."""
    
    def load_entity(self, name: str):
        """Load test entity."""
        # Try to find file
        test_file = self.data_path / "data" / f"{name}.test"
        if test_file.exists():
            # Auto-register if not already registered
            if not self.entity_exists(name):
                self.entity_registry.register_entity(
                    name=name,
                    format_type=self.processor_type,
                    file_path=str(test_file.relative_to(self.paths.data_root)),
                    metadata={"auto_discovered": True}
                )
            return test_file.read_text()
        
        # Check if registered
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info:
            file_path = Path(entity_info.file_path)
            if not file_path.is_absolute():
                file_path = Path(self.paths.data_root) / file_path
            if file_path.exists():
                return file_path.read_text()
        
        return None
    
    def save_entity(self, name: str, data: Any, metadata: Optional[dict] = None):
        """Save test entity."""
        # Create subdirectory
        data_dir = self.data_path / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        
        # Save file
        file_path = data_dir / f"{name}.test"
        file_path.write_text(str(data))
        
        # Register entity
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=str(file_path.relative_to(self.paths.data_root)),
            metadata=metadata or {}
        )


class TestBaseProcessorNew:
    """Test new BaseProcessor implementation."""
    
    def test_zero_configuration(self):
        """Test processor works with zero configuration."""
        # Should work without any setup
        processor = TestProcessor()
        
        # Should have created all components
        assert processor.paths is not None
        assert processor.entity_registry is not None
        assert processor.dataset_manager is not None
        assert processor.data_path.exists()
    
    def test_accepts_protospaths(self):
        """Test processor accepts ProtosPaths instance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = TestProcessor(paths=paths)
            
            # Should use provided paths
            assert processor.paths == paths
            assert str(processor.paths.data_root) == tmpdir
    
    def test_processor_type_detection(self):
        """Test automatic processor type detection."""
        # Test with our TestProcessor
        processor = TestProcessor()
        assert processor.processor_type == "test"
        
        # Test with mock classes that implement abstract methods
        class StructureProcessor(BaseProcessor):
            def load_entity(self, name: str):
                return None
            def save_entity(self, name: str, data: Any, metadata=None):
                pass
        
        class GRNProcessor(BaseProcessor):
            def load_entity(self, name: str):
                return None
            def save_entity(self, name: str, data: Any, metadata=None):
                pass
        
        class SequenceProcessor(BaseProcessor):
            def load_entity(self, name: str):
                return None
            def save_entity(self, name: str, data: Any, metadata=None):
                pass
        
        struct = StructureProcessor()
        assert struct.processor_type == "structure"
        
        grn = GRNProcessor()
        assert grn.processor_type == "grn"
        
        seq = SequenceProcessor()
        assert seq.processor_type == "sequence"
    
    def test_save_and_load_entity(self):
        """Test saving and loading entities with human names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Save entity
            test_data = {"key": "value", "number": 42}
            save_path = processor.save_entity("my_test", test_data)
            
            # Should be registered
            assert processor.entity_exists("my_test")
            
            # Load entity
            loaded = processor.load_entity("my_test")
            assert loaded == str(test_data)  # Base implementation saves as string
    
    def test_list_entities(self):
        """Test listing entities returns human names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Save some entities
            processor.save_entity("test1", "data1")
            processor.save_entity("test2", "data2")
            processor.save_entity("test3", "data3")
            
            # List should return human names
            entities = processor.list_entities()
            assert "test1" in entities
            assert "test2" in entities
            assert "test3" in entities
            assert len(entities) >= 3
    
    def test_delete_entity(self):
        """Test deleting entities."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Save and verify
            processor.save_entity("to_delete", "data")
            assert processor.entity_exists("to_delete")
            
            # Delete
            success = processor.delete_entity("to_delete")
            assert success
            
            # Should be gone
            assert not processor.entity_exists("to_delete")
            loaded = processor.load_entity("to_delete")
            assert loaded is None
    
    def test_drag_and_drop_workflow(self):
        """Test drag-and-drop file discovery."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Manually create a file (simulating drag-and-drop)
            data_dir = Path(tmpdir) / "test" / "data"
            data_dir.mkdir(parents=True, exist_ok=True)
            
            test_file = data_dir / "dropped_file.test"
            test_file.write_text("dropped content")
            
            # Load should find and register it
            loaded = processor.load_entity("dropped_file")
            assert loaded == "dropped content"
            
            # Should now be registered
            assert processor.entity_exists("dropped_file")
    
    def test_dataset_operations(self):
        """Test dataset management."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Create entities
            processor.save_entity("entity1", "data1")
            processor.save_entity("entity2", "data2")
            processor.save_entity("entity3", "data3")
            
            # Create dataset
            dataset_name = processor.create_dataset(
                "my_dataset",
                ["entity1", "entity2"],
                {"description": "Test dataset"}
            )
            
            assert dataset_name == "my_dataset"
            
            # List datasets
            datasets = processor.list_datasets()
            assert "my_dataset" in datasets
            
            # Get dataset info
            info = processor.get_dataset_info("my_dataset")
            assert info["name"] == "my_dataset"
            assert info["entity_count"] == 2
            
            # Add to dataset
            processor.add_to_dataset("my_dataset", ["entity3"])
            info = processor.get_dataset_info("my_dataset")
            assert info["entity_count"] == 3
            
            # Remove from dataset
            processor.remove_from_dataset("my_dataset", ["entity1"])
            info = processor.get_dataset_info("my_dataset")
            assert info["entity_count"] == 2
    
    def test_metadata_tracking(self):
        """Test metadata is properly initialized."""
        processor = TestProcessor(name="test_meta")
        
        assert processor.metadata["name"] == "test_meta"
        assert processor.metadata["processor_type"] == "TestProcessor"
        assert "created_at" in processor.metadata
        assert "data_path" in processor.metadata
    
    def test_subdirectory_paths(self):
        """Test getting subdirectory paths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Create a processor that would have subdirectories
            class StructureTestProcessor(BaseProcessor):
                def _get_processor_type(self):
                    return "structure"
                def load_entity(self, name: str):
                    return None
                def save_entity(self, name: str, data: Any, metadata=None):
                    pass
            
            processor = StructureTestProcessor(paths=paths)
            
            # Should be able to get subdirectory paths
            cache_path = processor.get_subdirectory_path("cache_dir")
            assert "cache" in str(cache_path)
    
    def test_no_path_parameters(self):
        """Test processor constructor has no path parameters."""
        import inspect
        
        # Get init signature
        sig = inspect.signature(BaseProcessor.__init__)
        params = sig.parameters
        
        # Should not have any path-related parameters except 'paths'
        param_names = list(params.keys())
        
        # These parameters should NOT exist
        assert "data_root" not in param_names
        assert "user_data_root" not in param_names  
        assert "ref_data_root" not in param_names
        assert "processor_data_dir" not in param_names
        assert "data_dir" not in param_names
        
        # Only 'paths' parameter is allowed
        assert "paths" in param_names
    
    def test_save_load_data_formats(self):
        """Test saving and loading data in different formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Test DataFrame
            df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
            processor.save_data("test_df", df, file_format='csv')
            loaded_df = processor.load_data("test_df", file_format='csv')
            # CSV adds index column, so we need to handle that
            loaded_df = loaded_df.drop(columns=['Unnamed: 0'], errors='ignore')
            pd.testing.assert_frame_equal(df, loaded_df)
            
            # Test JSON
            data_dict = {'key': 'value', 'number': 42}
            processor.save_data("test_json", data_dict, file_format='json')
            loaded_dict = processor.load_data("test_json", file_format='json')
            assert loaded_dict == data_dict
            
            # Test Pickle
            complex_data = {'array': np.array([1, 2, 3]), 'string': 'test'}
            processor.save_data("test_pkl", complex_data, file_format='pkl')
            loaded_pkl = processor.load_data("test_pkl", file_format='pkl')
            assert loaded_pkl['string'] == complex_data['string']
            np.testing.assert_array_equal(loaded_pkl['array'], complex_data['array'])
            
            # Test NumPy
            arr = np.array([[1, 2], [3, 4]])
            processor.save_data("test_npy", arr, file_format='npy')
            loaded_arr = processor.load_data("test_npy", file_format='npy')
            np.testing.assert_array_equal(arr, loaded_arr)
    
    def test_is_dataset_available(self):
        """Test dataset availability checking."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = TestProcessor(paths=paths)
            
            # Create some entities and a dataset
            processor.save_entity("entity1", "data1")
            processor.save_entity("entity2", "data2")
            
            processor.create_dataset("test_dataset", ["entity1", "entity2"])
            
            # Check availability
            assert processor.is_dataset_available("test_dataset")
            assert not processor.is_dataset_available("non_existent")