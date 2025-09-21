"""
Test BaseProcessor integration with new entity system.

Tests:
1. Processors correctly register entities
2. Dataset operations work through processors
3. Entity relationships are tracked
4. Human names used in all APIs
"""

import unittest
from unittest.mock import Mock, MagicMock
from pathlib import Path
from ...test_base import MigrationTestBase
from protos.io.core.base_processor import BaseProcessor


class MockProcessor(BaseProcessor):
    """Mock processor for testing base functionality."""
    
    def __init__(self, name="test_processor", paths=None):
        super().__init__(name=name, paths=paths)
        self._storage = {}  # Maps entity_id -> data
    
    def _get_processor_type(self) -> str:
        """Return test processor type."""
        return "test"
    
    def load_entity(self, name: str):
        """Load entity by name."""
        # Check if entity is registered
        entity_info = self.entity_registry.find_entity(name, format_type="test")
        if not entity_info:
            return None
        
        # Return from storage using entity ID
        return self._storage.get(entity_info.entity_id)
    
    def save_entity(self, name: str, data, metadata=None):
        """Save entity with name."""
        # Register with entity registry first
        file_path = f"test/{name}.test"
        self.entity_registry.register_entity(
            name=name,
            format_type="test",
            file_path=file_path,
            metadata=metadata or {}
        )
        
        # Get the entity info to get ID
        entity_info = self.entity_registry.find_entity(name, format_type="test")
        if entity_info:
            # Store data by entity ID
            self._storage[entity_info.entity_id] = data


class TestProcessorIntegration(MigrationTestBase):
    """Test processor integration with entity system."""
    
    def setUp(self):
        """Set up test environment."""
        super().setUp()
        
        # Create mock processor
        self.processor = MockProcessor(paths=self.paths)
    
    def test_processor_registers_entities(self):
        """Test that save_entity registers with entity registry."""
        # Save an entity
        test_data = {"value": 42, "name": "test"}
        self.processor.save_entity("TEST1", test_data, metadata={"type": "test"})
        
        # Check entity was registered
        self.assertTrue(self.processor.entity_registry.entity_exists("TEST1"))
        
        # Check format is correct
        formats = self.processor.entity_registry.get_entity_formats("TEST1")
        self.assertEqual(formats, ["test"])
        
        # Check metadata was stored
        metadata = self.processor.entity_registry.get_entity_metadata("TEST1", "test")
        self.assertEqual(metadata["type"], "test")
    
    def test_processor_loads_registered_entities(self):
        """Test that load_entity only works for registered entities."""
        # Try to load non-existent entity
        result = self.processor.load_entity("NONEXISTENT")
        self.assertIsNone(result)
        
        # Save an entity
        test_data = {"value": 123}
        self.processor.save_entity("TEST2", test_data)
        
        # Load it back
        loaded = self.processor.load_entity("TEST2")
        self.assertEqual(loaded, test_data)
    
    def test_dataset_operations_through_processor(self):
        """Test dataset operations work correctly."""
        # Save some entities
        for i in range(3):
            self.processor.save_entity(f"ENTITY{i}", {"id": i})
        
        # Create dataset
        self.processor.create_dataset("test_set", ["ENTITY0", "ENTITY1", "ENTITY2"],
                                    metadata={"purpose": "testing"})
        
        # List datasets
        datasets = self.processor.list_datasets()
        self.assertIn("test_set", datasets)
        
        # Load dataset
        loaded_data = self.processor.load_dataset("test_set")
        self.assertEqual(len(loaded_data), 3)
        self.assertIn("ENTITY0", loaded_data)
        self.assertEqual(loaded_data["ENTITY0"]["id"], 0)
        
        # Get dataset info
        info = self.processor.get_dataset_info("test_set")
        self.assertEqual(info["entity_count"], 3)
        self.assertEqual(info["metadata"]["purpose"], "testing")
    
    def test_dataset_handles_missing_entities(self):
        """Test dataset handles missing entities gracefully."""
        # Create entities
        self.processor.save_entity("EXISTS1", {"data": 1})
        self.processor.save_entity("EXISTS2", {"data": 2})
        
        # Create dataset including non-existent entity
        with self.assertWarns(UserWarning):
            self.processor.create_dataset("mixed_set", 
                                        ["EXISTS1", "MISSING", "EXISTS2"])
        
        # Load dataset should skip missing
        loaded = self.processor.load_dataset("mixed_set")
        self.assertEqual(len(loaded), 2)
        self.assertIn("EXISTS1", loaded)
        self.assertIn("EXISTS2", loaded)
        self.assertNotIn("MISSING", loaded)
    
    def test_entity_ids_not_exposed(self):
        """Test that entity IDs are never exposed in processor APIs."""
        # Save entity
        self.processor.save_entity("HUMAN_NAME", {"test": True})
        
        # Create dataset
        self.processor.create_dataset("test_dataset", ["HUMAN_NAME"])
        
        # All operations should return human names
        datasets = self.processor.list_datasets()
        self.assertOnlyHumanNamesReturned(datasets)
        
        loaded = self.processor.load_dataset("test_dataset")
        self.assertOnlyHumanNamesReturned(list(loaded.keys()))
        
        info = self.processor.get_dataset_info("test_dataset")
        entity_names = [e["name"] for e in info["entities"]]
        self.assertOnlyHumanNamesReturned(entity_names)
    
    def test_processor_type_specific_datasets(self):
        """Test datasets are processor-type specific."""
        # Create another mock processor with different type
        class OtherProcessor(MockProcessor):
            def _get_processor_type(self):
                return "other"
        
        other_processor = OtherProcessor(paths=self.paths)
        
        # Create datasets with same name in each processor
        self.processor.save_entity("ENTITY1", {"proc": "test"})
        self.processor.create_dataset("shared_name", ["ENTITY1"])
        
        other_processor.save_entity("ENTITY2", {"proc": "other"})
        other_processor.create_dataset("shared_name", ["ENTITY2"])
        
        # Each should see only their own dataset
        test_loaded = self.processor.load_dataset("shared_name")
        self.assertEqual(list(test_loaded.keys()), ["ENTITY1"])
        
        other_loaded = other_processor.load_dataset("shared_name")
        self.assertEqual(list(other_loaded.keys()), ["ENTITY2"])
    
    def test_processor_handles_entity_rename(self):
        """Test processor handles entity renames correctly."""
        # Save entity
        self.processor.save_entity("OLD_NAME", {"value": "test"})
        
        # Create dataset with entity
        self.processor.create_dataset("rename_dataset", ["OLD_NAME"])
        
        # Rename entity
        self.processor.entity_registry.rename_entity("OLD_NAME", "NEW_NAME")
        
        # Load through processor should work
        loaded = self.processor.load_entity("NEW_NAME")
        self.assertEqual(loaded, {"value": "test"})
        
        # Dataset should show new name
        dataset_data = self.processor.load_dataset("rename_dataset")
        self.assertIn("NEW_NAME", dataset_data)
        self.assertNotIn("OLD_NAME", dataset_data)


class TestStructureProcessorIntegration(MigrationTestBase):
    """Test StructureProcessor specific integration."""
    
    def setUp(self):
        """Set up test environment."""
        super().setUp()
        
        # Import StructureProcessor if available
        try:
            from protos.processing.structure.structure_processor import StructureProcessor
            self.StructureProcessor = StructureProcessor
        except ImportError:
            self.skipTest("StructureProcessor not available")
    
    def test_structure_processor_dual_format(self):
        """Test StructureProcessor handles CIF and PKL formats."""
        # This test would verify:
        # 1. Saving CIF registers with entity registry
        # 2. Creating PKL cache updates entity with second format
        # 3. Loading prefers PKL when available
        # 4. Both formats tracked in entity registry
        
        # Skip for now as it requires actual StructureProcessor implementation
        self.skipTest("StructureProcessor implementation needed")


if __name__ == '__main__':
    unittest.main()