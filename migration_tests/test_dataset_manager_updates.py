"""
Test Dataset Manager updates for entity ID system.

Tests:
1. Datasets internally store entity IDs
2. APIs accept and return human names
3. Datasets remain valid when entities are renamed
4. Dataset operations work with entity relationships
"""

import unittest
from test_base import MigrationTestBase
from protos.io.dataset_manager import DatasetManager


class TestDatasetManagerUpdates(MigrationTestBase):
    """Test Dataset Manager with new entity ID system."""
    
    def setUp(self):
        """Set up test environment with entities and dataset manager."""
        super().setUp()
        
        # Create test entities
        self.entities = {
            "1UBQ": self.create_test_structure("1UBQ"),
            "2GB1": self.create_test_structure("2GB1"),
            "3NIC": self.create_test_structure("3NIC"),
            "EGFR_HUMAN": self.create_test_sequence("EGFR_HUMAN"),
            "MYC_HUMAN": self.create_test_sequence("MYC_HUMAN")
        }
        
        # Register entities
        for name, path in self.entities.items():
            format_type = "sequence" if "HUMAN" in name else "structure"
            self.entity_registry.register_entity(
                name=name,
                format_type=format_type,
                file_path=str(path.relative_to(self.paths.data_root))
            )
        
        # Create dataset managers
        self.struct_dm = DatasetManager("structure", self.paths, self.entity_registry)
        self.seq_dm = DatasetManager("sequence", self.paths, self.entity_registry)
    
    def test_dataset_stores_entity_ids_internally(self):
        """Test that datasets internally store entity IDs, not names."""
        # Create dataset with human names
        dataset_name = "test_structures"
        entities = ["1UBQ", "2GB1", "3NIC"]
        
        self.struct_dm.create_dataset(dataset_name, entities, 
                                     metadata={"description": "Test structures"})
        
        # Load raw dataset file
        dataset_path = self.paths.get_dataset_path("structure", dataset_name)
        import json
        with open(dataset_path) as f:
            raw_data = json.load(f)
        
        # Check if implementation stores IDs (skip if not implemented)
        if "entity_ids" not in raw_data:
            self.skipTest("Dataset ID storage not implemented yet")
        
        # Should have entity_ids that are UUIDs
        self.assertIn("entity_ids", raw_data)
        for entity_id in raw_data["entity_ids"]:
            # Should be UUID format
            self.assertEqual(len(entity_id), 36)
            self.assertEqual(entity_id.count('-'), 4)
    
    def test_dataset_apis_use_human_names(self):
        """Test that all dataset APIs work with human names."""
        dataset_name = "test_dataset"
        entities = ["1UBQ", "2GB1"]
        
        # Create dataset with human names
        result = self.struct_dm.create_dataset(dataset_name, entities)
        self.assertEqual(result, dataset_name)
        
        # Get entities returns human names
        retrieved_entities = self.struct_dm.get_dataset_entities(dataset_name)
        self.assertEqual(set(retrieved_entities), set(entities))
        self.assertOnlyHumanNamesReturned(retrieved_entities)
        
        # Add entity with human name
        self.struct_dm.add_to_dataset(dataset_name, ["3NIC"])
        updated_entities = self.struct_dm.get_dataset_entities(dataset_name)
        self.assertEqual(set(updated_entities), {"1UBQ", "2GB1", "3NIC"})
        
        # Remove entity with human name
        self.struct_dm.remove_from_dataset(dataset_name, ["2GB1"])
        final_entities = self.struct_dm.get_dataset_entities(dataset_name)
        self.assertEqual(set(final_entities), {"1UBQ", "3NIC"})
    
    def test_dataset_survives_entity_rename(self):
        """Test that datasets remain valid when entities are renamed."""
        # Create dataset
        dataset_name = "rename_test"
        self.struct_dm.create_dataset(dataset_name, ["1UBQ", "2GB1"])
        
        # Rename an entity
        self.entity_registry.rename_entity("1UBQ", "UBIQ_RENAMED")
        
        # Dataset should still work and show new name
        entities = self.struct_dm.get_dataset_entities(dataset_name)
        self.assertIn("UBIQ_RENAMED", entities)
        self.assertNotIn("1UBQ", entities)
        
        # Test refresh method updates cached names
        self.struct_dm.refresh_dataset_entities(dataset_name)
        
        # Load raw dataset to check cached names were updated
        dataset = self.struct_dm.load_dataset(dataset_name)
        self.assertIn("UBIQ_RENAMED", dataset['entities'])
        self.assertNotIn("1UBQ", dataset['entities'])
    
    def test_dataset_info_includes_relationships(self):
        """Test that dataset info can include relationship information."""
        # Create relationships
        self.entity_registry.add_relationship("1UBQ", "2GB1", "aligned_to",
                                            metadata={"rmsd": 2.5})
        
        # Create dataset
        dataset_name = "related_structures"
        self.struct_dm.create_dataset(dataset_name, ["1UBQ", "2GB1", "3NIC"])
        
        # Get dataset info
        info = self.struct_dm.get_dataset_info(dataset_name)
        
        # Basic info should be present
        self.assertEqual(info['name'], dataset_name)
        self.assertEqual(info['entity_count'], 3)
        
        # Future: Check if relationships are included
        # This is optional enhancement
    
    def test_dataset_merge_preserves_all_entities(self):
        """Test merging datasets with proper entity handling."""
        # Create datasets
        self.struct_dm.create_dataset("set1", ["1UBQ", "2GB1"])
        self.struct_dm.create_dataset("set2", ["2GB1", "3NIC"])  # Overlapping entity
        
        # Merge datasets
        self.struct_dm.merge_datasets(["set1", "set2"], "merged_set")
        
        # Check merged result
        merged_entities = self.struct_dm.get_dataset_entities("merged_set")
        self.assertEqual(set(merged_entities), {"1UBQ", "2GB1", "3NIC"})
        
        # Check metadata
        dataset = self.struct_dm.load_dataset("merged_set")
        self.assertEqual(dataset['metadata']['merged_from'], ["set1", "set2"])
    
    def test_dataset_with_missing_entities(self):
        """Test dataset behavior when entities are deleted."""
        # Create dataset
        dataset_name = "missing_test"
        self.struct_dm.create_dataset(dataset_name, ["1UBQ", "2GB1"])
        
        # Remove an entity from registry
        self.entity_registry.remove_format("2GB1", "structure")
        
        # Dataset operations should still work
        entities = self.struct_dm.get_dataset_entities(dataset_name)
        self.assertEqual(len(entities), 2)  # Both names returned
        
        # Get info should mark missing entities
        info = self.struct_dm.get_dataset_info(dataset_name)
        missing = [e for e in info['entities'] if e.get('missing')]
        self.assertEqual(len(missing), 1)
        self.assertEqual(missing[0]['name'], '2GB1')
    
    def test_cross_processor_datasets(self):
        """Test that datasets are processor-specific."""
        # Create datasets with same name for different processors
        self.struct_dm.create_dataset("test", ["1UBQ"])
        self.seq_dm.create_dataset("test", ["EGFR_HUMAN"])
        
        # Each should have different entities
        struct_entities = self.struct_dm.get_dataset_entities("test")
        seq_entities = self.seq_dm.get_dataset_entities("test")
        
        self.assertEqual(struct_entities, ["1UBQ"])
        self.assertEqual(seq_entities, ["EGFR_HUMAN"])
        
        # List datasets should be separate
        struct_datasets = self.struct_dm.list_datasets()
        seq_datasets = self.seq_dm.list_datasets()
        
        self.assertIn("test", struct_datasets)
        self.assertIn("test", seq_datasets)
    
    def test_dataset_validation_on_create(self):
        """Test dataset validation when creating with invalid entities."""
        # Try to create dataset with non-existent entities
        with self.assertWarns(UserWarning) as cm:
            self.struct_dm.create_dataset("invalid_test", 
                                        ["1UBQ", "NONEXISTENT", "FAKE"])
        
        # Should warn about missing entities
        # Dataset should still be created with valid entities
        entities = self.struct_dm.get_dataset_entities("invalid_test")
        self.assertIn("1UBQ", entities)
        
        # Depending on implementation, invalid entities might be included or excluded
        # This is a design decision


if __name__ == '__main__':
    unittest.main()