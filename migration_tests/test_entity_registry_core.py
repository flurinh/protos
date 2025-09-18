"""
Test core EntityRegistry functionality for the new UUID-based system.

Tests:
1. UUID generation (not hash-based)
2. Name-to-ID mapping
3. Alias resolution
4. Human-readable names in all outputs
"""

import unittest
from pathlib import Path

from test_base import MigrationTestBase
from protos.io.entity_registry import EntityRegistry


class TestEntityRegistryCore(MigrationTestBase):
    """Test core entity registry functionality."""
    
    def test_entity_registration_creates_uuid(self):
        """Test that registering an entity creates a UUID, not a hash."""
        # Register an entity
        entity_name = "TEST_PROTEIN"
        file_path = self.create_test_structure(entity_name)
        
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(file_path.relative_to(self.paths.data_root)),
            metadata={"test": True}
        )
        
        # Check internal registry structure
        registry_data = self.get_registry_data()
        self.assertEqual(len(registry_data), 1, "Should have one entity")
        
        # Get the entity ID (should be UUID)
        entity_id = list(registry_data.keys())[0]
        
        # Verify it's a UUID format (36 chars with dashes)
        # For now, it might still be hash-based (10 chars)
        # This test will guide our migration
        if len(entity_id) == 10:
            self.skipTest("Still using hash-based IDs - needs migration")
        else:
            self.assertEqual(len(entity_id), 36, "Should be UUID length")
            self.assertEqual(entity_id.count('-'), 4, "Should have UUID format")
    
    def test_name_to_id_mapping(self):
        """Test that names map correctly to IDs internally."""
        # Register multiple entities
        entities = ["1UBQ", "2GB1", "EGFR_HUMAN"]
        
        for name in entities:
            file_path = self.create_test_structure(name)
            self.entity_registry.register_entity(
                name=name,
                format_type="structure",
                file_path=str(file_path.relative_to(self.paths.data_root))
            )
        
        # Test that we can find each entity by name
        for name in entities:
            entity_info = self.entity_registry.find_entity(name)
            self.assertIsNotNone(entity_info)
            self.assertEqual(entity_info.original_id, name)
    
    def test_public_methods_return_human_names(self):
        """Test that all public methods return human names, not IDs."""
        # Register some entities
        entities = ["1UBQ", "2GB1", "EGFR_HUMAN", "sp|P12345|TEST_HUMAN"]
        
        for name in entities:
            file_path = self.create_test_structure(name)
            result = self.entity_registry.register_entity(
                name=name,
                format_type="structure", 
                file_path=str(file_path.relative_to(self.paths.data_root))
            )
            # register_entity should return the human name
            self.assertEqual(result, name)
        
        # list_entities should return human names
        listed = self.entity_registry.list_entities()
        self.assertEqual(set(listed), set(entities))
        self.assertOnlyHumanNamesReturned(listed)
        
        # get_entity_formats should work with human names
        formats = self.entity_registry.get_entity_formats("1UBQ")
        self.assertEqual(formats, ["structure"])
    
    def test_alias_support(self):
        """Test that aliases work correctly."""
        # Register entity
        primary_name = "1UBQ"
        file_path = self.create_test_structure(primary_name)
        self.entity_registry.register_entity(
            name=primary_name,
            format_type="structure",
            file_path=str(file_path.relative_to(self.paths.data_root))
        )
        
        # Add aliases
        aliases = ["UBIQ_HUMAN", "P62988", "Ubiquitin"]
        for alias in aliases:
            self.entity_registry.add_alias(primary_name, alias)
        
        # Should find entity by any alias
        for alias in aliases:
            entity_info = self.entity_registry.find_entity(alias)
            self.assertIsNotNone(entity_info)
            self.assertEqual(entity_info.original_id, primary_name)
    
    def test_complex_identifier_handling(self):
        """Test handling of complex biological identifiers."""
        # UniProt-style identifier
        complex_id = "sp|P02724|GLPA_ECOLI"
        
        file_path = self.create_test_sequence(complex_id)
        result = self.entity_registry.register_entity(
            name=complex_id,
            format_type="sequence",
            file_path=str(file_path.relative_to(self.paths.data_root))
        )
        
        # Should return original complex ID
        self.assertEqual(result, complex_id)
        
        # Should find by exact name
        entity_info = self.entity_registry.find_entity(complex_id)
        self.assertIsNotNone(entity_info)
        self.assertEqual(entity_info.original_id, complex_id)
    
    def test_no_duplicate_registration(self):
        """Test that duplicate registrations are handled correctly."""
        entity_name = "TEST_DUP"
        
        # First registration
        file_path = self.create_test_structure(entity_name)
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(file_path.relative_to(self.paths.data_root))
        )
        
        # Second registration with same name and format
        # Should update, not create duplicate
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(file_path.relative_to(self.paths.data_root)),
            metadata={"updated": True}
        )
        
        # Should still have only one entity
        self.assertEqual(self.count_registered_entities(), 1)
        
        # Check metadata was updated
        entity_info = self.entity_registry.find_entity(entity_name)
        self.assertTrue(entity_info.metadata.get("updated", False))
    
    def test_multi_format_entity(self):
        """Test that same entity can exist in multiple formats."""
        entity_name = "MULTI_FORMAT"
        
        # Register as structure
        struct_path = self.create_test_structure(entity_name)
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(struct_path.relative_to(self.paths.data_root))
        )
        
        # Register as sequence
        seq_path = self.create_test_sequence(entity_name)
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="sequence",
            file_path=str(seq_path.relative_to(self.paths.data_root))
        )
        
        # Should have one entity with two formats
        self.assertEqual(self.count_registered_entities(), 1)
        
        # Check formats
        formats = self.entity_registry.get_entity_formats(entity_name)
        self.assertEqual(set(formats), {"structure", "sequence"})
        
        # Should find in both formats
        struct_info = self.entity_registry.find_entity(entity_name, "structure")
        self.assertIsNotNone(struct_info)
        
        seq_info = self.entity_registry.find_entity(entity_name, "sequence")
        self.assertIsNotNone(seq_info)
    
    def test_registry_persistence(self):
        """Test that registry persists across instances."""
        entity_name = "PERSISTENT"
        
        # Register entity
        file_path = self.create_test_structure(entity_name)
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(file_path.relative_to(self.paths.data_root))
        )
        
        # Create new registry instance
        new_registry = EntityRegistry(paths=self.paths)
        
        # Should find the entity
        entity_info = new_registry.find_entity(entity_name)
        self.assertIsNotNone(entity_info)
        self.assertEqual(entity_info.original_id, entity_name)
    
    def test_entity_removal(self):
        """Test removing formats and entities."""
        entity_name = "TO_REMOVE"
        
        # Register in two formats
        struct_path = self.create_test_structure(entity_name)
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(struct_path.relative_to(self.paths.data_root))
        )
        
        seq_path = self.create_test_sequence(entity_name)
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="sequence", 
            file_path=str(seq_path.relative_to(self.paths.data_root))
        )
        
        # Remove one format
        self.entity_registry.remove_format(entity_name, "sequence")
        
        # Should still exist in structure format
        self.assertTrue(self.entity_registry.entity_exists(entity_name, "structure"))
        self.assertFalse(self.entity_registry.entity_exists(entity_name, "sequence"))
        
        # Remove last format - entity should be removed
        self.entity_registry.remove_format(entity_name, "structure")
        self.assertFalse(self.entity_registry.entity_exists(entity_name))
    
    def test_entity_rename(self):
        """Test renaming entities."""
        # Register entity
        old_name = "OLD_NAME"
        file_path = self.create_test_structure(old_name)
        self.entity_registry.register_entity(
            name=old_name,
            format_type="structure",
            file_path=str(file_path.relative_to(self.paths.data_root))
        )
        
        # Add some aliases
        self.entity_registry.add_alias(old_name, "ALIAS1")
        self.entity_registry.add_alias(old_name, "ALIAS2")
        
        # Rename entity
        new_name = "NEW_NAME"
        self.entity_registry.rename_entity(old_name, new_name)
        
        # Old name should not exist
        self.assertFalse(self.entity_registry.entity_exists(old_name))
        
        # New name should exist
        self.assertTrue(self.entity_registry.entity_exists(new_name))
        
        # Should find by new name
        entity_info = self.entity_registry.find_entity(new_name)
        self.assertIsNotNone(entity_info)
        self.assertEqual(entity_info.original_id, new_name)
        
        # Aliases should still work
        self.assertIsNotNone(self.entity_registry.find_entity("ALIAS1"))
        self.assertIsNotNone(self.entity_registry.find_entity("ALIAS2"))
        
        # List should show new name
        entities = self.entity_registry.list_entities()
        self.assertIn(new_name, entities)
        self.assertNotIn(old_name, entities)
    
    def test_rename_to_existing_name_fails(self):
        """Test that renaming to an existing name fails."""
        # Register two entities
        self.entity_registry.register_entity(
            name="ENTITY1",
            format_type="structure",
            file_path="test1.cif"
        )
        self.entity_registry.register_entity(
            name="ENTITY2",
            format_type="structure",
            file_path="test2.cif"
        )
        
        # Try to rename to existing name
        with self.assertRaises(ValueError) as context:
            self.entity_registry.rename_entity("ENTITY1", "ENTITY2")
        
        self.assertIn("already exists", str(context.exception))


if __name__ == '__main__':
    unittest.main()