"""
Test entity relationship tracking functionality.

Tests:
1. Add relationships between entities
2. Query relationships by type and direction
3. Get related entities
4. Remove relationships
5. Handle symmetric relationships
"""

import unittest
from test_base import MigrationTestBase


class TestEntityRelationships(MigrationTestBase):
    """Test entity relationship functionality."""
    
    def setUp(self):
        """Set up test entities for relationship testing."""
        super().setUp()
        
        # Create some test entities
        self.entities = {
            "1UBQ": self.create_test_structure("1UBQ"),
            "1UBQ_PKL": self.create_test_structure("1UBQ_PKL"),
            "1UBQ_chainA": self.create_test_structure("1UBQ_chainA"),
            "1UBQ_aligned": self.create_test_structure("1UBQ_aligned"),
            "GRN_1UBQ": self.create_test_table("grn", "GRN_1UBQ")
        }
        
        # Register all entities
        for name, path in self.entities.items():
            format_type = "grn" if "GRN" in name else "structure"
            self.entity_registry.register_entity(
                name=name,
                format_type=format_type,
                file_path=str(path.relative_to(self.paths.data_root))
            )
    
    def test_add_relationship(self):
        """Test adding relationships between entities."""
        # Add derived_from relationship
        self.entity_registry.add_relationship(
            source_name="1UBQ",
            target_name="1UBQ_PKL",
            rel_type="derived_from",
            metadata={"method": "to_pkl", "compression": "gzip"}
        )
        
        # Check it was added
        relationships = self.entity_registry.get_relationships("1UBQ_PKL")
        self.assertEqual(len(relationships), 1)
        self.assertEqual(relationships[0]['type'], 'derived_from')
        self.assertEqual(relationships[0]['source_name'], '1UBQ')
        self.assertEqual(relationships[0]['metadata']['method'], 'to_pkl')
    
    def test_query_relationships_by_direction(self):
        """Test querying relationships by direction."""
        # Set up relationships
        self.entity_registry.add_relationship("1UBQ", "1UBQ_chainA", "subset_of")
        self.entity_registry.add_relationship("GRN_1UBQ", "1UBQ", "annotated_by")
        
        # Query incoming relationships for 1UBQ
        incoming = self.entity_registry.get_relationships("1UBQ", direction="incoming")
        self.assertEqual(len(incoming), 1)
        self.assertEqual(incoming[0]['type'], 'annotated_by')
        self.assertEqual(incoming[0]['source_name'], 'GRN_1UBQ')
        
        # Query outgoing relationships for 1UBQ  
        outgoing = self.entity_registry.get_relationships("1UBQ", direction="outgoing")
        self.assertEqual(len(outgoing), 1)
        self.assertEqual(outgoing[0]['type'], 'contains')  # Inverse of subset_of
        self.assertEqual(outgoing[0]['target_name'], '1UBQ_chainA')
    
    def test_query_relationships_by_type(self):
        """Test filtering relationships by type."""
        # Add multiple relationships
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from")
        self.entity_registry.add_relationship("1UBQ", "1UBQ_chainA", "subset_of")
        self.entity_registry.add_relationship("1UBQ", "1UBQ_aligned", "aligned_to")
        
        # Query only derived_from relationships
        derived = self.entity_registry.get_relationships("1UBQ_PKL", rel_type="derived_from")
        self.assertEqual(len(derived), 1)
        self.assertEqual(derived[0]['source_name'], '1UBQ')
        
        # Query non-existent type
        version = self.entity_registry.get_relationships("1UBQ", rel_type="version_of")
        self.assertEqual(len(version), 0)
    
    def test_get_related_entities(self):
        """Test getting related entity names."""
        # Create relationship network
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from")
        self.entity_registry.add_relationship("1UBQ", "1UBQ_chainA", "subset_of")
        self.entity_registry.add_relationship("GRN_1UBQ", "1UBQ", "annotated_by")
        
        # Get all related entities
        related = self.entity_registry.get_related_entities("1UBQ")
        self.assertEqual(set(related), {"1UBQ_PKL", "1UBQ_chainA", "GRN_1UBQ"})
        
        # Get only incoming related
        incoming = self.entity_registry.get_related_entities("1UBQ", direction="incoming")
        self.assertEqual(incoming, ["GRN_1UBQ"])
        
        # Get only specific type
        derived = self.entity_registry.get_related_entities("1UBQ_PKL", rel_type="derived_from")
        self.assertEqual(derived, ["1UBQ"])
    
    def test_symmetric_relationships(self):
        """Test symmetric relationships like aligned_to."""
        # Add symmetric relationship
        self.entity_registry.add_relationship("1UBQ", "1UBQ_aligned", "aligned_to")
        
        # Check both directions show aligned_to (not aligned_with)
        rel1 = self.entity_registry.get_relationships("1UBQ_aligned", direction="incoming")
        self.assertEqual(rel1[0]['type'], 'aligned_to')
        
        rel2 = self.entity_registry.get_relationships("1UBQ", direction="outgoing") 
        self.assertEqual(rel2[0]['type'], 'aligned_to')  # Symmetric
    
    def test_remove_relationship(self):
        """Test removing relationships."""
        # Add relationships
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from")
        self.entity_registry.add_relationship("1UBQ", "1UBQ_chainA", "subset_of")
        
        # Verify they exist
        self.assertEqual(len(self.entity_registry.get_relationships("1UBQ_PKL")), 1)
        
        # Remove one relationship
        self.entity_registry.remove_relationship("1UBQ", "1UBQ_PKL", "derived_from")
        
        # Check it's gone
        self.assertEqual(len(self.entity_registry.get_relationships("1UBQ_PKL")), 0)
        
        # Other relationship should still exist
        self.assertEqual(len(self.entity_registry.get_relationships("1UBQ_chainA")), 1)
    
    def test_duplicate_relationship_prevention(self):
        """Test that duplicate relationships are not added."""
        # Add relationship
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from",
                                            metadata={"version": 1})
        
        # Try to add same relationship again
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from",
                                            metadata={"version": 2})
        
        # Should still have only one relationship
        relationships = self.entity_registry.get_relationships("1UBQ_PKL")
        self.assertEqual(len(relationships), 1)
        # Original metadata should be preserved
        self.assertEqual(relationships[0]['metadata']['version'], 1)
    
    def test_invalid_relationship_type(self):
        """Test handling of invalid relationship types."""
        with self.assertRaises(ValueError) as context:
            self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "invalid_type")
        
        self.assertIn("Invalid relationship type", str(context.exception))
    
    def test_relationship_with_nonexistent_entities(self):
        """Test adding relationships to non-existent entities."""
        # Try with non-existent source
        with self.assertRaises(ValueError) as context:
            self.entity_registry.add_relationship("NONEXISTENT", "1UBQ", "derived_from")
        self.assertIn("Source entity", str(context.exception))
        
        # Try with non-existent target
        with self.assertRaises(ValueError) as context:
            self.entity_registry.add_relationship("1UBQ", "NONEXISTENT", "derived_from")
        self.assertIn("Target entity", str(context.exception))
    
    def test_relationships_persist(self):
        """Test that relationships persist across registry reloads."""
        # Add relationship
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from",
                                            metadata={"persistent": True})
        
        # Create new registry instance
        new_registry = self.entity_registry.__class__(paths=self.paths)
        
        # Check relationship still exists
        relationships = new_registry.get_relationships("1UBQ_PKL")
        self.assertEqual(len(relationships), 1)
        self.assertEqual(relationships[0]['metadata']['persistent'], True)
    
    def test_no_entity_ids_exposed(self):
        """Test that relationship APIs never expose entity IDs."""
        # Add relationships
        self.entity_registry.add_relationship("1UBQ", "1UBQ_PKL", "derived_from")
        
        # Get relationships - should only have human names
        relationships = self.entity_registry.get_relationships("1UBQ_PKL")
        
        # Check that no UUIDs are exposed
        for rel in relationships:
            self.assertNotIn('source', rel)  # Should be source_name
            self.assertNotIn('target', rel)  # Should be target_name
            self.assertIn('source_name', rel)
            self.assertIn('target_name', rel)
            
            # Verify names don't look like UUIDs
            self.assertOnlyHumanNamesReturned([rel['source_name'], rel['target_name']])


if __name__ == '__main__':
    unittest.main()