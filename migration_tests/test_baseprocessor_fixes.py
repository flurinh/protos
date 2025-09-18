"""
Test BaseProcessor fixes for identified weaknesses.

Tests:
1. _sanitize_filename is used consistently
2. Complex identifiers are handled correctly
3. Error handling improvements
4. Path resolution consistency
"""

import unittest
import logging
from pathlib import Path
from test_base import MigrationTestBase
from protos.core.base_processor import BaseProcessor


class FixedMockProcessor(BaseProcessor):
    """Mock processor with fixes applied."""
    
    def _get_processor_type(self) -> str:
        return "test"
    
    def load_entity(self, name: str):
        """Load entity with proper filename sanitization."""
        entity_info = self.entity_registry.find_entity(name, format_type="test")
        if not entity_info:
            return None
        
        # Get the actual file path from entity info
        file_path = Path(entity_info.file_path)
        if not file_path.is_absolute():
            file_path = Path(self.paths.data_root) / file_path
        
        # Mock implementation - return data if file exists
        if file_path.exists():
            return {"loaded": True, "name": name}
        
        # If file doesn't exist, raise exception
        raise FileNotFoundError(f"Entity file not found: {file_path}")
    
    def save_entity(self, name: str, data, metadata=None):
        """Save entity with proper filename sanitization."""
        # Sanitize filename
        safe_name = self._sanitize_filename(name)
        file_path = f"test/{safe_name}.test"
        
        # Register with original name
        self.entity_registry.register_entity(
            name=name,  # Original name preserved
            format_type="test",
            file_path=file_path,  # Sanitized path
            metadata=metadata or {}
        )
        
        # Simulate file creation
        actual_path = self.data_path / f"{safe_name}.test"
        actual_path.parent.mkdir(parents=True, exist_ok=True)
        actual_path.write_text(str(data))


class TestBaseProcessorFixes(MigrationTestBase):
    """Test fixes for BaseProcessor weaknesses."""
    
    def setUp(self):
        """Set up test environment."""
        super().setUp()
        self.processor = FixedMockProcessor(paths=self.paths)
    
    def test_sanitize_filename_with_complex_identifiers(self):
        """Test that complex identifiers are sanitized correctly."""
        # UniProt-style identifier
        complex_name = "sp|P12345|TEST_HUMAN"
        
        # Save entity
        self.processor.save_entity(complex_name, {"data": "test"})
        
        # Check that entity is registered with original name
        self.assertTrue(self.processor.entity_exists(complex_name))
        
        # Check that file was created with sanitized name
        safe_name = self.processor._sanitize_filename(complex_name)
        self.assertEqual(safe_name, "sp_P12345_TEST_HUMAN")
        
        file_path = self.processor.data_path / f"{safe_name}.test"
        self.assertTrue(file_path.exists())
        
        # Load should work with original name
        loaded = self.processor.load_entity(complex_name)
        self.assertIsNotNone(loaded)
        self.assertEqual(loaded["name"], complex_name)
    
    def test_sanitize_filename_with_windows_problematic_chars(self):
        """Test sanitization of Windows-problematic characters."""
        test_cases = {
            "file:name": "file_name",
            "file<>name": "file__name",
            "file*name": "file_name",
            "file?name": "file_name",
            "file|name": "file_name",
            "file\\name": "file_name",
            "file/name": "file_name",
            "file name": "file_name",
            '"filename"': "_filename_"
        }
        
        for original, expected in test_cases.items():
            sanitized = self.processor._sanitize_filename(original)
            self.assertEqual(sanitized, expected, 
                           f"Failed to sanitize '{original}' correctly")
    
    def test_delete_entity_with_sanitized_filename(self):
        """Test that delete_entity handles sanitized filenames correctly."""
        # Create entity with complex name
        complex_name = "test:entity|2025"
        self.processor.save_entity(complex_name, {"test": True})
        
        # Verify it exists
        self.assertTrue(self.processor.entity_exists(complex_name))
        
        # Delete should work with original name
        result = self.processor.delete_entity(complex_name)
        self.assertTrue(result)
        
        # Should no longer exist
        self.assertFalse(self.processor.entity_exists(complex_name))
        
        # File should be deleted
        safe_name = self.processor._sanitize_filename(complex_name)
        file_path = self.processor.data_path / f"{safe_name}.test"
        self.assertFalse(file_path.exists())
    
    def test_error_handling_uses_logger(self):
        """Test that errors are logged properly, not just printed."""
        # Create an entity that exists in registry but file is missing
        self.processor.save_entity("EXISTS", {"test": True})
        
        # Create dataset
        self.processor.create_dataset("test_dataset", ["EXISTS"])
        
        # Delete the physical file but keep registry entry
        file_path = self.processor.data_path / "EXISTS.test"
        file_path.unlink()
        
        # Capture logs
        with self.assertLogs(self.processor.logger, level=logging.WARNING) as cm:
            # This should log a warning when file not found
            result = self.processor.load_dataset("test_dataset")
            
        # Check that warning was logged
        self.assertTrue(any("EXISTS" in log for log in cm.output))
        
        # Result should be empty (failed loads excluded)
        self.assertEqual(result, {})
    
    def test_dataset_with_renamed_and_deleted_entities(self):
        """Test dataset handling with entity lifecycle changes."""
        # Create entities
        self.processor.save_entity("ORIGINAL", {"v": 1})
        self.processor.save_entity("KEEPER", {"v": 2})
        self.processor.save_entity("TO_DELETE", {"v": 3})
        
        # Create dataset
        self.processor.create_dataset("lifecycle_test", 
                                     ["ORIGINAL", "KEEPER", "TO_DELETE"])
        
        # Rename one entity
        self.processor.entity_registry.rename_entity("ORIGINAL", "RENAMED")
        
        # Delete another
        self.processor.delete_entity("TO_DELETE")
        
        # Load dataset should handle changes gracefully
        loaded = self.processor.load_dataset("lifecycle_test")
        
        # Should have renamed entity and keeper, but not deleted
        self.assertIn("RENAMED", loaded)
        self.assertIn("KEEPER", loaded)
        self.assertNotIn("ORIGINAL", loaded)
        self.assertNotIn("TO_DELETE", loaded)
        
        self.assertEqual(len(loaded), 2)


if __name__ == '__main__':
    unittest.main()