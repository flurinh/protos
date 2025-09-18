"""
Base test class for entity system migration tests.

This provides the foundation for all migration tests, ensuring:
1. ProtosPaths handles all data structure setup
2. Temporary test directories are used
3. Proper cleanup after tests
"""

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from typing import Optional

from protos.io.paths import ProtosPaths
from protos.io.entity_registry import EntityRegistry
from protos.io.dataset_manager import DatasetManager


class MigrationTestBase(unittest.TestCase):
    """Base class for migration tests with ProtosPaths setup."""
    
    def setUp(self):
        """Set up test environment with ProtosPaths managing all paths."""
        # Create temporary directory for test data
        self.temp_dir = tempfile.mkdtemp(prefix="protos_migration_test_")
        self.test_data_root = Path(self.temp_dir) / "data"
        self.test_data_root.mkdir(parents=True, exist_ok=True)
        
        # Set environment variable for ProtosPaths (following conftest.py pattern)
        os.environ["PROTOS_DATA_ROOT"] = str(self.test_data_root.absolute())
        
        # Initialize ProtosPaths - let it use the environment variable
        self.paths = ProtosPaths()
        
        # Ensure key directories are created by calling getters
        for processor_type in self.paths.processor_dirs:
            self.paths.get_processor_path(processor_type)
        
        # Ensure global registry path exists
        self.paths.get_global_registry_path()
        
        # Initialize core components
        self.entity_registry = EntityRegistry(paths=self.paths)
        
        # Ensure we're in a clean state
        self.assertEqual(len(self.entity_registry._registry), 0, 
                        "Registry should be empty at start")
    
    def tearDown(self):
        """Clean up test environment."""
        # Clean up environment variable
        if "PROTOS_DATA_ROOT" in os.environ:
            del os.environ["PROTOS_DATA_ROOT"]
        
        # Remove temporary directory
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)
    
    def create_test_file(self, subdir: str, filename: str, content: str = "test") -> Path:
        """
        Create a test file in the appropriate directory.
        
        Args:
            subdir: Subdirectory under processor path (e.g., "mmcif", "fasta")
            filename: Name of the file to create
            content: Content to write to the file
            
        Returns:
            Path to the created file
        """
        # Let ProtosPaths determine the correct path
        if "structure" in subdir or "mmcif" in subdir:
            base_dir = Path(self.paths.get_processor_path("structure"))
        elif "sequence" in subdir or "fasta" in subdir:
            base_dir = Path(self.paths.get_processor_path("sequence"))
        elif "grn" in subdir:
            base_dir = Path(self.paths.get_processor_path("grn"))
        elif "property" in subdir:
            base_dir = Path(self.paths.get_processor_path("property"))
        else:
            base_dir = Path(self.paths.get_processor_path("generic"))
        
        # Create subdirectory if specified
        if "/" in subdir:
            file_dir = base_dir / subdir
        else:
            file_dir = base_dir / subdir if subdir else base_dir
            
        file_dir.mkdir(parents=True, exist_ok=True)
        
        # Create file
        file_path = file_dir / filename
        file_path.write_text(content)
        
        return file_path
    
    def create_test_structure(self, pdb_id: str, content: Optional[str] = None) -> Path:
        """Create a test structure file."""
        if content is None:
            content = f"# Test structure for {pdb_id}\nATOM 1 CA ALA A 1 0.0 0.0 0.0"
        return self.create_test_file("mmcif", f"{pdb_id}.cif", content)
    
    def create_test_sequence(self, name: str, content: Optional[str] = None) -> Path:
        """Create a test sequence file."""
        if content is None:
            content = f">{name}\nMKFLVLALLGLLVFSVATVQA"
        return self.create_test_file("fasta", f"{name}.fasta", content)
    
    def create_test_table(self, processor_type: str, table_name: str, 
                         content: Optional[str] = None) -> Path:
        """Create a test table file for GRN or Property processor."""
        if content is None:
            if processor_type == "grn":
                content = "protein_id,1.50,2.50,3.50\nTEST1,A,L,V\nTEST2,A,I,V"
            else:  # property
                content = "entity_name,property1,property2\nTEST1,10.5,active\nTEST2,7.2,inactive"
        
        subdir = "tables"
        return self.create_test_file(f"{processor_type}/{subdir}", 
                                    f"{table_name}.csv", content)
    
    def assertEntityRegistered(self, entity_name: str, format_type: str):
        """Assert that an entity is registered with the given format."""
        entity_info = self.entity_registry.find_entity(entity_name, format_type)
        self.assertIsNotNone(entity_info, 
                           f"Entity '{entity_name}' should be registered for format '{format_type}'")
        return entity_info
    
    def assertEntityNotRegistered(self, entity_name: str):
        """Assert that an entity is not registered."""
        entity_info = self.entity_registry.find_entity(entity_name)
        self.assertIsNone(entity_info, 
                         f"Entity '{entity_name}' should not be registered")
    
    def assertOnlyHumanNamesReturned(self, result_list: list):
        """Assert that a list contains only human-readable names (no UUIDs/hashes)."""
        for item in result_list:
            # Check that item doesn't look like a hash or UUID
            self.assertFalse(
                len(item) == 10 and all(c in '0123456789abcdef' for c in item),
                f"Hash ID exposed: {item}"
            )
            self.assertFalse(
                len(item) == 36 and item.count('-') == 4,
                f"UUID exposed: {item}"
            )
    
    def get_registry_data(self) -> dict:
        """Get the raw registry data for inspection."""
        return self.entity_registry._registry
    
    def count_registered_entities(self) -> int:
        """Count the number of registered entities."""
        return len(self.entity_registry._registry)


if __name__ == '__main__':
    unittest.main()