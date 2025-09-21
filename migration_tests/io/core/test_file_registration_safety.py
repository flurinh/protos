"""
Test safe file registration features.

Tests:
1. Content hash detection for duplicates
2. Name conflict handling
3. Validation before registration
4. Input folder workflow
"""

import hashlib
from pathlib import Path

from ...test_base import MigrationTestBase


class TestFileRegistrationSafety(MigrationTestBase):
    """Test safe file registration features."""
    
    def compute_file_hash(self, file_path: Path) -> str:
        """Compute SHA256 hash of file content."""
        hasher = hashlib.sha256()
        with open(file_path, 'rb') as f:
            while chunk := f.read(8192):
                hasher.update(chunk)
        return hasher.hexdigest()
    
    def test_duplicate_content_detection(self):
        """Test that duplicate content is detected even with different names."""
        # Skip if not implemented yet
        if not hasattr(self.entity_registry, 'find_by_content_hash'):
            self.skipTest("Content hash detection not implemented yet")
        
        # Create two files with identical content
        content = "ATOM 1 CA ALA A 1 10.0 20.0 30.0"
        file1 = self.create_test_structure("protein1", content)
        file2 = self.create_test_structure("protein2", content)
        
        # Register first file
        self.entity_registry.register_entity(
            name="protein1",
            format_type="structure",
            file_path=str(file1.relative_to(self.paths.data_root))
        )
        
        # Attempt to register second file with same content
        # Should detect duplicate content
        result = self.entity_registry.register_file_safely(
            file_path=file2,
            entity_name="protein2",
            format_type="structure"
        )
        
        # Should not create new entity
        self.assertEqual(self.count_registered_entities(), 1)
        
        # Could either skip or add as alias
        # Implementation will determine behavior
    
    def test_name_conflict_handling(self):
        """Test handling when entity name already exists."""
        entity_name = "conflict_test"
        
        # Register first file
        file1 = self.create_test_structure(entity_name, "Content 1")
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(file1.relative_to(self.paths.data_root))
        )
        
        # Try to register different file with same name
        file2 = self.create_test_structure("temp_conflict", "Content 2")
        
        # Skip if safe registration not implemented
        if not hasattr(self.entity_registry, 'register_file_safely'):
            self.skipTest("Safe registration not implemented yet")
        
        # Should handle conflict (version, skip, or error)
        result = self.entity_registry.register_file_safely(
            file_path=file2,
            entity_name=entity_name,
            format_type="structure",
            conflict_strategy="version"  # Create versioned name
        )
        
        # Original should still exist
        self.assertTrue(self.entity_registry.entity_exists(entity_name))
    
    def test_structure_dual_format_handling(self):
        """Test that StructureProcessor correctly handles CIF and PKL formats."""
        entity_name = "1UBQ"
        
        # Register CIF file
        cif_path = self.create_test_file("structure/mmcif", f"{entity_name}.cif", 
                                        "# CIF content")
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure",
            file_path=str(cif_path.relative_to(self.paths.data_root)),
            metadata={"format": "cif"}
        )
        
        # Register PKL file for same entity
        pkl_path = self.create_test_file("structure/cache", f"{entity_name}.pkl",
                                        "# PKL content")
        self.entity_registry.register_entity(
            name=entity_name,
            format_type="structure", 
            file_path=str(pkl_path.relative_to(self.paths.data_root)),
            metadata={"format": "pkl", "priority": 1}
        )
        
        # Should have one entity with structure format
        self.assertEqual(self.count_registered_entities(), 1)
        
        # Should track both files under same entity
        entity_info = self.entity_registry.find_entity(entity_name, "structure")
        self.assertIsNotNone(entity_info)
        
        # Future: Check that PKL has higher priority than CIF
    
    def test_table_based_entity_registration(self):
        """Test registration of table-based entities (GRN/Property)."""
        table_name = "protein_properties"
        
        # Create property table with multiple entities
        table_content = """entity_name,lambda_max,activity
PROTEIN1,550,high
PROTEIN2,480,medium
PROTEIN3,620,low"""
        
        table_path = self.create_test_table("property", table_name, table_content)
        
        # Register table (each row is an entity)
        # This tests a different pattern than file-based
        entities_in_table = ["PROTEIN1", "PROTEIN2", "PROTEIN3"]
        
        for entity_name in entities_in_table:
            self.entity_registry.register_entity(
                name=entity_name,
                format_type="property",
                file_path=str(table_path.relative_to(self.paths.data_root)),
                metadata={"table": table_name, "row": entity_name}
            )
        
        # Should have three entities
        self.assertEqual(self.count_registered_entities(), 3)
        
        # Each should be findable
        for entity_name in entities_in_table:
            entity_info = self.entity_registry.find_entity(entity_name, "property")
            self.assertIsNotNone(entity_info)
            self.assertEqual(entity_info.metadata["table"], table_name)
    
    def test_input_folder_workflow(self):
        """Test the input folder registration workflow."""
        # Create input directories
        input_dir = Path(self.paths.data_root) / "input"
        processed_dir = input_dir / "processed"
        rejected_dir = input_dir / "rejected"
        
        for dir in [input_dir, processed_dir, rejected_dir]:
            dir.mkdir(parents=True, exist_ok=True)
        
        # Place files in input folder
        input_files = {
            "new_protein.cif": "# Valid CIF content",
            "sequences.fasta": ">SEQ1\nMKLLVTAA\n>SEQ2\nMFLVLLA",
            "bad_file.xyz": "Invalid format content"
        }
        
        for filename, content in input_files.items():
            (input_dir / filename).write_text(content)
        
        # Skip if InputManager not implemented
        try:
            from protos.io.input_manager import InputManager
            manager = InputManager(self.paths)
        except ImportError:
            self.skipTest("InputManager not implemented yet")
        
        # Process input files
        report = manager.process_input_files(dry_run=False)
        
        # Check results
        # Valid files should be in processed/
        # Invalid files should be in rejected/
        # Registry should have new entities
    
    def test_registry_health_check(self):
        """Test finding unregistered files."""
        # Create some files without registering them
        unregistered_files = [
            self.create_test_structure("unregistered1"),
            self.create_test_structure("unregistered2"),
            self.create_test_sequence("unregistered_seq")
        ]
        
        # Register one file properly
        registered_file = self.create_test_structure("registered")
        self.entity_registry.register_entity(
            name="registered",
            format_type="structure",
            file_path=str(registered_file.relative_to(self.paths.data_root))
        )
        
        # Skip if health checker not implemented
        try:
            from protos.io.registry_health import RegistryHealthChecker
            checker = RegistryHealthChecker(self.entity_registry, self.paths)
        except ImportError:
            self.skipTest("RegistryHealthChecker not implemented yet")
        
        # Find unregistered structure files
        unregistered = checker.find_unregistered_files("structure")
        
        # Should find the two unregistered structure files
        self.assertEqual(len(unregistered), 2)
        
        # Check names
        unregistered_names = {f.stem for f in unregistered}
        self.assertEqual(unregistered_names, {"unregistered1", "unregistered2"})


if __name__ == '__main__':
    unittest.main()