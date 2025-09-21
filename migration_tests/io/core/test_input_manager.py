"""
Test safe file registration with InputManager.

Tests:
1. Input folder setup and README creation
2. File type detection
3. Content hash duplicate detection  
4. Conflict resolution strategies
5. File validation
6. Processing workflow
"""

import unittest
import tempfile
import shutil
from pathlib import Path
from datetime import datetime

from ...test_base import MigrationTestBase
from protos.io.core.input_manager import InputManager, ConflictResolutionStrategy, InputFile


class TestInputManager(MigrationTestBase):
    """Test InputManager functionality."""
    
    def setUp(self):
        """Set up test environment."""
        super().setUp()
        self.manager = InputManager(paths=self.paths, entity_registry=self.entity_registry)
        
        # Create some test files
        self.test_files = {
            'test_structure.cif': self._create_test_cif(),
            'test_sequence.fasta': self._create_test_fasta(),
            'test_properties.csv': self._create_test_csv()
        }
    
    def _create_test_cif(self):
        """Create a minimal valid CIF file."""
        content = """data_TEST
_entry.id TEST
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA ALA A 1 1.000 2.000 3.000
"""
        return content
    
    def _create_test_fasta(self):
        """Create a valid FASTA file."""
        return """>TEST_PROTEIN Test protein sequence
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
"""
    
    def _create_test_csv(self):
        """Create a valid property table."""
        return """protein,property1,property2
TEST1,1.23,active
TEST2,4.56,inactive
"""
    
    def test_input_folder_setup(self):
        """Test that input folders are created with README."""
        # Folders should exist
        self.assertTrue(self.manager.input_dir.exists())
        self.assertTrue(self.manager.processed_dir.exists())
        self.assertTrue(self.manager.rejected_dir.exists())
        
        # README should exist
        readme = self.manager.input_dir / "README.txt"
        self.assertTrue(readme.exists())
        self.assertIn("Place your data files here", readme.read_text())
    
    def test_file_type_detection(self):
        """Test file type detection by extension and content."""
        # Test by extension
        self.assertEqual(self.manager._detect_file_type(Path("test.cif")), "structure")
        self.assertEqual(self.manager._detect_file_type(Path("test.fasta")), "sequence")
        self.assertEqual(self.manager._detect_file_type(Path("test.csv")), "property")
        self.assertEqual(self.manager._detect_file_type(Path("test.npy")), "embedding")
        
        # Test by content (create actual files)
        test_file = self.manager.input_dir / "no_extension"
        
        # Test FASTA content
        test_file.write_text(">TEST\nACGT")
        self.assertEqual(self.manager._detect_by_content(test_file), "sequence")
        
        # Test CIF content
        test_file.write_text("data_TEST\n_atom_site.id")
        self.assertEqual(self.manager._detect_by_content(test_file), "structure")
        
        # Clean up
        test_file.unlink()
    
    def test_scan_input_folder(self):
        """Test scanning input folder for files."""
        # Add test files to input folder
        for filename, content in self.test_files.items():
            file_path = self.manager.input_dir / filename
            file_path.write_text(content)
        
        # Scan folder
        input_files = self.manager.scan_input_folder()
        
        # Check results
        self.assertEqual(len(input_files), 3)
        
        # Check each file
        filenames = [f.path.name for f in input_files]
        self.assertIn('test_structure.cif', filenames)
        self.assertIn('test_sequence.fasta', filenames)
        self.assertIn('test_properties.csv', filenames)
        
        # Check file types detected correctly
        for f in input_files:
            if 'structure' in f.path.name:
                self.assertEqual(f.file_type, 'structure')
                self.assertEqual(f.entity_name, 'test_structure')
            elif 'sequence' in f.path.name:
                self.assertEqual(f.file_type, 'sequence')
                self.assertEqual(f.entity_name, 'test_sequence')
            elif 'properties' in f.path.name:
                self.assertEqual(f.file_type, 'property')
                self.assertEqual(f.entity_name, 'test_properties')
    
    def test_content_hash_duplicate_detection(self):
        """Test detection of duplicate files by content hash."""
        # Create a file and register it
        file1 = self.manager.input_dir / "protein1.cif"
        file1.write_text(self.test_files['test_structure.cif'])
        
        # Process it
        report = self.manager.process_input_files()
        self.assertEqual(len(report.processed), 1)
        
        # Create identical file with different name
        file2 = self.manager.input_dir / "protein2.cif"
        file2.write_text(self.test_files['test_structure.cif'])
        
        # Process again - should be skipped as duplicate
        report = self.manager.process_input_files()
        self.assertEqual(len(report.skipped), 1)
        self.assertEqual(report.skipped[0].details['reason'], 'duplicate_content')
    
    def test_name_conflict_resolution_skip(self):
        """Test SKIP strategy for name conflicts."""
        # Create and process a file
        file1 = self.manager.input_dir / "myprotein.cif"
        file1.write_text(self.test_files['test_structure.cif'])
        report = self.manager.process_input_files()
        self.assertEqual(len(report.processed), 1)
        
        # Create different file with same base name
        file2 = self.manager.input_dir / "myprotein.cif"
        file2.write_text(self.test_files['test_structure.cif'] + "\n# Modified")
        
        # Process with SKIP strategy
        report = self.manager.process_input_files(
            conflict_strategy=ConflictResolutionStrategy.SKIP
        )
        self.assertEqual(len(report.skipped), 1)
        self.assertEqual(report.skipped[0].details['reason'], 'name_exists')
    
    def test_name_conflict_resolution_version(self):
        """Test VERSION strategy for name conflicts."""
        # Create and process a file
        file1 = self.manager.input_dir / "myprotein.cif"
        file1.write_text(self.test_files['test_structure.cif'])
        report = self.manager.process_input_files()
        self.assertEqual(len(report.processed), 1)
        
        # Create different file with same base name
        file2 = self.manager.input_dir / "myprotein.cif"
        modified_content = self.test_files['test_structure.cif'].replace("1.000", "2.000")
        file2.write_text(modified_content)
        
        # Process with VERSION strategy
        report = self.manager.process_input_files(
            conflict_strategy=ConflictResolutionStrategy.VERSION
        )
        self.assertEqual(len(report.processed), 1)
        # Should create versioned name
        self.assertEqual(report.processed[0].entity_name, 'myprotein_v2')
    
    def test_file_validation(self):
        """Test file validation before registration."""
        # Test empty file
        empty_file = self.manager.input_dir / "empty.cif"
        empty_file.touch()
        
        report = self.manager.process_input_files()
        self.assertEqual(len(report.rejected), 1)
        self.assertIn("empty", report.rejected[0].error_message.lower())
        
        # Test invalid CIF
        bad_cif = self.manager.input_dir / "bad.cif"
        bad_cif.write_text("This is not a valid CIF file")
        
        report = self.manager.process_input_files()
        self.assertEqual(len(report.rejected), 1)
        self.assertIn("Invalid CIF format", report.rejected[0].error_message)
        
        # Test invalid FASTA
        bad_fasta = self.manager.input_dir / "bad.fasta"
        bad_fasta.write_text("No header line\nACGT")
        
        report = self.manager.process_input_files()
        self.assertEqual(len(report.rejected), 1)
        self.assertIn("Invalid FASTA format", report.rejected[0].error_message)
    
    def test_dry_run_mode(self):
        """Test dry run doesn't make changes."""
        # Add test file
        test_file = self.manager.input_dir / "dryrun.cif"
        test_file.write_text(self.test_files['test_structure.cif'])
        
        # Process in dry run mode
        report = self.manager.process_input_files(dry_run=True)
        
        # Check report shows what would happen
        self.assertEqual(len(report.processed), 1)
        self.assertEqual(report.processed[0].action_taken, 'would_register')
        
        # File should still be in input folder
        self.assertTrue(test_file.exists())
        
        # No files in processed folder
        processed_files = list(self.manager.processed_dir.iterdir())
        self.assertEqual(len(processed_files), 0)
        
        # Entity should not be registered
        self.assertFalse(self.entity_registry.entity_exists('dryrun'))
    
    def test_rejected_files_with_error_logs(self):
        """Test rejected files are moved with error logs."""
        # Create invalid file
        bad_file = self.manager.input_dir / "invalid.cif"
        bad_file.write_text("Invalid content")
        
        # Process
        report = self.manager.process_input_files()
        
        # Check file was rejected
        self.assertEqual(len(report.rejected), 1)
        
        # Original file should be gone
        self.assertFalse(bad_file.exists())
        
        # File should be in rejected folder
        rejected_files = list(self.manager.rejected_dir.glob("*invalid.cif"))
        self.assertEqual(len(rejected_files), 1)
        
        # Error log should exist
        error_logs = list(self.manager.rejected_dir.glob("*_errors.txt"))
        self.assertTrue(len(error_logs) > 0)
        
        # Check error log content
        error_log = error_logs[0]
        content = error_log.read_text()
        self.assertIn("Invalid CIF format", content)
    
    def test_processed_files_archival(self):
        """Test successfully processed files are archived."""
        # Add valid file
        test_file = self.manager.input_dir / "success.cif"
        test_file.write_text(self.test_files['test_structure.cif'])
        
        # Process
        report = self.manager.process_input_files()
        self.assertEqual(len(report.processed), 1)
        
        # Original file should be gone
        self.assertFalse(test_file.exists())
        
        # File should be in processed folder with timestamp
        processed_files = list(self.manager.processed_dir.glob("*success.cif"))
        self.assertEqual(len(processed_files), 1)
        
        # Filename should have timestamp prefix
        processed_name = processed_files[0].name
        self.assertRegex(processed_name, r'^\d{8}_\d{6}_success\.cif$')
    
    def test_entity_registration_with_metadata(self):
        """Test entities are registered with proper metadata."""
        # Add file
        test_file = self.manager.input_dir / "metadata_test.cif"
        test_file.write_text(self.test_files['test_structure.cif'])
        
        # Process
        report = self.manager.process_input_files()
        
        # Check entity was registered
        self.assertTrue(self.entity_registry.entity_exists('metadata_test'))
        
        # Get entity info
        entity_info = self.entity_registry.find_entity('metadata_test', 'structure')
        self.assertIsNotNone(entity_info)
        
        # Check metadata
        self.assertEqual(entity_info.metadata['source'], 'user_input')
        self.assertEqual(entity_info.metadata['original_name'], 'metadata_test.cif')
        self.assertIn('content_hash', entity_info.metadata)
        self.assertIn('registered', entity_info.metadata)
        
        # Content hash should be set
        self.assertTrue(len(entity_info.metadata['content_hash']) > 0)


if __name__ == '__main__':
    unittest.main()