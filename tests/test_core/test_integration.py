"""
Integration tests for core components: ProtosPaths, EntityRegistry, DatasetManager, BaseProcessor.

These tests verify that all components work together correctly before
we proceed to update specific processors.
"""

import pytest
import tempfile
from pathlib import Path
from typing import Any, Optional

from protos.io.paths import ProtosPaths
from protos.io.entity_registry import EntityRegistry
from protos.io.dataset_manager import DatasetManager
from protos.core.base_processor import BaseProcessor


class TestStructureProcessor(BaseProcessor):
    """Test implementation of a structure processor."""
    
    def load_entity(self, name: str) -> Any:
        """Load structure entity."""
        # Check mmcif directory
        cif_path = self.get_subdirectory_path('structure_dir') / f"{name}.cif"
        if cif_path.exists():
            # Auto-register if not already registered
            if not self.entity_exists(name):
                self.entity_registry.register_entity(
                    name=name,
                    format_type=self.processor_type,
                    file_path=str(cif_path.relative_to(self.paths.data_root)),
                    metadata={"auto_discovered": True}
                )
            return f"Structure data for {name}"
        
        # Check if registered
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info:
            file_path = Path(entity_info.file_path)
            if not file_path.is_absolute():
                file_path = Path(self.paths.data_root) / file_path
            if file_path.exists():
                return f"Structure data for {name}"
        
        return None
    
    def save_entity(self, name: str, data: Any, metadata: Optional[dict] = None):
        """Save structure entity."""
        # Save to mmcif directory
        mmcif_dir = self.get_subdirectory_path('structure_dir')
        mmcif_dir.mkdir(parents=True, exist_ok=True)
        
        file_path = mmcif_dir / f"{name}.cif"
        file_path.write_text(str(data))
        
        # Register entity
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=str(file_path.relative_to(self.paths.data_root)),
            metadata=metadata or {}
        )


class TestSequenceProcessor(BaseProcessor):
    """Test implementation of a sequence processor."""
    
    def _get_processor_type(self):
        return 'sequence'
    
    def load_entity(self, name: str) -> Any:
        """Load sequence entity."""
        # Check fasta directory
        fasta_path = self.get_subdirectory_path('fasta_dir') / f"{name}.fasta"
        if fasta_path.exists():
            if not self.entity_exists(name):
                self.entity_registry.register_entity(
                    name=name,
                    format_type=self.processor_type,
                    file_path=str(fasta_path.relative_to(self.paths.data_root)),
                    metadata={"auto_discovered": True}
                )
            return fasta_path.read_text()
        
        # Check registry
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info:
            file_path = Path(entity_info.file_path)
            if not file_path.is_absolute():
                file_path = Path(self.paths.data_root) / file_path
            if file_path.exists():
                return file_path.read_text()
        
        return None
    
    def save_entity(self, name: str, data: Any, metadata: Optional[dict] = None):
        """Save sequence entity."""
        # Save to fasta directory
        fasta_dir = self.get_subdirectory_path('fasta_dir')
        fasta_dir.mkdir(parents=True, exist_ok=True)
        
        file_path = fasta_dir / f"{name}.fasta"
        file_path.write_text(str(data))
        
        # Register entity
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=str(file_path.relative_to(self.paths.data_root)),
            metadata=metadata or {}
        )


class TestCoreIntegration:
    """Test all core components working together."""
    
    def test_shared_protospaths(self):
        """Test that all components can share the same ProtosPaths instance."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create single ProtosPaths instance
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            
            # Create components sharing the same paths
            entity_registry = EntityRegistry(paths=paths)
            
            struct_proc = TestStructureProcessor(name="struct", paths=paths)
            seq_proc = TestSequenceProcessor(name="seq", paths=paths)
            
            # Verify all use the same paths
            assert struct_proc.paths == paths
            assert seq_proc.paths == paths
            assert struct_proc.entity_registry.paths == paths
            assert seq_proc.entity_registry.paths == paths
    
    def test_cross_format_entity_tracking(self):
        """Test that entities can be tracked across multiple formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Create processors
            struct_proc = TestStructureProcessor(paths=paths)
            seq_proc = TestSequenceProcessor(paths=paths)
            
            # Save same entity in both formats
            struct_proc.save_entity("ubiquitin", "ATOM 1 CA MET A 1 ...", 
                                   {"resolution": 1.8})
            seq_proc.save_entity("ubiquitin", ">ubiquitin\nMQIFVKTLTG...",
                                {"length": 76})
            
            # Check entity exists in both formats
            assert struct_proc.entity_exists("ubiquitin")
            assert seq_proc.entity_exists("ubiquitin")
            
            # Refresh registry to see updates from other processor
            struct_proc.entity_registry.refresh()
            
            # Check shared registry sees both formats
            formats = struct_proc.entity_registry.get_entity_formats("ubiquitin")
            assert "structure" in formats
            assert "sequence" in formats
    
    def test_dataset_with_real_entities(self):
        """Test dataset creation with entities that actually exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            struct_proc = TestStructureProcessor(paths=paths)
            
            # Save some real entities
            struct_proc.save_entity("1ubq", "Structure of ubiquitin")
            struct_proc.save_entity("3sn6", "GPCR structure")
            struct_proc.save_entity("7zvl", "Another structure")
            
            # Create dataset
            struct_proc.create_dataset(
                "test_structures",
                ["1ubq", "3sn6", "7zvl"],
                {"description": "Test dataset"}
            )
            
            # Load dataset
            loaded = struct_proc.load_dataset("test_structures")
            assert len(loaded) == 3
            assert "1ubq" in loaded
            assert loaded["1ubq"] == "Structure data for 1ubq"
    
    def test_drag_and_drop_workflow(self):
        """Test that files can be dropped into directories and immediately used."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            struct_proc = TestStructureProcessor(paths=paths)
            
            # Simulate user dropping a file
            mmcif_dir = Path(tmpdir) / "structure" / "mmcif"
            mmcif_dir.mkdir(parents=True, exist_ok=True)
            
            dropped_file = mmcif_dir / "6xyz.cif"
            dropped_file.write_text("HEADER DROPPED STRUCTURE")
            
            # Should be immediately loadable
            data = struct_proc.load_entity("6xyz")
            assert data == "Structure data for 6xyz"
            
            # Should now be registered
            assert struct_proc.entity_exists("6xyz")
    
    def test_no_hardcoded_paths(self):
        """Verify that no hardcoded paths are used anywhere."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Set custom data root
            custom_root = Path(tmpdir) / "custom" / "location"
            paths = ProtosPaths(data_root=str(custom_root))
            
            struct_proc = TestStructureProcessor(paths=paths)
            
            # Save entity
            struct_proc.save_entity("test", "data")
            
            # Verify file is in custom location
            expected_file = custom_root / "structure" / "mmcif" / "test.cif"
            assert expected_file.exists()
            
            # Verify registry is in custom location
            registry_file = custom_root / "global_registry.json"
            assert registry_file.exists()
    
    def test_human_names_only_in_apis(self):
        """Verify that all public APIs use human-readable names only."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            struct_proc = TestStructureProcessor(paths=paths)
            
            # Save entities with various names
            names = ["1ubq", "EGFR_HUMAN", "my-protein_v2.1", "sp|P02724|GLPA_ECOLI"]
            for name in names:
                struct_proc.save_entity(name, f"Data for {name}")
            
            # List entities - should return human names
            entities = struct_proc.list_entities()
            for name in names:
                assert name in entities
            
            # Create dataset - accepts human names
            struct_proc.create_dataset("test_set", names[:2])
            
            # Get dataset info - returns human names
            info = struct_proc.get_dataset_info("test_set")
            for entity_info in info['entities']:
                assert entity_info['name'] in names
                # Should never expose hash IDs
                assert len(entity_info['name']) != 10 or not entity_info['name'].isalnum()
    
    def test_dataset_manager_integration(self):
        """Test DatasetManager works correctly with EntityRegistry."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            entity_registry = EntityRegistry(paths=paths)
            dataset_manager = DatasetManager(
                processor_type="structure",
                paths=paths,
                entity_registry=entity_registry
            )
            
            # Register some entities
            entity_registry.register_entity("1ubq", "structure", "structure/mmcif/1ubq.cif")
            entity_registry.register_entity("3sn6", "structure", "structure/mmcif/3sn6.cif")
            
            # Create dataset
            dataset_manager.create_dataset(
                "test_dataset",
                ["1ubq", "3sn6"],
                {"study": "integration test"}
            )
            
            # Check dataset file exists
            dataset_file = Path(tmpdir) / "structure" / "datasets" / "test_dataset.json"
            assert dataset_file.exists()
            
            # Load and verify
            dataset_info = dataset_manager.get_dataset_info("test_dataset")
            assert dataset_info['entity_count'] == 2
            assert len(dataset_info['entities']) == 2
    
    def test_abstract_methods_enforcement(self):
        """Test that processors must implement abstract methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Try to create processor without implementing abstract methods
            class BadProcessor(BaseProcessor):
                pass
            
            # Should fail to instantiate
            with pytest.raises(TypeError, match="abstract"):
                BadProcessor(paths=paths)
    
    def test_processor_type_detection(self):
        """Test automatic processor type detection from class name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Test various processor names
            struct_proc = TestStructureProcessor(paths=paths)
            assert struct_proc.processor_type == "structure"
            
            seq_proc = TestSequenceProcessor(paths=paths)
            assert seq_proc.processor_type == "sequence"
            
            # Test that data paths are correct
            assert "structure" in str(struct_proc.data_path)
            assert "sequence" in str(seq_proc.data_path)