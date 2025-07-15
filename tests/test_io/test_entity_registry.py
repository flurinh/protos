"""
Tests for EntityRegistry implementation.

Tests focus on:
- Human-readable names in all public APIs
- Hash IDs used internally only
- ProtosPaths integration
- Cross-format entity tracking
"""

import pytest
import json
from pathlib import Path
import tempfile

from protos.io.entity_registry import EntityRegistry, EntityInfo
from protos.io.paths import ProtosPaths


class TestEntityRegistry:
    """Test EntityRegistry functionality."""
    
    def test_init_with_protospaths(self):
        """Test EntityRegistry accepts ProtosPaths in __init__."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            registry = EntityRegistry(paths=paths)
            
            # Should use ProtosPaths for registry file
            expected_path = Path(tmpdir) / "global_registry.json"
            assert registry.registry_file == expected_path
    
    def test_register_entity_returns_human_name(self):
        """Test register_entity returns human-readable name, not hash."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register entity with human name
            result = registry.register_entity(
                name="1ubq",
                format_type="structure",
                file_path="structure/mmcif/1ubq.cif",
                metadata={"resolution": 1.8}
            )
            
            # Should return human name, not hash
            assert result == "1ubq"
            assert len(result) == 4  # Not a 10-char hash
    
    def test_list_entities_returns_human_names(self):
        """Test list_entities returns human-readable names only."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register multiple entities
            registry.register_entity("1ubq", "structure", "structure/mmcif/1ubq.cif")
            registry.register_entity("EGFR_HUMAN", "sequence", "sequence/fasta/EGFR_HUMAN.fasta")
            registry.register_entity("my_protein", "structure", "structure/mmcif/my_protein.cif")
            
            # List should contain human names only
            entities = registry.list_entities()
            assert "1ubq" in entities
            assert "EGFR_HUMAN" in entities
            assert "my_protein" in entities
            
            # Should not contain hash IDs (alphanumeric 10-char strings)
            for entity in entities:
                # Check if it's actually a hash ID format (10 alphanumeric chars)
                if len(entity) == 10 and entity.isalnum():
                    # This would be a hash ID - should not appear in public API
                    assert False, f"Hash ID {entity} exposed in public API"
    
    def test_find_entity_by_human_name(self):
        """Test finding entities by human-readable name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register entity
            registry.register_entity(
                "1ubq",
                "structure",
                "structure/mmcif/1ubq.cif",
                {"resolution": 1.8, "method": "X-ray"}
            )
            
            # Find by human name
            info = registry.find_entity("1ubq")
            assert info is not None
            assert info.original_id == "1ubq"
            assert info.format_type == "structure"
            assert info.file_path == "structure/mmcif/1ubq.cif"
            assert info.metadata["resolution"] == 1.8
    
    def test_file_paths_use_human_names(self):
        """Test that file paths always use human-readable names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register with human-readable path
            registry.register_entity(
                "sp_P02724_GLPA_ECOLI",
                "sequence",
                "sequence/fasta/sp_P02724_GLPA_ECOLI.fasta"
            )
            
            # Retrieve and check path
            info = registry.find_entity("sp_P02724_GLPA_ECOLI")
            assert "sp_P02724_GLPA_ECOLI" in info.file_path
            assert info.file_path.endswith(".fasta")
    
    def test_aliases_work_correctly(self):
        """Test entity aliases functionality."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register entity
            registry.register_entity("P62988", "sequence", "sequence/fasta/P62988.fasta")
            
            # Add aliases
            registry.add_alias("P62988", "UBIQ_HUMAN")
            registry.add_alias("P62988", "ubiquitin")
            
            # Should find by any alias
            assert registry.find_entity("P62988") is not None
            assert registry.find_entity("UBIQ_HUMAN") is not None
            assert registry.find_entity("ubiquitin") is not None
            
            # All should resolve to same entity
            info1 = registry.find_entity("P62988")
            info2 = registry.find_entity("UBIQ_HUMAN")
            info3 = registry.find_entity("ubiquitin")
            
            assert info1.hash_id == info2.hash_id == info3.hash_id
            assert info1.original_id == "P62988"
    
    def test_multi_format_entity(self):
        """Test entity with multiple formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register same entity in multiple formats
            registry.register_entity("1ubq", "structure", "structure/mmcif/1ubq.cif")
            registry.register_entity("1ubq", "sequence", "sequence/fasta/1ubq.fasta")
            
            # Check formats
            formats = registry.get_entity_formats("1ubq")
            assert "structure" in formats
            assert "sequence" in formats
            
            # Find specific format
            struct_info = registry.find_entity("1ubq", "structure")
            assert struct_info.file_path == "structure/mmcif/1ubq.cif"
            
            seq_info = registry.find_entity("1ubq", "sequence")
            assert seq_info.file_path == "sequence/fasta/1ubq.fasta"
    
    def test_entity_exists(self):
        """Test entity existence checking."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register entity
            registry.register_entity("test_protein", "structure", "structure/mmcif/test_protein.cif")
            
            # Check existence
            assert registry.entity_exists("test_protein")
            assert registry.entity_exists("test_protein", "structure")
            assert not registry.entity_exists("test_protein", "sequence")
            assert not registry.entity_exists("unknown_protein")
    
    def test_metadata_operations(self):
        """Test metadata get/update operations."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register with metadata
            registry.register_entity(
                "3sn6",
                "structure",
                "structure/mmcif/3sn6.cif",
                {"resolution": 3.2, "organism": "Homo sapiens"}
            )
            
            # Get metadata
            metadata = registry.get_entity_metadata("3sn6", "structure")
            assert metadata["resolution"] == 3.2
            assert metadata["organism"] == "Homo sapiens"
            
            # Update metadata
            registry.update_metadata("3sn6", "structure", {"chains": ["A", "B", "G"]})
            
            # Check updated
            metadata = registry.get_entity_metadata("3sn6", "structure")
            assert metadata["chains"] == ["A", "B", "G"]
            assert metadata["resolution"] == 3.2  # Original still there
    
    def test_remove_format(self):
        """Test removing a format from an entity."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register multi-format entity
            registry.register_entity("test", "structure", "structure/mmcif/test.cif")
            registry.register_entity("test", "sequence", "sequence/fasta/test.fasta")
            
            # Remove one format
            registry.remove_format("test", "sequence")
            
            # Should still exist in other format
            assert registry.entity_exists("test")
            assert registry.entity_exists("test", "structure")
            assert not registry.entity_exists("test", "sequence")
            
            # Remove last format - entity should be gone
            registry.remove_format("test", "structure")
            assert not registry.entity_exists("test")
    
    def test_persistence(self):
        """Test registry persistence across instances."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            
            # Create and populate registry
            registry1 = EntityRegistry(paths=paths)
            registry1.register_entity("persistent", "structure", "structure/mmcif/persistent.cif")
            registry1.add_alias("persistent", "persist_alias")
            
            # Create new instance - should load existing data
            registry2 = EntityRegistry(paths=paths)
            
            # Check data persisted
            assert registry2.entity_exists("persistent")
            assert registry2.entity_exists("persist_alias")
            
            info = registry2.find_entity("persistent")
            assert info.original_id == "persistent"
            assert info.file_path == "structure/mmcif/persistent.cif"
    
    def test_hash_ids_not_exposed(self):
        """Test that hash IDs are never exposed in public API."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Register entities
            registry.register_entity("test1", "structure", "test1.cif")
            registry.register_entity("test2", "structure", "test2.cif")
            
            # Check registry file to ensure hash IDs are used internally
            with open(registry.registry_file) as f:
                data = json.load(f)
            
            # Internal storage should use hash IDs as keys
            entity_keys = list(data['entities'].keys())
            assert all(len(key) == 10 for key in entity_keys)
            
            # But public API should never expose them
            entities = registry.list_entities()
            assert "test1" in entities
            assert "test2" in entities
            assert all(len(name) != 10 or name in ["test1", "test2"] for name in entities)
    
    def test_complex_identifiers(self):
        """Test handling of complex biological identifiers."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            # Complex identifiers
            complex_ids = [
                "sp|P02724|GLPA_ECOLI",
                "tr|Q9Y6K9|NEMO_HUMAN",
                "pdb|1ABC|A",
                "my-protein_v2.1"
            ]
            
            for complex_id in complex_ids:
                # Should handle without issues
                result = registry.register_entity(
                    complex_id,
                    "sequence",
                    f"sequence/fasta/{complex_id.replace('|', '_')}.fasta"
                )
                assert result == complex_id
                
                # Should find by exact name
                info = registry.find_entity(complex_id)
                assert info is not None
                assert info.original_id == complex_id