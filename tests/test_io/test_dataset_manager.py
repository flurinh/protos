"""
Tests for DatasetManager implementation.

Tests focus on:
- Dataset creation and management
- Human-readable names only
- ProtosPaths integration
- JSON file storage in datasets/ directories
"""

import pytest
import json
from pathlib import Path
import tempfile

from protos.io.dataset_manager import DatasetManager
from protos.io.entity_registry import EntityRegistry
from protos.io.paths import ProtosPaths


class TestDatasetManager:
    """Test DatasetManager functionality."""
    
    def test_init_with_protospaths(self):
        """Test DatasetManager accepts ProtosPaths in __init__."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            
            manager = DatasetManager("structure", paths=paths, entity_registry=registry)
            
            # Should create datasets directory
            datasets_dir = Path(tmpdir) / "structure" / "datasets"
            assert datasets_dir.exists()
    
    def test_create_dataset_with_human_names(self):
        """Test creating dataset with human-readable entity names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("structure", paths=paths, entity_registry=registry)
            
            # Register some entities
            registry.register_entity("1ubq", "structure", "structure/mmcif/1ubq.cif")
            registry.register_entity("3sn6", "structure", "structure/mmcif/3sn6.cif")
            
            # Create dataset with human names
            result = manager.create_dataset(
                "kinases",
                ["1ubq", "3sn6"],
                {"description": "Kinase structures"}
            )
            
            # Should return dataset name
            assert result == "kinases"
            
            # Check JSON file created
            dataset_path = Path(tmpdir) / "structure" / "datasets" / "kinases.json"
            assert dataset_path.exists()
            
            # Load and verify content
            with open(dataset_path) as f:
                data = json.load(f)
            
            assert data["name"] == "kinases"
            assert data["entities"] == ["1ubq", "3sn6"]
            assert data["metadata"]["description"] == "Kinase structures"
    
    def test_list_datasets_returns_names(self):
        """Test listing datasets returns human-readable names."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("sequence", paths=paths, entity_registry=registry)
            
            # Create multiple datasets
            manager.create_dataset("human_proteins", ["P62988", "P00533"])
            manager.create_dataset("viral_sequences", ["NC_045512", "NC_001802"])
            manager.create_dataset("test_set", ["TEST1", "TEST2"])
            
            # List should return names
            datasets = manager.list_datasets()
            assert "human_proteins" in datasets
            assert "viral_sequences" in datasets
            assert "test_set" in datasets
            assert len(datasets) == 3
    
    def test_load_dataset(self):
        """Test loading dataset returns correct information."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("grn", paths=paths, entity_registry=registry)
            
            # Create dataset
            manager.create_dataset(
                "opsin_study",
                ["BACR", "ChR2", "NpHR"],
                {"organism": "Various", "study": "Optogenetics"}
            )
            
            # Load dataset
            dataset = manager.load_dataset("opsin_study")
            
            assert dataset["name"] == "opsin_study"
            assert dataset["entities"] == ["BACR", "ChR2", "NpHR"]
            assert dataset["metadata"]["organism"] == "Various"
            assert dataset["processor_type"] == "grn"
    
    def test_add_to_dataset(self):
        """Test adding entities to existing dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("structure", paths=paths, entity_registry=registry)
            
            # Create initial dataset
            manager.create_dataset("my_structures", ["1ubq", "2gb1"])
            
            # Add more entities
            manager.add_to_dataset("my_structures", ["3sn6", "4mqa"])
            
            # Verify additions
            dataset = manager.load_dataset("my_structures")
            assert "3sn6" in dataset["entities"]
            assert "4mqa" in dataset["entities"]
            assert len(dataset["entities"]) == 4
            
            # Test no duplicates
            manager.add_to_dataset("my_structures", ["1ubq", "5new"])
            dataset = manager.load_dataset("my_structures")
            assert dataset["entities"].count("1ubq") == 1
            assert "5new" in dataset["entities"]
    
    def test_remove_from_dataset(self):
        """Test removing entities from dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("property", paths=paths, entity_registry=registry)
            
            # Create dataset
            manager.create_dataset("test_props", ["A", "B", "C", "D", "E"])
            
            # Remove some entities
            manager.remove_from_dataset("test_props", ["B", "D"])
            
            # Verify removals
            dataset = manager.load_dataset("test_props")
            assert "B" not in dataset["entities"]
            assert "D" not in dataset["entities"]
            assert dataset["entities"] == ["A", "C", "E"]
    
    def test_get_dataset_info(self):
        """Test getting detailed dataset information."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("embedding", paths=paths, entity_registry=registry)
            
            # Register entities
            registry.register_entity("seq1", "sequence", "sequence/fasta/seq1.fasta")
            registry.register_entity("seq1", "embedding", "embedding/embeddings/seq1.pkl")
            registry.register_entity("seq2", "sequence", "sequence/fasta/seq2.fasta")
            
            # Create dataset
            manager.create_dataset("embed_test", ["seq1", "seq2"])
            
            # Get info
            info = manager.get_dataset_info("embed_test")
            
            assert info["name"] == "embed_test"
            assert info["entity_count"] == 2
            assert len(info["entities"]) == 2
            
            # Check entity details
            seq1_info = next(e for e in info["entities"] if e["name"] == "seq1")
            assert "sequence" in seq1_info["formats"]
            assert "embedding" in seq1_info["formats"]
            
            seq2_info = next(e for e in info["entities"] if e["name"] == "seq2")
            assert "sequence" in seq2_info["formats"]
    
    def test_dataset_exists(self):
        """Test checking dataset existence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("ligand", paths=paths, entity_registry=registry)
            
            # Initially doesn't exist
            assert not manager.dataset_exists("my_ligands")
            
            # Create dataset
            manager.create_dataset("my_ligands", ["ATP", "GTP"])
            
            # Now exists
            assert manager.dataset_exists("my_ligands")
    
    def test_delete_dataset(self):
        """Test deleting a dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("structure", paths=paths, entity_registry=registry)
            
            # Create and verify
            manager.create_dataset("temp_dataset", ["A", "B"])
            assert manager.dataset_exists("temp_dataset")
            
            # Delete
            manager.delete_dataset("temp_dataset")
            
            # Should not exist
            assert not manager.dataset_exists("temp_dataset")
            
            # File should be gone
            dataset_path = Path(tmpdir) / "structure" / "datasets" / "temp_dataset.json"
            assert not dataset_path.exists()
    
    def test_update_metadata(self):
        """Test updating dataset metadata."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("sequence", paths=paths, entity_registry=registry)
            
            # Create dataset with initial metadata
            manager.create_dataset(
                "test_seqs",
                ["seq1", "seq2"],
                {"version": "1.0", "author": "Alice"}
            )
            
            # Update metadata
            manager.update_metadata("test_seqs", {
                "version": "1.1",
                "reviewer": "Bob",
                "validated": True
            })
            
            # Check updates
            dataset = manager.load_dataset("test_seqs")
            assert dataset["metadata"]["version"] == "1.1"
            assert dataset["metadata"]["author"] == "Alice"  # Original preserved
            assert dataset["metadata"]["reviewer"] == "Bob"
            assert dataset["metadata"]["validated"] is True
    
    def test_copy_dataset(self):
        """Test copying a dataset."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("grn", paths=paths, entity_registry=registry)
            
            # Create original
            manager.create_dataset(
                "original",
                ["A", "B", "C"],
                {"type": "test", "version": "1.0"}
            )
            
            # Copy
            manager.copy_dataset("original", "copy")
            
            # Verify copy
            copy = manager.load_dataset("copy")
            assert copy["entities"] == ["A", "B", "C"]
            assert copy["metadata"]["type"] == "test"
            assert copy["name"] == "copy"  # Name updated
    
    def test_merge_datasets(self):
        """Test merging multiple datasets."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("property", paths=paths, entity_registry=registry)
            
            # Create datasets to merge
            manager.create_dataset("set1", ["A", "B", "C"], {"source": "exp1"})
            manager.create_dataset("set2", ["C", "D", "E"], {"source": "exp2"})
            manager.create_dataset("set3", ["E", "F", "G"], {"source": "exp3"})
            
            # Merge
            manager.merge_datasets(["set1", "set2", "set3"], "merged")
            
            # Verify merge
            merged = manager.load_dataset("merged")
            
            # Should have all unique entities
            assert set(merged["entities"]) == {"A", "B", "C", "D", "E", "F", "G"}
            
            # Should track merge info
            assert merged["metadata"]["merged_from"] == ["set1", "set2", "set3"]
            assert "merge_date" in merged["metadata"]
    
    def test_warning_on_missing_entities(self):
        """Test warning when creating dataset with unregistered entities."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("structure", paths=paths, entity_registry=registry)
            
            # Register only some entities
            registry.register_entity("1ubq", "structure", "structure/mmcif/1ubq.cif")
            
            # Create dataset with mix of registered and unregistered
            # Should not fail, just warn
            result = manager.create_dataset(
                "mixed",
                ["1ubq", "unknown1", "unknown2"]
            )
            
            assert result == "mixed"
            
            # Dataset should still be created
            dataset = manager.load_dataset("mixed")
            assert dataset["entities"] == ["1ubq", "unknown1", "unknown2"]
    
    def test_json_file_format(self):
        """Test that dataset files have correct JSON structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            registry = EntityRegistry(paths=paths)
            manager = DatasetManager("sequence", paths=paths, entity_registry=registry)
            
            # Create dataset
            manager.create_dataset(
                "test_format",
                ["entity1", "entity2"],
                {"key": "value"}
            )
            
            # Read raw JSON
            dataset_path = Path(tmpdir) / "sequence" / "datasets" / "test_format.json"
            with open(dataset_path) as f:
                data = json.load(f)
            
            # Check structure
            assert "name" in data
            assert "processor_type" in data
            assert "entities" in data
            assert "metadata" in data
            assert "created" in data
            assert "modified" in data
            
            # Check formatting (should be indented)
            with open(dataset_path) as f:
                content = f.read()
            assert "  " in content  # Has indentation