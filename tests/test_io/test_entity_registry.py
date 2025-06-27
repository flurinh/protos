"""
Tests for the entity registry system.
"""

import os
import json
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from protos.io.data_access import EntityRegistry, generate_entity_id
from protos.core.base_processor import BaseProcessor
from protos.io.paths.path_config import ProtosPaths


class TestEntityIDGeneration:
    """Test entity ID generation functionality."""
    
    def test_generate_entity_id_basic(self):
        """Test basic entity ID generation."""
        # Generate ID from simple content
        id1 = generate_entity_id("1ABC")
        assert len(id1) == 10
        assert id1.isalnum()
        
        # Same content should give same ID
        id2 = generate_entity_id("1ABC")
        assert id1 == id2
        
        # Different content should give different ID
        id3 = generate_entity_id("2DEF")
        assert id1 != id3
    
    def test_generate_entity_id_with_prefix(self):
        """Test entity ID generation with prefixes."""
        # Different prefixes should give different IDs
        id1 = generate_entity_id("1ABC", prefix="structure")
        id2 = generate_entity_id("1ABC", prefix="sequence")
        assert id1 != id2
        
        # Same prefix should give same ID
        id3 = generate_entity_id("1ABC", prefix="structure")
        assert id1 == id3
    
    def test_entity_id_consistency(self):
        """Test that entity IDs are consistent across runs."""
        # Generate multiple IDs
        content_list = ["protein1", "sequence_ABC", "structure_123"]
        ids = {}
        
        # Generate IDs twice
        for content in content_list:
            ids[content] = [
                generate_entity_id(content),
                generate_entity_id(content)
            ]
        
        # Check consistency
        for content, id_list in ids.items():
            assert id_list[0] == id_list[1], f"IDs for {content} are not consistent"


class TestEntityRegistry:
    """Test EntityRegistry functionality."""
    
    @pytest.fixture
    def registry(self, tmp_path):
        """Create a temporary entity registry."""
        registry_file = tmp_path / "entity_registry.json"
        return EntityRegistry(str(registry_file))
    
    def test_register_entity(self, registry):
        """Test registering an entity."""
        entity_id = generate_entity_id("1ABC", prefix="structure")
        
        registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id="1ABC",
            file_path="/data/structures/1ABC.cif",
            metadata={"resolution": 2.5, "method": "X-RAY"},
            datasets=["test_structures"]
        )
        
        # Check entity was registered
        assert registry.entity_exists(entity_id)
        entity_info = registry.get_entity(entity_id)
        assert entity_info is not None
        assert "structure" in entity_info["formats"]
        assert entity_info["original_id"] == "1ABC"
    
    def test_find_entity_by_original_id(self, registry):
        """Test finding entities by original ID."""
        # Register entities
        entity1_id = generate_entity_id("1ABC", prefix="structure")
        registry.register_entity(
            entity_id=entity1_id,
            entity_type="structure",
            original_id="1ABC"
        )
        
        entity2_id = generate_entity_id("P12345", prefix="sequence")
        registry.register_entity(
            entity_id=entity2_id,
            entity_type="sequence",
            original_id="P12345"
        )
        
        # Find by original ID
        found_id = registry.find_entity_by_original_id("1ABC")
        assert found_id == entity1_id
        
        # Find with format type hint
        found_id = registry.find_entity_by_original_id("P12345", format_type="sequence")
        assert found_id == entity2_id
        
        # Not found
        found_id = registry.find_entity_by_original_id("NOTFOUND")
        assert found_id is None
    
    def test_list_entities(self, registry):
        """Test listing entities with filters."""
        # Register multiple entities
        entities = [
            ("struct1", "structure", "1ABC", ["dataset1"]),
            ("struct2", "structure", "2DEF", ["dataset1", "dataset2"]),
            ("seq1", "sequence", "P12345", ["dataset2"]),
            ("grn1", "grn", "table1", ["dataset3"])
        ]
        
        for entity_id, entity_type, original_id, datasets in entities:
            registry.register_entity(
                entity_id=entity_id,
                entity_type=entity_type,
                original_id=original_id,
                datasets=datasets
            )
        
        # List all entities
        all_entities = registry.list_entities()
        assert len(all_entities) == 4
        
        # List by format type
        structure_entities = registry.list_entities(format_type="structure")
        assert len(structure_entities) == 2
        assert "struct1" in structure_entities
        assert "struct2" in structure_entities
        
        # List by dataset
        dataset1_entities = registry.list_entities(dataset="dataset1")
        assert len(dataset1_entities) == 2
        assert "struct1" in dataset1_entities
        assert "struct2" in dataset1_entities
        
        dataset2_entities = registry.list_entities(dataset="dataset2")
        assert len(dataset2_entities) == 2
        assert "struct2" in dataset2_entities
        assert "seq1" in dataset2_entities
    
    def test_dataset_operations(self, registry):
        """Test dataset-related operations."""
        # Register entity
        entity_id = "test_entity"
        registry.register_entity(
            entity_id=entity_id,
            entity_type="test",
            datasets=["dataset1"]
        )
        
        # Add to dataset
        success = registry.add_entity_to_dataset(entity_id, "dataset2")
        assert success
        
        # Check datasets
        entity_info = registry.get_entity(entity_id)
        assert "dataset1" in entity_info["datasets"]
        assert "dataset2" in entity_info["datasets"]
        
        # Get dataset entities
        dataset1_entities = registry.get_dataset_entities("dataset1")
        assert entity_id in dataset1_entities
        
        dataset2_entities = registry.get_dataset_entities("dataset2")
        assert entity_id in dataset2_entities
        
        # Remove from dataset
        success = registry.remove_entity_from_dataset(entity_id, "dataset1")
        assert success
        
        entity_info = registry.get_entity(entity_id)
        assert "dataset1" not in entity_info["datasets"]
        assert "dataset2" in entity_info["datasets"]
    
    def test_entity_metadata(self, registry):
        """Test entity metadata operations."""
        entity_id = "test_entity"
        initial_metadata = {"key1": "value1", "key2": 2}
        
        registry.register_entity(
            entity_id=entity_id,
            entity_type="test",
            metadata=initial_metadata
        )
        
        # Update metadata
        new_metadata = {"key2": 3, "key3": "value3"}
        success = registry.update_entity_metadata(entity_id, new_metadata, format_type="test")
        assert success
        
        # Check merged metadata
        entity_info = registry.get_entity(entity_id)
        metadata = entity_info["formats"]["test"]["metadata"]
        assert metadata["key1"] == "value1"  # Original preserved
        assert metadata["key2"] == 3  # Updated
        assert metadata["key3"] == "value3"  # New added
    
    def test_entity_stats(self, registry):
        """Test entity statistics."""
        # Register entities
        registry.register_entity("e1", "structure", datasets=["d1"])
        registry.register_entity("e2", "structure", datasets=["d1", "d2"])
        registry.register_entity("e3", "sequence", datasets=["d2"])
        registry.register_entity("e4", "grn")  # No dataset (orphaned)
        
        stats = registry.get_entity_stats()
        
        assert stats["total_entities"] == 4
        assert stats["by_format"]["structure"] == 2
        assert stats["by_format"]["sequence"] == 1
        assert stats["by_format"]["grn"] == 1
        assert stats["by_dataset"]["d1"] == 2
        assert stats["by_dataset"]["d2"] == 2
        assert stats["orphaned"] == 1
    
    def test_registry_persistence(self, tmp_path):
        """Test that registry persists to disk."""
        registry_file = tmp_path / "entity_registry.json"
        
        # Create and populate registry
        registry1 = EntityRegistry(str(registry_file))
        entity_id = generate_entity_id("test_data")
        registry1.register_entity(
            entity_id=entity_id,
            entity_type="test",
            original_id="test123",
            metadata={"test": True}
        )
        
        # Create new registry instance from same file
        registry2 = EntityRegistry(str(registry_file))
        
        # Check data persisted
        assert registry2.entity_exists(entity_id)
        entity_info = registry2.get_entity(entity_id)
        assert entity_info["original_id"] == "test123"
        assert entity_info["formats"]["test"]["metadata"]["test"] is True
        
        # Check file format
        with open(registry_file) as f:
            data = json.load(f)
        
        assert "entities" in data
        assert "datasets" in data
        assert "version" in data
        assert entity_id in data["entities"]


class TestBaseProcessorEntityMethods:
    """Test BaseProcessor entity management methods."""
    
    @pytest.fixture
    def processor(self, tmp_path):
        """Create a test processor."""
        ProtosPaths.set_data_root(str(tmp_path))
        processor = BaseProcessor(
            name="test_processor",
            processor_data_dir="test"
        )
        yield processor
        ProtosPaths.set_data_root(None)
    
    def test_save_and_load_entity(self, processor):
        """Test saving and loading entities through BaseProcessor."""
        # Create test data
        test_data = pd.DataFrame({
            'id': ['A', 'B', 'C'],
            'value': [1, 2, 3]
        })
        
        # Save as entity
        entity_id = processor.save_entity(
            data=test_data,
            original_id="test_data_1",
            metadata={"description": "Test dataset"},
            datasets=["test_dataset"]
        )
        
        # Verify entity ID format
        assert len(entity_id) == 10
        assert entity_id.isalnum()
        
        # Check entity exists
        assert processor.entity_exists(entity_id)
        
        # Find by original ID
        found_id = processor.find_entity_by_original_id("test_data_1")
        assert found_id == entity_id
        
        # Get metadata
        metadata = processor.get_entity_metadata(entity_id)
        assert metadata["description"] == "Test dataset"
    
    def test_list_entities(self, processor):
        """Test listing entities through BaseProcessor."""
        # Save multiple entities
        entity_ids = []
        for i in range(3):
            data = pd.DataFrame({'value': [i]})
            entity_id = processor.save_entity(
                data=data,
                original_id=f"test_{i}",
                datasets=["dataset1"] if i < 2 else ["dataset2"]
            )
            entity_ids.append(entity_id)
        
        # List all entities
        all_entities = processor.list_entities()
        for entity_id in entity_ids:
            assert entity_id in all_entities
        
        # List by dataset
        dataset1_entities = processor.list_entities(dataset="dataset1")
        assert len(dataset1_entities) == 2
        assert entity_ids[0] in dataset1_entities
        assert entity_ids[1] in dataset1_entities
        assert entity_ids[2] not in dataset1_entities
    
    def test_entity_dataset_operations(self, processor):
        """Test entity-dataset operations through BaseProcessor."""
        # Save entity
        data = pd.DataFrame({'value': [1]})
        entity_id = processor.save_entity(
            data=data,
            original_id="test_entity",
            datasets=["dataset1"]
        )
        
        # Add to another dataset
        success = processor.add_entity_to_dataset(entity_id, "dataset2")
        assert success
        
        # Get dataset entities
        dataset1_entities = processor.get_dataset_entities("dataset1")
        assert entity_id in dataset1_entities
        
        dataset2_entities = processor.get_dataset_entities("dataset2")
        assert entity_id in dataset2_entities
    
    def test_entity_id_generation_for_different_types(self, processor):
        """Test entity ID generation for different data types."""
        # DataFrame
        df_data = pd.DataFrame({'a': [1, 2, 3]})
        df_id = processor.save_entity(df_data, original_id="df_test")
        
        # Dictionary
        dict_data = {'key': 'value', 'number': 42}
        dict_id = processor.save_entity(dict_data, original_id="dict_test")
        
        # NumPy array
        np_data = np.array([1, 2, 3, 4, 5])
        np_id = processor.save_entity(np_data, original_id="np_test")
        
        # All should have valid entity IDs
        for entity_id in [df_id, dict_id, np_id]:
            assert len(entity_id) == 10
            assert entity_id.isalnum()
            assert processor.entity_exists(entity_id)
        
        # All should be different
        assert df_id != dict_id
        assert dict_id != np_id
        assert df_id != np_id


class TestMultiFormatEntities:
    """Test multi-format entity support."""
    
    @pytest.fixture
    def registry(self, tmp_path):
        """Create a temporary entity registry."""
        registry_file = tmp_path / "entity_registry.json"
        return EntityRegistry(str(registry_file))
    
    def test_multi_format_entity(self, registry):
        """Test registering same entity in multiple formats."""
        # Use same original ID for all formats
        original_id = "P12345"
        entity_id = generate_entity_id(original_id)
        
        # Register as sequence
        registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=original_id,
            file_path="/data/sequences/P12345.fasta",
            metadata={"length": 350, "organism": "Human"},
            datasets=["human_proteins"]
        )
        
        # Register same entity as structure (after AlphaFold prediction)
        registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id=original_id,
            file_path="/data/structures/AF-P12345-F1.cif",
            metadata={"source": "alphafold", "confidence": 0.95},
            datasets=["human_proteins", "alphafold_structures"]
        )
        
        # Register same entity as GRN
        registry.register_entity(
            entity_id=entity_id,
            entity_type="grn",
            original_id=original_id,
            file_path="/data/grn/human_gpcr.csv",
            metadata={"grn_system": "gpcrdb", "positions": 127},
            datasets=["human_proteins", "gpcr_grns"]
        )
        
        # Register same entity as embedding
        registry.register_entity(
            entity_id=entity_id,
            entity_type="embedding",
            original_id=original_id,
            file_path="/data/embeddings/P12345_esm2.pkl",
            metadata={"model": "esm2_t33_650M", "dim": 1280},
            datasets=["human_proteins", "protein_embeddings"]
        )
        
        # Check entity has all formats
        formats = registry.get_entity_formats(entity_id)
        assert len(formats) == 4
        assert "sequence" in formats
        assert "structure" in formats
        assert "grn" in formats
        assert "embedding" in formats
        
        # Check entity info
        entity_info = registry.get_entity(entity_id)
        assert entity_info["original_id"] == original_id
        assert len(entity_info["datasets"]) == 4  # Union of all datasets
        
        # Get format-specific info
        seq_info = registry.get_entity_by_format(entity_id, "sequence")
        assert seq_info["metadata"]["length"] == 350
        
        struct_info = registry.get_entity_by_format(entity_id, "structure")
        assert struct_info["metadata"]["source"] == "alphafold"
        
        grn_info = registry.get_entity_by_format(entity_id, "grn")
        assert grn_info["metadata"]["grn_system"] == "gpcrdb"
        
        emb_info = registry.get_entity_by_format(entity_id, "embedding")
        assert emb_info["metadata"]["model"] == "esm2_t33_650M"
        
        # Find by original ID works for any format
        found_id = registry.find_entity_by_original_id(original_id)
        assert found_id == entity_id
        
        # Filter by format type
        found_id_seq = registry.find_entity_by_original_id(original_id, format_type="sequence")
        assert found_id_seq == entity_id
        
        found_id_struct = registry.find_entity_by_original_id(original_id, format_type="structure")
        assert found_id_struct == entity_id
    
    def test_register_entity_format_method(self, registry):
        """Test register_entity_format method."""
        # First register as sequence
        entity_id = generate_entity_id("Q67890")
        registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id="Q67890",
            file_path="/data/sequences/Q67890.fasta"
        )
        
        # Add structure format using register_entity_format
        success = registry.register_entity_format(
            entity_id=entity_id,
            format_type="structure",
            file_path="/data/structures/Q67890.cif",
            metadata={"resolution": 2.5}
        )
        assert success
        
        # Check both formats exist
        formats = registry.get_entity_formats(entity_id)
        assert len(formats) == 2
        assert "sequence" in formats
        assert "structure" in formats
        
        # Try to add format to non-existent entity
        success = registry.register_entity_format(
            entity_id="nonexistent",
            format_type="structure",
            file_path="/data/structures/fake.cif"
        )
        assert not success
    
    def test_entity_stats_multi_format(self, registry):
        """Test entity stats with multi-format entities."""
        # Register some single-format entities
        registry.register_entity("e1", "sequence", original_id="P111")
        registry.register_entity("e2", "structure", original_id="1ABC")
        
        # Register multi-format entity
        entity_id = generate_entity_id("P999")
        registry.register_entity(entity_id, "sequence", original_id="P999")
        registry.register_entity(entity_id, "structure", original_id="P999")
        registry.register_entity(entity_id, "grn", original_id="P999")
        
        stats = registry.get_entity_stats()
        
        assert stats["total_entities"] == 3  # e1, e2, and P999
        assert stats["multi_format_entities"] == 1  # Only P999
        assert stats["by_format"]["sequence"] == 2  # e1 and P999
        assert stats["by_format"]["structure"] == 2  # e2 and P999
        assert stats["by_format"]["grn"] == 1  # Only P999


class TestGRNEntityFormat:
    """Test how GRN tables work with the entity system."""
    
    def test_grn_table_with_entity_ids(self):
        """Test GRN table format with entity IDs."""
        # Create example GRN table
        grn_data = {
            'sequence_id': ['BR1_HUMAN', 'BR2_MOUSE', 'BR3_BOVIN'],
            '1.50': ['L45', 'L46', 'L47'],
            '2.50': ['V87', 'V88', 'V89'],
            '3.50': ['I123', 'I124', 'I125']
        }
        grn_df = pd.DataFrame(grn_data)
        
        # Generate entity IDs for each row (based on sequence ID only!)
        entity_ids = []
        for seq_id in grn_df['sequence_id']:
            # Universal entity ID based on sequence ID only
            entity_id = generate_entity_id(seq_id)
            entity_ids.append(entity_id)
        
        # Add entity_id column
        grn_df['entity_id'] = entity_ids
        
        # Reorder columns to put entity_id first
        cols = ['entity_id'] + [col for col in grn_df.columns if col != 'entity_id']
        grn_df = grn_df[cols]
        
        # Verify format
        assert grn_df.columns[0] == 'entity_id'
        assert len(grn_df['entity_id'].unique()) == 3
        
        # Each entity ID should be unique (in this case)
        assert grn_df['entity_id'].nunique() == len(grn_df)
        
        # Entity IDs should have correct format
        for entity_id in grn_df['entity_id']:
            assert len(entity_id) == 10
            assert entity_id.isalnum()
        
        # Test that same sequence ID always gives same entity ID
        test_id = generate_entity_id('BR1_HUMAN')
        assert test_id == entity_ids[0]