"""
Tests for multi-format entity support across processors.
"""

import pytest
import os
import shutil
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.sequence.seq_processor import SeqProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.io.data_access import GlobalRegistry, generate_entity_id


@pytest.fixture
def setup_test_environment(tmp_path):
    """Set up test environment with real data."""
    # Set up paths
    ProtosPaths.set_data_root(str(tmp_path))
    
    # Copy test data
    test_data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/tests/test-data")
    
    # Copy structure data
    struct_dir = tmp_path / "structure" / "mmcif"
    struct_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy a test structure
    source_cif = test_data_dir / "structure" / "mmcif" / "1ubq.cif"
    if source_cif.exists():
        shutil.copy(source_cif, struct_dir / "1ubq.cif")
    
    # Copy sequence data
    seq_dir = tmp_path / "sequence" / "fasta" 
    seq_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy GRN data
    grn_dir = tmp_path / "grn" / "ref"
    grn_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy a test GRN table
    source_grn = test_data_dir / "grn" / "ref" / "mo_ref.csv"
    if source_grn.exists():
        shutil.copy(source_grn, grn_dir / "mo_ref.csv")
    
    return tmp_path


class TestMultiFormatEntity:
    """Test that the same biological entity can exist in multiple formats."""
    
    def test_structure_to_sequence_same_entity(self, setup_test_environment):
        """Test that a structure and its extracted sequence share the same entity ID."""
        # Create processors
        struct_processor = CifBaseProcessor(name="struct_test")
        seq_processor = SeqProcessor(name="seq_test")
        
        # Load a real structure
        structure = struct_processor.load_structure("1ubq")
        assert structure is not None
        
        # Extract sequence from structure
        seq_dict = struct_processor.get_seq_dict()
        assert "1ubq" in seq_dict
        
        # The entity ID should be based on "1ubq"
        expected_entity_id = generate_entity_id("1ubq")
        
        # Check structure entity
        global_registry = GlobalRegistry()
        entity_info = global_registry.entity_registry.get_entity(expected_entity_id)
        assert entity_info is not None
        assert entity_info['original_id'] == "1ubq"
        assert 'structure' in entity_info['formats']
        
        # Save the sequence with the same entity ID
        # In real usage, the sequence processor would register the same entity
        # when it processes sequence data for "1ubq"
        sequence_data = seq_dict["1ubq"]["A"]  # Chain A sequence
        
        # Register sequence format for the same entity
        global_registry.entity_registry.register_entity(
            entity_id=expected_entity_id,
            entity_type="sequence",
            original_id="1ubq",
            file_path=str(setup_test_environment / "sequence" / "fasta" / "1ubq.fasta"),
            metadata={"length": len(sequence_data), "chain": "A"}
        )
        
        # Verify the entity now has both formats
        updated_entity = global_registry.entity_registry.get_entity(expected_entity_id)
        assert 'structure' in updated_entity['formats']
        assert 'sequence' in updated_entity['formats']
        
        # Both formats should reference the same original ID
        assert updated_entity['original_id'] == "1ubq"
    
    def test_resolve_identifier_finds_multi_format_entity(self, setup_test_environment):
        """Test that resolve_identifier works correctly for multi-format entities."""
        global_registry = GlobalRegistry()
        
        # Create an entity with multiple formats
        entity_id = generate_entity_id("P12345")
        
        # Register as sequence
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id="P12345",
            file_path="/test/P12345.fasta",
            metadata={"length": 350}
        )
        
        # Register as structure
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id="P12345",
            file_path="/test/P12345.cif",
            metadata={"chains": ["A", "B"]}
        )
        
        # Register as GRN
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="grn",
            original_id="P12345",
            file_path=None,  # GRN entries are in tables
            metadata={"protein_family": "test_family"}
        )
        
        # resolve_identifier should find the entity regardless of format
        resolved_seq = global_registry.entity_registry.resolve_identifier("P12345", format_type="sequence")
        resolved_struct = global_registry.entity_registry.resolve_identifier("P12345", format_type="structure")
        resolved_grn = global_registry.entity_registry.resolve_identifier("P12345", format_type="grn")
        
        # All should resolve to the same entity ID
        assert resolved_seq == entity_id
        assert resolved_struct == entity_id
        assert resolved_grn == entity_id
        
        # Without format type, should still resolve
        resolved_any = global_registry.entity_registry.resolve_identifier("P12345")
        assert resolved_any == entity_id
    
    def test_list_entities_by_format(self, setup_test_environment):
        """Test listing entities filtered by format type."""
        global_registry = GlobalRegistry()
        
        # Create entities with different format combinations
        # Entity 1: structure only
        entity1 = generate_entity_id("1ABC")
        global_registry.entity_registry.register_entity(
            entity_id=entity1,
            entity_type="structure",
            original_id="1ABC",
            file_path="/test/1ABC.cif"
        )
        
        # Entity 2: sequence only  
        entity2 = generate_entity_id("P54321")
        global_registry.entity_registry.register_entity(
            entity_id=entity2,
            entity_type="sequence",
            original_id="P54321",
            file_path="/test/P54321.fasta"
        )
        
        # Entity 3: both structure and sequence
        entity3 = generate_entity_id("2DEF")
        global_registry.entity_registry.register_entity(
            entity_id=entity3,
            entity_type="structure",
            original_id="2DEF",
            file_path="/test/2DEF.cif"
        )
        global_registry.entity_registry.register_entity(
            entity_id=entity3,
            entity_type="sequence",
            original_id="2DEF",
            file_path="/test/2DEF.fasta"
        )
        
        # List by format
        structures = global_registry.entity_registry.list_entities(format_type="structure")
        sequences = global_registry.entity_registry.list_entities(format_type="sequence")
        
        # Check results
        assert entity1 in structures
        assert entity3 in structures
        assert entity2 not in structures
        
        assert entity2 in sequences
        assert entity3 in sequences
        assert entity1 not in sequences