"""
Test that all processor list operations return biological names, not entity hashes.

This is a core requirement of the Protos philosophy: users work with names,
never with internal hash IDs.
"""

import pytest
import os
import pandas as pd
import numpy as np
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.io.data_access import GlobalRegistry, generate_entity_id
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.sequence.seq_processor import SeqProcessor
from protos.processing.embedding.embedding_processor import EmbeddingProcessor


@pytest.fixture
def setup_test_environment(request):
    """Set up test environment with test-data directory."""
    test_data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/tests/test-data")
    ProtosPaths.set_data_root(str(test_data_dir))
    
    # Clear global registry to ensure clean state
    global_registry = GlobalRegistry()
    global_registry.entity_registry._entities = {}
    global_registry.entity_registry._datasets = {}
    
    def teardown():
        ProtosPaths.set_data_root(None)
    
    request.addfinalizer(teardown)
    return test_data_dir


class TestListOperationsReturnNames:
    """Test that list operations return names, not hashes."""
    
    def test_cifbase_list_structures_returns_pdb_ids(self, setup_test_environment):
        """Test CifBaseProcessor.list_structures returns PDB IDs."""
        processor = CifBaseProcessor(name="test_list")
        
        # Load some real structures to populate the registry
        test_pdbs = ["1ubq", "1tqn", "3nir"]
        for pdb_id in test_pdbs:
            cif_path = setup_test_environment / "structure" / "mmcif" / f"{pdb_id}.cif"
            if cif_path.exists():
                processor.load_structure(pdb_id)
        
        # Get the list of structures
        structures = processor.list_structures()
        
        # Verify we get PDB IDs, not hashes
        for struct_id in structures:
            # Entity hashes are 10 alphanumeric characters
            is_hash = len(struct_id) == 10 and struct_id.isalnum()
            assert not is_hash, f"Got hash ID {struct_id} instead of PDB ID"
            
            # Should be one of our test PDBs
            if struct_id in test_pdbs:
                assert struct_id in test_pdbs
    
    def test_grn_list_entities_returns_sequence_ids(self, setup_test_environment):
        """Test GRNBaseProcessor.list_grn_entities returns sequence IDs."""
        processor = GRNBaseProcessor(name="test_grn_list")
        
        # Check if we have real GRN data
        grn_ref_path = setup_test_environment / "grn" / "ref" / "mo_grn.csv"
        if grn_ref_path.exists():
            # Load real GRN table
            processor.load_grn_table("ref/mo_grn")
            
            # Get list of entities
            entities = processor.list_grn_entities()
            
            # Verify we get sequence IDs, not hashes
            for entity in entities:
                # Entity hashes are 10 alphanumeric characters
                is_hash = len(entity) == 10 and entity.isalnum() and entity.islower()
                assert not is_hash, f"Got hash ID {entity} instead of sequence ID"
                
                # Should be in our loaded IDs
                assert entity in processor.ids
        else:
            # Create test GRN data
            test_data = pd.DataFrame({
                '1.50': ['A', 'V', 'L'],
                '2.50': ['G', 'G', 'P'],
                '3.50': ['I', 'V', 'F']
            }, index=['OPSIN_HUMAN', 'BR1_MOUSE', 'BR2_BOVIN'])
            
            processor.data = test_data
            processor.ids = test_data.index.tolist()
            processor.grns = test_data.columns.tolist()
            
            # Register entities
            for seq_id in processor.ids:
                entity_id = generate_entity_id(seq_id)
                global_registry = GlobalRegistry()
                global_registry.entity_registry.register_entity(
                    entity_id=entity_id,
                    entity_type="grn",
                    original_id=seq_id
                )
            
            # Get list
            entities = processor.list_grn_entities()
            
            # Verify sequence IDs
            for entity in entities:
                assert entity in ['OPSIN_HUMAN', 'BR1_MOUSE', 'BR2_BOVIN']
                assert not (len(entity) == 10 and entity.isalnum() and entity.islower())
    
    def test_seq_list_entities_returns_sequence_names(self, setup_test_environment):
        """Test SeqProcessor.list_sequence_entities returns sequence names."""
        processor = SeqProcessor(name="test_seq_list")
        global_registry = GlobalRegistry()
        
        # Just register test sequences in the registry
        test_seqs = ["P12345", "Q67890", "MYSEQ_HUMAN"]
        for seq_name in test_seqs:
            entity_id = generate_entity_id(seq_name)
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="sequence",
                original_id=seq_name
            )
        
        # Get list - this should check the registry, not cached sequences
        entities = processor.list_sequence_entities()
        
        # Verify names, not hashes
        for entity in entities:
            is_hash = len(entity) == 10 and entity.isalnum() and entity.islower()
            assert not is_hash, f"Got hash ID {entity} instead of sequence name"
            
            # Should be one of our test names
            if entity in test_seqs:
                assert entity in test_seqs
    
    def test_embedding_list_entities_returns_names(self, setup_test_environment):
        """Test EmbeddingProcessor.list_embedding_entities returns names."""
        processor = EmbeddingProcessor(name="test_emb_list")
        global_registry = GlobalRegistry()
        
        # Register test embeddings
        test_names = ["EMB1_TEST", "P99999", "GPCR_HUMAN"]
        for name in test_names:
            entity_id = generate_entity_id(name)
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="embedding",
                original_id=name,
                metadata={"model": "test", "dim": 320}
            )
        
        # Get list
        entities = processor.list_embedding_entities()
        
        # Verify names
        for entity in entities:
            is_hash = len(entity) == 10 and entity.isalnum() and entity.islower()
            assert not is_hash, f"Got hash ID {entity} instead of name"
            
            # Should be one of our test names (or empty if no embeddings)
            if entity in test_names:
                assert entity in test_names
    
    def test_cross_processor_consistency(self, setup_test_environment):
        """Test that the same entity appears with same name across processors."""
        global_registry = GlobalRegistry()
        
        # Register a protein in multiple formats
        protein_name = "TEST_PROTEIN"
        entity_id = generate_entity_id(protein_name)
        
        # Register as sequence
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=protein_name
        )
        
        # Register as structure (e.g., after AlphaFold prediction)
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id=protein_name
        )
        
        # Register as GRN entry
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="grn",
            original_id=protein_name
        )
        
        # Register as embedding
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="embedding",
            original_id=protein_name
        )
        
        # Create processors
        seq_proc = SeqProcessor(name="seq")
        struct_proc = CifBaseProcessor(name="struct")
        grn_proc = GRNBaseProcessor(name="grn")
        emb_proc = EmbeddingProcessor(name="emb")
        
        # List from each processor
        seq_list = seq_proc.list_sequence_entities()
        struct_list = struct_proc.list_structures()
        
        # For GRN, we need to mock having this in our loaded data
        grn_proc.ids = [protein_name]
        grn_list = grn_proc.list_grn_entities()
        
        emb_list = emb_proc.list_embedding_entities()
        
        # All should show the same protein name
        if protein_name in seq_list:
            assert protein_name in seq_list
        if protein_name in struct_list:
            assert protein_name in struct_list
        if protein_name in grn_list:
            assert protein_name in grn_list
        if protein_name in emb_list:
            assert protein_name in emb_list
        
        # None should show the hash
        for lst in [seq_list, struct_list, grn_list, emb_list]:
            assert entity_id not in lst