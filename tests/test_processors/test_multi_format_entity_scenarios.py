"""
Test multi-format entity scenarios.

This ensures that the same biological entity can be tracked across different
data formats (sequence, structure, GRN, embedding), which is core to the
Protos entity system design.
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


class TestMultiFormatEntityScenarios:
    """Test scenarios where entities exist in multiple formats."""
    
    def test_sequence_to_structure_workflow(self, setup_test_environment):
        """Test typical workflow: sequence → structure prediction."""
        global_registry = GlobalRegistry()
        
        # Start with a protein sequence
        protein_id = "P12345"
        entity_id = generate_entity_id(protein_id)
        
        # 1. Register as sequence
        seq_processor = SeqProcessor(name="seq_workflow")
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=protein_id,
            file_path=str(setup_test_environment / "sequence" / "fasta" / f"{protein_id}.fasta"),
            metadata={
                "length": 350,
                "organism": "Homo sapiens",
                "gene_name": "GPCR1"
            }
        )
        
        # 2. After AlphaFold prediction, register as structure
        struct_processor = CifBaseProcessor(name="struct_workflow")
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id=protein_id,
            file_path=str(setup_test_environment / "structure" / "mmcif" / f"AF-{protein_id}-F1.cif"),
            metadata={
                "source": "alphafold",
                "confidence": 0.95,
                "plddt_mean": 87.3
            }
        )
        
        # Verify both formats are tracked
        formats = global_registry.entity_registry.get_entity_formats(entity_id)
        assert "sequence" in formats
        assert "structure" in formats
        
        # Verify we can find by original ID in both formats
        seq_found = global_registry.entity_registry.find_entity_by_original_id(
            protein_id, format_type="sequence"
        )
        struct_found = global_registry.entity_registry.find_entity_by_original_id(
            protein_id, format_type="structure"
        )
        
        assert seq_found == entity_id
        assert struct_found == entity_id
        
        # Verify metadata is format-specific
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        assert entity_info["formats"]["sequence"]["metadata"]["organism"] == "Homo sapiens"
        assert entity_info["formats"]["structure"]["metadata"]["source"] == "alphafold"
    
    def test_complete_entity_lifecycle(self, setup_test_environment):
        """Test entity through complete analysis pipeline."""
        global_registry = GlobalRegistry()
        
        # UniProt ID as canonical identifier
        uniprot_id = "Q9Y5Y4"
        entity_id = generate_entity_id(uniprot_id)
        
        # 1. Start with sequence from UniProt
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=uniprot_id,
            metadata={
                "source": "uniprot",
                "accession": uniprot_id,
                "protein_name": "G-protein coupled receptor 139",
                "organism": "Homo sapiens",
                "length": 456
            },
            datasets=["human_gpcrs", "uniprot_2024"]
        )
        
        # 2. Structure solved/predicted
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id=uniprot_id,
            metadata={
                "pdb_id": "7XYZ",  # Hypothetical PDB entry
                "method": "X-RAY DIFFRACTION",
                "resolution": 2.8,
                "chains": ["A", "B"]
            },
            datasets=["human_gpcrs", "pdb_structures"]
        )
        
        # 3. GRN annotation applied
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="grn",
            original_id=uniprot_id,
            metadata={
                "grn_system": "gpcrdb",
                "family": "Class A",
                "positions_annotated": 127,
                "in_table": "human_gpcr_grn"
            },
            datasets=["human_gpcrs", "gpcrdb_grns"]
        )
        
        # 4. Embeddings generated
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="embedding",
            original_id=uniprot_id,
            metadata={
                "model": "esm2_t33_650M",
                "embedding_dim": 1280,
                "layers": 33
            },
            datasets=["human_gpcrs", "esm2_embeddings"]
        )
        
        # Verify all formats present
        formats = global_registry.entity_registry.get_entity_formats(entity_id)
        assert len(formats) == 4
        for fmt in ["sequence", "structure", "grn", "embedding"]:
            assert fmt in formats
        
        # Verify datasets accumulated correctly
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        assert "human_gpcrs" in entity_info["datasets"]
        assert len(entity_info["datasets"]) == 5  # All unique datasets
    
    def test_entity_id_consistency_across_processors(self, setup_test_environment):
        """Test that different processors generate same entity ID for same identifier."""
        # Test with various biological identifiers
        test_ids = [
            "P12345",      # UniProt
            "1ABC",        # PDB
            "BR1_HUMAN",   # Sequence name
            "GPCR_001"     # Custom ID
        ]
        
        for test_id in test_ids:
            # Each processor should generate the same entity ID
            seq_entity_id = generate_entity_id(test_id)
            struct_entity_id = generate_entity_id(test_id)
            grn_entity_id = generate_entity_id(test_id)
            emb_entity_id = generate_entity_id(test_id)
            
            # All should be identical
            assert seq_entity_id == struct_entity_id
            assert struct_entity_id == grn_entity_id
            assert grn_entity_id == emb_entity_id
            
            # Should be 10-char alphanumeric
            assert len(seq_entity_id) == 10
            assert seq_entity_id.isalnum()
    
    def test_grn_table_multi_entity_format(self, setup_test_environment):
        """Test GRN table with multiple entities in different formats."""
        global_registry = GlobalRegistry()
        grn_processor = GRNBaseProcessor(name="grn_multi")
        
        # Create test GRN table
        grn_data = pd.DataFrame({
            '1.50': ['A', 'V', 'L'],
            '2.50': ['G', 'G', 'P'],
            '3.50': ['I', 'V', 'F'],
            '7.50': ['Y', 'Y', 'W']
        }, index=['PROT1', 'PROT2', 'PROT3'])
        
        # Register each protein in multiple formats
        for protein_id in grn_data.index:
            entity_id = generate_entity_id(protein_id)
            
            # As sequence
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="sequence",
                original_id=protein_id,
                metadata={"source": "test"}
            )
            
            # As GRN entry
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="grn",
                original_id=protein_id,
                metadata={"in_table": "test_table"}
            )
            
            # PROT1 also has structure
            if protein_id == "PROT1":
                global_registry.entity_registry.register_entity(
                    entity_id=entity_id,
                    entity_type="structure",
                    original_id=protein_id,
                    metadata={"pdb_id": "TEST1"}
                )
        
        # Save GRN table with entity IDs
        grn_processor.data = grn_data
        grn_processor.ids = grn_data.index.tolist()
        grn_processor.grns = grn_data.columns.tolist()
        
        # Save to test-data structure
        saved_path = grn_processor.save_grn_table("test_multi_format", include_entity_ids=True)
        
        # Verify entity_id column
        saved_df = pd.read_csv(saved_path, index_col=0)
        assert 'entity_id' in saved_df.columns
        assert saved_df.columns[0] == 'entity_id'
        
        # Verify each entity
        for idx, protein_id in enumerate(grn_data.index):
            expected_entity_id = generate_entity_id(protein_id)
            assert saved_df.loc[protein_id, 'entity_id'] == expected_entity_id
            
            # Check formats
            formats = global_registry.entity_registry.get_entity_formats(expected_entity_id)
            assert "sequence" in formats
            assert "grn" in formats
            if protein_id == "PROT1":
                assert "structure" in formats
    
    def test_cross_format_metadata_independence(self, setup_test_environment):
        """Test that metadata is independent across formats."""
        global_registry = GlobalRegistry()
        
        protein_id = "TEST_PROTEIN"
        entity_id = generate_entity_id(protein_id)
        
        # Register with different metadata in each format
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=protein_id,
            metadata={
                "length": 300,
                "mw": 35000,
                "pi": 6.5
            }
        )
        
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id=protein_id,
            metadata={
                "resolution": 2.1,
                "r_factor": 0.22,
                "chains": 2
            }
        )
        
        # Update structure metadata
        global_registry.entity_registry.update_entity_metadata(
            entity_id,
            {"resolution": 1.8, "method": "cryo-EM"},
            format_type="structure"
        )
        
        # Verify sequence metadata unchanged
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        seq_meta = entity_info["formats"]["sequence"]["metadata"]
        assert seq_meta["length"] == 300
        assert "resolution" not in seq_meta
        
        # Verify structure metadata updated
        struct_meta = entity_info["formats"]["structure"]["metadata"]
        assert struct_meta["resolution"] == 1.8
        assert struct_meta["method"] == "cryo-EM"
        assert struct_meta["r_factor"] == 0.22  # Original preserved
    
    def test_dataset_tracking_across_formats(self, setup_test_environment):
        """Test that datasets are properly tracked across formats."""
        global_registry = GlobalRegistry()
        
        protein_id = "DATASET_TEST"
        entity_id = generate_entity_id(protein_id)
        
        # Register in different formats with different datasets
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="sequence",
            original_id=protein_id,
            datasets=["human_proteins", "kinases"]
        )
        
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="structure",
            original_id=protein_id,
            datasets=["pdb_2024", "kinases"]
        )
        
        global_registry.entity_registry.register_entity(
            entity_id=entity_id,
            entity_type="embedding",
            original_id=protein_id,
            datasets=["esm_embeddings", "human_proteins"]
        )
        
        # Check accumulated datasets
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        all_datasets = entity_info["datasets"]
        
        # Should have union of all datasets
        expected = {"human_proteins", "kinases", "pdb_2024", "esm_embeddings"}
        assert set(all_datasets) == expected
        
        # Test filtering by dataset
        kinase_entities = global_registry.entity_registry.list_entities(dataset="kinases")
        assert entity_id in kinase_entities
        
        human_entities = global_registry.entity_registry.list_entities(dataset="human_proteins")
        assert entity_id in human_entities