"""
Test cross-format workflows in Protos.

These tests demonstrate how entities can be tracked and transformed across
different data formats while maintaining consistent identity.
"""

import pytest
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


class TestSequenceToStructureWorkflow:
    """Test sequence to structure prediction workflow (e.g., AlphaFold)."""
    
    def test_sequence_to_alphafold_structure(self, setup_test_environment):
        """Test workflow where a sequence is used to predict a structure."""
        # Initialize processors
        seq_proc = SeqProcessor(name="seq_to_struct")
        struct_proc = CifBaseProcessor(name="seq_to_struct")
        
        # 1. Start with a sequence
        protein_id = "TEST_PROTEIN"
        sequence = "MKLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
        
        # Save sequence (automatically registers entity)
        entity_id = seq_proc.save_sequence_entity(protein_id, sequence)
        
        # 2. Simulate AlphaFold prediction (in reality, this would be external)
        # For testing, we'll create a mock structure and register it
        mock_structure_df = pd.DataFrame({
            'pdb_id': [f'AF-{protein_id}-F1'] * 10,
            'model_id': [1] * 10,
            'auth_chain_id': ['A'] * 10,
            'auth_seq_id': range(1, 11),
            'auth_comp_id': list(sequence[:10]),
            'atom_name': ['CA'] * 10,
            'x': np.random.rand(10) * 50,
            'y': np.random.rand(10) * 50,
            'z': np.random.rand(10) * 50
        })
        
        # Save structure with same entity ID
        struct_entity_id = struct_proc.save_structure_as_entity(
            structure_df=mock_structure_df,
            pdb_id=f'AF-{protein_id}-F1',
            metadata={
                'source': 'alphafold',
                'confidence': 0.95,
                'original_sequence_id': protein_id
            }
        )
        
        # 3. Verify both formats share the same entity
        assert entity_id == struct_entity_id
        
        # 4. Verify entity has both formats registered
        global_registry = GlobalRegistry()
        entity_info = global_registry.entity_registry.get_entity(entity_id)
        assert entity_info is not None
        assert 'sequence' in entity_info['formats']
        assert 'structure' in entity_info['formats']
        
        # 5. Verify metadata is preserved
        seq_meta = entity_info['formats']['sequence']['metadata']
        struct_meta = entity_info['formats']['structure']['metadata']
        assert struct_meta['source'] == 'alphafold'
        assert struct_meta['original_sequence_id'] == protein_id
    
    def test_batch_sequence_to_structure_workflow(self, setup_test_environment):
        """Test batch processing of sequences to structures."""
        seq_proc = SeqProcessor(name="batch_seq_to_struct")
        struct_proc = CifBaseProcessor(name="batch_seq_to_struct")
        
        # Multiple sequences
        sequences = {
            "PROT1": "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT",
            "PROT2": "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKH",
            "PROT3": "MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVA"
        }
        
        entity_mappings = {}
        
        # Process each sequence
        for protein_id, sequence in sequences.items():
            # Save sequence
            seq_entity_id = seq_proc.save_sequence_entity(protein_id, sequence)
            
            # Create mock structure
            mock_structure_df = pd.DataFrame({
                'pdb_id': [f'AF-{protein_id}-F1'] * len(sequence),
                'model_id': [1] * len(sequence),
                'auth_chain_id': ['A'] * len(sequence),
                'auth_seq_id': range(1, len(sequence) + 1),
                'auth_comp_id': list(sequence),
                'atom_name': ['CA'] * len(sequence),
                'x': np.random.rand(len(sequence)) * 50,
                'y': np.random.rand(len(sequence)) * 50,
                'z': np.random.rand(len(sequence)) * 50
            })
            
            # Save structure with same entity
            struct_entity_id = struct_proc.save_structure_as_entity(
                structure_df=mock_structure_df,
                pdb_id=f'AF-{protein_id}-F1',
                metadata={'predicted_from': protein_id}
            )
            
            # Track mappings
            entity_mappings[protein_id] = seq_entity_id
            assert seq_entity_id == struct_entity_id
        
        # Verify all entities have both formats
        global_registry = GlobalRegistry()
        for protein_id, entity_id in entity_mappings.items():
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            assert 'sequence' in entity_info['formats']
            assert 'structure' in entity_info['formats']


class TestStructureToSequenceWorkflow:
    """Test structure to sequence extraction workflow."""
    
    def test_extract_sequence_from_structure(self, setup_test_environment):
        """Test extracting sequence from a loaded structure."""
        struct_proc = CifBaseProcessor(name="struct_to_seq")
        seq_proc = SeqProcessor(name="struct_to_seq")
        
        # Load a real structure
        pdb_id = "1ubq"
        structure_df = struct_proc.load_structure(pdb_id)
        
        if structure_df is not None and len(structure_df) > 0:
            # Extract sequence from chain A
            chain_a_df = structure_df[structure_df['auth_chain_id'] == 'A']
            sequence = struct_proc.extract_sequence_from_dataframe(chain_a_df)
            
            # Save extracted sequence with entity linking
            entity_id = generate_entity_id(pdb_id)
            seq_entity_id = seq_proc.save_sequence_entity(
                f"{pdb_id}_chain_A",
                sequence,
                metadata={
                    'extracted_from_pdb': pdb_id,
                    'chain': 'A',
                    'extraction_method': 'ca_atoms'
                }
            )
            
            # Verify entity tracking
            global_registry = GlobalRegistry()
            entity_info = global_registry.entity_registry.get_entity(seq_entity_id)
            assert entity_info is not None
            assert entity_info['formats']['sequence']['metadata']['extracted_from_pdb'] == pdb_id
    
    def test_extract_all_chains_workflow(self, setup_test_environment):
        """Test extracting sequences from all chains in a structure."""
        struct_proc = CifBaseProcessor(name="all_chains")
        seq_proc = SeqProcessor(name="all_chains")
        
        # Load structure
        pdb_id = "1ubq"
        structure_df = struct_proc.load_structure(pdb_id)
        
        if structure_df is not None:
            # Get unique chains
            chains = structure_df['auth_chain_id'].unique()
            
            extracted_sequences = {}
            
            for chain in chains:
                # Extract sequence for each chain
                chain_df = structure_df[structure_df['auth_chain_id'] == chain]
                if len(chain_df) > 0:
                    sequence = struct_proc.extract_sequence_from_dataframe(chain_df)
                    
                    # Save with metadata
                    seq_entity_id = seq_proc.save_sequence_entity(
                        f"{pdb_id}_chain_{chain}",
                        sequence,
                        metadata={
                            'source_pdb': pdb_id,
                            'chain': chain,
                            'parent_entity': generate_entity_id(pdb_id)
                        }
                    )
                    
                    extracted_sequences[chain] = {
                        'sequence': sequence,
                        'entity_id': seq_entity_id
                    }
            
            # Verify at least one chain was extracted
            assert len(extracted_sequences) > 0


class TestSequenceToGRNWorkflow:
    """Test sequence to GRN assignment workflow."""
    
    def test_sequence_to_grn_assignment(self, setup_test_environment):
        """Test assigning GRNs to a sequence."""
        seq_proc = SeqProcessor(name="seq_to_grn")
        grn_proc = GRNBaseProcessor(name="seq_to_grn")
        
        # Start with a sequence
        protein_id = "TEST_GPCR"
        sequence = "MNGTEGPNFYVPFSNATGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVLGGFTSTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWVPYASVAFYIFTHQGSNFGPIFMTIPAFFAKSAAIYNPVIYIMMNKQFRNCMLTTICCGKNPLGDDEASATVSKTETSQVAPA"
        
        # Save sequence
        seq_entity_id = seq_proc.save_sequence_entity(protein_id, sequence)
        
        # Create test GRN reference
        grn_ref_data = pd.DataFrame({
            '1.50': ['N', 'V', 'L', 'I', 'A'],
            '2.50': ['L', 'L', 'P', 'A', 'G'],
            '3.50': ['V', 'I', 'F', 'L', 'M'],
            '7.50': ['Y', 'Y', 'W', 'F', 'H']
        }, index=['REF1', 'REF2', 'REF3', 'REF4', 'REF5'])
        
        grn_proc.data = grn_ref_data
        grn_proc.ids = grn_ref_data.index.tolist()
        grn_proc.grns = grn_ref_data.columns.tolist()
        
        # Save reference table
        grn_proc.save_grn_table("test_grn_ref")
        
        # Perform GRN assignment
        grn_assignment = grn_proc.assign_grns(
            sequence=sequence,
            protein_id=protein_id,
            reference_id="REF1",
            use_cached=False
        )
        
        # Save assignment result with entity tracking
        if grn_assignment is not None:
            # Create single-row GRN table for this sequence
            grn_row_data = pd.DataFrame([grn_assignment], index=[protein_id])
            grn_proc.data = grn_row_data
            grn_proc.ids = [protein_id]
            
            # Save with entity ID
            saved_path = grn_proc.save_grn_table(
                f"grn_assignment_{protein_id}",
                include_entity_ids=True
            )
            
            # Verify entity ID column
            saved_df = pd.read_csv(saved_path, index_col=0)
            assert 'entity_id' in saved_df.columns
            assert saved_df.loc[protein_id, 'entity_id'] == seq_entity_id


class TestAnyFormatToEmbeddingsWorkflow:
    """Test generating embeddings from any format."""
    
    def test_sequence_to_embeddings(self, setup_test_environment):
        """Test generating embeddings from sequences."""
        try:
            import torch
        except ImportError:
            pytest.skip("PyTorch not installed")
        
        seq_proc = SeqProcessor(name="seq_to_emb")
        emb_proc = EmbeddingProcessor(name="seq_to_emb")
        
        # Create test sequences
        sequences = {
            "PROT1": "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT",
            "PROT2": "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKH"
        }
        
        # Save sequences
        entity_ids = {}
        for protein_id, sequence in sequences.items():
            entity_id = seq_proc.save_sequence_entity(protein_id, sequence)
            entity_ids[protein_id] = entity_id
        
        # Mock embedding generation (in reality would use ESM model)
        for protein_id, entity_id in entity_ids.items():
            # Create mock embedding
            mock_embedding = torch.randn(len(sequences[protein_id]), 320)
            
            # Save embedding with entity tracking
            emb_path = setup_test_environment / "embedding" / f"{entity_id}.pkl"
            emb_path.parent.mkdir(parents=True, exist_ok=True)
            torch.save(mock_embedding, emb_path)
            
            # Register embedding format for entity
            global_registry = GlobalRegistry()
            global_registry.entity_registry.register_entity(
                entity_id=entity_id,
                entity_type="embedding",
                original_id=protein_id,
                file_path=str(emb_path),
                metadata={
                    'model': 'test_model',
                    'dim': 320,
                    'source_format': 'sequence'
                }
            )
        
        # Verify entities have both sequence and embedding formats
        for protein_id, entity_id in entity_ids.items():
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            assert 'sequence' in entity_info['formats']
            assert 'embedding' in entity_info['formats']
            assert entity_info['formats']['embedding']['metadata']['source_format'] == 'sequence'


class TestConversionLineageTracking:
    """Test tracking conversion lineage in metadata."""
    
    def test_full_workflow_lineage(self, setup_test_environment):
        """Test tracking lineage through a complete workflow."""
        seq_proc = SeqProcessor(name="lineage")
        struct_proc = CifBaseProcessor(name="lineage")
        grn_proc = GRNBaseProcessor(name="lineage")
        
        # Start with a sequence
        protein_id = "LINEAGE_TEST"
        sequence = "MKLSPADKTNVKAAWGKVGAHAGEYGAEAL"
        
        # Step 1: Save sequence with initial metadata
        seq_entity_id = seq_proc.save_sequence_entity(
            protein_id,
            sequence,
            metadata={
                'source': 'user_input',
                'date': '2025-06-27',
                'lineage': []
            }
        )
        
        # Step 2: Predict structure (mock)
        mock_structure_df = pd.DataFrame({
            'pdb_id': [f'AF-{protein_id}-F1'] * len(sequence),
            'model_id': [1] * len(sequence),
            'auth_chain_id': ['A'] * len(sequence),
            'auth_seq_id': range(1, len(sequence) + 1),
            'auth_comp_id': list(sequence),
            'atom_name': ['CA'] * len(sequence),
            'x': np.random.rand(len(sequence)) * 50,
            'y': np.random.rand(len(sequence)) * 50,
            'z': np.random.rand(len(sequence)) * 50
        })
        
        # Save with lineage tracking
        struct_entity_id = struct_proc.save_structure_as_entity(
            structure_df=mock_structure_df,
            pdb_id=f'AF-{protein_id}-F1',
            metadata={
                'prediction_method': 'alphafold',
                'source_sequence': protein_id,
                'lineage': [
                    {
                        'step': 1,
                        'operation': 'sequence_to_structure',
                        'source_format': 'sequence',
                        'source_id': protein_id,
                        'method': 'alphafold_v2'
                    }
                ]
            }
        )
        
        assert seq_entity_id == struct_entity_id
        
        # Step 3: Extract sequence back from structure
        extracted_seq = struct_proc.extract_sequence_from_dataframe(mock_structure_df)
        
        # Save extracted sequence with lineage
        extracted_entity_id = seq_proc.save_sequence_entity(
            f"{protein_id}_extracted",
            extracted_seq,
            metadata={
                'extraction_source': f'AF-{protein_id}-F1',
                'lineage': [
                    {
                        'step': 1,
                        'operation': 'sequence_to_structure',
                        'source_format': 'sequence',
                        'source_id': protein_id,
                        'method': 'alphafold_v2'
                    },
                    {
                        'step': 2,
                        'operation': 'structure_to_sequence',
                        'source_format': 'structure',
                        'source_id': f'AF-{protein_id}-F1',
                        'method': 'ca_extraction'
                    }
                ]
            }
        )
        
        # Verify lineage tracking
        global_registry = GlobalRegistry()
        entity_info = global_registry.entity_registry.get_entity(struct_entity_id)
        
        # Check structure metadata has lineage
        struct_lineage = entity_info['formats']['structure']['metadata']['lineage']
        assert len(struct_lineage) == 1
        assert struct_lineage[0]['operation'] == 'sequence_to_structure'
        assert struct_lineage[0]['source_id'] == protein_id