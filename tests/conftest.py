"""
Configuration for pytest.
This file contains fixtures that can be reused across multiple test files.

IMPORTANT: All tests use a single data root (tests/test-data) managed by ProtosPaths.
Tests should NEVER create their own directories or use tmp_path for data.
"""

import os
import sys
import pytest
import pandas as pd
import numpy as np
from pathlib import Path

# Add the project root to the path so we can import protos
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import ProtosPaths for global configuration
from protos.io.paths.path_config import ProtosPaths

@pytest.fixture(scope="session", autouse=True)
def configure_test_paths():
    """
    Configure ProtosPaths to use test-data directory for ALL tests.
    
    This is the ONLY place where data root should be set.
    All processors will automatically use this path.
    """
    # Set global data root to tests/test-data
    test_data_root = Path(__file__).parent / "test-data"
    ProtosPaths.set_data_root(str(test_data_root.absolute()))
    
    # Create the complete directory structure using ProtosPaths
    paths = ProtosPaths(create_dirs=True)
    
    yield paths
    
    # Reset after tests
    ProtosPaths.set_data_root(None)

@pytest.fixture
def sample_grn_data():
    """Create a sample GRN table for testing."""
    data = {
        '1x50': ['N42', 'N35', 'N30'],
        '2x50': ['D83', 'D74', 'D71'],
        '3x50': ['R135', 'R131', 'R128']
    }
    index = ['protein1', 'protein2', 'protein3']
    return pd.DataFrame(data, index=index)

@pytest.fixture
def sample_sequence_data():
    """Create sample sequence data for testing."""
    return {
        'protein1': 'MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA',
        'protein2': 'MNGTEGLNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLALTWVMALACAAPPLAGWSRYIPEGLQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIIIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSNFGPVFMTIPAFFAKSASIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASATVSKTETSQVAPA',
        'protein3': 'MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLAFTWVMALACAAPPLAGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSNFGPVFMTIPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASATVSKTETSQVAPA'
    }

@pytest.fixture
def sample_structure_data():
    """Create a minimal sample structure dataframe for testing."""
    data = {
        'pdb_id': ['1u19', '1u19', '1u19', '6oz2', '6oz2'],
        'auth_chain_id': ['A', 'A', 'A', 'A', 'A'],
        'group': ['ATOM', 'ATOM', 'ATOM', 'ATOM', 'ATOM'],
        'res_name3l': ['ALA', 'ARG', 'GLY', 'PHE', 'SER'],
        'res_name1l': ['A', 'R', 'G', 'F', 'S'],
        'res_atom_name': ['CA', 'CA', 'CA', 'CA', 'CA'],
        'atom_id': [1, 2, 3, 4, 5],
        'x': [10.0, 13.5, 17.0, 20.5, 24.0],
        'y': [5.0, 8.5, 12.0, 15.5, 19.0],
        'z': [2.0, 5.5, 9.0, 12.5, 16.0],
        'gen_seq_id': [1, 2, 3, 4, 5]
    }
    return pd.DataFrame(data)

@pytest.fixture
def sample_embedding_data():
    """Create sample embedding data for testing."""
    # Create mock embeddings for 3 proteins, each with 5 residues and 10 embedding dimensions
    embeddings = {}
    for prot_id in ['protein1', 'protein2', 'protein3']:
        # Create random embedding matrix (5 residues x 10 dimensions)
        embeddings[prot_id] = np.random.random((5, 10))
    return embeddings

@pytest.fixture(scope="session")
def test_structures():
    """
    Ensure test structures are available using ProtosPaths.
    
    This fixture uses the download functionality to get real structures
    if they're not already present in test-data.
    """
    from protos.processing.structure.struct_base_processor import CifBaseProcessor
    from protos.loaders.download_structures import download_protein_structures
    
    # Create processor (uses ProtosPaths automatically)
    processor = CifBaseProcessor(name="test_fixtures")
    
    # List of small PDB structures for testing
    test_pdbs = ["1ubq", "1tqn", "3nir"]
    available_structures = []
    
    # Check which structures we already have
    for pdb_id in test_pdbs:
        try:
            # Try to load - if it works, we have it
            processor.load_structure(pdb_id)
            available_structures.append(pdb_id)
        except:
            # Structure not available, try to download using processor
            try:
                from protos.loaders.download_structures import download_structures_with_processor
                successful, failed = download_structures_with_processor([pdb_id], processor)
                if successful:
                    available_structures.extend(successful)
            except Exception as e:
                print(f"Could not download {pdb_id}: {e}")
                continue
    
    return available_structures


def ensure_grn_reference_data():
    """
    Ensure GRN reference data is available for tests.
    
    Creates minimal reference data if not present, using ProtosPaths.
    """
    from protos.processing.grn.grn_base_processor import GRNBaseProcessor
    
    # Create GRN processor (uses ProtosPaths automatically)
    processor = GRNBaseProcessor(name="test_grn_refs")
    
    # Check if we have minimal reference data
    ref_files_needed = [
        "mo_ref",  # Microbial opsin reference
        "gpcr_ref",  # GPCR reference
    ]
    
    for ref_file in ref_files_needed:
        ref_path = f"ref/{ref_file}"
        if not processor.is_dataset_available(ref_path):
            # Create minimal reference data
            create_minimal_grn_reference(processor, ref_file)
    
    return ref_files_needed


def create_minimal_grn_reference(processor, ref_type):
    """Create minimal GRN reference data for testing."""
    import pandas as pd
    
    if ref_type == "mo_ref":
        # Minimal microbial opsin reference
        data = {
            '1.50': ['M1', 'M1', 'M1'],
            '2.50': ['D2', 'D2', 'D2'], 
            '3.50': ['W3', 'W3', 'W3'],
            '7.50': ['K7', 'K7', 'K7']  # Schiff base
        }
        index = ['TEST_MO_1', 'TEST_MO_2', 'TEST_MO_3']
    elif ref_type == "gpcr_ref":
        # Minimal GPCR reference
        data = {
            '1.50': ['N1', 'N1', 'N1'],
            '2.50': ['D2', 'D2', 'D2'],
            '3.50': ['R3', 'R3', 'R3'],
            '7.43': ['K7', 'K7', 'K7']  # GPCR Schiff base
        }
        index = ['TEST_GPCR_1', 'TEST_GPCR_2', 'TEST_GPCR_3']
    else:
        # Generic reference
        data = {
            '1.50': ['X1', 'X1', 'X1'],
            '2.50': ['X2', 'X2', 'X2'],
            '3.50': ['X3', 'X3', 'X3']
        }
        index = ['TEST_REF_1', 'TEST_REF_2', 'TEST_REF_3']
    
    ref_df = pd.DataFrame(data, index=index)
    processor.data = ref_df
    processor.save_grn_table(f"ref/{ref_type}")


@pytest.fixture(scope="session")
def test_grn_references():
    """
    Ensure GRN reference data is available for tests.
    """
    return ensure_grn_reference_data()


def ensure_sequence_test_data():
    """
    Ensure sequence test data is available for tests.
    """
    from protos.processing.sequence.seq_processor import SeqProcessor
    
    # Create sequence processor (uses ProtosPaths automatically)
    processor = SeqProcessor(name="test_seq_refs")
    
    # Test sequences for various purposes
    test_sequences = {
        'TEST_PROTEIN_1': 'MTNSDKWIFGLQELCVEMLAYGFQFGHILGWQYLVSPFHWQGTRTRDYPYGSMVTLQPYLQLVN',
        'TEST_PROTEIN_2': 'MSDKWIFGLQELCVEMLAYGFQFGHILGWQYLVSPFHWQGTRTRDYPYGSMVTLQPYLQLVNARK',
        'TEST_OPSIN_1': 'MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLR',
        'TEST_GPCR_1': 'MSGNTDTLTLVEMSYAQKFGYDVPQCLNDTVTNESSVAKFDKLTNYLVRNGLFRCTFHKFPSLEVVKLIV'
    }
    
    available_sequences = []
    for seq_id, sequence in test_sequences.items():
        try:
            # Try to load sequence to see if it exists
            entity_id = processor.get_entity_id_for_sequence(sequence)
            if entity_id:
                available_sequences.append(seq_id)
        except:
            # Save the sequence
            try:
                processor.save_sequence_entity(seq_id, sequence)
                available_sequences.append(seq_id)
            except Exception as e:
                print(f"Could not save test sequence {seq_id}: {e}")
    
    return available_sequences


@pytest.fixture(scope="session") 
def test_sequence_data():
    """
    Ensure sequence test data is available for tests.
    """
    return ensure_sequence_test_data()