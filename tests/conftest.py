import os
import sys
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
from typing import Optional

# Add the project root to the path so we can import protos
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

# Import ProtosPaths for global configuration
from protos.io.paths.path_config import ProtosPaths, get_default_data_root


@pytest.fixture(scope="session", autouse=True)
def configure_test_paths():
    """
    Configure ProtosPaths to use test-data directory for ALL tests.

    This is the ONLY place where data root should be set.
    All processors will automatically use this path.
    """
    test_data_root = Path(__file__).parent / "test-data"
    test_data_root.mkdir(exist_ok=True)

    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())

    # Create an instance
    paths = ProtosPaths()

    # Ensure key directories are created (since lazy, call getters)
    for processor_type in paths.processor_dirs:
        paths.get_processor_path(processor_type)
        if processor_type in paths.subdirs:
            for subdir_type in paths.subdirs[processor_type]:
                paths.get_subdir_path(processor_type, subdir_type)

    paths.get_global_registry_path()

    # Run setup script to copy reference data if it exists
    setup_script = Path(__file__).parent.parent / "setup_test_data_from_reference.py"
    if setup_script.exists():
        try:
            subprocess.run([sys.executable, str(setup_script)], check=True)
        except subprocess.CalledProcessError:
            print("Warning: Could not set up reference test data")

    yield paths

    # Cleanup
    if "PROTOS_DATA_ROOT" in os.environ:
        del os.environ["PROTOS_DATA_ROOT"]


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
        'protein3': 'MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIVLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHGYFVFGPTGGNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGLAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTLKPEVNNESFVIYMFVVHFTIPMIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSNFGPVFMTIPAFFAKSSSIYNPVIYIMMNKQFRNCMLTTLCCGKNPLGDDEASATVSKTETSQVAPA'
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
        'atom_name': ['CA', 'CA', 'CA', 'CA', 'CA'],
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
    embeddings = {}
    for prot_id in ['protein1', 'protein2', 'protein3']:
        embeddings[prot_id] = np.random.random((5, 10))
    return embeddings


@pytest.fixture(scope="session")
def test_structures():
    """
    Ensure test structures are available using ProtosPaths.
    """
    from protos.processing.structure import StructureProcessor

    processor = StructureProcessor(name="test_fixtures")

    test_pdbs = ["1ubq", "1tqn", "3nir"]
    available_structures = []

    for pdb_id in test_pdbs:
        try:
            processor.load_structure(pdb_id)
            available_structures.append(pdb_id)
        except:
            try:
                from protos.loaders.download_structures import download_structures_with_processor
                successful, failed = download_structures_with_processor([pdb_id], processor)
                if pdb_id in successful:
                    available_structures.append(pdb_id)
            except Exception as e:
                print(f"Could not download {pdb_id}: {e}")

    return available_structures


@pytest.fixture(scope="session")
def test_grn_references():
    """
    Ensure GRN reference data is available for tests.
    """
    return ensure_grn_reference_data()


def ensure_grn_reference_data():
    from protos.processing.grn import GRNProcessor

    processor = GRNProcessor(name="test_grn_refs")

    ref_files_needed = ["mo_ref", "gpcr_ref"]

    for ref_file in ref_files_needed:
        ref_path = f"ref/{ref_file}"
        if not processor.is_dataset_available(ref_path):
            create_minimal_grn_reference(processor, ref_file)

    return ref_files_needed


def create_minimal_grn_reference(processor, ref_type):
    import pandas as pd

    if ref_type == "mo_ref":
        data = {
            '1.50': ['M1', 'M1', 'M1'],
            '2.50': ['D2', 'D2', 'D2'],
            '3.50': ['W3', 'W3', 'W3'],
            '7.50': ['K7', 'K7', 'K7']
        }
        index = ['TEST_MO_1', 'TEST_MO_2', 'TEST_MO_3']
    elif ref_type == "gpcr_ref":
        data = {
            '1.50': ['N1', 'N1', 'N1'],
            '2.50': ['D2', 'D2', 'D2'],
            '3.50': ['R3', 'R3', 'R3'],
            '7.43': ['K7', 'K7', 'K7']
        }
        index = ['TEST_GPCR_1', 'TEST_GPCR_2', 'TEST_GPCR_3']
    else:
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
def test_sequence_data():
    """
    Ensure sequence test data is available for tests.
    """
    return ensure_sequence_test_data()


def ensure_sequence_test_data():
    from protos.processing.sequence import SequenceProcessor

    processor = SequenceProcessor(name="test_seq_refs")

    test_sequences = {
        'TEST_PROTEIN_1': 'MTNSDKWIFGLQELCVEMLAYGFQFGHILGWQYLVSPFHWQGTRTRDYPYGSMVTLQPYLQLVN',
        'TEST_PROTEIN_2': 'MSDKWIFGLQELCVEMLAYGFQFGHILGWQYLVSPFHWQGTRTRDYPYGSMVTLQPYLQLVNARK',
        'TEST_OPSIN_1': 'MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLR',
        'TEST_GPCR_1': 'MSGNTDTLTLVEMSYAQKFGYDVPQCLNDTVTNESSVAKFDKLTNYLVRNGLFRCTFHKFPSLEVVKLIV'
    }

    available_sequences = []
    for seq_id, sequence in test_sequences.items():
        try:
            entity_id = processor.get_entity_id_for_sequence(sequence)
            if entity_id:
                available_sequences.append(seq_id)
        except:
            try:
                processor.save_sequence_entity(seq_id, sequence)
                available_sequences.append(seq_id)
            except Exception as e:
                print(f"Could not save test sequence {seq_id}: {e}")

    return available_sequences