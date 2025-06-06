"""
Configuration for pytest.
This file contains fixtures that can be reused across multiple test files.
"""

import os
import sys
import pytest
import pandas as pd
import numpy as np

# Add the project root to the path so we can import protos
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

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
@pytest.fixture
def data_path(tmp_path):
    """Base path for test data."""
    path = tmp_path / "data" / "structure"
    path.mkdir(parents=True, exist_ok=True)
    return str(path)
