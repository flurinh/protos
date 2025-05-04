"""
Tests for the UniProt loader functionality in the protos package.

This module tests downloading and processing protein sequences from UniProt.
"""

import os
import tempfile
# Set RUN_NETWORK_TESTS=1 by default
os.environ["RUN_NETWORK_TESTS"] = "1"
import pytest
import pandas as pd
from pathlib import Path

from protos.io.fasta_utils import read_fasta, write_fasta
from protos.loaders.uniprot_utils import get_uniprot, map_uniprot_to_pdb
from protos.loaders.uniprot_loader import UniprotDL
from protos.io.paths.path_config import ProtosPaths, DataSource, join_path


@pytest.fixture
def test_uniprot_ids():
    """Define test UniProt IDs to use for real data tests."""
    # Return list of UniProt IDs for proteins with diverse features
    return ["P00533", "P01308", "P05067", "P02751"]


@pytest.fixture
def test_paths():
    """Create a temporary directory with ProtosPaths structure for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize ProtosPaths with temp directory
        paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True)
        
        yield paths, temp_dir


@pytest.fixture
def prepared_test_environment(test_paths, test_uniprot_ids):
    """Create a properly structured test environment with dataset file."""
    # Get test paths
    paths, temp_dir = test_paths
    
    # Create dataset name
    dataset_name = "test_dataset"
    
    # Get metadata directory where dataset files should be stored
    metadata_dir = paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    
    # Create dataset file with test IDs
    dataset_file = join_path(metadata_dir, f"{dataset_name}.txt")
    
    # Write test IDs to dataset file
    os.makedirs(os.path.dirname(dataset_file), exist_ok=True)
    with open(dataset_file, 'w') as f:
        f.write(" ".join(test_uniprot_ids[:2]))  # Use only 2 IDs for speed
    
    # Return paths object and dataset name
    yield paths, dataset_name


@pytest.fixture
def mock_uniprot_data():
    """Create mock data for testing without network calls."""
    def _create_mock_data(uniprot_id, sequence="TESTSEQUENCE"):
        return {
            "uniprot": uniprot_id,
            "seq": sequence,
            "gene": f"{uniprot_id}_GENE",
            "species": "TEST_SPECIES",
            "organism": "Test Organism",
            "info": f"Test info for {uniprot_id}",
            "dataset": "test_dataset"
        }
    return _create_mock_data


# Tests for utility functions - they don't need changes

@pytest.mark.skipif(not os.environ.get("RUN_NETWORK_TESTS"), 
                   reason="Network-dependent tests are disabled")
def test_get_uniprot(test_uniprot_ids):
    """Test downloading a single protein from UniProt."""
    # Only test the first ID to avoid too many requests
    test_id = test_uniprot_ids[0]  # P00533 - EGFR
    
    # Download the information
    uniprot_data = get_uniprot(test_id, reviewed=True)
    
    # Verify data format
    assert isinstance(uniprot_data, pd.DataFrame)
    assert not uniprot_data.empty
    assert "Sequence" in uniprot_data.columns
    assert "Entry" in uniprot_data.columns
    
    # Verify content
    assert len(uniprot_data) >= 1
    assert uniprot_data["Entry"].iloc[0] == test_id
    
    # Verify sequence data
    sequence = uniprot_data["Sequence"].iloc[0]
    assert len(sequence) > 0
    assert all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence[:20])


@pytest.mark.skipif(not os.environ.get("RUN_NETWORK_TESTS"), 
                   reason="Network-dependent tests are disabled")
def test_map_uniprot_to_pdb(test_uniprot_ids):
    """Test mapping UniProt IDs to PDB structures."""
    # Test mapping a couple of IDs
    test_ids = test_uniprot_ids[:2]
    
    # Perform the mapping
    mapping_results = map_uniprot_to_pdb(test_ids)
    
    # Verify the results
    assert isinstance(mapping_results, pd.DataFrame)
    assert "uid" in mapping_results.columns
    assert "pdb_id" in mapping_results.columns
    
    # Check that at least some of our test IDs were mapped
    assert not mapping_results.empty
    uids = set(mapping_results["uid"].values)
    assert any(uid in uids for uid in test_ids)


# Tests using the ProtosPaths approach

def test_uniprot_loader_initialization(prepared_test_environment):
    """Test initializing the UniProt loader with ProtosPaths."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader with data_root
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=2
    )
    
    # Verify initialization
    assert loader.dataset == dataset_name
    assert loader.limit == 2
    
    # Verify paths were set up correctly
    assert loader.fasta_dir == paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert loader.metadata_dir == paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    assert os.path.exists(loader.dataset_file)


def test_uniprot_loader_load_dataset(prepared_test_environment, test_uniprot_ids):
    """Test loading a dataset using the loader."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name
    )
    
    # Load the dataset
    loaded_ids = loader.load_dataset()
    
    # Verify loaded IDs
    assert len(loaded_ids) == 2
    assert loaded_ids[0] in test_uniprot_ids
    assert loaded_ids[1] in test_uniprot_ids


@pytest.mark.skipif(not os.environ.get("RUN_NETWORK_TESTS"), 
                   reason="Network-dependent tests are disabled")
def test_uniprot_loader_download_gene(prepared_test_environment, test_uniprot_ids):
    """Test downloading a single gene from UniProt."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name
    )
    
    # Download a single gene
    test_id = test_uniprot_ids[0]
    result = loader.download_gene_single_query(test_id)
    
    # Verify result format
    assert isinstance(result, list)
    assert len(result) == 7  # 7 columns as defined in COLS
    
    # Verify result content
    assert result[0] == test_id  # uniprot ID
    assert len(result[1]) > 0    # sequence
    assert result[6] == dataset_name  # dataset


def test_uniprot_loader_save_fasta(prepared_test_environment, test_uniprot_ids, mock_uniprot_data, monkeypatch):
    """Test saving sequences as FASTA files in the standard location."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=1
    )
    
    # Mock download function to avoid network calls
    def mock_download(uniprot_id):
        data = mock_uniprot_data(uniprot_id)
        return [
            data["uniprot"], 
            data["seq"], 
            data["gene"], 
            data["species"], 
            data["organism"], 
            data["info"], 
            data["dataset"]
        ]
    
    monkeypatch.setattr(loader, "download_gene_single_query", mock_download)
    
    # Load dataset and mock download a gene
    loader.load_dataset()
    test_id = loader.genes[0]
    
    # Populate the data
    loader.data_list.append(mock_download(test_id))
    loader.gene_df = pd.DataFrame(loader.data_list, columns=['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset'])
    
    # Save as FASTA file (should go to standard fasta directory)
    saved_files = loader.save_uniprot_fasta(uniprot=test_id, mode='entry')
    
    # Verify file was created in the standard location
    assert len(saved_files) == 1
    
    # Verify correct location
    standard_fasta_dir = paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert os.path.dirname(saved_files[0]) == standard_fasta_dir
    
    # Verify filename and content
    assert test_id in os.path.basename(saved_files[0])
    
    # Read the file
    loaded_sequences = read_fasta(saved_files[0])
    assert test_id in loaded_sequences
    assert loaded_sequences[test_id] == "TESTSEQUENCE"


@pytest.mark.skipif(not os.environ.get("RUN_NETWORK_TESTS"), 
                   reason="Network-dependent tests are disabled")
def test_uniprot_loader_integrated(prepared_test_environment, test_uniprot_ids, mock_uniprot_data, monkeypatch):
    """Integration test for the UniProt loader with standardized paths."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=2
    )
    
    # Mock download functions to avoid network calls
    def mock_download_single(uniprot_id):
        data = mock_uniprot_data(uniprot_id)
        return [
            data["uniprot"], 
            data["seq"], 
            data["gene"], 
            data["species"], 
            data["organism"], 
            data["info"], 
            data["dataset"]
        ]
    
    def mock_download_batch(*args, **kwargs):
        # Populate loader.data_list and loader.gene_df with test data
        for uid in loader.genes:
            loader.data_list.append(mock_download_single(uid))
        
        loader.gene_df = pd.DataFrame(
            loader.data_list, 
            columns=['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset']
        )
        return loader.gene_df
    
    monkeypatch.setattr(loader, "download_gene_single_query", mock_download_single)
    monkeypatch.setattr(loader, "download_genes_single_query", mock_download_batch)
    
    # 1. Load dataset
    loader.load_dataset()
    assert len(loader.genes) == 2
    
    # 2. Download genes
    loader.download_genes_single_query(batchsize=2, save=True)
    assert len(loader.gene_df) == 2
    
    # 3. Save as FASTA in both modes
    entry_files = loader.save_uniprot_fasta(mode='entry')
    db_file = loader.save_uniprot_fasta(mode='database')
    
    # 4. Verify individual files
    for uid in loader.genes:
        # Find the corresponding file
        uid_file = [f for f in entry_files if uid in os.path.basename(f)]
        assert len(uid_file) == 1
        
        # Verify the file is in the standard fasta directory
        standard_fasta_dir = paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
        assert os.path.dirname(uid_file[0]) == standard_fasta_dir
        
        # Verify content
        loaded_sequences = read_fasta(uid_file[0])
        assert uid in loaded_sequences
    
    # 5. Verify database file
    assert len(db_file) == 1
    db_path = db_file[0]
    
    # Verify database file location (should be in metadata directory)
    metadata_dir = paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    assert os.path.dirname(db_path) == metadata_dir
    
    # Verify database content
    db_sequences = read_fasta(db_path)
    assert len(db_sequences) == 2
    for uid in loader.genes:
        assert uid in db_sequences


def test_save_to_standard_location(prepared_test_environment, test_uniprot_ids, mock_uniprot_data, monkeypatch):
    """Test saving a UniProt sequence to standard location."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=1
    )
    
    # Mock download function
    def mock_download(uniprot_id):
        data = mock_uniprot_data(uniprot_id)
        return [
            data["uniprot"], 
            data["seq"], 
            data["gene"], 
            data["species"], 
            data["organism"], 
            data["info"], 
            data["dataset"]
        ]
    
    monkeypatch.setattr(loader, "download_gene_single_query", mock_download)
    
    # Load dataset and add mock data
    loader.load_dataset()
    test_id = loader.genes[0]
    
    loader.data_list.append(mock_download(test_id))
    loader.gene_df = pd.DataFrame(loader.data_list, columns=loader.gene_df.columns if not loader.gene_df.empty else ['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset'])
    
    # Save to standard location
    file_path = loader.save_to_standard_location(test_id)
    
    # Verify file exists in standard location
    assert os.path.exists(file_path)
    
    # Verify it's in the standard fasta directory
    standard_fasta_dir = paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert os.path.dirname(file_path) == standard_fasta_dir
    
    # Verify content
    loaded_sequences = read_fasta(file_path)
    assert test_id in loaded_sequences
    assert loaded_sequences[test_id] == "TESTSEQUENCE"
