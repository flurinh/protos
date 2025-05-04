"""
Tests for the FASTA utilities in the protos package.

This module contains tests for reading, writing, and manipulating FASTA files,
including actual downloads from UniProt for testing with real protein sequences.
"""

import os
import tempfile
import pytest
import pandas as pd
from pathlib import Path

from protos.io.fasta_utils import (
    read_fasta, 
    write_fasta, 
    clean_sequence, 
    validate_fasta_format,
    generate_protein_id,
    create_fasta_id_map
)
from protos.loaders.uniprot_utils import get_uniprot
from protos.loaders.uniprot_loader import UniprotDL


@pytest.fixture
def sample_fasta_content():
    """Create sample FASTA content for testing."""
    return """>sp|P00533|EGFR_HUMAN Epidermal growth factor receptor
MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLTQLGTFEDHFLSLQRMFNNCEVVLGNLEITYVQRNYDL
SFLKTIQEVAGYVLIALNTVERIPLENLQIIRGNMYYENSYALAVLSNYDANKTGLKELPMRNLQEILHGAVRFS
NNPALCNVESIQWRDIVSSDFLSNMSMDFQNHLGSCQKCDPSCPNGSCWGAGEENCQKLTKIICAQQCSGRCRGK
SPSDCCHNQCAAGCTGPRESDCLVCRKFRDEATCKDTCPPLMLYNPTTYQMDVNPEGKYSFGATCVKKCPRNYV
>sp|P01308|INS_HUMAN Insulin
MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGA
GSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"""


@pytest.fixture
def temp_fasta_file(sample_fasta_content):
    """Create a temporary FASTA file for testing."""
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp_file:
        temp_file.write(sample_fasta_content)
        temp_path = temp_file.name
    
    yield temp_path
    
    # Clean up the temporary file
    if os.path.exists(temp_path):
        os.remove(temp_path)


@pytest.fixture
def test_uniprot_ids():
    """Define test UniProt IDs to use for real data tests."""
    # Return list of UniProt IDs for proteins with diverse features
    return ["P00533", "P01308", "P05067", "P02751"]


def test_read_fasta(temp_fasta_file):
    """Test reading FASTA files."""
    sequences = read_fasta(temp_fasta_file)
    
    # Verify the correct sequences were read
    assert len(sequences) == 2
    assert "sp|P00533|EGFR_HUMAN" in sequences
    assert "sp|P01308|INS_HUMAN" in sequences
    
    # Verify the sequences are correct
    assert sequences["sp|P00533|EGFR_HUMAN"].startswith("MRPSGTAGAALLALLAALCPASR")
    assert sequences["sp|P01308|INS_HUMAN"].startswith("MALWMRLLPLLALLALWGPDPAA")


def test_write_fasta(temp_fasta_file):
    """Test writing FASTA files."""
    # Read the original sequences
    original_sequences = read_fasta(temp_fasta_file)
    
    # Create a new temporary file for writing
    output_path = temp_fasta_file + ".out"
    
    # Write the sequences to the new file
    write_fasta(original_sequences, output_path)
    
    # Read back the written sequences
    written_sequences = read_fasta(output_path)
    
    # Verify the sequences match
    assert set(original_sequences.keys()) == set(written_sequences.keys())
    for key in original_sequences:
        assert original_sequences[key] == written_sequences[key]
    
    # Clean up
    os.remove(output_path)


def test_clean_sequence():
    """Test cleaning of sequence data."""
    # Test with whitespace and lowercase
    test_seq = "  ACGT ACGT  acgt acgt  "
    cleaned = clean_sequence(test_seq)
    assert cleaned == "ACGTACGTACGTACGT"
    
    # Test with invalid characters
    test_seq = "ACGT1234*-ACGT"
    cleaned = clean_sequence(test_seq)
    assert cleaned == "ACGTACGT"
    
    # Test with only invalid characters
    test_seq = "12345*-"
    cleaned = clean_sequence(test_seq)
    assert cleaned == ""


def test_validate_fasta_format():
    """Test validation of FASTA format."""
    # Valid FASTA
    valid_fasta = ">test_seq\nACGTACGT"
    assert validate_fasta_format(valid_fasta) is True
    
    # Invalid FASTA - no header
    invalid_fasta1 = "ACGTACGT"
    assert validate_fasta_format(invalid_fasta1) is False
    
    # Invalid FASTA - empty sequence
    invalid_fasta2 = ">test_seq\n"
    assert validate_fasta_format(invalid_fasta2) is False
    
    # Invalid FASTA - empty content
    invalid_fasta3 = ""
    assert validate_fasta_format(invalid_fasta3) is False


def test_generate_protein_id():
    """Test generating protein IDs from descriptors."""
    descriptor1 = "EGFR_HUMAN"
    descriptor2 = "INS_HUMAN"
    
    # Test with default parameters
    id1 = generate_protein_id(descriptor1)
    id2 = generate_protein_id(descriptor2)
    
    # Verify the IDs are correctly formatted
    assert id1.startswith("O")
    assert id2.startswith("O")
    assert len(id1) == 10
    assert len(id2) == 10
    assert id1 != id2
    
    # Test with custom parameters
    id1_custom = generate_protein_id(descriptor1, protein_family="P", descriptor_length=8)
    assert id1_custom.startswith("P")
    assert len(id1_custom) == 8


def test_create_fasta_id_map():
    """Test creating a mapping of FASTA IDs."""
    # Create a sample FASTA dictionary
    sample_dict = {
        "sp|P00533|EGFR_HUMAN": "MRPSGTAGAALLALLAALCPASR",
        "sp|P01308|INS_HUMAN": "MALWMRLLPLLALLALWGPDPAA"
    }
    
    # Create the ID map
    id_map = create_fasta_id_map(sample_dict)
    
    # Verify the map format
    assert isinstance(id_map, pd.DataFrame)
    assert "protein_id" in id_map.columns
    assert "description" in id_map.columns
    assert len(id_map) == 2
    
    # Verify the map content
    descriptions = list(id_map["description"])
    assert "sp|P00533|EGFR_HUMAN" in descriptions
    assert "sp|P01308|INS_HUMAN" in descriptions


@pytest.mark.skipif(not os.environ.get("RUN_NETWORK_TESTS"), 
                   reason="Network-dependent tests are disabled")
def test_uniprot_download(test_uniprot_ids):
    """Test downloading sequences from UniProt."""
    # Only test the first ID to avoid too many requests
    test_id = test_uniprot_ids[0]  # P00533 - EGFR
    
    # Download the sequence
    uniprot_data = get_uniprot(test_id, reviewed=True)
    
    # Verify we got data back
    assert not uniprot_data.empty
    assert "Sequence" in uniprot_data.columns
    
    # Extract the protein sequence
    sequence = uniprot_data["Sequence"].values[0]
    
    # Verify the sequence is not empty and looks like a protein sequence
    assert len(sequence) > 0
    assert all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence[:20])  # Check first 20 amino acids


@pytest.mark.skipif(not os.environ.get("RUN_NETWORK_TESTS"), 
                   reason="Network-dependent tests are disabled")
def test_uniprot_loader(test_uniprot_ids, tmp_path):
    """Test the UniprotDL class."""
    # Create a temporary directory structure for the test
    tmp_uniprot_dir = tmp_path / "uniprot"
    tmp_uniprot_dir.mkdir()
    (tmp_uniprot_dir / "datasets").mkdir()
    (tmp_uniprot_dir / "table").mkdir()
    (tmp_uniprot_dir / "fastas").mkdir()
    
    # Create a file with test IDs
    dataset_name = "test_dataset"
    dataset_file = tmp_uniprot_dir / "datasets" / f"{dataset_name}.txt"
    with open(dataset_file, 'w') as f:
        f.write(" ".join(test_uniprot_ids[:2]))  # Use only 2 IDs for speed
    
    # Initialize the loader
    loader = UniprotDL(
        path=str(tmp_uniprot_dir) + "/",
        dataset=dataset_name,
        limit=2  # Limit to 2 sequences for speed
    )
    
    # Load dataset IDs
    loaded_ids = loader.load_dataset()
    assert len(loaded_ids) == 2
    assert loaded_ids[0] in test_uniprot_ids
    assert loaded_ids[1] in test_uniprot_ids
    
    # Download the sequences (this might take a bit longer)
    # Skip actual download to avoid slow tests
    # loader.download_genes_single_query(batchsize=2, save=True)
    # assert len(loader.gene_df) == 2


def test_integration_read_process_write(temp_fasta_file):
    """Integration test for reading, processing, and writing FASTA files."""
    # Read the original sequences
    original_sequences = read_fasta(temp_fasta_file)
    
    # Create a processed version with cleaned sequences
    processed_sequences = {
        seq_id: clean_sequence(sequence)
        for seq_id, sequence in original_sequences.items()
    }
    
    # Write the processed sequences to a new file
    output_path = temp_fasta_file + ".processed"
    write_fasta(processed_sequences, output_path)
    
    # Read back the processed sequences
    processed_sequences_read = read_fasta(output_path)
    
    # Verify the processing worked correctly
    assert set(processed_sequences.keys()) == set(processed_sequences_read.keys())
    for key in processed_sequences:
        assert processed_sequences[key] == processed_sequences_read[key]
    
    # Clean up
    os.remove(output_path)