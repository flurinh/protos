"""
Tests for the FASTA file utilities in protos.io.fasta_utils
"""

import pytest
import os
import tempfile
from protos.io.fasta_utils import (
    read_fasta,
    write_fasta,
    generate_protein_id,
    create_fasta_id_map
)

class TestFastaUtils:
    """Test suite for FASTA file utilities."""
    
    @pytest.fixture
    def sample_fasta_content(self):
        """Sample FASTA content for testing."""
        return (
            ">protein1 description for protein 1\n"
            "MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIMLGFPIN\n"
            "FLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGC\n"
            ">protein2 description for protein 2\n"
            "MNGTEGLNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIVLGFPIN\n"
            "FLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHGYFVFGPTGC\n"
            ">protein3 description for protein 3\n"
            "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIVLGFPIN\n"
            "FLTLYVTVQHKKLRTPLNYILLNLAVANLFMVFGGFTTTLYTSLHGYFVFGPTGC\n"
        )
    
    @pytest.fixture
    def sample_fasta_file(self, sample_fasta_content):
        """Create a temporary FASTA file with sample content."""
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta') as temp:
            temp.write(sample_fasta_content)
            temp_path = temp.name
        
        yield temp_path
        
        # Cleanup
        os.unlink(temp_path)
    
    def test_read_fasta(self, sample_fasta_file):
        """Test reading a FASTA file."""
        # Read the test FASTA file
        sequences = read_fasta(sample_fasta_file)
        
        # Check that we have the expected sequences
        assert len(sequences) == 3
        assert 'protein1' in sequences
        assert 'protein2' in sequences
        assert 'protein3' in sequences
        
        # Check the sequence contents
        assert sequences['protein1'].startswith('MNGTEGPNFYVPFSNKTGVVRSPFEYPQ')
        assert sequences['protein2'].startswith('MNGTEGLNFYVPFSNKTGVVRSPFEYPQ')
        assert sequences['protein3'].startswith('MNGTEGPNFYVPFSNKTGVVRSPFEAPQ')
        
        # Check that there are no newlines in the sequences
        for seq_id, seq in sequences.items():
            assert '\n' not in seq
    
    def test_write_fasta(self):
        """Test writing sequences to a FASTA file."""
        # Create test sequences
        sequences = {
            'protein1': 'MNGTEGPNFYVPFSNKTGVVRSPFEYPQ',
            'protein2': 'MNGTEGLNFYVPFSNKTGVVRSPFEYPQ',
            'protein3': 'MNGTEGPNFYVPFSNKTGVVRSPFEAPQ'
        }
        
        # Create a temporary file
        with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as temp:
            temp_path = temp.name
        
        try:
            # Write sequences to the file
            write_fasta(sequences, temp_path)
            
            # Read back the sequences
            read_sequences = read_fasta(temp_path)
            
            # Check that the sequences match what we wrote
            assert read_sequences == sequences
            
            # Check the file content manually
            with open(temp_path, 'r') as f:
                content = f.read()
            
            # Check that each sequence ID is in the file
            for seq_id in sequences:
                assert f">{seq_id}" in content
            
            # Check that each sequence is in the file
            for seq in sequences.values():
                assert seq in content.replace('\n', '')
                
        finally:
            # Cleanup
            os.unlink(temp_path)
    
    def test_generate_protein_id(self):
        """Test generation of protein IDs from descriptions."""
        # Test standard case
        protein_id = generate_protein_id("test description")
        assert protein_id.startswith('O')
        assert len(protein_id) == 10
        
        # Test with different protein family
        protein_id = generate_protein_id("test description", protein_family='G')
        assert protein_id.startswith('G')
        assert len(protein_id) == 10
        
        # Test with different length
        protein_id = generate_protein_id("test description", descriptor_length=15)
        assert len(protein_id) == 15
        
        # Test that same input gives same output
        assert generate_protein_id("test description") == generate_protein_id("test description")
    
    def test_create_fasta_id_map(self):
        """Test creation of ID mapping dataframe."""
        # Create test dictionary
        fasta_dict = {
            "description1": "ACGTACGT",
            "description2": "GTACGTAC"
        }
        
        # Generate ID mapping
        id_map = create_fasta_id_map(fasta_dict)
        
        # Check that the dataframe has the right structure
        assert "protein_id" in id_map.columns
        assert "description" in id_map.columns
        assert len(id_map) == 2
        
        # Check that the descriptions are preserved
        assert "description1" in id_map["description"].values
        assert "description2" in id_map["description"].values
        
        # Check that the protein IDs follow the expected format (start with 'O')
        for pid in id_map["protein_id"]:
            assert pid.startswith('O')
            assert len(pid) == 10