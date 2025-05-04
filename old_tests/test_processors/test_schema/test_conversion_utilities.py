"""
Tests for the conversion_utilities module.

These old_tests verify that the conversion utility functions correctly transform data
between different formats and representations.
"""

import pytest
import pandas as pd
import numpy as np
import re

from protos.processing.schema.conversion_utilities import (
    # Structure conversions
    cif_to_structure_df,
    structure_df_to_cif,
    dict_to_structure_df,
    structure_df_to_dict,
    three_to_one_letter_code,
    one_to_three_letter_code,
    
    # GRN conversions
    grn_str_to_float,
    grn_float_to_str,
    sort_grns,
    grn_mapping_to_df,
    df_to_grn_mapping,
    
    # Sequence conversions
    fasta_to_dict,
    dict_to_fasta,
    alignment_result_to_df,
    df_to_alignment_result
)

class TestStructureConversions:
    """Tests for structure data conversion functions."""
    
    def test_cif_to_structure_df(self):
        """Test converting CIF DataFrame to standard structure format."""
        # Create a sample CIF DataFrame
        cif_df = pd.DataFrame({
            'group_PDB': ['ATOM', 'ATOM', 'ATOM'],
            'label_asym_id': ['A', 'A', 'A'],
            'auth_asym_id': ['A', 'A', 'A'],
            'label_seq_id': [1, 1, 2],
            'auth_seq_id': [1, 1, 2],
            'label_comp_id': ['ALA', 'ALA', 'LEU'],
            'label_atom_id': ['CA', 'CB', 'CA'],
            'Cartn_x': [10.0, 10.5, 14.0],
            'Cartn_y': [20.0, 20.5, 24.0],
            'Cartn_z': [30.0, 30.5, 34.0],
            'id': [1, 2, 3]
        })
        
        # Convert to structure DataFrame
        structure_df = cif_to_structure_df(cif_df, structure_id='1abc')
        
        # Verify conversion
        assert 'pdb_id' in structure_df.columns
        assert 'group' in structure_df.columns
        assert 'auth_chain_id' in structure_df.columns
        assert 'gen_chain_id' in structure_df.columns
        assert 'auth_seq_id' in structure_df.columns
        assert 'gen_seq_id' in structure_df.columns
        assert 'atom_name' in structure_df.columns
        assert 'x' in structure_df.columns
        assert 'y' in structure_df.columns
        assert 'z' in structure_df.columns
        
        # Check data conversion
        assert structure_df['pdb_id'].iloc[0] == '1abc'
        assert structure_df['group'].iloc[0] == 'ATOM'
        assert structure_df['auth_chain_id'].iloc[0] == 'A'
        assert structure_df['gen_chain_id'].iloc[0] == 'A'
        assert structure_df['auth_seq_id'].iloc[0] == 1
        assert structure_df['gen_seq_id'].iloc[0] == 1
        assert structure_df['atom_name'].iloc[0] == 'CA'
        assert structure_df['res_name3l'].iloc[0] == 'ALA'
        assert structure_df['res_name1l'].iloc[0] == 'A'
        assert structure_df['x'].iloc[0] == 10.0
        assert structure_df['y'].iloc[0] == 20.0
        assert structure_df['z'].iloc[0] == 30.0
    
    def test_structure_df_to_cif(self):
        """Test converting standard structure format to CIF DataFrame."""
        # Create a sample structure DataFrame
        structure_df = pd.DataFrame({
            'pdb_id': ['1abc', '1abc', '1abc'],
            'group': ['ATOM', 'ATOM', 'ATOM'],
            'gen_chain_id': ['A', 'A', 'A'],
            'auth_chain_id': ['A', 'A', 'A'],
            'gen_seq_id': [1, 1, 2],
            'auth_seq_id': [1, 1, 2],
            'res_name3l': ['ALA', 'ALA', 'LEU'],
            'atom_name': ['CA', 'CB', 'CA'],
            'x': [10.0, 10.5, 14.0],
            'y': [20.0, 20.5, 24.0],
            'z': [30.0, 30.5, 34.0],
            'atom_id': [1, 2, 3],
            'res_name1l': ['A', 'A', 'L'],
            'res_atom_name': ['ALA.CA', 'ALA.CB', 'LEU.CA']
        })
        
        # Convert to CIF DataFrame
        cif_df = structure_df_to_cif(structure_df)
        
        # Verify conversion
        assert 'group_PDB' in cif_df.columns
        assert 'label_asym_id' in cif_df.columns
        assert 'auth_asym_id' in cif_df.columns
        assert 'label_seq_id' in cif_df.columns
        assert 'auth_seq_id' in cif_df.columns
        assert 'label_comp_id' in cif_df.columns
        assert 'label_atom_id' in cif_df.columns
        assert 'Cartn_x' in cif_df.columns
        assert 'Cartn_y' in cif_df.columns
        assert 'Cartn_z' in cif_df.columns
        assert 'id' in cif_df.columns
        assert 'type_symbol' in cif_df.columns
        
        # Check data conversion
        assert cif_df['group_PDB'].iloc[0] == 'ATOM'
        assert cif_df['label_asym_id'].iloc[0] == 'A'
        assert cif_df['auth_asym_id'].iloc[0] == 'A'
        assert cif_df['label_seq_id'].iloc[0] == 1
        assert cif_df['auth_seq_id'].iloc[0] == 1
        assert cif_df['label_comp_id'].iloc[0] == 'ALA'
        assert cif_df['label_atom_id'].iloc[0] == 'CA'
        assert cif_df['Cartn_x'].iloc[0] == 10.0
        assert cif_df['Cartn_y'].iloc[0] == 20.0
        assert cif_df['Cartn_z'].iloc[0] == 30.0
        assert cif_df['id'].iloc[0] == 1
        assert cif_df['type_symbol'].iloc[0] == 'C'  # First letter of 'CA'
    
    def test_dict_to_structure_df(self):
        """Test converting dictionary to structure DataFrame."""
        # Create a sample structure dictionary
        structure_dict = {
            'pdb_id': ['1abc', '1abc', '1abc'],
            'group': ['ATOM', 'ATOM', 'ATOM'],
            'gen_chain_id': ['A', 'A', 'A'],
            'auth_chain_id': ['A', 'A', 'A'],
            'gen_seq_id': [1, 1, 2],
            'auth_seq_id': [1, 1, 2],
            'res_name3l': ['ALA', 'ALA', 'LEU'],
            'atom_name': ['CA', 'CB', 'CA'],
            'x': [10.0, 10.5, 14.0],
            'y': [20.0, 20.5, 24.0],
            'z': [30.0, 30.5, 34.0],
            'atom_id': [1, 2, 3],
            'res_name1l': ['A', 'A', 'L'],
            'res_atom_name': ['ALA.CA', 'ALA.CB', 'LEU.CA']
        }
        
        # Convert to structure DataFrame
        structure_df = dict_to_structure_df(structure_dict)
        
        # Verify conversion
        assert isinstance(structure_df, pd.DataFrame)
        assert len(structure_df) == 3
        
        # Check columns
        for col in structure_dict:
            assert col in structure_df.columns
            
        # Check data
        assert structure_df['pdb_id'].iloc[0] == '1abc'
        assert structure_df['atom_name'].iloc[0] == 'CA'
        assert structure_df['x'].iloc[0] == 10.0
    
    def test_structure_df_to_dict(self):
        """Test converting structure DataFrame to dictionary."""
        # Create a sample structure DataFrame
        structure_df = pd.DataFrame({
            'pdb_id': ['1abc', '1abc', '1abc'],
            'group': ['ATOM', 'ATOM', 'ATOM'],
            'gen_chain_id': ['A', 'A', 'A'],
            'auth_chain_id': ['A', 'A', 'A'],
            'gen_seq_id': [1, 1, 2],
            'auth_seq_id': [1, 1, 2],
            'res_name3l': ['ALA', 'ALA', 'LEU'],
            'atom_name': ['CA', 'CB', 'CA'],
            'x': [10.0, 10.5, 14.0],
            'y': [20.0, 20.5, 24.0],
            'z': [30.0, 30.5, 34.0],
            'atom_id': [1, 2, 3],
            'res_name1l': ['A', 'A', 'L'],
            'res_atom_name': ['ALA.CA', 'ALA.CB', 'LEU.CA']
        })
        
        # Convert to dictionary
        structure_dict = structure_df_to_dict(structure_df)
        
        # Verify conversion
        assert isinstance(structure_dict, dict)
        
        # Check keys and values
        for col in structure_df.columns:
            assert col in structure_dict
            assert len(structure_dict[col]) == 3
            
        # Check data
        assert structure_dict['pdb_id'][0] == '1abc'
        assert structure_dict['atom_name'][0] == 'CA'
        assert structure_dict['x'][0] == 10.0
    
    def test_three_to_one_letter_code(self):
        """Test converting three-letter amino acid codes to one-letter codes."""
        # Test common amino acids
        assert three_to_one_letter_code('ALA') == 'A'
        assert three_to_one_letter_code('CYS') == 'C'
        assert three_to_one_letter_code('ASP') == 'D'
        assert three_to_one_letter_code('GLU') == 'E'
        assert three_to_one_letter_code('PHE') == 'F'
        
        # Test case insensitivity
        assert three_to_one_letter_code('ala') == 'A'
        assert three_to_one_letter_code('Cys') == 'C'
        
        # Test special cases
        assert three_to_one_letter_code('UNK') == 'X'
        assert three_to_one_letter_code('') == '-'
        
        # Test invalid codes
        assert three_to_one_letter_code('ZZZ') == 'X'
    
    def test_one_to_three_letter_code(self):
        """Test converting one-letter amino acid codes to three-letter codes."""
        # Test common amino acids
        assert one_to_three_letter_code('A') == 'ALA'
        assert one_to_three_letter_code('C') == 'CYS'
        assert one_to_three_letter_code('D') == 'ASP'
        assert one_to_three_letter_code('E') == 'GLU'
        assert one_to_three_letter_code('F') == 'PHE'
        
        # Test case insensitivity
        assert one_to_three_letter_code('a') == 'ALA'
        assert one_to_three_letter_code('c') == 'CYS'
        
        # Test special cases
        assert one_to_three_letter_code('X') == 'UNK'
        assert one_to_three_letter_code('-') == ''
        
        # Test invalid codes
        assert one_to_three_letter_code('Z') == 'GLX'
        assert one_to_three_letter_code('B') == 'ASX'
        assert one_to_three_letter_code('8') == 'UNK'


class TestGRNConversions:
    """Tests for GRN data conversion functions."""
    
    def test_grn_str_to_float(self):
        """Test converting GRN strings to float representations."""
        # Test standard GRN format (e.g., 1x50 -> 1.5)
        assert grn_str_to_float('1x50') == 1.5
        assert grn_str_to_float('2x40') == 2.4
        assert grn_str_to_float('7x65') == 7.65
        
        # Test N-terminal format (e.g., n.10 -> -0.1)
        assert grn_str_to_float('n.10') == -0.1
        assert grn_str_to_float('n.5') == -0.05
        
        # Test C-terminal format (e.g., c.5 -> 8.05)
        assert grn_str_to_float('c.5') == 8.05
        assert grn_str_to_float('c.10') == 8.1
        
        # Test loop format (e.g., 2.45 -> 2.045)
        assert grn_str_to_float('2.45') == 2.045
        assert grn_str_to_float('3.60') == 3.06
        
        # Test invalid formats
        assert grn_str_to_float('invalid') == 0.0
    
    def test_grn_float_to_str(self):
        """Test converting GRN float representations to strings."""
        # Test standard GRN format
        assert grn_float_to_str(1.5) == '1x50'
        assert grn_float_to_str(2.4) == '2x40'
        assert grn_float_to_str(7.65) == '7x65'
        
        # Test N-terminal format
        assert grn_float_to_str(-0.1) == 'n.10'
        assert grn_float_to_str(-0.05) == 'n.5'
        
        # Test C-terminal format
        assert grn_float_to_str(8.05) == 'c.5'
        assert grn_float_to_str(8.1) == 'c.10'
        
        # Test loop format
        assert grn_float_to_str(2.045) == '2.45'
        assert grn_float_to_str(3.06) == '3.60'
    
    def test_sort_grns(self):
        """Test sorting GRNs in standard order."""
        # Create an unsorted list of GRNs
        unsorted_grns = ['3x50', '1x50', '7x50', 'n.10', 'c.5', '2.45']
        
        # Sort the GRNs
        sorted_grns = sort_grns(unsorted_grns)
        
        # Verify sorting
        assert sorted_grns[0] == 'n.10'     # N-terminal should be first
        assert sorted_grns[1] == '1x50'     # Then standard GRNs in order
        assert sorted_grns[2] == '2.45'     # Loop notation should be between helix numbers
        assert sorted_grns[3] == '3x50'
        assert sorted_grns[4] == '7x50'
        assert sorted_grns[5] == 'c.5'      # C-terminal should be last
    
    def test_grn_mapping_to_df(self):
        """Test converting GRN mapping dictionary to DataFrame."""
        # Create a sample GRN mapping
        grn_mapping = {
            1: '1x50',
            2: '1x51',
            5: '2x40',
            10: '3x50'
        }
        
        # Provide a sequence for residue information
        sequence = 'ACDEFGHIKLM'  # 11 amino acids (positions 0-10)
        
        # Convert to DataFrame
        df = grn_mapping_to_df(grn_mapping, sequence)
        
        # Verify conversion
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 4  # 4 GRN mappings
        
        # Check columns
        assert 'position' in df.columns
        assert 'grn' in df.columns
        assert 'residue' in df.columns
        
        # Check data
        assert df['position'].tolist() == [1, 2, 5, 10]
        assert df['grn'].tolist() == ['1x50', '1x51', '2x40', '3x50']
        assert df['residue'].tolist() == ['A', 'C', 'F', 'L']
        
        # Test without sequence
        df_no_seq = grn_mapping_to_df(grn_mapping)
        assert 'residue' not in df_no_seq.columns
    
    def test_df_to_grn_mapping(self):
        """Test converting DataFrame with GRN mapping to dictionary."""
        # Create a sample DataFrame
        df = pd.DataFrame({
            'position': [1, 2, 5, 10],
            'grn': ['1x50', '1x51', '2x40', '3x50'],
            'residue': ['A', 'C', 'F', 'L']
        })
        
        # Convert to GRN mapping dictionary
        grn_mapping = df_to_grn_mapping(df)
        
        # Verify conversion
        assert isinstance(grn_mapping, dict)
        assert len(grn_mapping) == 4
        
        # Check data
        assert grn_mapping[1] == '1x50'
        assert grn_mapping[2] == '1x51'
        assert grn_mapping[5] == '2x40'
        assert grn_mapping[10] == '3x50'
        
        # Test with missing columns
        df_missing = pd.DataFrame({
            'position': [1, 2, 5, 10]
            # Missing 'grn' column
        })
        
        with pytest.raises(ValueError):
            df_to_grn_mapping(df_missing)


class TestSequenceConversions:
    """Tests for sequence data conversion functions."""
    
    def test_fasta_to_dict(self):
        """Test converting FASTA content to dictionary."""
        # Create sample FASTA content
        fasta_content = ">seq1 Description of sequence 1\n"
        fasta_content += "ACDEFGHIKLMNPQRSTVWY\n"
        fasta_content += ">seq2 Description of sequence 2\n"
        fasta_content += "ACDEFGHIKLM\n"
        fasta_content += "NPQRSTVWY\n"  # Line break in sequence
        
        # Convert to dictionary
        sequences = fasta_to_dict(fasta_content)
        
        # Verify conversion
        assert isinstance(sequences, dict)
        assert len(sequences) == 2
        
        # Check sequence IDs
        assert 'seq1' in sequences
        assert 'seq2' in sequences
        
        # Check sequence content
        assert sequences['seq1'] == 'ACDEFGHIKLMNPQRSTVWY'
        assert sequences['seq2'] == 'ACDEFGHIKLMNPQRSTVWY'  # Line break should be removed
    
    def test_dict_to_fasta(self):
        """Test converting dictionary to FASTA content."""
        # Create sample sequence dictionary
        sequences = {
            'seq1': 'ACDEFGHIKLMNPQRSTVWY',
            'seq2': 'ACDEFGHIKLMNPQRSTVWY'
        }
        
        # Convert to FASTA content
        fasta_content = dict_to_fasta(sequences, width=10)
        
        # Verify conversion
        assert isinstance(fasta_content, str)
        assert '>seq1' in fasta_content
        assert '>seq2' in fasta_content
        
        # Check sequence formatting
        assert 'ACDEFGHIKL' in fasta_content  # First line of seq1 (width=10)
        assert 'MNPQRSTVWY' in fasta_content  # Second line of seq1 (width=10)
        
        # Count lines
        lines = fasta_content.split('\n')
        assert len(lines) == 6  # 2 sequences * (1 header + 2 sequence lines) = 6
    
    def test_alignment_result_to_df(self):
        """Test converting alignment results to DataFrame."""
        # Create sample alignment results
        query_id = 'seq1'
        target_id = 'seq2'
        aligned_query = 'ACDEFG-HI'
        aligned_target = 'ACDE-GYHI'
        score = 0.8
        
        # Convert to DataFrame
        df = alignment_result_to_df(query_id, target_id, aligned_query, aligned_target, score)
        
        # Verify conversion
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        
        # Check columns
        assert 'query_id' in df.columns
        assert 'target_id' in df.columns
        assert 'sequence_identity' in df.columns
        assert 'alignment_length' in df.columns
        assert 'mismatches' in df.columns
        assert 'gap_openings' in df.columns
        assert 'query_start' in df.columns
        assert 'query_end' in df.columns
        assert 'target_start' in df.columns
        assert 'target_end' in df.columns
        assert 'e_value' in df.columns
        assert 'bit_score' in df.columns
        
        # Check data
        assert df['query_id'].iloc[0] == 'seq1'
        assert df['target_id'].iloc[0] == 'seq2'
        assert df['alignment_length'].iloc[0] == 9
        assert df['bit_score'].iloc[0] == 0.8
        
        # Check calculated values
        assert df['sequence_identity'].iloc[0] > 0  # Should have some identity
        assert df['gap_openings'].iloc[0] == 2  # Two gaps in the alignment
    
    def test_df_to_alignment_result(self):
        """Test converting alignment DataFrame to dictionary."""
        # Create sample alignment DataFrame
        df = pd.DataFrame({
            'query_id': ['seq1'],
            'target_id': ['seq2'],
            'sequence_identity': [77.78],  # 7/9 matches
            'alignment_length': [9],
            'mismatches': [0],
            'gap_openings': [2],
            'query_start': [1],
            'query_end': [9],
            'target_start': [1],
            'target_end': [9],
            'e_value': [1e-10],
            'bit_score': [0.8]
        })
        
        # Convert to alignment result dictionary
        result = df_to_alignment_result(df)
        
        # Verify conversion
        assert isinstance(result, dict)
        
        # Check fields
        assert result['query_id'] == 'seq1'
        assert result['target_id'] == 'seq2'
        assert result['sequence_identity'] == 77.78
        assert result['alignment_length'] == 9
        assert result['query_start'] == 1
        assert result['query_end'] == 9
        assert result['target_start'] == 1
        assert result['target_end'] == 9
        assert result['score'] == 0.8
        
        # Test with empty DataFrame
        empty_df = pd.DataFrame()
        empty_result = df_to_alignment_result(empty_df)
        assert empty_result == {}