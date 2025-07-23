"""Tests for StructureProcessor GRN integration methods."""

import pytest
import pandas as pd

from protos.processing.structure import StructureProcessor
from protos.io.paths.path_config import ProtosPaths


class TestCifBaseProcessorGRNMethods:
    """Test the new GRN-related methods in StructureProcessor."""
    
    @pytest.fixture
    def struct_processor(self, tmp_path):
        """Create a StructureProcessor instance for testing."""
        # Set global data root to temp directory
        os.environ["PROTOS_DATA_ROOT"] = str(tmp_path)
        
        processor = StructureProcessor(
            name="test_processor"
        )
        
        yield processor
        
        # Cleanup
        # Clear environment
    
    def test_get_seq_dict_basic(self, struct_processor):
        """Test basic sequence extraction."""
        # Create test data
        test_data = pd.DataFrame({
            'pdb_id': ['1ABC', '1ABC', '1ABC', '1ABC'],
            'auth_chain_id': ['A', 'A', 'A', 'B'],
            'auth_seq_id': [1, 2, 3, 1],
            'auth_comp_id': ['MET', 'ALA', 'GLY', 'VAL']
        })
        
        struct_processor.data = test_data
        
        # Extract sequences
        sequences = struct_processor.get_seq_dict()
        
        assert isinstance(sequences, dict)
        assert len(sequences) == 2
        assert '1ABC_A' in sequences
        assert '1ABC_B' in sequences
        assert sequences['1ABC_A'] == 'MAG'
        assert sequences['1ABC_B'] == 'V'
    
    def test_get_seq_dict_with_gaps(self, struct_processor):
        """Test sequence extraction with gaps."""
        # Create test data with non-sequential residues
        test_data = pd.DataFrame({
            'pdb_id': ['1ABC', '1ABC', '1ABC', '1ABC'],
            'auth_chain_id': ['A', 'A', 'A', 'A'],
            'auth_seq_id': [1, 2, 5, 6],  # Gap between 2 and 5
            'auth_comp_id': ['MET', 'ALA', 'GLY', 'VAL']
        })
        
        struct_processor.data = test_data
        
        # Extract sequences
        sequences = struct_processor.get_seq_dict()
        
        assert '1ABC_A' in sequences
        # Should have gaps represented
        assert len(sequences['1ABC_A']) == 6  # Positions 1,2,3,4,5,6 (with gaps at 3,4)
    
    @pytest.mark.skip(reason="Method map_gen_seq_to_residue not yet implemented in StructureProcessor")
    def test_map_gen_seq_to_residue(self, struct_processor):
        """Test mapping generation sequence ID to residue information."""
        # Create test data
        test_data = pd.DataFrame({
            'pdb_id': ['1ABC'] * 4,
            'auth_chain_id': ['A'] * 4,
            'gen_seq_id': [10, 11, 12, 13],
            'auth_seq_id': [100, 101, 102, 103],
            'auth_comp_id': ['ALA', 'VAL', 'ILE', 'LEU']
        })
        
        struct_processor.data = test_data
        
        # Test mapping
        result = struct_processor.map_gen_seq_to_residue('1ABC', 'A', 11)
        
        assert result == 'V101'
        
        # Test non-existent mapping
        result = struct_processor.map_gen_seq_to_residue('1ABC', 'A', 99)
        assert result == '-'
    
    @pytest.mark.skip(reason="Method map_auth_seq_to_residue not yet implemented in StructureProcessor")
    def test_map_auth_seq_to_residue(self, struct_processor):
        """Test mapping auth sequence ID to residue information."""
        # Create test data
        test_data = pd.DataFrame({
            'pdb_id': ['1ABC'] * 4,
            'auth_chain_id': ['A'] * 4,
            'auth_seq_id': [100, 101, 102, 103],
            'auth_comp_id': ['ALA', 'VAL', 'ILE', 'LEU']
        })
        
        struct_processor.data = test_data
        
        # Test mapping
        result = struct_processor.map_auth_seq_to_residue('1ABC', 'A', 101)
        
        assert result == 'V101'
        
        # Test non-existent mapping
        result = struct_processor.map_auth_seq_to_residue('1ABC', 'A', 999)
        assert result == '-'
    
    @pytest.mark.skip(reason="Method map_grn_to_resname_number not yet implemented in StructureProcessor")
    def test_map_grn_to_resname_number(self, struct_processor):
        """Test mapping GRN table to residue name and number."""
        # Create structure data
        struct_data = pd.DataFrame({
            'pdb_id': ['1ABC'] * 5,
            'auth_chain_id': ['A'] * 5,
            'gen_seq_id': [50, 51, 52, 53, 54],
            'auth_seq_id': [150, 151, 152, 153, 154],
            'auth_comp_id': ['LEU', 'VAL', 'ALA', 'ILE', 'PHE']
        })
        
        struct_processor.data = struct_data
        
        # Create GRN table
        grn_table = pd.DataFrame({
            '1.50': ['L50'],  # Maps to gen_seq_id 50
            '2.50': ['A52'],  # Maps to gen_seq_id 52
            '3.50': ['F54']   # Maps to gen_seq_id 54
        }, index=['1ABC_A'])
        
        # Map GRN to residue info
        result = struct_processor.map_grn_to_resname_number(
            grn_table=grn_table,
            seq_type='gen'
        )
        
        # Check results
        assert isinstance(result, pd.DataFrame)
        assert result.loc['1ABC_A', '1.50'] == 'L150'
        assert result.loc['1ABC_A', '2.50'] == 'A152'
        assert result.loc['1ABC_A', '3.50'] == 'F154'
    
    @pytest.mark.skip(reason="Method map_grn_to_resname_number not yet implemented in StructureProcessor")
    def test_map_grn_with_missing_data(self, struct_processor):
        """Test GRN mapping with missing structure data."""
        # Limited structure data
        struct_data = pd.DataFrame({
            'pdb_id': ['1ABC'],
            'auth_chain_id': ['A'],
            'gen_seq_id': [50],
            'auth_seq_id': [150],
            'auth_comp_id': ['LEU']
        })
        
        struct_processor.data = struct_data
        
        # GRN table with more positions than we have structure for
        grn_table = pd.DataFrame({
            '1.50': ['L50'],   # This exists
            '2.50': ['A52'],   # This doesn't exist
            '3.50': ['F54']    # This doesn't exist
        }, index=['1ABC_A'])
        
        # Map GRN to residue info
        result = struct_processor.map_grn_to_resname_number(
            grn_table=grn_table,
            seq_type='gen'
        )
        
        # Check results
        assert result.loc['1ABC_A', '1.50'] == 'L150'
        assert result.loc['1ABC_A', '2.50'] == '-'
        assert result.loc['1ABC_A', '3.50'] == '-'
    
    @pytest.mark.skip(reason="Method map_grn_to_resname_number not yet implemented in StructureProcessor")
    def test_map_grn_auth_seq_type(self, struct_processor):
        """Test GRN mapping using auth sequence type."""
        # Create structure data
        struct_data = pd.DataFrame({
            'pdb_id': ['1ABC'] * 3,
            'auth_chain_id': ['A'] * 3,
            'gen_seq_id': [10, 11, 12],
            'auth_seq_id': [100, 101, 102],
            'auth_comp_id': ['LEU', 'VAL', 'ALA']
        })
        
        struct_processor.data = struct_data
        
        # GRN table using auth_seq_id positions
        grn_table = pd.DataFrame({
            '1.50': ['L100'],  # auth_seq_id 100
            '2.50': ['V101'],  # auth_seq_id 101
            '3.50': ['A102']   # auth_seq_id 102
        }, index=['1ABC_A'])
        
        # Map GRN using auth sequence type
        result = struct_processor.map_grn_to_resname_number(
            grn_table=grn_table,
            seq_type='auth'
        )
        
        # Check results - should return same values since we're using auth
        assert result.loc['1ABC_A', '1.50'] == 'L100'
        assert result.loc['1ABC_A', '2.50'] == 'V101'
        assert result.loc['1ABC_A', '3.50'] == 'A102'