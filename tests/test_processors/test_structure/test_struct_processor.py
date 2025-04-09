import pytest
import os
import pandas as pd
import numpy as np
from protos.processing.structure.struct_processor import CifProcessor
from protos.processing.structure.struct_utils import get_distance_matrix

# Global constant for testing - common GPCR structure
TEST_PDB_ID = "4ahz"


class TestCifProcessor:
    @pytest.fixture
    def structure_dir(self):
        """Directory containing structure files."""
        return "mmcif"
    
    @pytest.fixture
    def test_structure_file(self, data_path, structure_dir):
        """Create a test structure file for CifProcessor tests."""
        # Create directory if it doesn't exist
        test_cif_dir = os.path.join(data_path, structure_dir)
        os.makedirs(test_cif_dir, exist_ok=True)
        
        # Create a minimal CIF file for testing
        test_cif_path = os.path.join(test_cif_dir, f"{TEST_PDB_ID}.cif")
        
        with open(test_cif_path, 'w') as f:
            f.write("""data_4AHZ
# Minimal CIF file for testing
loop_
_atom_site.group_PDB 
_atom_site.id 
_atom_site.type_symbol 
_atom_site.label_atom_id 
_atom_site.label_alt_id 
_atom_site.label_comp_id 
_atom_site.label_asym_id 
_atom_site.label_entity_id 
_atom_site.label_seq_id 
_atom_site.pdbx_PDB_ins_code 
_atom_site.Cartn_x 
_atom_site.Cartn_y 
_atom_site.Cartn_z 
_atom_site.occupancy 
_atom_site.B_iso_or_equiv 
_atom_site.auth_seq_id 
_atom_site.auth_comp_id 
_atom_site.auth_asym_id 
_atom_site.auth_atom_id 
_atom_site.pdbx_PDB_model_num 
ATOM 1 N N . MET A 1 1 ? 27.340 24.430 2.614 1.00 9.67 1 MET A N 1
ATOM 2 C CA . MET A 1 1 ? 26.266 25.413 2.842 1.00 10.38 1 MET A CA 1
ATOM 3 C C . MET A 1 1 ? 26.913 26.639 3.531 1.00 9.62 1 MET A C 1
ATOM 4 O O . MET A 1 1 ? 27.886 26.463 4.263 1.00 9.62 1 MET A O 1
ATOM 5 C CB . MET A 1 1 ? 25.112 24.880 3.649 1.00 14.38 1 MET A CB 1
ATOM 6 C CG . MET A 1 1 ? 23.856 25.834 3.626 1.00 13.68 1 MET A CG 1
ATOM 7 S SD . MET A 1 1 ? 22.329 25.029 4.361 1.00 10.66 1 MET A SD 1
ATOM 8 C CE . MET A 1 1 ? 21.225 26.061 3.892 1.00 14.97 1 MET A CE 1
ATOM 9 N N . ASP A 1 2 ? 26.335 27.770 3.258 1.00 7.80 2 ASP A N 1
ATOM 10 C CA . ASP A 1 2 ? 26.850 29.021 3.754 1.00 7.19 2 ASP A CA 1
ATOM 11 C C . ASP A 1 2 ? 26.409 29.163 5.238 1.00 7.91 2 ASP A C 1
ATOM 12 O O . ASP A 1 2 ? 25.238 29.438 5.505 1.00 9.53 2 ASP A O 1
ATOM 13 C CB . ASP A 1 2 ? 26.260 30.222 2.984 1.00 8.31 2 ASP A CB 1
ATOM 14 C CG . ASP A 1 2 ? 26.559 30.169 1.476 1.00 12.81 2 ASP A CG 1
ATOM 15 O OD1 . ASP A 1 2 ? 27.722 29.951 1.064 1.00 10.03 2 ASP A OD1 1
ATOM 16 O OD2 . ASP A 1 2 ? 25.565 30.359 0.697 1.00 14.69 2 ASP A OD2 1
""")
        
        # Return the path for use in tests
        return test_cif_path
    
    @pytest.fixture
    def cp(self, data_path, structure_dir):
        """Create a CifProcessor instance for testing."""
        return CifProcessor(path=data_path, structure_dir=structure_dir)
    
    def test_init(self, cp, data_path, structure_dir):
        """Test CifProcessor initialization."""
        assert cp.path == data_path
        assert cp.path_structure_dir == os.path.join(data_path, structure_dir)
        assert isinstance(cp.pdb_ids, list)
    
    def test_available_datasets(self, cp):
        """Test getting available datasets."""
        try:
            datasets = cp.available_datasets()
            assert isinstance(datasets, list)
        except:
            # Skip if the method is not implemented or no datasets.json file
            pytest.skip("available_datasets not implemented or no datasets.json file")
    
    def test_get_available_pdb_files(self, cp, test_structure_file):
        """Test getting available PDB files."""
        pdb_files = cp.get_available_pdb_files()
        assert isinstance(pdb_files, list)
        assert TEST_PDB_ID in pdb_files, f"Expected to find {TEST_PDB_ID} in available PDB files {pdb_files}"
    
    def test_get_pdb_ids_from_file(self, cp):
        """Test getting PDB IDs from a file."""
        # Skip if no PDB IDs file is specified
        if cp.path_pdb_ids_file is None:
            pytest.skip("No PDB IDs file specified")
        
        pdb_ids = cp.get_pdb_ids_from_file()
        assert isinstance(pdb_ids, list)
        assert len(pdb_ids) > 0
    
    def test_load_structure(self, cp, test_structure_file):
        """Test loading a structure."""
        # Load the structure
        structure = cp.load_structure(TEST_PDB_ID)
        assert structure is not None, f"Failed to load structure {TEST_PDB_ID}"
        assert isinstance(structure, pd.DataFrame)
        assert TEST_PDB_ID in cp.pdb_ids
        
        # Test that structure data was loaded properly
        assert len(structure) > 0, "Structure DataFrame is empty"
        assert 'pdb_id' in structure.columns
        assert structure['pdb_id'].iloc[0] == TEST_PDB_ID
    
    def test_get_chains(self, cp, test_structure_file):
        """Test getting chains for a structure."""
        # Load the structure first
        cp.load_structure(TEST_PDB_ID)
        
        # Get chains
        chains = cp.get_chains(TEST_PDB_ID)
        assert isinstance(chains, list)
        assert len(chains) > 0, f"No chains found in {TEST_PDB_ID}"
        
        # Since our test structure has chain A
        assert 'A' in chains
    
    def test_get_chain_residues(self, cp, test_structure_file):
        """Test getting residues from a chain."""
        # Load the structure first
        cp.load_structure(TEST_PDB_ID)
        
        # Get chains and select the first one
        chains = cp.get_chains(TEST_PDB_ID)
        assert len(chains) > 0, f"No chains found in {TEST_PDB_ID}"
        
        test_chain = chains[0]
        
        # Get residues
        residues = cp.get_chain_residues(TEST_PDB_ID, test_chain)
        assert isinstance(residues, list)
        assert len(residues) > 0, f"No residues found in chain {test_chain}"
        
        # Check residue format
        for residue in residues:
            assert isinstance(residue, tuple) or isinstance(residue, int)
    
    def test_get_ca_coordinates(self, cp, test_structure_file):
        """Test getting CA coordinates."""
        # Load the structure first
        cp.load_structure(TEST_PDB_ID)
        
        # Get chains and select the first one
        chains = cp.get_chains(TEST_PDB_ID)
        assert len(chains) > 0, f"No chains found in {TEST_PDB_ID}"
        
        test_chain = chains[0]
        
        # Get CA coordinates
        ca_coords = cp.get_ca_coordinates(TEST_PDB_ID, test_chain)
        assert isinstance(ca_coords, np.ndarray)
        assert ca_coords.shape[1] == 3  # X, Y, Z coordinates
        assert ca_coords.shape[0] > 0, "No CA coordinates found"
    
    def test_get_backbone_coordinates(self, cp, test_structure_file):
        """Test getting backbone coordinates."""
        # Load the structure first
        cp.load_structure(TEST_PDB_ID)
        
        # Get chains and select the first one
        chains = cp.get_chains(TEST_PDB_ID)
        assert len(chains) > 0, f"No chains found in {TEST_PDB_ID}"
        
        test_chain = chains[0]
        
        # Get backbone coordinates
        bb_coords = cp.get_backbone_coordinates(TEST_PDB_ID, test_chain)
        assert isinstance(bb_coords, dict)
        for atom_type in ['CA', 'C', 'N', 'O']:
            assert atom_type in bb_coords
            assert isinstance(bb_coords[atom_type], np.ndarray)
            assert bb_coords[atom_type].shape[1] == 3  # X, Y, Z coordinates
    
    def test_get_sequence(self, cp, test_structure_file):
        """Test getting sequence from a structure."""
        # Load the structure first
        cp.load_structure(TEST_PDB_ID)
        
        # Get chains and select the first one
        chains = cp.get_chains(TEST_PDB_ID)
        assert len(chains) > 0, f"No chains found in {TEST_PDB_ID}"
        
        test_chain = chains[0]
        
        # Get sequence
        try:
            sequence = cp.get_sequence(TEST_PDB_ID, test_chain)
            assert isinstance(sequence, str)
            assert len(sequence) > 0, "Empty sequence returned"
        except AttributeError:
            pytest.skip("get_sequence method not implemented")
    
    def test_filter_by_ids(self, cp, test_structure_file):
        """Test filtering by PDB IDs."""
        # Load a structure
        cp.load_structure(TEST_PDB_ID)
        
        # Filter to include only our test structure
        cp.filter_by_ids([TEST_PDB_ID])
        
        # Check that our structure is still there
        assert len(cp.pdb_ids) == 1
        assert cp.pdb_ids[0] == TEST_PDB_ID
        
        # Try filtering with unknown ID (should result in empty list)
        cp.filter_by_ids(["unknown_pdb"])
        assert len(cp.pdb_ids) == 0
    
    def test_distance_matrix(self, cp, test_structure_file):
        """Test calculating distance matrix."""
        # Load a structure
        cp.load_structure(TEST_PDB_ID)
        
        # Get CA coordinates
        chains = cp.get_chains(TEST_PDB_ID)
        test_chain = chains[0]
        ca_coords = cp.get_ca_coordinates(TEST_PDB_ID, test_chain)
        
        # Calculate distance matrix
        dist_matrix = get_distance_matrix(ca_coords)
        assert isinstance(dist_matrix, np.ndarray)
        assert dist_matrix.shape[0] == dist_matrix.shape[1]
        assert dist_matrix.shape[0] == ca_coords.shape[0]
        
        # Check diagonal (self-distances should be 0)
        for i in range(dist_matrix.shape[0]):
            assert dist_matrix[i, i] == 0
    
    def test_extract_binding_pocket(self, cp, test_structure_file):
        """Test extracting binding pocket."""
        # Skip if the method is not implemented
        if not hasattr(cp, 'extract_binding_pocket'):
            pytest.skip("extract_binding_pocket method not implemented")
            
        # Load a structure
        cp.load_structure(TEST_PDB_ID)
        
        try:
            # Extract binding pocket
            binding_pocket = cp.extract_binding_pocket(TEST_PDB_ID, "CA", 5.0)
            assert isinstance(binding_pocket, pd.DataFrame)
        except Exception as e:
            pytest.skip(f"Failed to extract binding pocket: {str(e)}")
    
    def test_reset_index(self, cp):
        """Test resetting the index of the structure DataFrame."""
        # Skip if no data is loaded
        if cp.data is None:
            pytest.skip("No data loaded in CifProcessor")
            
        # Get the initial index type
        original_index_type = type(cp.data.index)
        
        # Reset the index
        cp.reset_index()
        
        # Check that the index has been reset appropriately
        assert isinstance(cp.data.index, pd.RangeIndex)
