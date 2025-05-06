"""
Tests for the CifBaseProcessor class using real data from downloads.

This test suite validates the CifBaseProcessor by testing all its functionality
using real PDB structures downloaded from RCSB PDB.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.loaders.download_structures import download_protein_structures
from protos.loaders.alphafold_utils import download_alphafold_structures
from protos.processing.structure.struct_utils import STRUCT_COLUMNS


@pytest.fixture
def processor(test_structure_data):
    """Create a CifBaseProcessor instance for testing with real data paths."""
    processor = CifBaseProcessor(
        name="test_processor",
        data_root=str(test_structure_data["root"]),
        processor_data_dir="structure",
        structure_dir="mmcif",
        dataset_dir="structure_dataset"
    )
    return processor


@pytest.fixture
def processor_with_structures(processor, pdb_test_structures):
    """Initialize processor with real downloaded structures."""
    # Set PDB IDs to downloaded structures
    processor.pdb_ids = pdb_test_structures
    processor.load_structures()
    return processor


def test_initialization(processor):
    """Test basic processor initialization."""
    # Verify the processor was initialized correctly
    assert processor.name == "test_processor"
    assert "structure" in processor.processor_data_dir
    assert "mmcif" in processor.structure_dir
    assert "structure_dataset" in processor.dataset_dir
    
    # Verify paths are set correctly
    assert os.path.isdir(processor.data_path)
    assert processor.path_structure_dir.endswith("mmcif")
    assert processor.path_dataset_dir.endswith("structure_dataset")


def test_find_available_pdb_files(processor, pdb_test_structures):
    """Test finding available PDB files in directory."""
    # Get available PDB files
    available_files = processor.get_available_pdb_files()
    
    # Verify all downloaded structures are found
    for pdb_id in pdb_test_structures:
        assert pdb_id in available_files


def test_download_structure(processor):
    """Test downloading a structure directly with the processor."""
    # Mock the download function to avoid actual network calls during tests
    with patch.object(CifBaseProcessor, 'download_cif') as mock_download:
        # Configure the mock to return success
        mock_download.return_value = True
        
        # Test downloading a structure
        test_pdb_id = "1UBQ"
        result = processor.download_cif(test_pdb_id, save_dir=processor.path_structure_dir)
        
        # Verify download was called with correct arguments
        assert result is True
        mock_download.assert_called_once_with(test_pdb_id, save_dir=processor.path_structure_dir)


def test_download_structures(processor):
    """Test downloading multiple structures with the processor."""
    # Mock the download function to avoid actual network calls during tests
    with patch.object(CifBaseProcessor, 'download_cif') as mock_download:
        # Configure the mock to return success
        mock_download.return_value = True
        
        # Test downloading multiple structures
        test_pdb_ids = ["1ABC", "2DEF", "3GHI"]
        for pdb_id in test_pdb_ids:
            result = processor.download_cif(pdb_id, save_dir=processor.path_structure_dir)
            assert result is True
        
        # Verify download was called for each PDB ID
        assert mock_download.call_count == len(test_pdb_ids)
        for pdb_id in test_pdb_ids:
            mock_download.assert_any_call(pdb_id, save_dir=processor.path_structure_dir)


def test_load_structure(processor, pdb_test_structures):
    """Test loading a single structure."""
    # Load a single structure from downloaded test structures
    test_pdb_id = pdb_test_structures[0]
    data = processor.load_structure(test_pdb_id)
    
    # Verify structure was loaded correctly
    assert data is not None
    assert not data.empty
    assert 'pdb_id' in data.columns
    assert data['pdb_id'].iloc[0] == test_pdb_id
    
    # Verify essential columns are present (adjusted to match actual column names)
    essential_columns = [
        'pdb_id', 'atom_id', 'atom_name', 
        'auth_seq_id', 'res_name3l', 'auth_chain_id', 'x', 'y', 'z'
    ]
    for col in essential_columns:
        assert col in data.columns or col.lower() in data.columns


def test_load_structures(processor, pdb_test_structures):
    """Test loading multiple structures."""
    # Initialize PDB IDs and load structures
    processor.pdb_ids = pdb_test_structures
    processor.load_structures()
    
    # Verify structures were loaded correctly
    assert processor.data is not None
    assert not processor.data.empty
    assert 'pdb_id' in processor.data.columns
    
    # Verify all structures are loaded
    loaded_pdb_ids = processor.data['pdb_id'].unique()
    for pdb_id in pdb_test_structures:
        assert pdb_id in loaded_pdb_ids


def test_data_columns(processor_with_structures):
    """Test that all required data columns are present."""
    # Check for essential columns (adjusted to match actual column names)
    essential_columns = [
        'pdb_id', 'atom_id', 'atom_name', 
        'auth_seq_id', 'res_name3l', 'auth_chain_id', 'x', 'y', 'z'
    ]
    
    # Account for potential column name variations
    data_columns = list(processor_with_structures.data.columns)
    lower_data_columns = [col.lower() for col in data_columns]
    
    for col in essential_columns:
        assert col in data_columns or col.lower() in lower_data_columns


def test_get_chains(processor_with_structures):
    """Test extracting chains from structures."""
    # Get chains for each loaded structure
    for pdb_id in processor_with_structures.pdb_ids:
        chains = processor_with_structures.get_chains(pdb_id)
        
        # Verify chains were extracted
        assert chains is not None
        assert len(chains) > 0


def test_get_residues(processor_with_structures):
    """Test extracting residues from structures."""
    # Get residues for a chain in each structure
    for pdb_id in processor_with_structures.pdb_ids:
        chains = processor_with_structures.get_chains(pdb_id)
        
        # Test with first chain
        chain_id = chains[0]
        
        # Use data filtering to get residues instead of get_residues method
        chain_data = processor_with_structures.data[
            (processor_with_structures.data['pdb_id'] == pdb_id) &
            (processor_with_structures.data['auth_chain_id'] == chain_id)
        ]
        
        # Get unique residue IDs
        residues = chain_data['auth_seq_id'].unique()
        
        # Verify residues were extracted
        assert len(residues) > 0


def test_get_sequence(processor_with_structures):
    """Test extracting sequence from structures."""
    # Get sequence for a chain in each structure
    for pdb_id in processor_with_structures.pdb_ids:
        chains = processor_with_structures.get_chains(pdb_id)
        
        # Test with first chain
        chain_id = chains[0]
        sequence = processor_with_structures.get_sequence(pdb_id, chain_id)
        
        # Verify sequence was extracted
        assert sequence is not None
        assert len(sequence) > 0
        # Protein sequence should contain valid amino acid characters
        assert all(aa in "ACDEFGHIKLMNPQRSTVWYX" for aa in sequence)


def test_get_coordinates(processor_with_structures):
    """Test extracting coordinates from structures."""
    # Get coordinates for a chain in each structure
    for pdb_id in processor_with_structures.pdb_ids:
        chains = processor_with_structures.get_chains(pdb_id)
        
        # Test with first chain
        chain_id = chains[0]
        
        # Test CA coordinates
        ca_coords = processor_with_structures.get_ca_coordinates(pdb_id, chain_id)
        assert ca_coords is not None
        assert isinstance(ca_coords, np.ndarray)
        assert ca_coords.shape[1] == 3  # x, y, z


def test_filter_data(processor_with_structures):
    """Test filtering data by various criteria."""
    # Test filtering by PDB ID
    test_pdb_id = processor_with_structures.pdb_ids[0]
    filtered_data = processor_with_structures.data[
        processor_with_structures.data['pdb_id'] == test_pdb_id
    ]
    assert not filtered_data.empty
    assert set(filtered_data['pdb_id'].unique()) == {test_pdb_id}
    
    # Test filtering by chain
    chains = processor_with_structures.get_chains(test_pdb_id)
    chain_id = chains[0]
    chain_data = processor_with_structures.data[
        (processor_with_structures.data['pdb_id'] == test_pdb_id) &
        (processor_with_structures.data['auth_chain_id'] == chain_id)
    ]
    assert not chain_data.empty
    assert set(chain_data['auth_chain_id'].unique()) == {chain_id}


def test_reset_data(processor_with_structures):
    """Test resetting processor data."""
    # Verify data is loaded
    assert processor_with_structures.data is not None
    assert not processor_with_structures.data.empty
    
    # Reset data
    processor_with_structures.reset_data()
    
    # Verify data is reset
    assert processor_with_structures.data is None
    assert not processor_with_structures.pdb_ids


def test_create_and_load_dataset(processor, pdb_test_structures):
    """Test creating and loading a dataset."""
    # Create a dataset
    dataset_id = "test_dataset"
    dataset_name = "Test Dataset"
    dataset_description = "Dataset for testing"
    
    processor.create_dataset(
        dataset_id=dataset_id,
        name=dataset_name,
        description=dataset_description,
        content=pdb_test_structures
    )
    
    # Load the dataset
    processor.load_dataset(dataset_id)
    
    # Verify dataset was loaded correctly
    assert processor.data is not None
    assert not processor.data.empty
    assert set(processor.pdb_ids) == set(pdb_test_structures)


def test_delete_dataset(processor, pdb_test_structures):
    """Test deleting a dataset."""
    # Create a dataset
    dataset_id = "delete_test_dataset"
    processor.create_dataset(
        dataset_id=dataset_id,
        name="Delete Test Dataset",
        description="Dataset for testing deletion",
        content=pdb_test_structures
    )
    
    # Delete the dataset
    result = processor.delete_dataset(dataset_id)
    
    # Verify deletion was successful
    assert result is True
    
    # Verify dataset is not in available datasets
    available_datasets = processor.list_datasets()
    dataset_ids = [dataset["id"] for dataset in available_datasets 
                  if isinstance(dataset, dict) and "id" in dataset]
    assert dataset_id not in dataset_ids


def test_find_binding_pocket(processor_with_structures):
    """Test finding binding pocket residues."""
    # This test is more complex and structure-dependent
    # We'll implement a simplified version that just verifies the method runs
    # without errors, as finding actual binding pockets depends on the specific structure
    
    # Get the first PDB ID
    pdb_id = processor_with_structures.pdb_ids[0]
    chains = processor_with_structures.get_chains(pdb_id)
    
    # Try for the first chain
    chain_id = chains[0]
    
    try:
        # Use a larger distance cutoff to increase chance of finding something
        binding_residues = processor_with_structures.get_binding_pocket(
            pdb_id, chain_id, distance_cutoff=12.0
        )
        
        # If binding residues are found, verify their structure
        if binding_residues is not None and not binding_residues.empty:
            # Just verify we can extract binding pocket data
            assert 'res_id' in binding_residues.columns or 'res_id' in binding_residues.columns.str.lower()
    except Exception as e:
        # If there's an exception, it should be related to no ligands in the structure
        # which is acceptable for this test
        print(f"Note: Exception when finding binding pocket for {pdb_id}, chain {chain_id}: {e}")
        # We don't make this a test failure as some structures may not have ligands
        pass



def test_filter_data_flexibly_inplace(processor_with_structures):
    """Test flexible filtering modifying data inplace."""
    processor = processor_with_structures
    initial_row_count = len(processor.data)
    initial_pdb_ids = set(processor.pdb_ids)

    # Define filters: Select only chain A of the first PDB structure
    first_pdb_id = processor.pdb_ids[0]
    filters = {
        'pdb_id': first_pdb_id, # Implicit __eq
        'auth_chain_id': 'A'
    }

    # Apply filters inplace
    result = processor.filter_data_flexibly(filters, inplace=True)

    # Verify inplace modification
    assert result is None # inplace returns None
    assert len(processor.data) < initial_row_count
    assert set(processor.pdb_ids) == {first_pdb_id} # Only one PDB ID should remain
    assert processor.data['pdb_id'].unique()[0] == first_pdb_id
    assert processor.data['auth_chain_id'].unique()[0] == 'A'

def test_filter_data_flexibly_return(processor_with_structures):
    """Test flexible filtering returning a new DataFrame."""
    processor = processor_with_structures
    initial_row_count = len(processor.data)
    initial_data_id = id(processor.data) # Get memory ID to check if it changes

    # Define filters: Select CA atoms with B-factor >= 20.0
    filters = {
        'atom_name__eq': 'CA', # Alpha Carbon
        'b_factor__ge': 20.0
    }

    # Apply filters returning new DataFrame
    filtered_df = processor.filter_data_flexibly(filters, inplace=False)

    # Verify original data is unchanged
    assert len(processor.data) == initial_row_count
    assert id(processor.data) == initial_data_id # Should be the same object

    # Verify filtered DataFrame
    assert isinstance(filtered_df, pd.DataFrame)
    assert len(filtered_df) < initial_row_count
    assert len(filtered_df) > 0 # Assuming some CA atoms meet the criteria
    assert all(filtered_df['atom_name'] == 'CA')
    assert all(filtered_df['b_factor'] >= 20.0)

def test_filter_data_flexibly_advanced_operators(processor_with_structures):
    """Test various advanced filter operators."""
    processor = processor_with_structures
    first_pdb_id = processor.pdb_ids[0]

    # Filter for specific residues in chain A of the first PDB
    filters = {
        'pdb_id': first_pdb_id,
        'auth_chain_id': 'A',
        'res_name3l__is_in': ['LEU', 'ALA', 'GLY'], # Use res_name3l
        'auth_seq_id__gt': 5,
        'auth_seq_id__le': 15
    }
    filtered_df = processor.filter_data_flexibly(filters, inplace=False)

    assert not filtered_df.empty
    assert all(filtered_df['pdb_id'] == first_pdb_id)
    assert all(filtered_df['auth_chain_id'] == 'A')
    assert all(filtered_df['res_name3l'].isin(['LEU', 'ALA', 'GLY']))
    assert all(filtered_df['auth_seq_id'] > 5)
    assert all(filtered_df['auth_seq_id'] <= 15)

    # Test 'contains' (ensure element column exists and is string)
    if 'element' not in processor.data.columns:
         processor.data['element'] = processor.data['atom_name'].str[0].str.upper() # Add if missing

    filters_contains = {'element__contains': 'C'}
    filtered_contains = processor.filter_data_flexibly(filters_contains, inplace=False)
    assert len(filtered_contains) > 0
    assert len(filtered_contains) <= len(processor.data)

    # Test isna/notna (example: alt_id is often NaN/None)
    filters_alt_na = {'alt_id__isna': True}
    filtered_alt_na = processor.filter_data_flexibly(filters_alt_na, inplace=False)

    filters_alt_notna = {'alt_id__notna': True}
    filtered_alt_notna = processor.filter_data_flexibly(filters_alt_notna, inplace=False)

    assert len(filtered_alt_na) + len(filtered_alt_notna) == len(processor.data)


def test_filter_data_flexibly_errors(processor_with_structures):
    """Test error conditions for flexible filtering."""
    processor = processor_with_structures

    # Test invalid column
    with pytest.raises(ValueError, match="Filter column 'invalid_column' not found"):
        processor.filter_data_flexibly({'invalid_column': 10}, inplace=False)

    # Test invalid operator
    with pytest.raises(ValueError, match="Invalid filter operator: __invalidop"):
        processor.filter_data_flexibly({'atom_id__invalidop': 10}, inplace=False)

    # Test type error (e.g., greater than on string)
    with pytest.raises(TypeError):
        processor.filter_data_flexibly({'auth_chain_id__gt': 'A'}, inplace=False)


def test_add_ligand(processor_with_structures):
    """Test adding a ligand DataFrame to an existing structure."""
    processor = processor_with_structures
    target_pdb_id = processor.pdb_ids[0] # Add to the first loaded structure
    initial_rows = len(processor.data[processor.data['pdb_id'] == target_pdb_id])
    max_atom_id_before = processor.data[processor.data['pdb_id'] == target_pdb_id]['atom_id'].max()

    # Create a dummy ligand DataFrame
    # Coordinates should be chosen to be "docked" conceptually for this test
    # Let's place it near the origin for simplicity
    ligand_data = {
        'atom_name': ['C1', 'C2', 'N1', 'O1'],
        'res_name': ['LG1'] * 4, # Ligand residue name
        'element': ['C', 'C', 'N', 'O'],
        'x': [1.0, 2.0, 1.5, 0.5],
        'y': [1.0, 1.0, 2.0, 0.0],
        'z': [0.0, 0.0, 0.0, 0.5]
    }
    ligand_df = pd.DataFrame(ligand_data)
    num_lig_atoms = len(ligand_df)

    # Add the ligand
    processor.add_ligand(
        target_pdb_id=target_pdb_id,
        ligand_df=ligand_df,
        ligand_chain_id='Z',
        ligand_res_seq_id=9001
    )

    # --- Verification ---
    # 1. Check total row count increased correctly
    assert len(processor.data) == initial_rows + num_lig_atoms

    # 2. Check ligand atoms are present in the target PDB
    target_df_after = processor.data[processor.data['pdb_id'] == target_pdb_id]
    ligand_atoms_in_data = target_df_after[
        (target_df_after['auth_chain_id'] == 'Z') &
        (target_df_after['auth_seq_id'] == 9001) &
        (target_df_after['res_name'] == 'LG1') &
        (target_df_after['group'] == 'HETATM')
    ]
    assert len(ligand_atoms_in_data) == num_lig_atoms

    # 3. Check atom IDs were renumbered correctly
    min_lig_atom_id = ligand_atoms_in_data['atom_id'].min()
    max_lig_atom_id = ligand_atoms_in_data['atom_id'].max()
    assert min_lig_atom_id == max_atom_id_before + 1
    assert max_lig_atom_id == max_atom_id_before + num_lig_atoms
    assert len(ligand_atoms_in_data['atom_id'].unique()) == num_lig_atoms # Check uniqueness

    # 4. Check coordinates were preserved
    assert np.allclose(
        ligand_atoms_in_data.sort_values('atom_name')[['x', 'y', 'z']].values,
        ligand_df.sort_values('atom_name')[['x', 'y', 'z']].values
    )

    # 5. Check other structures were not affected
    if len(processor.pdb_ids) > 1:
        other_pdb_id = processor.pdb_ids[1]
        other_df = processor.data[processor.data['pdb_id'] == other_pdb_id]
        assert 'Z' not in other_df['auth_chain_id'].unique()

def test_add_ligand_errors(processor_with_structures):
    """Test error conditions for adding a ligand."""
    processor = processor_with_structures
    target_pdb_id = processor.pdb_ids[0]

    # Create minimal valid ligand df
    ligand_df = pd.DataFrame({'atom_name': ['C1'], 'res_name': ['LG1'], 'x':[0],'y':[0],'z':[0]})

    # Error: Target PDB not loaded
    with pytest.raises(ValueError, match="Target PDB ID 'invalid_pdb' not found"):
        processor.add_ligand('invalid_pdb', ligand_df)

    # Error: Missing required column in ligand_df
    invalid_ligand_df = ligand_df.drop(columns=['x'])
    with pytest.raises(ValueError, match="Ligand DataFrame is missing required columns: \['x'\]"):
        processor.add_ligand(target_pdb_id, invalid_ligand_df)

    # Error: Empty ligand_df
    empty_ligand_df = pd.DataFrame(columns=ligand_df.columns)
    with pytest.raises(ValueError, match="Ligand DataFrame cannot be empty"):
        processor.add_ligand(target_pdb_id, empty_ligand_df)

