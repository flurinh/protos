"""
Tests for the refactored structure utilities.

This module tests the functionality of the refactored structure utilities
that use standardized schemas and the path management system.
"""

import os
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from protos.processing.schema.structure.structure_utils import (
    load_structure,
    cif_block_to_df,
    extract_torsion_angles
)
from protos.processing.schema.schema_definitions import (
    STRUCTURE_CORE_COLUMNS,
    STRUCTURE_COLUMN_ORDER
)
from protos.io.paths import get_structure_path


# Path to test data directory
TEST_DIR = Path(__file__).parent.parent.parent.parent / 'test_data'
TEST_STRUCTURE_DIR = TEST_DIR / 'structure_test' / 'mmcif'


class TestLoadStructure:
    """Tests for the load_structure function."""
    
    def test_load_structure_file_not_found(self):
        """Test that load_structure raises FileNotFoundError when file not found."""
        with pytest.raises(FileNotFoundError):
            load_structure('nonexistent_file.cif')
    
    @pytest.mark.parametrize("pdb_id", ["6ott", "4zw9"])
    def test_load_structure_pdb_id(self, pdb_id):
        """Test loading structure by PDB ID."""
        # Mock the get_structure_path function to return a test file path
        test_file = TEST_STRUCTURE_DIR / f"{pdb_id}.cif"
        
        # Skip test if file doesn't exist
        if not test_file.exists():
            pytest.skip(f"Test file {test_file} not found")
        
        with patch('protos.processing.schema.structure.structure_utils.get_structure_path', return_value=str(test_file)):
            # Load structure
            structure_df = load_structure(pdb_id)
            
            # Check that the structure was loaded
            assert structure_df is not None
            assert not structure_df.empty
            
            # Check that the structure has the correct columns
            for col in STRUCTURE_CORE_COLUMNS:
                assert col in structure_df.columns
            
            # Check that the structure has the correct PDB ID
            assert structure_df['pdb_id'].iloc[0] == pdb_id.lower()
    
    def test_load_structure_with_file_path(self):
        """Test loading structure from a file path."""
        # Find a test file
        test_files = list(TEST_STRUCTURE_DIR.glob('*.cif'))
        if not test_files:
            pytest.skip("No test CIF files found")
        
        test_file = test_files[0]
        
        # Load structure
        structure_df = load_structure(test_file)
        
        # Check that the structure was loaded
        assert structure_df is not None
        assert not structure_df.empty
        
        # Check that the structure has the correct columns
        for col in STRUCTURE_CORE_COLUMNS:
            assert col in structure_df.columns
        
        # Check that the columns are in the correct order
        assert list(structure_df.columns) == STRUCTURE_COLUMN_ORDER
    
    def test_load_structure_with_structure_id(self):
        """Test loading structure with a custom structure ID."""
        # Find a test file
        test_files = list(TEST_STRUCTURE_DIR.glob('*.cif'))
        if not test_files:
            pytest.skip("No test CIF files found")
        
        test_file = test_files[0]
        custom_id = "TEST_STRUCTURE"
        
        # Load structure
        structure_df = load_structure(test_file, structure_id=custom_id)
        
        # Check that the structure has the correct PDB ID
        assert structure_df['pdb_id'].iloc[0] == custom_id
    
    def test_torsion_angle_extraction(self):
        """Test that torsion angles are correctly extracted."""
        # Find a test file
        test_files = list(TEST_STRUCTURE_DIR.glob('*.cif'))
        if not test_files:
            pytest.skip("No test CIF files found")
        
        test_file = test_files[0]
        
        # Load structure
        structure_df = load_structure(test_file)
        
        # Check that torsion angles were extracted
        assert 'phi' in structure_df.columns
        assert 'psi' in structure_df.columns
        assert 'omega' in structure_df.columns
        
        # Check for valid values (not all NaN)
        assert not structure_df['phi'].isna().all()
        assert not structure_df['psi'].isna().all()
        assert not structure_df['omega'].isna().all()
    
    def test_column_data_types(self):
        """Test that column data types are correctly set."""
        # Find a test file
        test_files = list(TEST_STRUCTURE_DIR.glob('*.cif'))
        if not test_files:
            pytest.skip("No test CIF files found")
        
        test_file = test_files[0]
        
        # Load structure
        structure_df = load_structure(test_file)
        
        # Check data types
        assert structure_df['pdb_id'].dtype == np.dtype('O')  # string
        assert structure_df['gen_seq_id'].dtype == np.dtype('int64')
        assert structure_df['auth_seq_id'].dtype == np.dtype('int64')
        assert structure_df['atom_id'].dtype == np.dtype('int64')
        assert structure_df['x'].dtype == np.dtype('float64')
        assert structure_df['y'].dtype == np.dtype('float64')
        assert structure_df['z'].dtype == np.dtype('float64')
        
    def test_res_atom_name_construction(self):
        """Test that res_atom_name is correctly constructed."""
        # Find a test file
        test_files = list(TEST_STRUCTURE_DIR.glob('*.cif'))
        if not test_files:
            pytest.skip("No test CIF files found")
        
        test_file = test_files[0]
        
        # Load structure
        structure_df = load_structure(test_file)
        
        # Check res_atom_name construction
        for i, row in structure_df.iterrows():
            if pd.notna(row['res_name3l']) and pd.notna(row['atom_name']):
                expected = f"{row['res_name3l']}.{row['atom_name']}"
                assert row['res_atom_name'] == expected