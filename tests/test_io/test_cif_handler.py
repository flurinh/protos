"""
Tests for CIF file handling functionality.
"""

import os
import tempfile
import pytest
import pandas as pd
import numpy as np

# Import the CifHandler class
from protos.io.cif_handler import CifHandler


@pytest.fixture
def sample_cif_df():
    """Create a sample DataFrame representing CIF data with all required columns."""
    data = {
        'pdb_id': ['test_cif'] * 5,
        'auth_chain_id': ['A'] * 5,
        'group': ['ATOM'] * 5,
        'res_name': ['ALA', 'ARG', 'GLY', 'PHE', 'SER'],  # Required column
        'res_name3l': ['ALA', 'ARG', 'GLY', 'PHE', 'SER'],
        'res_name1l': ['A', 'R', 'G', 'F', 'S'],
        'atom_name': ['CA', 'CA', 'CA', 'CA', 'CA'],  # Required column
        'res_atom_name': ['CA', 'CA', 'CA', 'CA', 'CA'],
        'atom_id': [1, 2, 3, 4, 5],
        'x': [10.0, 13.5, 17.0, 20.5, 24.0],
        'y': [5.0, 8.5, 12.0, 15.5, 19.0],
        'z': [2.0, 5.5, 9.0, 12.5, 16.0],
        'gen_seq_id': [1, 2, 3, 4, 5],
        'auth_seq_id': [1, 2, 3, 4, 5],
        'element': ['C', 'C', 'C', 'C', 'C'],
        'b_factor': [10.0, 12.0, 14.0, 16.0, 18.0],
        'occupancy': [1.0, 1.0, 1.0, 1.0, 1.0]
    }
    return pd.DataFrame(data)


@pytest.fixture
def cif_handler():
    """Create a CifHandler instance."""
    return CifHandler()

@pytest.fixture
def temp_cif_file(sample_cif_df, cif_handler):
    """Create a temporary CIF file from sample data."""
    with tempfile.NamedTemporaryFile(suffix='.cif', delete=False) as tmp:
        # Use the handler to write the file
        tmp_name = cif_handler.write(tmp.name, sample_cif_df, force_overwrite=True)
    
    yield tmp_name
    
    # Clean up after test
    if os.path.exists(tmp_name):
        os.unlink(tmp_name)


def test_df_to_cif(sample_cif_df, cif_handler, tmp_path):
    """Test converting DataFrame to CIF format string."""
    # Use the handler to write to a temporary file
    output_file = os.path.join(tmp_path, 'test_df_to_cif.cif')
    cif_handler.write(output_file, sample_cif_df)
    
    # Verify file exists
    assert os.path.exists(output_file)
    
    # Read the file content
    with open(output_file, 'r') as f:
        cif_content = f.read()
    
    # Check basic CIF format elements
    assert cif_content.startswith('data_')
    assert '_atom_site.group_PDB' in cif_content
    assert '_atom_site.Cartn_x' in cif_content
    
    # Check that all atoms are included
    for i in range(1, 6):
        assert f"ATOM   {i}" in cif_content


def test_write_cif_file(sample_cif_df, cif_handler, tmp_path):
    """Test writing DataFrame to CIF file."""
    # Define output path
    output_file = os.path.join(tmp_path, 'output.cif')
    
    # Write the file using the handler
    cif_handler.write(
        output_file,
        sample_cif_df,
        force_overwrite=True
    )
    
    # Check file was created
    assert os.path.exists(output_file)
    assert os.path.getsize(output_file) > 0
    
    # Read content and verify
    with open(output_file, 'r') as f:
        content = f.read()
        assert 'data_' in content
        assert 'ATOM   1' in content


def test_write_cif_file_versioned(sample_cif_df, cif_handler, tmp_path):
    """Test writing DataFrame to CIF file with versioning."""
    # Define output path
    output_file = os.path.join(tmp_path, 'versioned_output.cif')
    
    # Write file with versioning using the handler's write_with_versioning method
    result_path = cif_handler.write_with_versioning(
        file_path=output_file,
        data=sample_cif_df,
        versioned=True
    )
    
    # Check file was created with version in filename
    assert result_path is not None
    assert os.path.exists(result_path)
    assert result_path != output_file  # Should have version added
    assert '_v1' in os.path.basename(result_path)
    
    # Write another version
    result_path_2 = cif_handler.write_with_versioning(
        file_path=output_file,
        data=sample_cif_df,
        versioned=True
    )
    
    # Check new version was created
    assert result_path_2 is not None
    assert os.path.exists(result_path_2)
    assert '_v2' in os.path.basename(result_path_2)


def test_cif_file_overwrite_protection(sample_cif_df, cif_handler, tmp_path):
    """Test CIF file overwrite protection."""
    # Define output path
    output_file = os.path.join(tmp_path, 'protected.cif')
    
    # Write initial file
    cif_handler.write(
        output_file,
        sample_cif_df
    )
    
    # Verify file exists
    assert os.path.exists(output_file)
    
    # Try to overwrite without force flag
    try:
        cif_handler.write(
            output_file,
            sample_cif_df,
            force_overwrite=False
        )
        # If we get here, the test failed
        assert False, "Should have raised an error trying to overwrite file"
    except Exception:
        # Expected behavior
        pass
    
    # Now try with force_overwrite=True
    try:
        cif_handler.write(
            output_file,
            sample_cif_df,
            force_overwrite=True
        )
        # Should succeed
        assert os.path.exists(output_file)
    except Exception as e:
        assert False, f"Should not have raised an error with force_overwrite=True: {str(e)}"