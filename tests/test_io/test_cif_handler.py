"""
Tests for CIF file handling functionality.

Uses pytest's tmp_path fixture for clean test isolation.
"""

import pytest
import pandas as pd
import numpy as np

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


def test_df_to_cif(sample_cif_df, cif_handler, tmp_path):
    """Test converting DataFrame to CIF format."""
    # Use the handler to write to a temporary file
    output_file = tmp_path / 'test_df_to_cif.cif'
    cif_handler.write(str(output_file), sample_cif_df)
    
    # Verify file exists and has content
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    
    # Read the file content
    cif_content = output_file.read_text()
    
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
    output_file = tmp_path / 'output.cif'
    
    # Write the file using the handler
    cif_handler.write(
        str(output_file),
        sample_cif_df,
        force_overwrite=True
    )
    
    # Check file was created
    assert output_file.exists()
    assert output_file.stat().st_size > 0
    
    # Read content and verify
    content = output_file.read_text()
    assert 'data_' in content
    assert 'ATOM   1' in content


def test_write_cif_file_versioned(sample_cif_df, cif_handler, tmp_path):
    """Test writing DataFrame to CIF file with versioning."""
    # Define output path
    output_file = tmp_path / 'versioned_output.cif'
    
    # Write file with versioning
    result_path = cif_handler.write_with_versioning(
        file_path=str(output_file),
        data=sample_cif_df,
        versioned=True
    )
    
    # Check file was created with version in filename
    assert result_path is not None
    # Convert result_path to Path for easier checking
    result_path_obj = tmp_path / result_path.split('/')[-1]
    assert '_v1' in result_path_obj.name
    
    # Write another version
    result_path_2 = cif_handler.write_with_versioning(
        file_path=str(output_file),
        data=sample_cif_df,
        versioned=True
    )
    
    # Check new version was created
    assert result_path_2 is not None
    result_path_2_obj = tmp_path / result_path_2.split('/')[-1]
    assert '_v2' in result_path_2_obj.name


def test_cif_file_overwrite_protection(sample_cif_df, cif_handler, tmp_path):
    """Test CIF file overwrite protection."""
    # Define output path
    output_file = tmp_path / 'protected.cif'
    
    # Write initial file
    cif_handler.write(
        str(output_file),
        sample_cif_df
    )
    
    # Verify file exists
    assert output_file.exists()
    
    # Try to overwrite without force flag - should raise an exception
    with pytest.raises(Exception):
        cif_handler.write(
            str(output_file),
            sample_cif_df,
            force_overwrite=False
        )
    
    # Now try with force_overwrite=True - should succeed
    cif_handler.write(
        str(output_file),
        sample_cif_df,
        force_overwrite=True
    )
    assert output_file.exists()


def test_cif_format_requirements(sample_cif_df, cif_handler, tmp_path):
    """Test that CIF files meet format requirements."""
    output_file = tmp_path / 'format_test.cif'
    cif_handler.write(str(output_file), sample_cif_df)
    
    content = output_file.read_text()
    
    # Check required CIF sections
    required_sections = [
        'data_',
        'loop_',
        '_atom_site.',
        'ATOM'
    ]
    
    for section in required_sections:
        assert section in content, f"Missing required section: {section}"
    
    # Check coordinate format (should have decimal places)
    assert '10.000' in content or '10.0' in content
    assert '13.500' in content or '13.5' in content