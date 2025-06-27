"""
Functional tests for download functionality that test error handling and edge cases.
These tests don't require network access.
"""

import os
import pytest
import tempfile
from unittest.mock import patch, MagicMock
import requests

from protos.loaders.download_structures import download_protein_structures
from protos.loaders.alphafold_utils import download_alphafold_structures
from protos.io.paths.path_config import ProtosPaths


@pytest.fixture
def temp_download_dir():
    """Create a temporary directory for downloads."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Save current global setting
        original_global = ProtosPaths.get_global_data_root()
        
        try:
            # Set global data root for this test
            ProtosPaths.set_data_root(temp_dir)
            yield temp_dir
        finally:
            # Restore original global setting
            if original_global:
                ProtosPaths.set_data_root(original_global)
            else:
                ProtosPaths._global_data_root = None


def test_download_creates_directory(temp_download_dir):
    """Test that download functions create directories if they don't exist."""
    # Test with a non-existent directory
    new_dir = os.path.join(temp_download_dir, "new_downloads")
    assert not os.path.exists(new_dir)
    
    # Mock the actual download to avoid network calls
    from Bio.PDB import PDBList
    with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
        mock_retrieve.return_value = os.path.join(new_dir, "test.cif")
        
        # Call the download function
        download_protein_structures(["test"], target_folder=new_dir)
        
        # Directory should have been created
        assert os.path.exists(new_dir)


def test_alphafold_download_creates_directory():
    """Test that AlphaFold download creates directory if needed."""
    with tempfile.TemporaryDirectory() as temp_dir:
        new_dir = os.path.join(temp_dir, "alphafold", "structures")
        assert not os.path.exists(new_dir)
        
        # Mock the requests to avoid network calls
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.content = b"mock protein structure data"
        
        with patch('requests.get', return_value=mock_response):
            # Call the download function
            download_alphafold_structures("P12345", max_models=1, output_dir=new_dir)
            
            # Directory should have been created
            assert os.path.exists(new_dir)
            
            # File should have been created
            expected_file = os.path.join(new_dir, "AF-P12345-F1-model_v1.cif")
            assert os.path.exists(expected_file)
            
            # Verify content was written
            with open(expected_file, 'rb') as f:
                assert f.read() == b"mock protein structure data"


def test_alphafold_download_handles_failure():
    """Test that AlphaFold download handles HTTP errors gracefully."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Mock failed request
        mock_response = MagicMock()
        mock_response.status_code = 404
        
        # Capture print output
        import io
        import sys
        captured_output = io.StringIO()
        sys.stdout = captured_output
        
        try:
            with patch('requests.get', return_value=mock_response):
                # Call the download function
                download_alphafold_structures("INVALID", max_models=1, output_dir=temp_dir)
                
                # Should not crash, check output
                output = captured_output.getvalue()
                assert "Failed to download" in output
                assert "404" in output
                
                # No file should be created
                expected_file = os.path.join(temp_dir, "AF-INVALID-F1-model_v1.cif")
                assert not os.path.exists(expected_file)
        finally:
            sys.stdout = sys.__stdout__


def test_alphafold_multiple_models():
    """Test downloading multiple AlphaFold models."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Mock successful responses for 3 models
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.content = b"mock structure"
        
        with patch('requests.get', return_value=mock_response) as mock_get:
            # Download 3 models
            download_alphafold_structures("P12345", max_models=3, output_dir=temp_dir)
            
            # Should have made 3 requests
            assert mock_get.call_count == 3
            
            # Check all 3 files were created
            for i in range(1, 4):
                expected_file = os.path.join(temp_dir, f"AF-P12345-F1-model_v{i}.cif")
                assert os.path.exists(expected_file)


def test_download_empty_pdb_list():
    """Test handling of empty PDB ID list."""
    with tempfile.TemporaryDirectory() as temp_dir:
        from Bio.PDB import PDBList
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            # Download with empty list
            download_protein_structures([], target_folder=temp_dir)
            
            # Should not call retrieve_pdb_file
            mock_retrieve.assert_not_called()


def test_download_with_paths_integration(temp_download_dir):
    """Test integration with ProtosPaths system."""
    paths = ProtosPaths(create_dirs=True)
    
    # Get standard directories
    struct_dir = paths.get_structure_subdir_path("structure_dir")
    assert os.path.exists(struct_dir)
    
    # Test that download functions work with these paths
    from Bio.PDB import PDBList
    with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
        mock_retrieve.return_value = os.path.join(struct_dir, "test.cif")
        
        # Should not raise any errors
        download_protein_structures(["test"], target_folder=struct_dir)
        mock_retrieve.assert_called_once()


def test_download_with_absolute_path():
    """Test download with absolute path."""
    with tempfile.TemporaryDirectory() as temp_dir:
        abs_path = os.path.abspath(os.path.join(temp_dir, "downloads"))
        
        from Bio.PDB import PDBList
        with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
            mock_retrieve.return_value = os.path.join(abs_path, "test.cif")
            
            # Should create directory and work correctly
            download_protein_structures(["test"], target_folder=abs_path)
            assert os.path.exists(abs_path)


def test_download_with_relative_path():
    """Test download with relative path."""
    original_cwd = os.getcwd()
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            os.chdir(temp_dir)
            rel_path = "downloads"
            
            from Bio.PDB import PDBList
            with patch.object(PDBList, 'retrieve_pdb_file') as mock_retrieve:
                mock_retrieve.return_value = os.path.join(rel_path, "test.cif")
                
                # Should create directory relative to current dir
                download_protein_structures(["test"], target_folder=rel_path)
                assert os.path.exists(rel_path)
                assert os.path.abspath(rel_path) == os.path.join(temp_dir, rel_path)
    finally:
        os.chdir(original_cwd)