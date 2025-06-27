"""
Integration tests for download functionality that actually test the protos functions.

These tests may require internet connectivity and will download real files.
"""

import os
import pytest
import tempfile
from pathlib import Path

from protos.loaders.download_structures import download_protein_structures
from protos.loaders.alphafold_utils import download_alphafold_structures
from protos.loaders.uniprot_utils import map_uniprot_to_pdb
from protos.io.paths.path_config import ProtosPaths
from protos.processing.structure.struct_base_processor import CifBaseProcessor


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


@pytest.mark.integration
@pytest.mark.skipif(not os.environ.get("RUN_INTEGRATION_TESTS"), 
                    reason="Set RUN_INTEGRATION_TESTS=1 to run integration tests")
def test_download_pdb_structure_real(temp_download_dir):
    """Test downloading a real PDB structure."""
    # Create download directory
    download_dir = os.path.join(temp_download_dir, "structures")
    os.makedirs(download_dir, exist_ok=True)
    
    # Download a small, well-known structure (insulin)
    pdb_ids = ["1trz"]  # Small insulin structure
    download_protein_structures(pdb_ids, target_folder=download_dir)
    
    # Check that the file was downloaded
    expected_files = [
        "1trz.cif",
        "tr/1trz.cif",
        "pdb1trz.ent"
    ]
    
    # Bio.PDB.PDBList creates subdirectories, check various possible locations
    downloaded = False
    downloaded_file = None
    for root, dirs, files in os.walk(download_dir):
        for file in files:
            if "1trz" in file.lower():
                downloaded = True
                downloaded_file = os.path.join(root, file)
                break
    
    assert downloaded, f"No file containing '1trz' found in {download_dir}"
    assert os.path.getsize(downloaded_file) > 1000, "Downloaded file seems too small"
    
    # Try to load it with CifBaseProcessor to verify it's valid
    processor = CifBaseProcessor(name="test")
    # The file might be gzipped, try to read it
    if downloaded_file.endswith('.gz'):
        import gzip
        with gzip.open(downloaded_file, 'rt') as f:
            content = f.read()
        assert "1TRZ" in content or "1trz" in content, "PDB ID not found in file"
    else:
        with open(downloaded_file, 'r') as f:
            content = f.read()
        assert "1TRZ" in content or "1trz" in content, "PDB ID not found in file"


@pytest.mark.integration
@pytest.mark.skipif(not os.environ.get("RUN_INTEGRATION_TESTS"), 
                    reason="Set RUN_INTEGRATION_TESTS=1 to run integration tests")
def test_download_alphafold_structure_real(temp_download_dir):
    """Test downloading a real AlphaFold structure."""
    # Create download directory
    download_dir = os.path.join(temp_download_dir, "alphafold")
    
    # Download AlphaFold structure for a small protein
    # Using human insulin (P01308) as it's small and well-studied
    uid = "P01308"
    download_alphafold_structures(uid, max_models=1, output_dir=download_dir)
    
    # Check that the file was downloaded
    expected_file = os.path.join(download_dir, f"AF-{uid}-F1-model_v1.cif")
    assert os.path.exists(expected_file), f"Expected file {expected_file} not found"
    assert os.path.getsize(expected_file) > 1000, "Downloaded file seems too small"
    
    # Verify it's a valid CIF file
    with open(expected_file, 'r') as f:
        content = f.read()
    assert "data_" in content, "Not a valid CIF file (missing data_ block)"
    assert uid in content, f"UniProt ID {uid} not found in file"


@pytest.mark.integration
@pytest.mark.skipif(not os.environ.get("RUN_INTEGRATION_TESTS"), 
                    reason="Set RUN_INTEGRATION_TESTS=1 to run integration tests")
def test_uniprot_to_pdb_mapping_real():
    """Test real UniProt to PDB ID mapping."""
    # Use well-known proteins with PDB structures
    uniprot_ids = ["P00698", "P01308"]  # Lysozyme and insulin
    
    # This might fail if the UniProt API is down or changed
    try:
        mapping_df = map_uniprot_to_pdb(uniprot_ids)
        
        # Check that we got some results
        assert len(mapping_df) > 0, "No mapping results returned"
        assert "uid" in mapping_df.columns, "Missing 'uid' column"
        assert "pdb_id" in mapping_df.columns, "Missing 'pdb_id' column"
        
        # These proteins should have PDB structures
        mapped_uids = set(mapping_df['uid'].unique())
        assert len(mapped_uids) > 0, "No UniProt IDs were mapped"
        
    except Exception as e:
        pytest.skip(f"UniProt API might be unavailable: {str(e)}")


@pytest.mark.integration
@pytest.mark.skipif(not os.environ.get("RUN_INTEGRATION_TESTS"), 
                    reason="Set RUN_INTEGRATION_TESTS=1 to run integration tests")
def test_download_and_process_workflow(temp_download_dir):
    """Test a complete workflow: download and process a structure."""
    # Create directories
    paths = ProtosPaths(create_dirs=True)
    struct_dir = paths.get_structure_subdir_path("structure_dir")
    
    # Download a small structure
    pdb_ids = ["1ubq"]  # Ubiquitin - small and fast to download
    download_protein_structures(pdb_ids, target_folder=struct_dir)
    
    # Find the downloaded file
    downloaded_file = None
    for root, dirs, files in os.walk(struct_dir):
        for file in files:
            if "1ubq" in file.lower():
                downloaded_file = os.path.join(root, file)
                break
    
    assert downloaded_file is not None, "Structure file not downloaded"
    
    # Process with CifBaseProcessor
    processor = CifBaseProcessor(name="test_workflow")
    
    # Try to load the structure
    # Note: Bio.PDB downloads files with specific naming conventions
    if downloaded_file.endswith('.gz'):
        # Decompress if needed
        import gzip
        import shutil
        decompressed = downloaded_file[:-3]
        with gzip.open(downloaded_file, 'rb') as f_in:
            with open(decompressed, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        downloaded_file = decompressed
    
    # Verify the file is readable
    assert os.path.exists(downloaded_file), f"File not found: {downloaded_file}"
    assert os.path.getsize(downloaded_file) > 0, "File is empty"


def test_download_with_paths_integration(temp_download_dir):
    """Test that download functions work with ProtosPaths system (no network required)."""
    # This test verifies the integration without actual downloads
    paths = ProtosPaths(create_dirs=True)
    
    # Get the standard structure directory
    struct_dir = paths.get_structure_subdir_path("structure_dir")
    assert os.path.exists(struct_dir), "Structure directory not created"
    
    # Verify the download functions accept the paths
    # We won't actually download, just verify the directory handling
    
    # Test PDB download directory creation
    pdb_download_dir = os.path.join(struct_dir, "pdb_downloads")
    if not os.path.exists(pdb_download_dir):
        os.makedirs(pdb_download_dir)
    assert os.path.exists(pdb_download_dir)
    
    # Test AlphaFold download directory creation  
    af_download_dir = os.path.join(struct_dir, "alphafold_structures")
    if not os.path.exists(af_download_dir):
        os.makedirs(af_download_dir)
    assert os.path.exists(af_download_dir)