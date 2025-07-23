"""
Fixtures for structure processor tests.

This module provides common fixtures for testing the StructureProcessor,
ensuring all tests use real downloadable PDB structures.
"""

import pytest
from pathlib import Path
from protos.processing.structure import StructureProcessor
from protos.loaders.download_structures import download_protein_structures


# List of real, small PDB structures that are good for testing
REAL_TEST_PDB_IDS = [
    "1ubq",  # Ubiquitin - small, well-studied protein
    "1crn",  # Crambin - one of the smallest proteins
    "1l2y",  # Trp-cage miniprotein - very small
    "2gb1",  # Protein G B1 domain - small, stable
    "1tqn",  # Small peptide
    "1uaz",  # Small protein
    "3nir",  # Small opsin fragment
]


@pytest.fixture
def real_pdb_ids():
    """Provide a list of real PDB IDs for testing."""
    return REAL_TEST_PDB_IDS[:3]  # Use first 3 for faster tests


@pytest.fixture
def structure_processor(configure_test_paths):
    """Create a StructureProcessor instance for testing."""
    processor = StructureProcessor(name="test_processor")
    return processor


@pytest.fixture
def processor_with_real_structures(structure_processor, real_pdb_ids):
    """
    Create a StructureProcessor with real downloaded structures.
    
    This fixture ensures that real PDB structures are available for testing,
    downloading them if necessary.
    """
    processor = structure_processor
    
    # Check which structures are already available
    available = processor.get_available_pdb_files()
    
    # Download any missing structures
    missing = [pdb_id for pdb_id in real_pdb_ids if pdb_id not in available]
    if missing:
        try:
            # Use the download functionality
            download_protein_structures(
                missing,
                output_dir=processor.paths.get_subdir_path("structure", "structure_dir"),
                file_format='cif'
            )
        except Exception as e:
            pytest.skip(f"Could not download test structures: {e}")
    
    # Load the structures
    processor.pdb_ids = real_pdb_ids
    processor.load_structures()
    
    return processor


@pytest.fixture
def real_structure_data():
    """
    Provide a sample of real structure data for testing.
    
    This is based on the actual structure of ubiquitin (1UBQ).
    """
    import pandas as pd
    
    # Minimal but realistic structure data
    data = pd.DataFrame({
        'group': ['ATOM', 'ATOM', 'ATOM', 'ATOM', 'ATOM'],
        'pdb_id': ['1ubq', '1ubq', '1ubq', '1ubq', '1ubq'],
        'auth_chain_id': ['A', 'A', 'A', 'A', 'A'],
        'auth_seq_id': [1, 1, 1, 1, 2],
        'auth_comp_id': ['MET', 'MET', 'MET', 'MET', 'GLN'],
        'atom_name': ['N', 'CA', 'C', 'O', 'N'],
        'x': [27.340, 27.130, 25.760, 25.430, 24.860],
        'y': [24.430, 22.930, 22.540, 21.390, 23.460],
        'z': [2.614, 2.749, 2.199, 1.916, 2.078],
        'element': ['N', 'C', 'C', 'O', 'N'],
        'b_factor': [25.00, 15.00, 20.00, 25.00, 30.00],
        'model_num': [0, 0, 0, 0, 0]
    })
    
    return data


def ensure_test_structures_available():
    """
    Ensure that test structures are available in the test data directory.
    
    This function is called during test setup to make sure we have
    real structures to test with.
    """
    from protos.io.paths import ProtosPaths
    
    paths = ProtosPaths()
    processor = StructureProcessor(name="test_setup")
    
    # Get available structures
    available = processor.get_available_pdb_files()
    
    # We need at least 3 structures for testing
    if len(available) < 3:
        # Download some small, reliable structures
        to_download = []
        for pdb_id in REAL_TEST_PDB_IDS:
            if pdb_id not in available:
                to_download.append(pdb_id)
            if len(to_download) + len(available) >= 3:
                break
        
        if to_download:
            try:
                download_protein_structures(
                    to_download,
                    output_dir=processor.paths.get_subdir_path("structure", "structure_dir"),
                    file_format='cif'
                )
            except Exception as e:
                print(f"Warning: Could not download test structures: {e}")
    
    return available