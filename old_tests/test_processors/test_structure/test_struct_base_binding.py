"""
Tests for the binding pocket extraction functionality of CifBaseProcessor.

These old_tests validate the binding pocket extraction functionality
using real experimental mmCIF data.
"""

import os
import pytest
import pandas as pd
import numpy as np

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.structure.struct_utils import ALPHA_CARBON

# Example data with known ligands
# These are real structures with retinal ligand
TEST_PDB_IDS = [
    "6qyj",  # Bacteriorhodopsin - has RET (retinal) ligand
    "1f88"   # Sensory rhodopsin - has RET (retinal) ligand
]

# Ligand code for retinal
RET_LIGAND = 'RET'

@pytest.fixture
def data_path():
    """Path to the test data directory."""
    test_data_path = os.path.join("protos", "data")
    os.makedirs(test_data_path, exist_ok=True)
    return test_data_path

@pytest.fixture
def structure_dir():
    """Directory containing structure files."""
    return "structure/mmcif"

@pytest.fixture
def test_mmcif_path(data_path, structure_dir):
    """Full path to the test mmCIF directory."""
    mmcif_path = os.path.join(data_path, structure_dir)
    os.makedirs(mmcif_path, exist_ok=True)
    return mmcif_path

@pytest.fixture
def processor(data_path, structure_dir):
    """Create a CifBaseProcessor instance for testing."""
    return CifBaseProcessor(
        name="test_processor",
        data_root=data_path,
        structure_dir=structure_dir
    )


class TestCifBaseProcessorBinding:
    """Tests focused on binding pocket extraction functionality."""
    
    def test_extract_binding_pocket(self, processor):
        """Test extracting binding pocket around a ligand."""
        test_pdb_id = TEST_PDB_IDS[0]
        
        # Check if structure file exists
        structure_path = os.path.join("protos", "data", "structure", "mmcif", f"{test_pdb_id}.cif")
        if not os.path.exists(structure_path):
            pytest.skip(f"Test structure {test_pdb_id}.cif not found at {structure_path}")
        
        # Load the structure
        processor.load_structure(test_pdb_id)
        
        # Extract binding pocket around RET ligand with default distance (5Å)
        binding_pocket = processor.extract_binding_pocket(
            pdb_id=test_pdb_id,
            ligand=RET_LIGAND,
            distance=5.0
        )
        
        # Verify the binding pocket was extracted
        assert binding_pocket is not None
        assert isinstance(binding_pocket, pd.DataFrame)
        assert len(binding_pocket) > 0
        
        # Verify binding pocket atoms are from the correct structure
        assert binding_pocket['pdb_id'].eq(test_pdb_id.lower()).all()
        
        # Check that all atoms are from protein (ATOM group)
        assert binding_pocket['group'].eq('ATOM').all()
        
        # Verify it contains multiple residues
        assert binding_pocket['auth_seq_id'].nunique() > 3
        
        # Try with a smaller distance to ensure it's responsive to distance parameter
        small_pocket = processor.extract_binding_pocket(
            pdb_id=test_pdb_id,
            ligand=RET_LIGAND,
            distance=3.0
        )
        
        # Verify the smaller pocket has fewer atoms
        assert len(small_pocket) < len(binding_pocket)
    
    def test_binding_pocket_with_nonexistent_ligand(self, processor):
        """Test extracting binding pocket with a nonexistent ligand."""
        test_pdb_id = TEST_PDB_IDS[0]
        
        # Check if structure file exists
        structure_path = os.path.join("protos", "data", "structure", "mmcif", f"{test_pdb_id}.cif")
        if not os.path.exists(structure_path):
            pytest.skip(f"Test structure {test_pdb_id}.cif not found at {structure_path}")
        
        # Load the structure
        processor.load_structure(test_pdb_id)
        
        # Try to extract binding pocket with a nonexistent ligand
        binding_pocket = processor.extract_binding_pocket(
            pdb_id=test_pdb_id,
            ligand="XYZ",  # Nonexistent ligand code
            distance=5.0
        )
        
        # Verify we get an empty dataframe, not an error
        assert binding_pocket is not None
        assert isinstance(binding_pocket, pd.DataFrame)
        assert len(binding_pocket) == 0
    
    def test_compare_binding_pockets(self, processor):
        """Test comparing binding pockets from two different structures."""
        # Check if structure files exist
        structure_paths = []
        for pdb_id in TEST_PDB_IDS:
            path = os.path.join("protos", "data", "structure", "mmcif", f"{pdb_id}.cif")
            structure_paths.append(path)
        
        if not all(os.path.exists(path) for path in structure_paths):
            pytest.skip(f"Not all test structures found at {structure_paths}")
            
        # Load both structures
        processor.load_structures(TEST_PDB_IDS)
        
        # Extract binding pockets with same parameters
        pocket1 = processor.extract_binding_pocket(
            pdb_id=TEST_PDB_IDS[0],
            ligand=RET_LIGAND,
            distance=5.0
        )
        
        pocket2 = processor.extract_binding_pocket(
            pdb_id=TEST_PDB_IDS[1],
            ligand=RET_LIGAND,
            distance=5.0
        )
        
        # Verify both pockets were extracted
        assert pocket1 is not None and len(pocket1) > 0
        assert pocket2 is not None and len(pocket2) > 0
        
        # Check that they have the same structure (columns)
        assert set(pocket1.columns) == set(pocket2.columns)
        
        # Each pocket should be from its own structure
        assert pocket1['pdb_id'].eq(TEST_PDB_IDS[0].lower()).all()
        assert pocket2['pdb_id'].eq(TEST_PDB_IDS[1].lower()).all()
        
        # Calculate number of residues in each pocket
        residues1 = pocket1[['auth_chain_id', 'auth_seq_id']].drop_duplicates()
        residues2 = pocket2[['auth_chain_id', 'auth_seq_id']].drop_duplicates()
        
        # Print information for debugging
        print(f"Binding pocket info for {TEST_PDB_IDS[0]}: {len(pocket1)} atoms, {len(residues1)} residues")
        print(f"Binding pocket info for {TEST_PDB_IDS[1]}: {len(pocket2)} atoms, {len(residues2)} residues")
        
        # Since both are rhodopsins with retinal, they should have similar numbers of binding residues
        # This is somewhat fragile, but should work for similar proteins
        min_expected_residues = 5  # Minimum reasonable number
        assert len(residues1) >= min_expected_residues
        assert len(residues2) >= min_expected_residues