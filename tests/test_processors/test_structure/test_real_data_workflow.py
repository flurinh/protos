"""
Test complete workflows with real PDB data to ensure everything works end-to-end.

This test module validates that we can:
1. Download real PDB structures
2. Create datasets from real PDB IDs
3. Process and analyze the real data
4. Save and reload datasets
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.loaders.download_structures import download_protein_structures
from protos.io.paths.path_config import ProtosPaths


class TestRealDataWorkflow:
    """Test complete workflows using real PDB data."""
    
    @pytest.fixture
    def real_pdb_ids(self):
        """Return a list of real, small PDB structures for testing."""
        return [
            "1ubq",  # Ubiquitin - small, well-studied protein
            "1l2y",  # Trp-cage miniprotein - very small
            "1crn",  # Crambin - small protein
            "2gb1",  # Protein G B1 domain - small, stable
            "1a3n",  # Apo-myoglobin - medium size
        ]
    
    @pytest.fixture
    def workflow_processor(self, test_data_root):
        """Create a processor for workflow testing."""
        # Use the test-data directory from conftest
        processor = CifBaseProcessor(
            name="workflow_test",
            processor_data_dir="structure"
        )
        
        yield processor
    
    def test_download_and_create_dataset(self, workflow_processor, real_pdb_ids):
        """Test downloading real structures and creating a dataset."""
        # Step 1: Download real structures
        # Use the full path to the structure directory which already includes mmcif
        target_dir = workflow_processor.path_structure_dir
        successful, failed = download_protein_structures(
            real_pdb_ids,
            target_folder=str(target_dir)
        )
        
        # Verify downloads
        assert len(successful) >= 3, f"Should download at least 3 structures, got {len(successful)}"
        assert len(failed) <= 2, f"Should have at most 2 failures, got {len(failed)}"
        
        # Step 2: Create dataset from downloaded structures
        dataset_id = "real_proteins"
        workflow_processor.create_dataset(
            dataset_id=dataset_id,
            name="Real Protein Structures",
            description="Small real protein structures for testing",
            content=successful  # Use only successfully downloaded
        )
        
        # Step 3: Load the dataset
        workflow_processor.load_dataset(dataset_id)
        
        # Step 4: Verify real data was loaded
        assert workflow_processor.data is not None
        assert not workflow_processor.data.empty
        assert len(workflow_processor.pdb_ids) == len(successful)
        
        # Verify we have real protein data
        for pdb_id in workflow_processor.pdb_ids:
            pdb_data = workflow_processor.data[workflow_processor.data['pdb_id'] == pdb_id]
            
            # Check we have atoms
            assert len(pdb_data) > 0, f"No atoms found for {pdb_id}"
            
            # Check we have expected columns
            expected_cols = ['pdb_id', 'res_atom_name', 'res_name3l', 'auth_chain_id', 
                           'auth_seq_id', 'x', 'y', 'z']
            for col in expected_cols:
                assert col in pdb_data.columns, f"Missing column {col} for {pdb_id}"
            
            # Check coordinates are numeric
            assert pdb_data['x'].dtype in [np.float64, np.float32], "X coordinates should be numeric"
            assert pdb_data['y'].dtype in [np.float64, np.float32], "Y coordinates should be numeric"
            assert pdb_data['z'].dtype in [np.float64, np.float32], "Z coordinates should be numeric"
            
            # Check we have different atom types
            atom_types = pdb_data['res_atom_name'].unique()
            assert 'CA' in atom_types, f"Should have CA atoms in {pdb_id}"
            assert len(atom_types) > 5, f"Should have multiple atom types in {pdb_id}"
            
            # Check we have residues
            residues = pdb_data['res_name3l'].unique()
            assert len(residues) > 10, f"Should have multiple residues in {pdb_id}"
    
    def test_analyze_real_protein_structure(self, workflow_processor, real_pdb_ids):
        """Test analyzing real protein structures."""
        # Download a specific well-studied protein
        pdb_id = "1ubq"  # Ubiquitin
        target_dir = Path(workflow_processor.structure_dir) / "mmcif"
        successful, failed = download_protein_structures(
            [pdb_id],
            target_folder=str(target_dir)
        )
        
        assert len(successful) == 1, "Should successfully download ubiquitin"
        
        # Load the structure
        workflow_processor.load_structure(pdb_id)
        
        # Analyze the structure
        chains = workflow_processor.get_chains(pdb_id)
        assert 'A' in chains, "Ubiquitin should have chain A"
        
        # Get sequence
        sequences = workflow_processor.get_seq_dict()
        assert f"{pdb_id}_A" in sequences
        
        # Ubiquitin sequence should start with MQI...
        seq = sequences[f"{pdb_id}_A"]
        assert seq.startswith("MQI"), f"Ubiquitin sequence should start with MQI, got {seq[:10]}"
        assert len(seq) == 76, f"Ubiquitin should have 76 residues, got {len(seq)}"
        
        # Get CA coordinates
        ca_coords = workflow_processor.get_ca_coordinates(pdb_id, 'A')
        assert ca_coords is not None, "Should get CA coordinates"
        assert ca_coords.shape == (76, 3), "Should have 76 CA atoms with xyz coordinates"
        
        # Calculate center of mass
        com = np.mean(ca_coords, axis=0)
        assert len(com) == 3, "Center of mass should have 3 coordinates"
        assert all(np.isfinite(com)), "Center of mass should be finite"
    
    def test_dataset_persistence(self, workflow_processor, real_pdb_ids):
        """Test that datasets can be saved and reloaded correctly."""
        # Download and create dataset
        target_dir = Path(workflow_processor.structure_dir) / "mmcif"
        successful, _ = download_protein_structures(
            real_pdb_ids[:2],  # Just first 2 for speed
            target_folder=str(target_dir)
        )
        
        dataset_id = "persistence_test"
        workflow_processor.create_dataset(
            dataset_id=dataset_id,
            name="Persistence Test",
            description="Test dataset persistence",
            content=successful
        )
        
        # Load and get statistics
        workflow_processor.load_dataset(dataset_id)
        original_count = len(workflow_processor.data)
        original_pdbs = set(workflow_processor.pdb_ids)
        
        # Create new processor instance
        new_processor = CifBaseProcessor(
            name="workflow_test",
            processor_data_dir="structure"
        )
        
        # Load the same dataset
        new_processor.load_dataset(dataset_id)
        
        # Verify data matches
        assert len(new_processor.data) == original_count
        assert set(new_processor.pdb_ids) == original_pdbs
        
        # Verify data integrity
        for pdb_id in original_pdbs:
            orig_data = workflow_processor.data[workflow_processor.data['pdb_id'] == pdb_id]
            new_data = new_processor.data[new_processor.data['pdb_id'] == pdb_id]
            
            # Check row counts match
            assert len(orig_data) == len(new_data), f"Row count mismatch for {pdb_id}"
            
            # Check key values match
            assert orig_data['atom_name'].nunique() == new_data['atom_name'].nunique()
            assert orig_data['auth_seq_id'].nunique() == new_data['auth_seq_id'].nunique()
    
    def test_multi_chain_protein_analysis(self, workflow_processor):
        """Test analyzing a multi-chain protein complex."""
        # Download a protein complex with multiple chains
        pdb_id = "1a3n"  # Has multiple chains
        successful, _ = download_protein_structures(
            [pdb_id],
            target_folder=workflow_processor.structure_dir
        )
        
        if not successful:
            pytest.skip(f"Could not download {pdb_id}")
        
        workflow_processor.load_structure(pdb_id)
        
        # Get all chains
        chains = workflow_processor.get_chains(pdb_id)
        assert len(chains) >= 1, "Should have at least one chain"
        
        # Analyze each chain
        for chain in chains:
            chain_data = workflow_processor.data[workflow_processor.data['auth_chain_id'] == chain]
            assert len(chain_data) > 0, f"Chain {chain} should have atoms"
            
            # Get chain sequence
            sequences = workflow_processor.get_seq_dict()
            chain_key = f"{pdb_id}_{chain}"
            if chain_key in sequences:
                seq = sequences[chain_key]
                assert len(seq) > 0, f"Chain {chain} should have a sequence"
    
    def test_filter_and_analyze_binding_site(self, workflow_processor):
        """Test filtering data to analyze a binding site."""
        # Use a protein known to have ligands
        pdb_id = "1a3n"
        successful, _ = download_protein_structures(
            [pdb_id],
            target_folder=workflow_processor.structure_dir
        )
        
        if not successful:
            pytest.skip(f"Could not download {pdb_id}")
        
        workflow_processor.load_structure(pdb_id)
        
        # Find potential binding pocket around a point
        # This is a simplified version - real binding site detection is more complex
        center = [10.0, 10.0, 10.0]  # Arbitrary point
        radius = 8.0
        
        # Filter atoms within radius of center
        data = workflow_processor.data.copy()
        distances = np.sqrt(
            (data['x'] - center[0])**2 + 
            (data['y'] - center[1])**2 + 
            (data['z'] - center[2])**2
        )
        
        pocket_atoms = data[distances <= radius]
        
        if len(pocket_atoms) > 0:
            # Analyze pocket composition
            pocket_residues = pocket_atoms['res_name3l'].unique()
            pocket_chains = pocket_atoms['auth_chain_id'].unique()
            
            assert len(pocket_residues) > 0, "Should have residues in pocket"
            assert len(pocket_chains) > 0, "Should have chains in pocket"
            
            # Check atom types in pocket
            atom_types = pocket_atoms['atom_name'].unique()
            assert len(atom_types) > 0, "Should have various atom types in pocket"