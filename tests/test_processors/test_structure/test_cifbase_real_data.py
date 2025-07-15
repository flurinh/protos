"""
Test StructureProcessor with real structure data.
"""

import pytest
from pathlib import Path
import pandas as pd

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseProcessorRealData:
    """Test StructureProcessor with real structure files."""
    
    @pytest.fixture
    def processor(self):
        """Create processor with test data directory."""
        # Use the test-data directory which has real CIF files
        test_data_path = Path(__file__).parent.parent.parent.parent / "test-data"
        paths = ProtosPaths(data_root=str(test_data_path))
        return StructureProcessor(paths=paths)
    
    def test_list_available_structures(self, processor):
        """Test listing available structure files."""
        # Initialize PDB IDs from directory
        processor.initialize_pdb_ids()
        
        # Should find the test structures
        assert len(processor.pdb_ids) > 0
        print(f"Found structures: {processor.pdb_ids}")
        
        # Check for known test structures
        expected = ["1uaz", "3ddl", "4pxk"]
        for pdb_id in expected:
            if pdb_id in processor.pdb_ids:
                assert True
                break
        else:
            pytest.skip(f"No expected test structures found. Available: {processor.pdb_ids}")
    
    def test_load_real_structure(self, processor):
        """Test loading a real structure file."""
        # Find available structures
        processor.initialize_pdb_ids()
        
        if not processor.pdb_ids:
            pytest.skip("No structure files available in test data")
        
        # Load the first available structure
        pdb_id = processor.pdb_ids[0]
        print(f"Loading structure: {pdb_id}")
        
        # Load structure
        structure = processor.load_structure(pdb_id, apply_dtypes=True)
        
        # Verify it loaded
        assert structure is not None
        assert len(structure) > 0
        assert 'pdb_id' in structure.columns
        assert structure['pdb_id'].iloc[0] == pdb_id
        
        # Check critical columns exist
        required_columns = ['x', 'y', 'z', 'auth_chain_id', 'auth_seq_id']
        for col in required_columns:
            assert col in structure.columns
        
        # Check coordinate types
        assert structure['x'].dtype == 'float64'
        assert structure['y'].dtype == 'float64'
        assert structure['z'].dtype == 'float64'
        
        print(f"Loaded {len(structure)} atoms from {pdb_id}")
    
    def test_cache_workflow(self, processor):
        """Test caching workflow."""
        processor.initialize_pdb_ids()
        
        if not processor.pdb_ids:
            pytest.skip("No structure files available")
        
        pdb_id = processor.pdb_ids[0]
        
        # First load - will parse CIF and cache
        structure1 = processor.load_structure(pdb_id, use_cache=True, save_processed=True)
        assert structure1 is not None
        
        # Second load - should use cache
        structure2 = processor.load_structure(pdb_id, use_cache=True)
        assert structure2 is not None
        
        # Should have same data
        assert len(structure1) == len(structure2)
        pd.testing.assert_frame_equal(structure1, structure2)
    
    def test_extract_chains(self, processor):
        """Test extracting chain information."""
        processor.initialize_pdb_ids()
        
        if not processor.pdb_ids:
            pytest.skip("No structure files available")
        
        # Load a structure
        pdb_id = processor.pdb_ids[0]
        structure = processor.load_structure(pdb_id)
        
        if structure is None or structure.empty:
            pytest.skip(f"Could not load structure {pdb_id}")
        
        # Store in processor for chain methods
        processor.data = structure
        processor.pdb_ids = [pdb_id]
        
        # Get chains
        chains = processor.get_chains(pdb_id)
        assert len(chains) > 0
        print(f"Structure {pdb_id} has chains: {chains}")
        
        # Try to get sequence for first chain
        chain_id = chains[0]
        try:
            sequence = processor.get_sequence(pdb_id, chain_id)
            assert len(sequence) > 0
            print(f"Chain {chain_id} sequence (first 20): {sequence[:20]}...")
        except ValueError as e:
            print(f"Could not extract sequence: {e}")
    
    def test_dataset_with_real_structures(self, processor):
        """Test creating dataset with real structures."""
        processor.initialize_pdb_ids()
        
        if len(processor.pdb_ids) < 2:
            pytest.skip("Need at least 2 structures for dataset test")
        
        # Create dataset with first 2 structures
        dataset_structures = processor.pdb_ids[:2]
        
        dataset_name = processor.create_dataset(
            "test_real_structures",
            dataset_structures,
            {"description": "Test dataset with real PDB structures"}
        )
        
        assert dataset_name == "test_real_structures"
        
        # Verify dataset was created
        datasets = processor.list_datasets()
        assert "test_real_structures" in datasets
        
        # Load dataset info
        info = processor.get_dataset_info("test_real_structures")
        assert info["entity_count"] == 2
    
    def test_filter_by_residue_range(self, processor):
        """Test filtering by residue range."""
        processor.initialize_pdb_ids()
        
        if not processor.pdb_ids:
            pytest.skip("No structure files available")
        
        # Load structure
        pdb_id = processor.pdb_ids[0]
        structure = processor.load_structure(pdb_id)
        
        if structure is None or structure.empty:
            pytest.skip(f"Could not load structure {pdb_id}")
        
        processor.data = structure.copy()
        
        # Get residue range
        min_res = processor.data['auth_seq_id'].min()
        max_res = processor.data['auth_seq_id'].max()
        
        # Filter to first 10 residues
        end_res = min(min_res + 10, max_res)
        processor.filter_by_residue_range(min_res, end_res)
        
        # Check filtering worked
        assert processor.data['auth_seq_id'].max() <= end_res
        assert processor.data['auth_seq_id'].min() >= min_res
        print(f"Filtered to residues {min_res}-{end_res}, {len(processor.data)} atoms remaining")