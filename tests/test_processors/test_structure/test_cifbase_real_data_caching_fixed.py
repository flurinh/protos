"""
Test StructureProcessor caching with real CIF files using test-data directory.

This test ensures that the caching functionality works with actual CIF files,
and saves cache files to the proper test-data directory.
"""

import pytest
import pandas as pd
from pathlib import Path
import time

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseRealDataCachingFixed:
    """Test caching functionality with real CIF files in test-data directory."""
    
    @pytest.fixture
    def processor(self):
        """Create a processor using the global test data directory."""
        # ProtosPaths will use the test-data directory configured in conftest
        processor = StructureProcessor()
        return processor
    
    def test_load_real_cif_and_cache(self, processor):
        """Test loading real CIF files and caching them."""
        print("\n=== Testing Real CIF Loading and Caching ===")
        
        # Test loading each CIF file that exists in test-data
        cif_files = ["1uaz", "3ddl", "4pxk"]
        
        for cif_name in cif_files:
            cif_path = Path(processor.paths.get_subdir_path("structure", "structure_dir")) / f"{cif_name}.cif"
            if not cif_path.exists():
                print(f"Skipping {cif_name} - CIF file not found")
                continue
                
            print(f"\nLoading {cif_name}...")
            
            # Clear any existing cache to test fresh
            cache_file = Path(processor.paths.get_subdir_path("structure", "cache_dir")) / f"{cif_name}.pkl"
            if cache_file.exists():
                cache_file.unlink()
            
            # First load - should parse CIF and create cache
            start_time = time.time()
            structure = processor.load_structure(cif_name, use_cache=False, save_processed=True)
            parse_time = time.time() - start_time
            
            assert structure is not None, f"Failed to load {cif_name}"
            assert len(structure) > 0, f"Empty structure for {cif_name}"
            print(f"  - Loaded {len(structure)} atoms in {parse_time:.3f}s")
            print(f"  - Chains: {structure['auth_chain_id'].unique()}")
            print(f"  - Residues: {structure['auth_seq_id'].nunique()}")
            
            # Check cache was created
            assert cache_file.exists(), f"Cache file not created for {cif_name}"
            print(f"  - Cache created: {cache_file}")
            
            # Second load - should use cache (much faster)
            start_time = time.time()
            structure2 = processor.load_structure(cif_name, use_cache=True)
            cache_time = time.time() - start_time
            
            assert structure2 is not None
            pd.testing.assert_frame_equal(structure, structure2)
            print(f"  - Cache load time: {cache_time:.3f}s (speedup: {parse_time/cache_time:.1f}x)")
            
            # Verify entity is registered
            assert processor.entity_exists(cif_name)
    
    def test_save_and_load_dataset_pkl(self, processor):
        """Test saving and loading dataset as PKL."""
        print("\n=== Testing Dataset PKL Save/Load ===")
        
        # Load structures using load_structures
        pdb_ids = ["1uaz", "3ddl", "4pxk"]
        available_ids = [pid for pid in pdb_ids 
                        if (Path(processor.paths.get_subdir_path("structure", "structure_dir")) / f"{pid}.cif").exists()]
        
        if not available_ids:
            pytest.skip("No CIF files available for testing")
        
        processor.load_structures(available_ids)
        
        # Check that data was loaded
        assert processor.data is not None
        total_atoms = len(processor.data)
        print(f"  - Loaded {total_atoms} total atoms from {len(available_ids)} structures")
        
        # Save as dataset PKL
        dataset_name = "test_combined_structures"
        processor.create_dataset(dataset_name, available_ids)
        
        save_start = time.time()
        pkl_path = processor.save_data(dataset_name)
        save_time = time.time() - save_start
        
        assert Path(pkl_path).exists()
        file_size_mb = Path(pkl_path).stat().st_size / (1024 * 1024)
        print(f"  - Saved dataset PKL: {pkl_path}")
        print(f"  - File size: {file_size_mb:.2f} MB")
        print(f"  - Save time: {save_time:.3f}s")
        
        # Clear processor state
        processor.data = None
        processor.pdb_ids = []
        
        # Load from PKL
        load_start = time.time()
        loaded_data = processor.load_data(dataset_name)
        load_time = time.time() - load_start
        
        assert loaded_data is not None
        assert len(loaded_data) == total_atoms
        print(f"  - Loaded {len(loaded_data)} atoms from PKL in {load_time:.3f}s")
        
        # Verify all PDB IDs are present
        loaded_pdb_ids = set(loaded_data['pdb_id'].unique())
        assert loaded_pdb_ids == set(available_ids)
    
    def test_verify_cache_persistence(self, processor):
        """Verify that cache files persist in test-data directory."""
        print("\n=== Testing Cache Persistence ===")
        
        # List current cache contents
        cache_files = list(processor.paths.get_subdir_path("structure", "cache_dir").glob("*.pkl"))
        print(f"\nCache directory: {processor.paths.get_subdir_path('structure', 'cache_dir')}")
        print(f"Cache files found: {len(cache_files)}")
        for cache_file in cache_files:
            print(f"  - {cache_file.name} ({cache_file.stat().st_size} bytes)")
        
        # List dataset contents
        dataset_files = list(processor.paths.get_subdir_path("structure", "dataset_dir").glob("*.pkl"))
        print(f"\nDataset directory: {processor.paths.get_subdir_path('structure', 'dataset_dir')}")
        print(f"Dataset files found: {len(dataset_files)}")
        for dataset_file in dataset_files:
            print(f"  - {dataset_file.name} ({dataset_file.stat().st_size} bytes)")