"""
Test StructureProcessor caching with real CIF files.

This test ensures that the caching functionality works with actual CIF files,
not just mock data.
"""

import pytest
import tempfile
import shutil
import pandas as pd
from pathlib import Path
import time

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseRealDataCaching:
    """Test caching functionality with real CIF files."""
    
    @pytest.fixture
    def real_cif_processor(self):
        """Create a processor with real CIF files."""
        # Create a temporary directory
        tmpdir = tempfile.mkdtemp()
        
        # Set up ProtosPaths
        paths = ProtosPaths(data_root=tmpdir)
        
        # Copy real CIF files to the temp directory
        source_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/tests/test-data/structure/mmcif")
        dest_dir = Path(tmpdir) / "structure" / "mmcif"
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy the CIF files
        cif_files = ["1uaz.cif", "3ddl.cif", "4pxk.cif"]
        for cif_file in cif_files:
            shutil.copy(source_dir / cif_file, dest_dir / cif_file)
        
        # Create processor
        processor = StructureProcessor(paths=paths)
        
        yield processor, tmpdir, cif_files
        
        # Cleanup
        shutil.rmtree(tmpdir)
    
    def test_load_real_cif_and_cache(self, real_cif_processor):
        """Test loading real CIF files and caching them."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Testing Real CIF Loading and Caching ===")
        
        # Test loading each CIF file
        for cif_name in ["1uaz", "3ddl", "4pxk"]:
            print(f"\nLoading {cif_name}...")
            
            # First load - should parse CIF and create cache
            # Use use_cache=False on first load to ensure we parse from CIF
            start_time = time.time()
            structure = processor.load_structure(cif_name, use_cache=False, save_processed=True)
            parse_time = time.time() - start_time
            
            assert structure is not None
            assert len(structure) > 0
            print(f"  - Loaded {len(structure)} atoms in {parse_time:.3f}s")
            print(f"  - Chains: {structure['auth_chain_id'].unique()}")
            print(f"  - Residues: {structure['auth_seq_id'].nunique()}")
            
            # Check cache was created
            cache_file = processor.paths.get_subdir_path("structure", "cache_dir") / f"{cif_name}.pkl"
            assert cache_file.exists()
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
    
    def test_create_dataset_with_real_structures(self, real_cif_processor):
        """Test creating a dataset from real structures."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Testing Dataset Creation with Real Structures ===")
        
        # Load all structures
        pdb_ids = ["1uaz", "3ddl", "4pxk"]
        for pdb_id in pdb_ids:
            processor.load_structure(pdb_id, save_processed=True)
        
        # Create dataset
        processor.create_dataset("test_structures", pdb_ids, metadata={
            "description": "Test dataset with real PDB structures",
            "organism": "various",
            "method": "X-ray crystallography"
        })
        
        # Verify dataset was created
        datasets = processor.list_datasets()
        assert "test_structures" in datasets
        print(f"  - Dataset created: test_structures")
        
        # Get dataset info
        info = processor.get_dataset_info("test_structures")
        print(f"  - Dataset contains {info['entity_count']} structures")
        entity_names = [e['name'] for e in info['entities']]
        assert set(entity_names) == set(pdb_ids)
    
    def test_save_and_load_dataset_pkl(self, real_cif_processor):
        """Test saving and loading dataset as PKL."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Testing Dataset PKL Save/Load ===")
        
        # Load structures using load_structures
        pdb_ids = ["1uaz", "3ddl", "4pxk"]
        processor.load_structures(pdb_ids)
        
        # Check that data was loaded
        assert processor.data is not None
        total_atoms = len(processor.data)
        print(f"  - Loaded {total_atoms} total atoms from {len(pdb_ids)} structures")
        
        # Save as dataset PKL
        dataset_name = "combined_structures"
        processor.create_dataset(dataset_name, pdb_ids)
        
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
        assert loaded_pdb_ids == set(pdb_ids)
    
    def test_fallback_loading_without_pkl(self, real_cif_processor):
        """Test that load_data falls back to individual files when PKL is missing."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Testing Fallback Loading ===")
        
        # Load and cache individual structures
        pdb_ids = ["1uaz", "3ddl"]
        for pdb_id in pdb_ids:
            processor.load_structure(pdb_id, save_processed=True)
        
        # Create dataset but don't save PKL
        processor.create_dataset("fallback_dataset", pdb_ids)
        
        # Ensure no PKL exists
        pkl_path = Path(processor.paths.get_subdir_path("structure", "dataset_dir")) / "fallback_dataset.pkl"
        assert not pkl_path.exists()
        
        # Load dataset - should fall back to individual files
        processor.data = None
        start_time = time.time()
        loaded = processor.load_data("fallback_dataset")
        fallback_time = time.time() - start_time
        
        assert loaded is not None
        assert set(loaded['pdb_id'].unique()) == set(pdb_ids)
        print(f"  - Loaded {len(loaded)} atoms via fallback in {fallback_time:.3f}s")
        
        # Now save as PKL
        processor.save_data("fallback_dataset")
        assert pkl_path.exists()
        
        # Load again - should use PKL (faster)
        processor.data = None
        start_time = time.time()
        loaded2 = processor.load_data("fallback_dataset")
        pkl_time = time.time() - start_time
        
        assert loaded2 is not None
        print(f"  - Loaded from PKL in {pkl_time:.3f}s (speedup: {fallback_time/pkl_time:.1f}x)")
    
    def test_mixed_format_dataset(self, real_cif_processor):
        """Test dataset with both cached and uncached structures."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Testing Mixed Format Dataset ===")
        
        # Cache only some structures
        processor.load_structure("1uaz", save_processed=True)  # Cached
        processor.load_structure("3ddl", save_processed=True)  # Cached
        # 4pxk will not be cached initially
        
        # Create dataset with all three
        all_pdb_ids = ["1uaz", "3ddl", "4pxk"]
        processor.create_dataset("mixed_dataset", all_pdb_ids)
        
        # Load dataset - should handle mixed sources
        loaded = processor.load_data("mixed_dataset")
        assert loaded is not None
        assert set(loaded['pdb_id'].unique()) == set(all_pdb_ids)
        print(f"  - Successfully loaded mixed dataset with {len(loaded)} atoms")
        
        # Check that 4pxk was loaded but not necessarily cached
        # (load_structures doesn't cache by default)
        cache_file = Path(processor.paths.get_subdir_path("structure", "cache_dir")) / "4pxk.pkl"
        if cache_file.exists():
            print("  - Previously uncached structure (4pxk) is now cached")
        else:
            print("  - Structure 4pxk was loaded from CIF (not cached during dataset load)")
    
    def test_verify_data_integrity(self, real_cif_processor):
        """Verify that cached data matches original CIF data."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Testing Data Integrity ===")
        
        pdb_id = "1uaz"
        
        # Load from CIF (no cache)
        structure_cif = processor.load_structure(pdb_id, use_cache=False, save_processed=False)
        
        # Save to cache
        processor.save_entity(pdb_id, structure_cif)
        
        # Load from cache
        structure_cache = processor.load_entity(pdb_id)
        
        # Compare key columns
        assert len(structure_cif) == len(structure_cache)
        assert set(structure_cif.columns) == set(structure_cache.columns)
        
        # Check numeric precision
        for coord in ['x', 'y', 'z']:
            if coord in structure_cif.columns:
                max_diff = abs(structure_cif[coord] - structure_cache[coord]).max()
                assert max_diff < 1e-6
                print(f"  - Maximum {coord} coordinate difference: {max_diff}")
        
        print("  - Data integrity verified: cache matches CIF")
    
    def test_performance_comparison(self, real_cif_processor):
        """Compare performance of different loading methods."""
        processor, tmpdir, cif_files = real_cif_processor
        
        print("\n=== Performance Comparison ===")
        
        pdb_ids = ["1uaz", "3ddl", "4pxk"]
        
        # 1. Load from CIF files (no cache)
        processor.data = None
        start = time.time()
        for pdb_id in pdb_ids:
            processor.load_structure(pdb_id, use_cache=False, save_processed=True)
        cif_time = time.time() - start
        print(f"  - CIF parsing time: {cif_time:.3f}s")
        
        # 2. Load from individual cache files
        processor.data = None
        start = time.time()
        for pdb_id in pdb_ids:
            processor.load_entity(pdb_id)
        cache_time = time.time() - start
        print(f"  - Individual cache load time: {cache_time:.3f}s (speedup: {cif_time/cache_time:.1f}x)")
        
        # 3. Create and load dataset PKL
        processor.load_structures(pdb_ids)
        processor.create_dataset("perf_test", pdb_ids)
        processor.save_data("perf_test")
        
        processor.data = None
        start = time.time()
        processor.load_data("perf_test")
        dataset_time = time.time() - start
        print(f"  - Dataset PKL load time: {dataset_time:.3f}s (speedup: {cif_time/dataset_time:.1f}x)")
        
        print("\n  Summary:")
        print(f"  - Cache is {cif_time/cache_time:.1f}x faster than CIF parsing")
        print(f"  - Dataset PKL is {cif_time/dataset_time:.1f}x faster than CIF parsing")
        print(f"  - Dataset PKL is {cache_time/dataset_time:.1f}x faster than individual cache")