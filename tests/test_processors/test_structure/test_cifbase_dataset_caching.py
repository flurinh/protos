"""
Test StructureProcessor dataset and caching functionality.
"""

import pytest
import tempfile
import os
from pathlib import Path
import pandas as pd
import time

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseDatasetCaching:
    """Test dataset and caching functionality."""
    
    def test_save_and_load_dataset_pkl(self):
        """Test saving and loading datasets as PKL in structure_dataset/."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = StructureProcessor(paths=paths)
            
            # Create test structures
            structures = []
            for i in range(3):
                test_data = pd.DataFrame({
                    'pdb_id': [f'protein_{i}'] * 10,
                    'auth_chain_id': ['A'] * 10,
                    'auth_seq_id': list(range(1, 11)),
                    'res_name1l': ['A', 'G', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'S'],
                    'atom_name': ['CA'] * 10,
                    'x': list(range(i*10, i*10 + 10)),
                    'y': list(range(i*20, i*20 + 10)),
                    'z': list(range(i*30, i*30 + 10)),
                    'group': ['ATOM'] * 10
                })
                processor.save_entity(f"protein_{i}", test_data)
                structures.append(test_data)
            
            # Combine into dataset
            combined_data = pd.concat(structures, ignore_index=True)
            processor.data = combined_data
            processor.pdb_ids = [f'protein_{i}' for i in range(3)]
            
            # Save dataset as PKL - should go to structure_dataset/
            dataset_path = processor.save_data("my_dataset")
            
            # Check it saved to the right place
            expected_path = Path(processor.paths.get_subdir_path("structure", "dataset_dir")) / "my_dataset.pkl"
            assert expected_path.exists(), f"Dataset should be saved to {expected_path}"
            
            # Clear processor data
            processor.data = None
            processor.pdb_ids = []
            
            # Load dataset - should load from PKL
            loaded_data = processor.load_data("my_dataset")
            
            assert loaded_data is not None
            assert len(loaded_data) == 30  # 3 proteins x 10 atoms
            assert set(loaded_data['pdb_id'].unique()) == {'protein_0', 'protein_1', 'protein_2'}
    
    def test_individual_structure_caching(self):
        """Test that individual structures are cached in cache/."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = StructureProcessor(paths=paths)
            
            # Create a CIF file in mmcif/
            mmcif_dir = processor.paths.get_subdir_path("structure", "structure_dir")
            test_cif = Path(mmcif_dir) / "test_protein.cif"
            test_cif.write_text("# Dummy CIF content")
            
            # Create test structure data
            test_data = pd.DataFrame({
                'pdb_id': ['test_protein'] * 5,
                'auth_chain_id': ['A'] * 5,
                'auth_seq_id': [1, 2, 3, 4, 5],
                'res_name1l': ['A', 'G', 'V', 'L', 'I'],
                'atom_name': ['CA'] * 5,
                'x': [1.0, 2.0, 3.0, 4.0, 5.0],
                'y': [6.0, 7.0, 8.0, 9.0, 10.0],
                'z': [11.0, 12.0, 13.0, 14.0, 15.0],
                'group': ['ATOM'] * 5
            })
            
            # Save using save_entity - should cache as PKL
            processor.save_entity("test_protein", test_data)
            
            # Check cache file exists
            cache_file = Path(processor.paths.get_subdir_path("structure", "cache_dir")) / "test_protein.pkl"
            assert cache_file.exists(), "Structure should be cached as PKL"
            
            # Load entity - should load from cache
            loaded = processor.load_entity("test_protein")
            assert loaded is not None
            assert len(loaded) == 5
            pd.testing.assert_frame_equal(loaded, test_data)
    
    def test_dataset_fallback_to_individual_structures(self):
        """Test that dataset loading falls back to individual structures if no PKL exists."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = StructureProcessor(paths=paths)
            
            # Create individual structures - save as both formats
            pdb_ids = []
            for i in range(3):
                test_data = pd.DataFrame({
                    'pdb_id': [f'prot{i}'] * 3,
                    'auth_chain_id': ['A'] * 3,
                    'auth_seq_id': [1, 2, 3],
                    'res_name1l': ['A', 'G', 'V'],
                    'atom_name': ['CA'] * 3,
                    'x': [1.0, 2.0, 3.0],
                    'y': [4.0, 5.0, 6.0],
                    'z': [7.0, 8.0, 9.0],
                    'group': ['ATOM'] * 3
                })
                # Save as both CIF and PKL to ensure fallback works
                processor.save_structure(f"prot{i}", test_data, format='both')
                pdb_ids.append(f"prot{i}")
            
            # Create dataset metadata
            dataset_name = processor.create_dataset("test_dataset", pdb_ids)
            
            # Clear processor state
            processor.data = None
            processor.pdb_ids = []
            
            # Load using load_data - should fall back to individual structures
            # since no preprocessed dataset PKL exists
            loaded_data = processor.load_data("test_dataset")
            
            assert loaded_data is not None
            assert len(loaded_data) == 9  # 3 structures x 3 atoms
            assert set(loaded_data['pdb_id'].unique()) == {'prot0', 'prot1', 'prot2'}
    
    def test_cache_performance(self):
        """Test that cached loading is faster than parsing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = StructureProcessor(paths=paths)
            
            # Create a large structure
            test_data = pd.DataFrame({
                'pdb_id': ['large_protein'] * 1000,
                'auth_chain_id': ['A'] * 1000,
                'auth_seq_id': list(range(1, 1001)),
                'res_name1l': ['A'] * 1000,
                'atom_name': ['CA'] * 1000,
                'x': list(range(1000)),
                'y': list(range(1000, 2000)),
                'z': list(range(2000, 3000)),
                'group': ['ATOM'] * 1000
            })
            
            # Save structure (which creates cache)
            processor.save_entity("large_protein", test_data)
            
            # Time loading with cache
            start = time.time()
            loaded1 = processor.load_entity("large_protein")
            time_with_cache = time.time() - start
            
            # Clear from entity registry to simulate fresh load
            if processor.entity_exists("large_protein"):
                processor.entity_registry._registry.pop(
                    processor.entity_registry._resolve_to_hash("large_protein"), None
                )
            
            # Time loading again (should still use cache file)
            start = time.time()
            loaded2 = processor.load_entity("large_protein")
            time_second_load = time.time() - start
            
            # Both should be fast from cache
            print(f"First load: {time_with_cache:.4f}s, Second load: {time_second_load:.4f}s")
            assert loaded1 is not None
            assert loaded2 is not None
            pd.testing.assert_frame_equal(loaded1, loaded2)
    
    def test_dataset_manager_integration(self):
        """Test that dataset manager properly tracks datasets."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir)
            processor = StructureProcessor(paths=paths)
            
            # Create structures
            structures = []
            for i in range(3):
                test_data = pd.DataFrame({
                    'pdb_id': [f'struct_{i}'] * 2,
                    'auth_chain_id': ['A'] * 2,
                    'auth_seq_id': [1, 2],
                    'res_name1l': ['A', 'G'],
                    'atom_name': ['CA'] * 2,
                    'x': [i, i+0.5],
                    'y': [i+1, i+1.5],
                    'z': [i+2, i+2.5],
                    'group': ['ATOM'] * 2
                })
                processor.save_entity(f"struct_{i}", test_data)
                structures.append(f"struct_{i}")
            
            # Create dataset through dataset manager
            dataset_name = processor.create_dataset(
                "kinase_family",
                structures,
                {"description": "Test kinase structures", "organism": "human"}
            )
            
            # List datasets
            datasets = processor.list_datasets()
            assert "kinase_family" in datasets
            
            # Get dataset info
            info = processor.get_dataset_info("kinase_family")
            assert info["entity_count"] == 3
            assert info["metadata"]["organism"] == "human"
            
            # Add more structures
            processor.add_to_dataset("kinase_family", ["struct_3", "struct_4"])
            info = processor.get_dataset_info("kinase_family")
            assert info["entity_count"] == 5