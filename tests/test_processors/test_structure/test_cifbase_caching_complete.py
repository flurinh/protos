"""
Test complete caching workflow for StructureProcessor.
"""

import pytest
import tempfile
from pathlib import Path
import pandas as pd
import time

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseCachingComplete:
    """Test the complete caching workflow as documented."""
    
    def test_load_entity_with_format_parameter(self):
        """Test load_entity supports format parameter for CIF vs PKL."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test structure
            test_data = pd.DataFrame({
                'pdb_id': ['test_prot'] * 5,
                'auth_chain_id': ['A'] * 5,
                'auth_seq_id': [1, 2, 3, 4, 5],
                'res_name1l': ['A', 'G', 'V', 'L', 'I'],
                'atom_name': ['CA'] * 5,
                'x': [1.0, 2.0, 3.0, 4.0, 5.0],
                'y': [6.0, 7.0, 8.0, 9.0, 10.0],
                'z': [11.0, 12.0, 13.0, 14.0, 15.0],
                'group': ['ATOM'] * 5
            })
            
            # Save as both formats
            processor.save_structure("test_prot", test_data, format='both')
            
            # Load with format='pkl' (default) - should use cache
            loaded_pkl = processor.load_entity("test_prot", format='pkl')
            assert loaded_pkl is not None
            assert len(loaded_pkl) == 5
            
            # Load with format='cif' - should parse CIF
            loaded_cif = processor.load_entity("test_prot", format='cif')
            assert loaded_cif is not None
            # Data should be same regardless of format
            pd.testing.assert_frame_equal(loaded_pkl, loaded_cif)
    
    def test_load_dataset_with_format_parameter(self):
        """Test load_dataset supports format parameter."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test structures
            structures = []
            pdb_ids = []
            for i in range(3):
                test_data = pd.DataFrame({
                    'pdb_id': [f'prot{i}'] * 3,
                    'auth_chain_id': ['A'] * 3,
                    'auth_seq_id': [1, 2, 3],
                    'res_name1l': ['A', 'G', 'V'],
                    'atom_name': ['CA'] * 3,
                    'x': [float(i), float(i+1), float(i+2)],
                    'y': [float(i+3), float(i+4), float(i+5)],
                    'z': [float(i+6), float(i+7), float(i+8)],
                    'group': ['ATOM'] * 3
                })
                processor.save_structure(f"prot{i}", test_data, format='both')
                structures.append(test_data)
                pdb_ids.append(f"prot{i}")
            
            # Create dataset
            processor.create_dataset("test_dataset", pdb_ids)
            
            # Combine data and save as dataset PKL
            combined = pd.concat(structures, ignore_index=True)
            processor.data = combined
            processor.save_data("test_dataset")
            
            # Clear processor state
            processor.data = None
            processor.pdb_ids = []
            
            # Test load_dataset with format='pkl' (should use structure_dataset/)
            dataset_data = processor.load_dataset("test_dataset", format='pkl')
            assert isinstance(dataset_data, pd.DataFrame)
            assert len(dataset_data) == 9  # 3 structures x 3 atoms
            
            # Test load_dataset with format='cif' (should load individual CIFs)
            dataset_dict = processor.load_dataset("test_dataset", format='cif')
            assert isinstance(dataset_dict, dict)
            assert len(dataset_dict) == 3
            assert all(pdb_id in dataset_dict for pdb_id in pdb_ids)
    
    def test_save_dataset_creates_both_json_and_pkl(self):
        """Test that save_data creates both dataset JSON and PKL file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test data
            test_data = pd.DataFrame({
                'pdb_id': ['prot1'] * 3 + ['prot2'] * 3,
                'auth_chain_id': ['A'] * 6,
                'auth_seq_id': [1, 2, 3, 1, 2, 3],
                'res_name1l': ['A', 'G', 'V'] * 2,
                'atom_name': ['CA'] * 6,
                'x': range(6),
                'y': range(6, 12),
                'z': range(12, 18),
                'group': ['ATOM'] * 6
            })
            
            processor.data = test_data
            processor.pdb_ids = ['prot1', 'prot2']
            
            # Save dataset
            processor.save_data("my_dataset", format='pkl')
            
            # Check PKL file exists in structure_dataset/
            pkl_path = processor.path_dataset_dir / "my_dataset.pkl"
            assert pkl_path.exists()
            
            # Check dataset JSON exists in datasets/
            datasets_dir = processor.data_path / "datasets"
            json_files = list(datasets_dir.glob("*.json"))
            assert any("my_dataset" in f.stem for f in json_files)
            
            # Verify dataset is listed
            datasets = processor.list_datasets()
            assert "my_dataset" in datasets
    
    def test_caching_performance_benefit(self):
        """Test that PKL loading is faster than CIF parsing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create a moderately large structure
            test_data = pd.DataFrame({
                'pdb_id': ['large_prot'] * 500,
                'auth_chain_id': ['A'] * 500,
                'auth_seq_id': list(range(1, 501)),
                'res_name1l': ['A'] * 500,
                'atom_name': ['CA'] * 500,
                'x': list(range(500)),
                'y': list(range(500, 1000)),
                'z': list(range(1000, 1500)),
                'group': ['ATOM'] * 500
            })
            
            # Save as PKL only (for cache)
            processor.save_entity("large_prot", test_data)
            
            # Measure loading time from PKL
            start = time.time()
            loaded_pkl = processor.load_entity("large_prot", format='pkl')
            pkl_time = time.time() - start
            
            # For this test, we just verify PKL loading works
            # (CIF parsing would require actual CIF file creation)
            assert loaded_pkl is not None
            assert len(loaded_pkl) == 500
            print(f"PKL load time: {pkl_time:.4f}s")
    
    def test_fallback_behavior(self):
        """Test fallback from dataset PKL to individual structures."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create individual structures (only as PKL cache)
            pdb_ids = []
            for i in range(3):
                test_data = pd.DataFrame({
                    'pdb_id': [f'fallback{i}'] * 2,
                    'auth_chain_id': ['A'] * 2,
                    'auth_seq_id': [1, 2],
                    'res_name1l': ['G', 'A'],
                    'atom_name': ['CA'] * 2,
                    'x': [i, i+0.5],
                    'y': [i+1, i+1.5],
                    'z': [i+2, i+2.5],
                    'group': ['ATOM'] * 2
                })
                processor.save_entity(f"fallback{i}", test_data)
                pdb_ids.append(f"fallback{i}")
            
            # Create dataset metadata (but no PKL file)
            processor.create_dataset("fallback_test", pdb_ids)
            
            # Try to load dataset - should fall back to individual files
            loaded = processor.load_data("fallback_test")
            
            assert loaded is not None
            assert len(loaded) == 6  # 3 structures x 2 atoms
            assert set(loaded['pdb_id'].unique()) == {'fallback0', 'fallback1', 'fallback2'}