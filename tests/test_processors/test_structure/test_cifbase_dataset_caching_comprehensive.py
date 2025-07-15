"""
Comprehensive tests for StructureProcessor dataset caching functionality.

Tests the two-tier caching system:
1. Individual structure cache in cache/
2. Dataset cache in structure_dataset/
"""

import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
import time

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


class TestCifBaseDatasetCaching:
    """Test the complete dataset caching workflow for StructureProcessor."""
    
    def test_individual_structure_caching(self):
        """Test that individual structures are cached as PKL files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test structure
            test_data = pd.DataFrame({
                'pdb_id': ['test1'] * 10,
                'group': ['ATOM'] * 10,
                'auth_chain_id': ['A'] * 10,
                'auth_seq_id': list(range(1, 11)),
                'res_name3l': ['ALA'] * 10,
                'res_name1l': ['A'] * 10,
                'atom_name': ['CA'] * 10,
                'x': np.arange(10, dtype=float),
                'y': np.arange(10, 20, dtype=float),
                'z': np.arange(20, 30, dtype=float)
            })
            
            # Save structure using save_entity
            processor.save_entity("test1", test_data, metadata={
                "test": True,
                "atom_count": len(test_data)
            })
            
            # Check cache file exists
            cache_file = processor.path_cache_dir / "test1.pkl"
            assert cache_file.exists()
            
            # Load through load_entity - should use cache
            loaded = processor.load_entity("test1")
            assert loaded is not None
            assert len(loaded) == 10
            pd.testing.assert_frame_equal(loaded, test_data)
            
            # Verify entity is registered
            assert processor.entity_exists("test1")
    
    def test_dataset_pkl_saving(self):
        """Test saving entire datasets as PKL files in structure_dataset/."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test dataset with multiple structures
            structures = []
            pdb_ids = ['prot1', 'prot2', 'prot3']
            
            for i, pdb_id in enumerate(pdb_ids):
                test_data = pd.DataFrame({
                    'pdb_id': [pdb_id] * 5,
                    'group': ['ATOM'] * 5,
                    'auth_chain_id': ['A'] * 5,
                    'auth_seq_id': list(range(1, 6)),
                    'res_name3l': ['ALA', 'GLY', 'VAL', 'LEU', 'ILE'],
                    'res_name1l': ['A', 'G', 'V', 'L', 'I'],
                    'atom_name': ['CA'] * 5,
                    'x': np.arange(i*5, (i+1)*5, dtype=float),
                    'y': np.arange(i*5+10, (i+1)*5+10, dtype=float),
                    'z': np.arange(i*5+20, (i+1)*5+20, dtype=float)
                })
                structures.append(test_data)
                # Save individual structure to cache
                processor.save_entity(pdb_id, test_data)
            
            # Combine into dataset
            combined_data = pd.concat(structures, ignore_index=True)
            processor.data = combined_data
            processor.pdb_ids = pdb_ids
            
            # Save dataset as PKL
            processor.save_data("test_dataset")
            
            # Check PKL file exists in structure_dataset/
            dataset_pkl = processor.path_dataset_dir / "test_dataset.pkl"
            assert dataset_pkl.exists()
            
            # Also create dataset metadata
            processor.create_dataset("test_dataset", pdb_ids)
            
            # Verify dataset is listed
            datasets = processor.list_datasets()
            assert "test_dataset" in datasets
    
    def test_load_data_priority(self):
        """Test that load_data prioritizes PKL files over individual structures."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create individual structures
            pdb_ids = ['struct1', 'struct2']
            structures = []
            for i, pdb_id in enumerate(pdb_ids):
                test_data = pd.DataFrame({
                    'pdb_id': [pdb_id] * 3,
                    'group': ['ATOM'] * 3,
                    'auth_chain_id': ['A'] * 3,
                    'auth_seq_id': [1, 2, 3],
                    'res_name3l': ['ALA', 'GLY', 'VAL'],
                    'res_name1l': ['A', 'G', 'V'],
                    'atom_name': ['CA'] * 3,
                    'x': [float(i), float(i+1), float(i+2)],
                    'y': [float(i+3), float(i+4), float(i+5)],
                    'z': [float(i+6), float(i+7), float(i+8)]
                })
                processor.save_entity(pdb_id, test_data)
                structures.append(test_data)
            
            # Create dataset metadata
            processor.create_dataset("priority_test", pdb_ids)
            
            # Verify dataset was created
            datasets = processor.list_datasets()
            assert "priority_test" in datasets
            
            # First, manually load and combine structures to establish the baseline
            combined = pd.concat(structures, ignore_index=True)
            processor.data = combined
            processor.pdb_ids = pdb_ids
            
            # Load without PKL - should load individual files through load_structures
            processor.data = None  # Clear data
            loaded1 = processor.load_data("priority_test")
            assert loaded1 is not None
            assert len(loaded1) == 6  # 2 structures x 3 atoms
            
            # Now save as dataset PKL
            processor.save_data("priority_test")
            
            # Verify PKL was created
            pkl_path = processor.path_dataset_dir / "priority_test.pkl"
            assert pkl_path.exists()
            
            # Clear processor state
            processor.data = None
            processor.pdb_ids = []
            
            # Load again - should use PKL (faster)
            loaded2 = processor.load_data("priority_test")
            assert loaded2 is not None
            assert len(loaded2) == 6
            
            # Data should be identical
            pd.testing.assert_frame_equal(loaded1.sort_values(['pdb_id', 'auth_seq_id']).reset_index(drop=True),
                                        loaded2.sort_values(['pdb_id', 'auth_seq_id']).reset_index(drop=True))
    
    def test_save_data_formats(self):
        """Test save_data supports both PKL and CSV formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test data
            test_data = pd.DataFrame({
                'pdb_id': ['format_test'] * 4,
                'group': ['ATOM'] * 4,
                'auth_chain_id': ['A'] * 4,
                'auth_seq_id': [1, 2, 3, 4],
                'res_name3l': ['ALA', 'GLY', 'VAL', 'LEU'],
                'res_name1l': ['A', 'G', 'V', 'L'],
                'atom_name': ['CA'] * 4,
                'x': [1.0, 2.0, 3.0, 4.0],
                'y': [5.0, 6.0, 7.0, 8.0],
                'z': [9.0, 10.0, 11.0, 12.0]
            })
            
            processor.data = test_data
            
            # Save as PKL
            pkl_path = processor.save_data("format_test_pkl", file_format="pkl")
            assert Path(pkl_path).exists()
            assert Path(pkl_path).suffix == '.pkl'
            
            # Save as CSV
            csv_path = processor.save_data("format_test_csv", file_format="csv")
            assert Path(csv_path).exists()
            assert Path(csv_path).suffix == '.csv'
            
            # Load both formats
            loaded_pkl = processor.load_data("format_test_pkl", file_format="pkl")
            loaded_csv = processor.load_data("format_test_csv", file_format="csv")
            
            assert loaded_pkl is not None
            assert loaded_csv is not None
            assert len(loaded_pkl) == len(loaded_csv) == 4
    
    def test_fallback_when_pkl_missing(self):
        """Test fallback to individual structures when dataset PKL is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create individual structures
            pdb_ids = ['fall1', 'fall2', 'fall3']
            for pdb_id in pdb_ids:
                test_data = pd.DataFrame({
                    'pdb_id': [pdb_id] * 2,
                    'group': ['ATOM'] * 2,
                    'auth_chain_id': ['A'] * 2,
                    'auth_seq_id': [1, 2],
                    'res_name3l': ['ALA', 'GLY'],
                    'res_name1l': ['A', 'G'],
                    'atom_name': ['CA'] * 2,
                    'x': [1.0, 2.0],
                    'y': [3.0, 4.0],
                    'z': [5.0, 6.0]
                })
                processor.save_entity(pdb_id, test_data)
            
            # Create dataset metadata WITHOUT creating PKL
            processor.create_dataset("fallback_test", pdb_ids)
            
            # Ensure no PKL exists
            pkl_path = processor.path_dataset_dir / "fallback_test.pkl"
            assert not pkl_path.exists()
            
            # Load dataset - should fall back to individual files
            loaded = processor.load_data("fallback_test")
            assert loaded is not None
            assert len(loaded) == 6  # 3 structures x 2 atoms
            assert set(loaded['pdb_id'].unique()) == set(pdb_ids)
    
    def test_load_structure_caching_workflow(self):
        """Test the complete load_structure workflow with caching."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Ensure cache doesn't exist
            cache_file = processor.path_cache_dir / "test_cif.pkl"
            if cache_file.exists():
                cache_file.unlink()
            
            # Create a mock CIF file
            cif_content = """data_test_cif
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1  N  N   . ALA A 1 1 ? 1.000 2.000 3.000 1.00 10.00 ? 1 ALA A N   1
ATOM   2  CA CA  . ALA A 1 1 ? 2.000 3.000 4.000 1.00 10.00 ? 1 ALA A CA  1
ATOM   3  C  C   . ALA A 1 1 ? 3.000 4.000 5.000 1.00 10.00 ? 1 ALA A C   1
"""
            cif_file = processor.path_structure_dir / "test_cif.cif"
            cif_file.write_text(cif_content)
            
            # First load - should parse CIF and save to cache
            # Force use_cache=False to ensure we're parsing from CIF
            loaded1 = processor.load_structure("test_cif", use_cache=False, save_processed=True, debug=True)
            print(f"Loaded1 type: {type(loaded1)}, shape: {loaded1.shape if loaded1 is not None else 'None'}")
            assert loaded1 is not None
            assert len(loaded1) > 0
            
            # Check cache was created
            cache_file = processor.path_cache_dir / "test_cif.pkl"
            print(f"Cache file path: {cache_file}")
            print(f"Cache file exists: {cache_file.exists()}")
            
            # Also check if it was saved anywhere else
            if not cache_file.exists():
                # List all files in cache dir
                if processor.path_cache_dir.exists():
                    print(f"Files in cache dir: {list(processor.path_cache_dir.iterdir())}")
                else:
                    print("Cache dir doesn't exist")
            
            assert cache_file.exists()
            
            # Delete the CIF file
            cif_file.unlink()
            
            # Second load - should use cache even though CIF is gone
            loaded2 = processor.load_structure("test_cif", use_cache=True)
            assert loaded2 is not None
            assert len(loaded2) == len(loaded1)
    
    def test_dataset_metadata_and_pkl_coordination(self):
        """Test that dataset metadata (JSON) and data (PKL) are coordinated."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test structures
            pdb_ids = ['coord1', 'coord2']
            all_data = []
            
            for pdb_id in pdb_ids:
                test_data = pd.DataFrame({
                    'pdb_id': [pdb_id] * 3,
                    'group': ['ATOM'] * 3,
                    'auth_chain_id': ['A'] * 3,
                    'auth_seq_id': [1, 2, 3],
                    'res_name3l': ['ALA', 'GLY', 'VAL'],
                    'res_name1l': ['A', 'G', 'V'],
                    'atom_name': ['CA'] * 3,
                    'x': [1.0, 2.0, 3.0],
                    'y': [4.0, 5.0, 6.0],
                    'z': [7.0, 8.0, 9.0]
                })
                processor.save_entity(pdb_id, test_data)
                all_data.append(test_data)
            
            # Create dataset and save PKL
            processor.create_dataset("coord_test", pdb_ids, metadata={"test": True})
            processor.data = pd.concat(all_data, ignore_index=True)
            processor.save_data("coord_test")
            
            # Verify both exist
            json_path = processor.data_path / "datasets" / "coord_test.json"
            pkl_path = processor.path_dataset_dir / "coord_test.pkl"
            assert json_path.exists()
            assert pkl_path.exists()
            
            # Load dataset info
            info = processor.get_dataset_info("coord_test")
            assert info is not None
            # Extract entity names from the list of dicts
            entity_names = [e['name'] for e in info['entities']]
            assert set(entity_names) == set(pdb_ids)
            assert info['metadata']['test'] is True
            
            # Load dataset data
            loaded = processor.load_data("coord_test")
            assert loaded is not None
            assert len(loaded) == 6  # 2 structures x 3 atoms
    
    def test_no_data_error_handling(self):
        """Test proper error handling when no data is available to save."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Try to save without data
            with pytest.raises(ValueError, match="No data to save"):
                processor.save_data("empty_dataset")
    
    def test_unsupported_format_error(self):
        """Test error handling for unsupported file formats."""
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
            processor = StructureProcessor(paths=paths)
            
            # Create test data
            test_data = pd.DataFrame({
                'pdb_id': ['test'] * 3,
                'group': ['ATOM'] * 3,
                'auth_chain_id': ['A'] * 3,
                'auth_seq_id': [1, 2, 3],
                'res_name3l': ['ALA', 'GLY', 'VAL'],
                'res_name1l': ['A', 'G', 'V'],
                'atom_name': ['CA'] * 3,
                'x': [1.0, 2.0, 3.0],
                'y': [4.0, 5.0, 6.0],
                'z': [7.0, 8.0, 9.0]
            })
            processor.data = test_data
            
            # Try to save with unsupported format
            with pytest.raises(ValueError, match="Unsupported format"):
                processor.save_data("bad_format", file_format="xyz")
            
            # Create a dataset so load_data can find it
            processor.create_dataset("format_test_dataset", ["test"])
            
            # First save in a valid format
            processor.save_data("format_test_dataset", file_format="pkl")
            
            # Try to load with unsupported format (the dataset exists, but format is wrong)
            with pytest.raises(ValueError, match="Unsupported format"):
                processor.load_data("format_test_dataset", file_format="xyz")