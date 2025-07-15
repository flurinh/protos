"""
Debug test for StructureProcessor load/save functionality.
"""

import pytest
import tempfile
import pandas as pd
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


def test_simple_load_save_workflow():
    """Test the simplest possible load/save workflow."""
    with tempfile.TemporaryDirectory() as tmpdir:
        paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
        processor = StructureProcessor(paths=paths)
        
        # Create simple test data
        test_data = pd.DataFrame({
            'pdb_id': ['test1'] * 3,
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
        
        # Save entity
        processor.save_entity("test1", test_data)
        
        # Create dataset
        processor.create_dataset("simple_test", ["test1"])
        
        # Check dataset exists
        datasets = processor.list_datasets()
        print(f"Available datasets: {datasets}")
        assert "simple_test" in datasets
        
        # Check dataset info
        info = processor.get_dataset_info("simple_test")
        print(f"Dataset info: {info}")
        assert info is not None
        # info['entities'] is a list of dicts with 'name' and 'formats'
        entity_names = [e['name'] for e in info['entities']]
        assert "test1" in entity_names
        
        # Try to load the dataset using load_data
        processor.data = test_data  # Set data first
        processor.pdb_ids = ["test1"]
        
        # Save as PKL
        pkl_path = processor.save_data("simple_test")
        print(f"Saved PKL to: {pkl_path}")
        assert Path(pkl_path).exists()
        
        # Clear and reload
        processor.data = None
        loaded = processor.load_data("simple_test")
        print(f"Loaded data shape: {loaded.shape if loaded is not None else 'None'}")
        assert loaded is not None
        assert len(loaded) == 3