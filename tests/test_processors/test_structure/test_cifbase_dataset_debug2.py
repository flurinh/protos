"""
Debug dataset loading issue.
"""

import pytest
import tempfile
import os
import json
import pandas as pd
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


def test_dataset_loading_steps():
    """Test dataset loading step by step."""
    with tempfile.TemporaryDirectory() as tmpdir:
        paths = ProtosPaths(data_root=tmpdir)
        processor = StructureProcessor(paths=paths)
        
        # Create and save a structure
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
        
        processor.save_entity("test1", test_data)
        
        # Create dataset
        processor.create_dataset("debug_test", ["test1"])
        
        # Step 1: Check if dataset exists in list
        datasets = processor.list_datasets()
        print(f"1. Available datasets: {datasets}")
        assert "debug_test" in datasets
        
        # Step 2: Try get_dataset
        dataset_info = processor.get_dataset("debug_test")
        print(f"2. get_dataset result: {dataset_info}")
        
        # Step 3: Check dataset manager directly
        if processor.dataset_manager:
            try:
                dm_result = processor.dataset_manager.load_dataset("debug_test")
                print(f"3. dataset_manager.load_dataset result: {dm_result}")
            except Exception as e:
                print(f"3. dataset_manager.load_dataset error: {e}")
        
        # Step 4: Check the actual JSON file
        json_path = Path(tmpdir) / "structure" / "datasets" / "debug_test.json"
        print(f"4. JSON path exists: {json_path.exists()}")
        if json_path.exists():
            with open(json_path) as f:
                json_content = json.load(f)
            print(f"4. JSON content: {json_content}")
        
        # Step 5: Try loading data without PKL
        loaded = processor.load_data("debug_test")
        print(f"5. load_data result: {loaded is not None}")
        
        # Step 6: If dataset_info is found, check its structure
        if dataset_info:
            pdb_ids = dataset_info.get('pdb_ids', dataset_info.get('content', dataset_info.get('entities', [])))
            print(f"6. PDB IDs from dataset: {pdb_ids}")
            
            # Try load_structures directly
            if pdb_ids:
                processor.load_structures(pdb_ids)
                print(f"7. After load_structures, data shape: {processor.data.shape if processor.data is not None else 'None'}")
        
        # Step 7: Check if the structure cache files exist
        cache_file = Path(processor.paths.get_subdir_path("structure", "cache_dir")) / "test1.pkl"
        print(f"8. Cache file exists: {cache_file.exists()}")