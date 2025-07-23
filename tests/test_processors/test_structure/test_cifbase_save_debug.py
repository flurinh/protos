"""
Debug why files aren't being saved to test-data.
"""

import tempfile
import shutil
from pathlib import Path
import pandas as pd
import os

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


def test_save_debug():
    """Debug save issue."""
    # Create temp directory
    tmpdir = tempfile.mkdtemp()
    print(f"\nWorking directory: {tmpdir}")
    
    try:
        # Set up ProtosPaths
        paths = ProtosPaths(data_root=tmpdir)
        
        # Create processor
        processor = StructureProcessor(paths=paths)
        
        # Create test data
        test_data = pd.DataFrame({
            'pdb_id': ['test1'] * 5,
            'auth_chain_id': ['A'] * 5,
            'auth_seq_id': [1, 2, 3, 4, 5],
            'atom_name': ['CA'] * 5,
            'x': [1.0, 2.0, 3.0, 4.0, 5.0],
            'y': [1.0, 2.0, 3.0, 4.0, 5.0],
            'z': [1.0, 2.0, 3.0, 4.0, 5.0]
        })
        
        print(f"\nCache directory: {processor.paths.get_subdir_path('structure', 'cache_dir')}")
        print(f"Cache directory exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'cache_dir'))}")
        
        # Save using save_entity
        print("\nSaving entity 'test1'...")
        processor.save_entity("test1", test_data, metadata={"test": True})
        
        # Check if file was created
        cache_file = Path(processor.paths.get_subdir_path('structure', 'cache_dir')) / 'test1.pkl'
        print(f'\nCache file path: {cache_file}')
        print(f'Cache file exists: {os.path.exists(cache_file)}')
        
        if os.path.exists(cache_file):
            print(f'Cache file size: {cache_file.stat().st_size} bytes')
        
        # List contents of cache directory
        print(f'\nContents of cache directory:')
        for item in Path(processor.paths.get_subdir_path('structure', 'cache_dir')).iterdir():
            print(f"  - {item.name}")
        
        # Test dataset save
        processor.data = test_data
        print(f"\nDataset directory: {processor.paths.get_subdir_path('structure', 'dataset_dir')}")
        print(f"Dataset directory exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'dataset_dir'))}")
        
        pkl_path = processor.save_data("test_dataset")
        print(f'\nSaved dataset to: {pkl_path}')
        print(f'Dataset file exists: {Path(pkl_path).exists()}')
        
        # List contents of dataset directory
        print(f"\nContents of dataset directory:")
        for item in Path(processor.paths.get_subdir_path('structure', 'dataset_dir')).iterdir():
            print(f"  - {item.name}")
            
    finally:
        # Cleanup
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_save_debug()