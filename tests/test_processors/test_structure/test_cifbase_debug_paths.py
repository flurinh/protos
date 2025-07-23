"""
Debug why cache directories aren't being created properly.
"""

import tempfile
import shutil
from pathlib import Path
import os

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


def test_debug_paths():
    """Debug path creation issue."""
    # Create temp directory
    tmpdir = tempfile.mkdtemp()
    print(f'\nWorking directory: {tmpdir}')
    
    try:
        # Set up ProtosPaths with create_dirs=True
        paths = ProtosPaths(data_root=tmpdir)
        print(f'\nProtosPaths data_root: {paths.data_root}')
        
        # Create processor
        processor = StructureProcessor(paths=paths)
        print(f'\nProcessor data_path: {processor.data_path}')
        print(f'Processor type: {processor.processor_type}')
        
        # Check paths
        print(f"\npath_resolver.get_subdir_path('structure', 'structure_dir'): {processor.paths.get_subdir_path('structure', 'structure_dir')}")
        print(f"path_resolver.get_subdir_path('structure', 'cache_dir'): {processor.paths.get_subdir_path('structure', 'cache_dir')}")
        print(f"path_resolver.get_subdir_path('structure', 'dataset_dir'): {processor.paths.get_subdir_path('structure', 'dataset_dir')}")
        
        # Check if directories exist
        print(f"\npath_resolver.get_subdir_path('structure', 'structure_dir') exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'structure_dir'))}")
        print(f"path_resolver.get_subdir_path('structure', 'cache_dir') exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'cache_dir'))}")
        print(f"path_resolver.get_subdir_path('structure', 'dataset_dir') exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'dataset_dir'))}")
        
        # Try to create directories manually through ProtosPaths
        print("\nManually creating directories through ProtosPaths...")
        # Directory creation handled automatically
        
        # Check again
        print(f"\nAfter ensure_processor_dirs:")
        print(f"path_resolver.get_subdir_path('structure', 'cache_dir') exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'cache_dir'))}")
        print(f"path_resolver.get_subdir_path('structure', 'dataset_dir') exists: {os.path.exists(processor.paths.get_subdir_path('structure', 'dataset_dir'))}")
        
        # List all directories under structure
        structure_path = Path(tmpdir) / "structure"
        if os.path.exists(structure_path):
            print(f'\nContents of {structure_path}:')
            for item in structure_path.iterdir():
                print(f'  - {item.name} (dir: {item.is_dir()})')
        
        # Test saving a file
        import pandas as pd
        test_data = pd.DataFrame({'a': [1, 2, 3]})
        
        # Try saving to cache
        cache_file = Path(processor.paths.get_subdir_path("structure", "cache_dir")) / "test.pkl"
        print(f'\nTrying to save to: {cache_file}')
        
        # Ensure parent directory exists
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        test_data.to_pickle(cache_file)
        
        print(f'File saved successfully: {os.path.exists(cache_file)}')
        
    finally:
        # Cleanup
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_debug_paths()