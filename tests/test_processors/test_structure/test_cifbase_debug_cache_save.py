"""
Debug why cache files aren't being saved.
"""

import tempfile
import shutil
from pathlib import Path

from protos.io.paths import ProtosPaths
from protos.processing.structure import StructureProcessor


def test_debug_cache_save():
    """Debug cache saving issue."""
    # Create temp directory
    tmpdir = tempfile.mkdtemp()
    print(f"\nWorking directory: {tmpdir}")
    
    try:
        # Set up ProtosPaths
        paths = ProtosPaths(data_root=tmpdir, create_dirs=True)
        
        # Copy real CIF file
        source_file = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/tests/test-data/structure/mmcif/1uaz.cif")
        dest_dir = Path(tmpdir) / "structure" / "mmcif"
        dest_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(source_file, dest_dir / "1uaz.cif")
        
        # Create processor
        processor = StructureProcessor(paths=paths)
        
        # Check paths
        print(f"Cache dir: {processor.path_cache_dir}")
        print(f"Cache dir exists: {processor.path_cache_dir.exists()}")
        
        # Load structure with debug
        print("\nLoading structure with save_processed=True...")
        structure = processor.load_structure("1uaz", use_cache=False, save_processed=True, debug=True)
        
        if structure is not None:
            print(f"Structure loaded: {len(structure)} atoms")
            print(f"Structure type: {type(structure)}")
            print(f"Columns: {list(structure.columns)}")
            
            # Check if save_entity was called
            print("\nChecking cache file...")
            cache_file = processor.path_cache_dir / "1uaz.pkl"
            print(f"Cache file path: {cache_file}")
            print(f"Cache file exists: {cache_file.exists()}")
            
            if not cache_file.exists():
                # Try to save manually
                print("\nTrying manual save_entity...")
                try:
                    processor.save_entity("1uaz", structure)
                    print(f"Manual save complete")
                    print(f"Cache file exists now: {cache_file.exists()}")
                except Exception as e:
                    print(f"Error during manual save: {e}")
                    import traceback
                    traceback.print_exc()
            
            # Check if entity is registered
            print(f"\nEntity registered: {processor.entity_exists('1uaz')}")
        else:
            print("Failed to load structure!")
    
    finally:
        # Cleanup
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    test_debug_cache_save()