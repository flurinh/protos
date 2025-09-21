#!/usr/bin/env python3
"""
Rebuild CCD database from scratch with proper SMILES extraction.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.io.ingest import ccd_loader_unified as ccd_loader

def rebuild_ccd():
    """Rebuild CCD database from scratch."""
    print("=== Rebuilding CCD Database from Scratch ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize ProtosPaths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Get CCD directory
    ccd_dir = Path(paths.get_subdir_path("ligand", "cache_dir")) / "databases" / "ccd"
    ccd_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"CCD directory: {ccd_dir}")
    
    # Check current status
    print("\n1. Checking current status...")
    if (ccd_dir / "components.cif.gz").exists():
        print("   ✓ CCD file exists")
        print(f"   Size: {(ccd_dir / 'components.cif.gz').stat().st_size / 1024 / 1024:.1f} MB")
    else:
        print("   ✗ CCD file not found")
    
    if (ccd_dir / "ccd_index.json").exists():
        print("   ✓ Index exists")
    else:
        print("   ✗ Index not found")
    
    # Download CCD if needed
    print("\n2. Downloading CCD...")
    if not (ccd_dir / "components.cif.gz").exists():
        success = ccd_loader.download_ccd(ccd_dir, force=False)
        if success:
            print("   ✓ Download successful")
        else:
            print("   ✗ Download failed")
            return False
    else:
        print("   ℹ CCD already downloaded")
    
    # Build index
    print("\n3. Building index (this may take a few minutes)...")
    success = ccd_loader.build_ccd_index(ccd_dir, force=True)
    
    if success:
        print("   ✓ Index build successful")
        
        # Get statistics
        stats = ccd_loader.get_ccd_statistics(ccd_dir)
        print(f"\n4. CCD Statistics:")
        print(f"   Total components: {stats['total_components']:,}")
        print(f"   With SMILES: {stats['with_smiles']:,}")
        print(f"   With InChI: {stats['with_inchi']:,}")
        print(f"   Component types: {len(stats['types'])}")
        
        # Test some components
        print("\n5. Testing component access...")
        test_components = ['ATP', 'NAD', 'HEM', 'FAD', 'COA']
        
        for comp_id in test_components:
            comp = ccd_loader.get_ccd_component(ccd_dir, comp_id)
            if comp:
                smiles = comp.get('smiles') or comp.get('smiles_canonical')
                if smiles:
                    print(f"   ✓ {comp_id}: {smiles[:50]}...")
                else:
                    print(f"   ✗ {comp_id}: No SMILES found")
            else:
                print(f"   ✗ {comp_id}: Not found")
        
        return True
    else:
        print("   ✗ Index build failed")
        return False


if __name__ == "__main__":
    success = rebuild_ccd()
    if success:
        print("\n✅ CCD database rebuilt successfully!")
    else:
        print("\n✗ Failed to rebuild CCD database")
        sys.exit(1)