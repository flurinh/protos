#!/usr/bin/env python3
"""
Test unified CCD loader functionality.
Tests download, indexing, and fast access.
"""

import os
import sys
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.io.ingest import ccd_loader_unified as ccd_loader

def test_ccd_loader():
    """Test the unified CCD loader."""
    print("=== Testing Unified CCD Loader ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize ProtosPaths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Get CCD directory - use cache_dir for databases
    ccd_dir = Path(paths.get_subdir_path("molecule", "cache_dir")) / "databases" / "ccd"
    ccd_dir.mkdir(parents=True, exist_ok=True)
    print(f"CCD directory: {ccd_dir}")
    print(f"Directory exists: {ccd_dir.exists()}")
    
    # 1. Test ensure_ccd_ready (download + index if needed)
    print("\n1. Testing ensure_ccd_ready()...")
    start = time.time()
    ready = ccd_loader.ensure_ccd_ready(ccd_dir)
    elapsed = time.time() - start
    print(f"   Result: {'✓ Success' if ready else '✗ Failed'}")
    print(f"   Time: {elapsed:.2f} seconds")
    
    if not ready:
        print("Failed to prepare CCD. Exiting.")
        return False
    
    # 2. Test statistics
    print("\n2. Getting CCD statistics...")
    stats = ccd_loader.get_ccd_statistics(ccd_dir)
    print(f"   Total components: {stats['total_components']:,}")
    print(f"   With SMILES: {stats['with_smiles']:,}")
    print(f"   With InChI: {stats['with_inchi']:,}")
    print(f"   Indexed: {stats['indexed']}")
    print(f"   Component types: {len(stats['types'])}")
    
    # 3. Test fast component lookup
    print("\n3. Testing fast component lookup...")
    test_components = ['ATP', 'NAD', 'HEM', 'SAM', 'COA', 'FAD']
    
    for comp_id in test_components:
        start = time.time()
        comp_data = ccd_loader.get_ccd_component(ccd_dir, comp_id)
        elapsed = time.time() - start
        
        if comp_data:
            print(f"\n   {comp_id}: ✓ Found in {elapsed*1000:.1f} ms")
            print(f"      Name: {comp_data.get('name', 'N/A')}")
            print(f"      Formula: {comp_data.get('formula', 'N/A')}")
            print(f"      Type: {comp_data.get('type', 'N/A')}")
            if comp_data.get('smiles'):
                print(f"      SMILES: {comp_data['smiles'][:50]}...")
        else:
            print(f"\n   {comp_id}: ✗ Not found ({elapsed*1000:.1f} ms)")
    
    # 4. Test search functionality
    print("\n4. Testing search functionality...")
    
    # Search by type
    print("\n   Searching for non-polymer ligands...")
    start = time.time()
    non_polymers = ccd_loader.search_ccd(ccd_dir, 'type', 'NON-POLYMER')
    elapsed = time.time() - start
    print(f"      Found {len(non_polymers)} non-polymer components in {elapsed*1000:.1f} ms")
    if non_polymers:
        print(f"      Examples: {', '.join(non_polymers[:5])}")
    
    # Search for components with SMILES
    print("\n   Searching for components with SMILES...")
    start = time.time()
    with_smiles = ccd_loader.search_ccd(ccd_dir, 'has_smiles', True)
    elapsed = time.time() - start
    print(f"      Found {len(with_smiles)} components with SMILES in {elapsed*1000:.1f} ms")
    
    # Search by formula
    print("\n   Searching by formula (C10H16N5O13P3)...")
    start = time.time()
    atp_like = ccd_loader.search_ccd(ccd_dir, 'formula', 'C10H16N5O13P3')
    elapsed = time.time() - start
    print(f"      Found {len(atp_like)} components in {elapsed*1000:.1f} ms")
    if atp_like:
        print(f"      Components: {', '.join(atp_like)}")
    
    # 5. Performance test - load multiple components
    print("\n5. Performance test - loading 100 random components...")
    import random
    
    # Get list of all components
    all_components = list(stats.get('components', {}).keys()) if stats.get('total_components', 0) > 0 else []
    
    if len(all_components) >= 100:
        random_components = random.sample(all_components[:1000], 100)  # Sample from first 1000
        
        start = time.time()
        found = 0
        for comp_id in random_components:
            comp_data = ccd_loader.get_ccd_component(ccd_dir, comp_id)
            if comp_data:
                found += 1
        elapsed = time.time() - start
        
        print(f"   Loaded {found}/100 components in {elapsed:.2f} seconds")
        print(f"   Average time per component: {elapsed/100*1000:.1f} ms")
    else:
        print("   Not enough components for performance test")
    
    # 6. Test common ligands availability
    print("\n6. Testing common ligands...")
    for ligand_id in ccd_loader.COMMON_LIGANDS[:5]:
        comp = ccd_loader.get_ccd_component(ccd_dir, ligand_id)
        if comp:
            print(f"   {ligand_id}: ✓ {comp.get('name', 'Unknown')}")
        else:
            print(f"   {ligand_id}: ✗ Not found")
    
    print("\n✅ All tests completed!")
    return True


if __name__ == "__main__":
    success = test_ccd_loader()
    sys.exit(0 if success else 1)