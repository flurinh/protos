#!/usr/bin/env python3
"""
Test QM9 loader functionality with compressed format.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.loaders import qm9_loader

def test_qm9_loader():
    """Test QM9 loader with compressed format."""
    print("=== Testing QM9 Loader ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize ProtosPaths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Get QM9 directory
    qm9_dir = Path(paths.get_subdir_path("ligand", "cache_dir")) / "databases" / "qm9"
    qm9_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"QM9 directory: {qm9_dir}")
    
    # 1. Check if QM9 is downloaded
    print("\n1. Checking QM9 status...")
    if qm9_loader.is_qm9_downloaded(qm9_dir):
        print("   ✓ QM9 is downloaded")
        # Check if extracted
        xyz_files = list(qm9_dir.glob("*.xyz"))
        print(f"   XYZ files: {len(xyz_files)}")
        if len(xyz_files) == 0:
            print("   ℹ QM9 archive exists but not extracted")
    else:
        print("   ✗ QM9 not downloaded")
    
    # 2. Test downloading (if not present)
    if not qm9_loader.is_qm9_downloaded(qm9_dir):
        print("\n2. Downloading QM9...")
        success = qm9_loader.download_qm9_dataset(qm9_dir)
        if success:
            print("   ✓ Download successful")
        else:
            print("   ✗ Download failed")
            return False
    
    # 3. Test extraction with the convenience function
    print("\n3. Testing auto-extraction...")
    success = qm9_loader.ensure_qm9_ready(qm9_dir)
    if success:
        print("   ✓ QM9 is ready")
        xyz_files = list(qm9_dir.glob("*.xyz"))
        print(f"   Extracted files: {len(xyz_files)}")
    else:
        print("   ✗ Failed to prepare QM9")
        return False
    
    # 4. Test loading a molecule with auto-extraction
    print("\n4. Testing molecule loading...")
    test_molecules = [1, 100, 1000, 10000]
    
    for mol_id in test_molecules:
        mol_data = qm9_loader.get_qm9_molecule_with_extraction(qm9_dir, mol_id)
        if mol_data:
            print(f"\n   Molecule {mol_id}:")
            print(f"     SMILES: {mol_data.get('smiles', 'N/A')}")
            if 'properties' in mol_data:
                props = mol_data['properties']
                print(f"     HOMO: {props.get('homo', 'N/A')} eV")
                print(f"     LUMO: {props.get('lumo', 'N/A')} eV")
                print(f"     Gap: {props.get('gap', 'N/A')} eV")
            if 'xyz' in mol_data:
                atoms = len(mol_data['xyz'])
                print(f"     Atoms: {atoms}")
        else:
            print(f"\n   ✗ Failed to load molecule {mol_id}")
    
    # 5. Test searching by property
    print("\n5. Testing property search...")
    # Search for molecules with small HOMO-LUMO gap
    results = qm9_loader.search_qm9_with_extraction(qm9_dir, 'gap', 0.1, 0.3)
    print(f"   Found {len(results)} molecules with gap between 0.1-0.3 eV")
    if results:
        print(f"   First result: ID={results[0]['id']}, gap={results[0]['properties']['gap']}")
    
    # 6. Test the extraction function directly
    print("\n6. Testing direct extraction...")
    extracted = qm9_loader.extract_qm9_dataset(qm9_dir)
    if extracted:
        print("   ✓ Extraction successful")
    else:
        print("   ℹ Already extracted or failed")
    
    print("\n✅ QM9 loader tests completed!")
    return True


if __name__ == "__main__":
    success = test_qm9_loader()
    sys.exit(0 if success else 1)