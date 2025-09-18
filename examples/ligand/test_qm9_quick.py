#!/usr/bin/env python3
"""
Quick test of QM9 loader with already extracted files.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.loaders import qm9_loader

def test_qm9_quick():
    """Test QM9 loader with existing extracted files."""
    print("=== Quick QM9 Loader Test ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize ProtosPaths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Get QM9 directory
    qm9_dir = Path(paths.get_subdir_path("ligand", "cache_dir")) / "databases" / "qm9"
    
    print(f"QM9 directory: {qm9_dir}")
    
    # 1. Check extracted files
    print("\n1. Checking extracted files...")
    mol_dir = qm9_dir / "qm9_molecules"
    if mol_dir.exists():
        xyz_files = list(mol_dir.glob("*.xyz"))
        print(f"   Found {len(xyz_files)} extracted molecules")
    else:
        print("   No extracted molecules found")
        return False
    
    # 2. Test loading specific molecules
    print("\n2. Testing molecule loading...")
    test_molecules = [1, 100, 1000, 10000]
    
    for mol_id in test_molecules:
        mol_data = qm9_loader.load_qm9_molecule(qm9_dir, f"dsgdb9nsd_{mol_id:06d}")
        if mol_data:
            print(f"\n   Molecule {mol_id}:")
            print(f"     ID: {mol_data['id']}")
            print(f"     SMILES: {mol_data.get('smiles', 'N/A')}")
            print(f"     Atoms: {mol_data['n_atoms']}")
            if 'properties' in mol_data:
                props = mol_data['properties']
                print(f"     HOMO: {props.get('homo', 'N/A'):.4f} eV")
                print(f"     LUMO: {props.get('lumo', 'N/A'):.4f} eV")
                print(f"     Gap: {props.get('gap', 'N/A'):.4f} eV")
                print(f"     Dipole: {props.get('mu', 'N/A'):.4f} Debye")
        else:
            print(f"\n   ✗ Failed to load molecule {mol_id}")
    
    # 3. Test property search on subset
    print("\n3. Testing property search...")
    results = qm9_loader.search_qm9_by_property(qm9_dir, 'gap', 0.1, 0.3, limit=10)
    print(f"   Found {len(results)} molecules with gap between 0.1-0.3 eV")
    if results:
        for i, mol in enumerate(results[:3]):
            print(f"   Result {i+1}: ID={mol['id']}, gap={mol['properties']['gap']:.4f} eV, SMILES={mol.get('smiles', 'N/A')}")
    
    # 4. Parse a single XYZ file to verify format
    print("\n4. Testing XYZ parsing...")
    sample_file = mol_dir / "dsgdb9nsd_000001.xyz"
    if sample_file.exists():
        mol_data = qm9_loader.parse_qm9_xyz(sample_file)
        if mol_data:
            print(f"   ✓ Successfully parsed {sample_file.name}")
            print(f"     Elements: {[atom['element'] for atom in mol_data['atoms']]}")
            print(f"     First atom coords: ({mol_data['atoms'][0]['x']:.3f}, {mol_data['atoms'][0]['y']:.3f}, {mol_data['atoms'][0]['z']:.3f})")
        else:
            print(f"   ✗ Failed to parse {sample_file.name}")
    
    # 5. Test index creation
    print("\n5. Testing index creation...")
    success = qm9_loader.create_qm9_index(qm9_dir)
    if success:
        print("   ✓ Index created successfully")
        
        # Load and check index
        index = qm9_loader.load_qm9_index(qm9_dir)
        if index:
            print(f"   Index contains {len(index['molecules'])} molecules")
            print(f"   Unique SMILES: {len(index['by_smiles'])}")
            print(f"   Atom count groups: {len(index['by_n_atoms'])}")
    
    print("\n✅ Quick QM9 tests completed!")
    return True


if __name__ == "__main__":
    success = test_qm9_quick()
    sys.exit(0 if success else 1)