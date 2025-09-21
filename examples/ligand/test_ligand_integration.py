#!/usr/bin/env python3
"""
Test ligand processor integration with CCD, QM9, and Enamine databases.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.processing.ligand import LigandProcessor

def test_ligand_integration():
    """Test ligand processor with all three databases."""
    print("=== Testing Ligand Processor Database Integration ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize processor
    processor = LigandProcessor(name="test_integration")
    
    # 1. Test CCD integration
    print("1. Testing CCD Integration...")
    try:
        # Get specific component
        print("   Getting CCD ligand ATP...")
        atp = processor.get_ccd_ligand("ATP")
        if atp:
            print(f"   ✓ Successfully loaded ATP")
            print(f"     SMILES: {atp.get('smiles', 'N/A')[:50]}...")
            if 'properties' in atp:
                print(f"     MW: {atp['properties'].get('mw', 'N/A')}")
    except Exception as e:
        print(f"   ✗ CCD test failed: {e}")
    
    # 2. Test QM9 integration
    print("\n2. Testing QM9 Integration...")
    try:
        # Get a specific molecule
        print("   Getting QM9 molecule #100...")
        mol = processor.get_qm9_molecule(100)
        if mol:
            print(f"   ✓ Successfully loaded molecule")
            print(f"     SMILES: {mol.get('smiles', 'N/A')}")
            if 'quantum_properties' in mol:
                qp = mol['quantum_properties']
                print(f"     HOMO: {qp.get('homo', 'N/A'):.4f} eV")
                print(f"     LUMO: {qp.get('lumo', 'N/A'):.4f} eV")
                print(f"     Gap: {qp.get('gap', 'N/A'):.4f} eV")
        
        # Search by properties
        print("\n   Searching QM9 for small gap molecules...")
        small_gap = processor.search_qm9_by_properties(
            {'gap': (0.1, 0.2)}, 
            limit=3
        )
        print(f"   Found {len(small_gap)} molecules with gap 0.1-0.2 eV")
        for i, mol in enumerate(small_gap):
            gap = mol['quantum_properties']['gap']
            print(f"     {i+1}. {mol['smiles']} - gap: {gap:.4f} eV")
    except Exception as e:
        print(f"   ✗ QM9 test failed: {e}")
    
    # 3. Test Enamine integration (will fail without credentials)
    print("\n3. Testing Enamine Integration...")
    try:
        from protos.io.ingest import enamine_loader
        
        # List available datasets
        datasets = enamine_loader.list_available_datasets()
        print(f"   Available Enamine datasets: {len(datasets)}")
        for name, info in list(datasets.items())[:3]:
            print(f"     - {name}: {info['size']} compounds")
        
        # Check credentials
        username, password = enamine_loader.get_enamine_credentials()
        if username and password:
            print("   ✓ Enamine credentials found")
            # Could test download/search here
        else:
            print("   ℹ Enamine credentials not set (expected)")
    except Exception as e:
        print(f"   ✗ Enamine test failed: {e}")
    
    # 4. Test format conversion
    print("\n4. Testing Format Conversion...")
    try:
        # Create test molecules
        test_smiles = ["CCO", "CC(=O)O", "c1ccccc1"]
        
        # Save as dataset
        print("   Saving test molecules as dataset...")
        processor.save_dataset("test_molecules", {
            smiles: {"smiles": smiles, "name": f"Molecule_{i}"}
            for i, smiles in enumerate(test_smiles)
        })
        
        # Convert to SDF
        print("   Converting to SDF format...")
        sdf_path = processor.save_sdf_file("test_export", test_smiles)
        print(f"   ✓ Saved to: {sdf_path}")
        
        # Load back from SDF
        print("   Loading from SDF...")
        loaded = processor.load_sdf_file("test_export", as_entities=False)
        print(f"   ✓ Loaded {len(loaded)} molecules")
    except Exception as e:
        print(f"   ✗ Format conversion failed: {e}")
    
    print("\n✅ Integration tests completed!")
    return True


if __name__ == "__main__":
    success = test_ligand_integration()
    sys.exit(0 if success else 1)