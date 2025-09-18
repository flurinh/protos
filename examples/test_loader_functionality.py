#!/usr/bin/env python3
"""
Test the current functionality of CCD and QM9 loaders.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.loaders import qm9_loader, ccd_loader

# Test directories
test_dir = Path(__file__).parent / "test_data" / "databases"
test_dir.mkdir(parents=True, exist_ok=True)

print("=== Testing Database Loaders ===")
print(f"Test directory: {test_dir}")

# Test 1: CCD Loader
print("\n--- Testing CCD Loader ---")
ccd_dir = test_dir / "ccd"
ccd_dir.mkdir(exist_ok=True)

# Check if downloaded
is_downloaded = ccd_loader.is_ccd_downloaded(ccd_dir)
print(f"CCD downloaded: {is_downloaded}")

if is_downloaded:
    # Try to load a component
    print("\nTrying to load ATP from CCD...")
    atp = ccd_loader.load_ccd_component(ccd_dir, "ATP")
    if atp:
        print(f"✓ ATP loaded successfully")
        print(f"  ID: {atp.get('id')}")
        print(f"  Name: {atp.get('name')}")
        print(f"  Formula: {atp.get('formula')}")
        
        # Get SMILES
        smiles = ccd_loader.get_ccd_smiles(ccd_dir, "ATP")
        if smiles:
            print(f"  SMILES: {smiles[:50]}...")
    else:
        print("✗ Failed to load ATP")
        
    # Check what files exist
    print(f"\nFiles in CCD directory:")
    for f in ccd_dir.iterdir():
        print(f"  {f.name} ({f.stat().st_size / 1024 / 1024:.1f} MB)")

# Test 2: QM9 Loader
print("\n--- Testing QM9 Loader ---")
qm9_dir = test_dir / "qm9"
qm9_dir.mkdir(exist_ok=True)

# Check if downloaded
is_downloaded = qm9_loader.is_qm9_downloaded(qm9_dir)
print(f"QM9 downloaded: {is_downloaded}")

if is_downloaded:
    # Check what files exist
    print(f"\nFiles in QM9 directory:")
    for f in qm9_dir.iterdir():
        print(f"  {f.name} ({f.stat().st_size / 1024 / 1024:.1f} MB)")
    
    # Check if extracted
    extract_dir = qm9_dir / "qm9_molecules"
    if extract_dir.exists():
        mol_count = len(list(extract_dir.glob("*.xyz")))
        print(f"\n✓ QM9 extracted: {mol_count} molecules found")
        
        # Try to load a molecule
        if mol_count > 0:
            first_mol = next(extract_dir.glob("*.xyz"))
            mol_data = qm9_loader.parse_qm9_xyz(first_mol)
            if mol_data:
                print(f"\nSample molecule: {mol_data['id']}")
                print(f"  Atoms: {mol_data['n_atoms']}")
                print(f"  SMILES: {mol_data.get('smiles')}")
                print(f"  HOMO-LUMO gap: {mol_data['properties'].get('gap')} eV")
    else:
        print("\n✗ QM9 not extracted yet")
        print("  Archive exists but molecules not extracted")
        
        # Check if we need to extract
        archive = qm9_dir / qm9_loader.QM9_FILENAME
        if archive.exists():
            print(f"\n  Archive found: {archive.name}")
            print("  Need to call extract_qm9_dataset() to extract molecules")

print("\n=== Loader Test Complete ===")