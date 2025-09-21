#!/usr/bin/env python3
"""
Test database loaders with proper extraction workflow.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.ingest import qm9_loader, ccd_loader

# Test directories
test_dir = Path(__file__).parent / "test_data" / "databases"
test_dir.mkdir(parents=True, exist_ok=True)

print("=== Testing Database Download and Extraction ===")
print(f"Test directory: {test_dir}")

# Test QM9 workflow
print("\n--- Testing QM9 Workflow ---")
qm9_dir = test_dir / "qm9"
qm9_dir.mkdir(exist_ok=True)

# Step 1: Check if downloaded
if not qm9_loader.is_qm9_downloaded(qm9_dir):
    print("QM9 not downloaded. Would download here...")
    print(f"URL: {qm9_loader.QM9_URL}")
    print(f"File: {qm9_loader.QM9_FILENAME}")
else:
    print("QM9 archive exists")
    
    # Step 2: Check if extracted
    extract_dir = qm9_dir / "qm9_molecules"
    if not extract_dir.exists() or not any(extract_dir.iterdir()):
        print("QM9 not extracted. Extracting...")
        success = qm9_loader.extract_qm9_dataset(qm9_dir)
        if success:
            print("✓ QM9 extracted successfully")
        else:
            print("✗ QM9 extraction failed")
    else:
        print("QM9 already extracted")
    
    # Step 3: Test search functionality
    if extract_dir.exists():
        print("\nSearching for molecules with specific properties...")
        results = qm9_loader.search_qm9_by_property(
            qm9_dir, 
            'gap',  # HOMO-LUMO gap
            min_value=0.1,
            max_value=0.3,
            limit=5
        )
        
        if results:
            print(f"✓ Found {len(results)} molecules:")
            for mol in results[:3]:
                print(f"  {mol['id']}: gap={mol['properties']['gap']:.3f} eV, {mol['n_atoms']} atoms")

# Test CCD workflow
print("\n--- Testing CCD Workflow ---")
ccd_dir = test_dir / "ccd"
ccd_dir.mkdir(exist_ok=True)

# CCD doesn't need extraction - it's read directly from gzip
if ccd_loader.is_ccd_downloaded(ccd_dir):
    print("CCD components.cif.gz exists")
    
    # Test loading
    print("\nLoading some common ligands...")
    ligands = ['ATP', 'NAD', 'HEM', 'FAD']
    
    for ligand_id in ligands:
        component = ccd_loader.load_ccd_component(ccd_dir, ligand_id)
        if component:
            smiles = ccd_loader.get_ccd_smiles(ccd_dir, ligand_id)
            print(f"✓ {ligand_id}: {component.get('name', 'N/A')}")
            if smiles:
                print(f"  SMILES: {smiles[:40]}...")
        else:
            print(f"✗ {ligand_id}: Not found")
else:
    print("CCD not downloaded")
    print(f"Would download from: {ccd_loader.CCD_COMPONENTS_URL}")

print("\n=== Test Complete ===")