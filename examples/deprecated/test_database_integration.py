#!/usr/bin/env python3
"""
Test database integration with proper data loading.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.processing.molecule import MoleculeProcessor
from protos.io.ingest import qm9_loader, ccd_loader

# Set up test environment
test_data_root = Path(__file__).parent / "test_data"
os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())

print("=== Testing Database Integration ===")
print(f"Data root: {test_data_root}")

# Initialize processor
lig_proc = MoleculeProcessor()

# Test 1: Direct loader test for CCD
print("\n--- Test 1: CCD Direct Loader Test ---")
ccd_dir = test_data_root / "ligand" / "databases" / "ccd"
ccd_dir.mkdir(parents=True, exist_ok=True)

# Check if we can load a component
if ccd_loader.ensure_ccd_ready(ccd_dir):
    print("CCD is ready")
    atp = ccd_loader.get_ccd_ligand_safe(ccd_dir, "ATP")
    if atp:
        print(f"✓ ATP loaded: {atp['name']}")
        print(f"  Formula: {atp['formula']}")
        print(f"  SMILES: {atp['smiles'][:40]}...")
else:
    print("CCD not downloaded. Would need to download first.")

# Test 2: Direct loader test for QM9  
print("\n--- Test 2: QM9 Direct Loader Test ---")
qm9_dir = test_data_root / "ligand" / "databases" / "qm9"
qm9_dir.mkdir(parents=True, exist_ok=True)

# Check if we can search molecules
if qm9_loader.ensure_qm9_ready(qm9_dir):
    print("QM9 is ready")
    # Try to load molecule #1
    mol = qm9_loader.get_qm9_molecule_with_extraction(qm9_dir, 1)
    if mol:
        print(f"✓ Molecule 1 loaded")
        print(f"  Atoms: {mol['n_atoms']}")
        print(f"  SMILES: {mol.get('smiles', 'N/A')}")
        print(f"  Gap: {mol['properties'].get('gap', 'N/A')} eV")
else:
    print("QM9 not ready. Archive exists but needs extraction or download.")

# Test 3: Ligand processor integration
print("\n--- Test 3: Ligand Processor Integration ---")

# Try CCD through processor
print("\nTesting CCD through processor:")
atp_data = lig_proc.get_ccd_ligand('ATP', download_if_missing=False)
if atp_data:
    print(f"✓ ATP loaded via processor")
    print(f"  Source: {atp_data.get('source')}")
    print(f"  CCD ID: {atp_data.get('ccd_id')}")
else:
    print("✗ Could not load ATP via processor")

# Try QM9 through processor
print("\nTesting QM9 through processor:")
mol_data = lig_proc.get_qm9_molecule(1, download_if_missing=False)
if mol_data:
    print(f"✓ QM9 molecule loaded via processor")
    print(f"  Source: {mol_data.get('source')}")
    print(f"  QM9 ID: {mol_data.get('qm9_id')}")
    print(f"  Quantum properties: {len(mol_data.get('quantum_properties', {}))} properties")
else:
    print("✗ Could not load QM9 molecule via processor")

# Test 4: Check format conversions
print("\n--- Test 4: Format Conversion Test ---")
if atp_data and atp_data.get('smiles'):
    print("Testing SMILES to CIF conversion...")
    cif_df = lig_proc.convert_to_cif_dataframe(
        atp_data['smiles'], 
        chain_id='L', 
        res_name='ATP'
    )
    if cif_df is not None:
        print(f"✓ Converted to CIF DataFrame: {len(cif_df)} atoms")
        print(f"  Columns: {list(cif_df.columns)[:5]}...")
    else:
        print("✗ CIF conversion failed (likely missing RDKit)")

print("\n=== Integration Test Complete ===")