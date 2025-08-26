"""
Quick test to verify CIF is the default format for ligand structure conversion.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.processing.ligand import LigandProcessor

# Set up environment
test_data_root = Path(__file__).parent / "data"
os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())

# Initialize processor
lig_proc = LigandProcessor()

# Test molecule
caffeine = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"

# Save entity
lig_proc.save_entity(caffeine, {'smiles': caffeine, 'name': 'Caffeine'})

# Test 1: Default format should be CIF
print("Test 1: Default format")
path = lig_proc.convert_to_structure_format(caffeine)
if path and path.endswith('.cif'):
    print("✓ Default format is CIF:", path)
else:
    print("✗ Default format is not CIF:", path)

# Test 2: Explicit CIF format
print("\nTest 2: Explicit CIF format")
path = lig_proc.convert_to_structure_format(caffeine, output_format='cif', res_name='CAF')
if path and path.endswith('.cif'):
    print("✓ Explicit CIF works:", path)
else:
    print("✗ Explicit CIF failed:", path)

# Test 3: PDB format still works
print("\nTest 3: PDB format still available")
path = lig_proc.convert_to_structure_format(caffeine, output_format='pdb')
if path and path.endswith('.pdb'):
    print("✓ PDB format works:", path)
else:
    print("✗ PDB format failed:", path)

print("\n✅ All tests completed!")