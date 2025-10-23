"""
Test script for SDF file handling in LigandProcessor.

This script demonstrates:
1. Loading and saving SDF files
2. Converting between formats
3. Integration with structure processors
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

# Activate conda environment
os.system("source /home/hidberf/miniconda/bin/activate protos")

from protos.processing.molecule import MoleculeProcessor


def test_sdf_operations():
    """Test SDF file operations."""
    print("=== Testing SDF File Operations ===\n")
    
    # Initialize processor
    lig_proc = MoleculeProcessor()
    print(f"✓ Initialized LigandProcessor")
    
    # Test data - some drug molecules
    test_molecules = [
        {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'name': 'Aspirin', 'chembl_id': 'CHEMBL25'},
        {'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'name': 'Ibuprofen', 'chembl_id': 'CHEMBL521'},
        {'smiles': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'name': 'Caffeine', 'chembl_id': 'CHEMBL113'},
        {'smiles': 'CC(C)NCC(COC1=CC=CC=C1OC)O', 'name': 'Propranolol', 'chembl_id': 'CHEMBL27'}
    ]
    
    # 1. Save molecules to SDF file
    print("\n1. Saving molecules to SDF file...")
    try:
        sdf_path = lig_proc.save_sdf_file('test_drugs', test_molecules)
        print(f"   ✓ Saved to: {sdf_path}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
        return
    
    # 2. Load SDF file
    print("\n2. Loading SDF file...")
    try:
        molecules = lig_proc.load_sdf_file('test_drugs', as_entities=True)
        print(f"   ✓ Loaded {len(molecules)} molecules")
        
        # Display first molecule
        if molecules:
            mol = molecules[0]
            print(f"   First molecule:")
            print(f"     SMILES: {mol.get('smiles', 'N/A')}")
            if 'sdf_properties' in mol:
                print(f"     Properties: {list(mol['sdf_properties'].keys())}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # 3. Test conversion to structure format
    print("\n3. Converting ligand to PDB format...")
    aspirin_smiles = 'CC(=O)OC1=CC=CC=C1C(=O)O'
    
    # First save as entity
    try:
        lig_proc.save_entity(aspirin_smiles, {'smiles': aspirin_smiles, 'name': 'Aspirin'})
        print("   ✓ Saved aspirin as entity")
    except:
        pass
    
    try:
        pdb_path = lig_proc.convert_to_structure_format(aspirin_smiles, output_format='pdb')
        if pdb_path:
            print(f"   ✓ Converted to PDB: {pdb_path}")
        else:
            print("   ✗ Conversion failed")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # 4. Test DataFrame operations
    print("\n4. Testing DataFrame operations...")
    import pandas as pd
    
    # Create DataFrame
    df = pd.DataFrame(test_molecules)
    
    try:
        sdf_path = lig_proc.save_sdf_file('test_drugs_df', df)
        print(f"   ✓ Saved DataFrame to SDF: {sdf_path}")
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    # 5. Test loading existing SDF files
    print("\n5. Testing SDF file validation...")
    from protos.io.sdf_utils import validate_sdf_file
    
    if Path(sdf_path).exists():
        is_valid, errors = validate_sdf_file(sdf_path)
        if is_valid:
            print("   ✓ SDF file is valid")
        else:
            print(f"   ✗ Validation errors: {errors}")
    
    print("\n✅ SDF operations test completed!")


def test_chembl_to_sdf():
    """Test downloading ChEMBL compounds and saving to SDF."""
    print("\n\n=== Testing ChEMBL to SDF Workflow ===\n")
    
    lig_proc = MoleculeProcessor()
    
    if not lig_proc.chembl_available:
        print("   ⚠️ ChEMBL not available, skipping test")
        return
    
    print("1. Downloading EGFR inhibitors...")
    try:
        compounds = lig_proc.get_protein_ligands('EGFR', limit=5, min_pchembl=7.0)
        print(f"   ✓ Downloaded {len(compounds)} compounds")
        
        if compounds:
            # Save to SDF
            print("\n2. Saving to SDF file...")
            sdf_path = lig_proc.save_sdf_file('egfr_inhibitors', compounds)
            print(f"   ✓ Saved to: {sdf_path}")
            
            # Extract properties
            from protos.io.sdf_utils import extract_unique_properties
            props = extract_unique_properties(sdf_path)
            print(f"   Properties in SDF: {props}")
            
    except Exception as e:
        print(f"   ✗ Error: {e}")
    
    print("\n✅ ChEMBL to SDF test completed!")


def main():
    """Run all tests."""
    # Set up environment
    test_data_root = Path(__file__).parent / "data"
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Run tests
    test_sdf_operations()
    test_chembl_to_sdf()
    
    print("\n✨ All tests completed!")


if __name__ == "__main__":
    main()