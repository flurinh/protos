"""
Test script for LigandProcessor functionality.

This script demonstrates the basic functionality of the ligand handling system
following Protos principles - all operations through the processor.
"""

import os
import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.processing.ligand.ligand_processor import LigandProcessor
from protos.io.paths.path_config import ProtosPaths

def test_ligand_processor():
    """Test basic LigandProcessor functionality."""
    
    print("=== Testing LigandProcessor ===\n")
    
    # Initialize processor
    print("1. Initializing LigandProcessor...")
    lig_proc = LigandProcessor()
    print(f"   Processor type: {lig_proc.processor_type}")
    print(f"   ChEMBL available: {lig_proc.chembl_available}")
    
    # Test SMILES validation and property calculation
    print("\n2. Testing SMILES validation and properties...")
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    
    # Load as entity (will create on-the-fly)
    aspirin_data = lig_proc.load_entity(aspirin_smiles)
    if aspirin_data:
        print(f"   SMILES: {aspirin_data['smiles']}")
        if aspirin_data.get('properties'):
            props = aspirin_data['properties']
            mw = props.get('mw', 'N/A')
            logp = props.get('logp', 'N/A')
            print(f"   MW: {mw:.2f}" if isinstance(mw, (int, float)) else f"   MW: {mw}")
            print(f"   LogP: {logp:.2f}" if isinstance(logp, (int, float)) else f"   LogP: {logp}")
            print(f"   HBA: {props.get('hba', 'N/A')}")
            print(f"   HBD: {props.get('hbd', 'N/A')}")
    
    # Test saving entity
    print("\n3. Testing entity save...")
    ligand_data = {
        'smiles': aspirin_smiles,
        'chembl_id': 'CHEMBL25',
        'properties': aspirin_data.get('properties', {})
    }
    
    try:
        lig_proc.save_entity(aspirin_smiles, ligand_data)
        print("   Entity saved successfully")
        
        # Check if it exists
        if lig_proc.entity_exists(aspirin_smiles):
            print("   Entity exists in registry")
    except Exception as e:
        print(f"   Error saving entity: {e}")
    
    # Test dataset creation
    print("\n4. Testing dataset creation...")
    test_smiles = [
        "CC(=O)OC1=CC=CC=C1C(=O)O",  # Aspirin
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",  # Ibuprofen
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"  # Caffeine
    ]
    
    # Save entities first
    for smiles in test_smiles:
        try:
            lig_proc.save_entity(smiles, {'smiles': smiles})
        except:
            pass
    
    try:
        lig_proc.create_dataset("test_drugs", test_smiles, 
                               metadata={"description": "Common drug molecules"})
        print("   Dataset created successfully")
        
        # List datasets
        datasets = lig_proc.list_datasets()
        print(f"   Available datasets: {datasets}")
    except Exception as e:
        print(f"   Error creating dataset: {e}")
    
    # Test similarity search
    print("\n5. Testing similarity search...")
    similar = lig_proc.search_similar_ligands(aspirin_smiles, similarity=0.5)
    print(f"   Found {len(similar)} similar ligands")
    
    print("\n=== Basic tests completed ===")


def test_chembl_integration():
    """Test ChEMBL integration through the processor."""
    
    print("\n=== Testing ChEMBL Integration ===\n")
    
    # Initialize processor
    lig_proc = LigandProcessor()
    
    if not lig_proc.chembl_available:
        print("   ChEMBL functionality not available")
        print("   Install with: pip install chembl_webresource_client")
        return
    
    print("1. ChEMBL functionality available")
    
    # Test protein mapping through processor
    print("\n2. Testing protein ligand download...")
    test_proteins = ["EGFR", "P00533", "1M17"]
    
    for protein in test_proteins:
        print(f"\n   Testing {protein}:")
        try:
            # Use processor method, not direct ChEMBL
            ligands = lig_proc.get_protein_ligands(protein, limit=2)
            print(f"   Found {len(ligands)} ligands")
            
            if ligands:
                first = ligands[0]
                print(f"   Example: {first.get('chembl_id', 'N/A')}")
                print(f"   Activity: {first.get('activity_type', 'N/A')} = {first.get('value_nm', 'N/A')} nM")
        except Exception as e:
            print(f"   Error: {e}")
    
    print("\n=== ChEMBL tests completed ===")


def main():
    """Run all tests."""
    test_ligand_processor()
    test_chembl_integration()
    print("\n✓ All tests completed!")


if __name__ == "__main__":
    main()