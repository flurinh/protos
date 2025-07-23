"""
Test script for LigandProcessor and ChEMBL loader.

This script demonstrates the basic functionality of the new ligand handling system.
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
    print(f"   Data path: {lig_proc.data_path}")
    print(f"   Processor type: {lig_proc.processor_type}")
    
    # Test SMILES validation and property calculation
    print("\n2. Testing SMILES validation and properties...")
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    
    # Load as entity (will create on-the-fly)
    aspirin_data = lig_proc.load_entity(aspirin_smiles)
    if aspirin_data:
        print(f"   SMILES: {aspirin_data['smiles']}")
        if aspirin_data.get('properties'):
            props = aspirin_data['properties']
            print(f"   MW: {props.get('mw', 'N/A'):.2f}")
            print(f"   LogP: {props.get('logp', 'N/A'):.2f}")
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


def test_chembl_loader():
    """Test ChEMBL loader functionality."""
    
    print("\n=== Testing ChEMBL Loader ===\n")
    
    from protos.loaders.chembl_loader import ChEMBLDL
    
    # Initialize loader
    print("1. Initializing ChEMBL loader...")
    chembl_loader = ChEMBLDL()
    print(f"   ChEMBL client available: {chembl_loader.chembl is not None}")
    
    if not chembl_loader.chembl:
        print("   ChEMBL client not available. Install chembl_webresource_client to test.")
        return
    
    # Test protein mapping
    print("\n2. Testing protein ID mapping...")
    test_proteins = ["EGFR", "P00533", "1M17"]  # Gene name, UniProt, PDB
    
    for protein in test_proteins:
        chembl_id = chembl_loader.map_protein_to_chembl_target(protein)
        print(f"   {protein} -> {chembl_id}")
    
    print("\n=== ChEMBL tests completed ===")


if __name__ == "__main__":
    # Test basic functionality
    test_ligand_processor()
    
    # Test ChEMBL integration (requires chembl_webresource_client)
    test_chembl_loader()
    
    print("\n✓ All tests completed!")