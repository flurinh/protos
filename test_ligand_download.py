"""
Minimal test script for ligand download and processing - FIXED version.

This script demonstrates proper usage of the LigandProcessor following Protos principles:
1. NO direct path manipulation
2. NO direct loader usage
3. ALL operations through the processor
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.processing.ligand import LigandProcessor


def main():
    """Test ligand processor functionality."""
    # Setup - only set environment variable, no path manipulation
    test_data_root = Path(__file__).parent / "data"
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Initialize processor - that's all!
    print("Initializing LigandProcessor...")
    lig_proc = LigandProcessor()
    print("✓ LigandProcessor initialized")
    
    # Test ChEMBL functionality
    print("\n" + "="*60)
    print("Testing ChEMBL Integration")
    print("="*60)
    
    if not lig_proc.chembl_available:
        print("\n⚠️  ChEMBL functionality not available!")
        print("To test ChEMBL functionality, install: pip install chembl_webresource_client")
        print("\nUsing synthetic test data instead...")
        
        # Create synthetic test data
        test_compounds = [
            {
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',  # Aspirin
                'chembl_id': 'CHEMBL25',
                'activity_type': 'IC50',
                'value_nm': 100,
                'protein_id': 'COX2'
            },
            {
                'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # Ibuprofen
                'chembl_id': 'CHEMBL521',
                'activity_type': 'IC50',
                'value_nm': 250,
                'protein_id': 'COX2'
            }
        ]
    else:
        print("✓ ChEMBL functionality available")
        
        # Download ligands using processor method
        print("\nDownloading ligands for EGFR (limit=5)...")
        try:
            test_compounds = lig_proc.get_protein_ligands(
                "EGFR",
                min_pchembl=6.0,
                limit=5
            )
            print(f"✓ Downloaded {len(test_compounds)} compounds")
            
            # If no compounds, try another target
            if not test_compounds:
                print("\nTrying COX2...")
                test_compounds = lig_proc.get_protein_ligands(
                    "COX2",
                    min_pchembl=5.0,
                    limit=5
                )
                print(f"✓ Downloaded {len(test_compounds)} compounds")
                
        except Exception as e:
            print(f"⚠️  Download failed: {e}")
            test_compounds = []
    
    # If still no compounds, use synthetic data
    if not test_compounds:
        print("\nUsing synthetic test data...")
        test_compounds = [
            {
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'chembl_id': 'CHEMBL25',
                'protein_id': 'TEST'
            }
        ]
    
    # Test entity operations
    print("\n" + "="*60)
    print("Testing Entity Operations")
    print("="*60)
    
    # Save entities (they're already saved if from ChEMBL, but let's ensure)
    registered_smiles = []
    for compound in test_compounds[:3]:  # Limit to 3 for testing
        try:
            smiles = compound['smiles']
            lig_proc.save_entity(smiles, compound)
            registered_smiles.append(smiles)
            print(f"✓ Saved: {compound.get('chembl_id', 'Unknown')}")
        except Exception as e:
            print(f"⚠️  Failed to save: {e}")
    
    # List entities
    print("\nListing entities...")
    entities = lig_proc.list_entities()
    print(f"Total entities: {len(entities)}")
    if entities:
        print("First 3:")
        for entity in entities[:3]:
            print(f"  - {entity[:50]}..." if len(entity) > 50 else f"  - {entity}")
    
    # Load entity
    if registered_smiles:
        print(f"\nLoading entity...")
        entity_data = lig_proc.load_entity(registered_smiles[0])
        if entity_data:
            print("✓ Entity loaded")
            print(f"  ChEMBL ID: {entity_data.get('chembl_id', 'N/A')}")
            props = entity_data.get('properties', {})
            if props:
                print(f"  MW: {props.get('mw', 'N/A')}")
    
    # Test datasets
    print("\n" + "="*60)
    print("Testing Dataset Operations")
    print("="*60)
    
    if registered_smiles:
        dataset_name = "test_ligands"
        
        # Create dataset
        print(f"\nCreating dataset '{dataset_name}'...")
        try:
            lig_proc.create_dataset(
                dataset_name,
                registered_smiles,
                metadata={"test": True, "source": "ChEMBL/synthetic"}
            )
            print("✓ Dataset created")
        except Exception as e:
            print(f"⚠️  Failed: {e}")
        
        # List datasets
        datasets = lig_proc.list_datasets()
        print(f"\nAvailable datasets: {datasets}")
        
        # Load dataset
        if dataset_name in datasets:
            print(f"\nLoading dataset '{dataset_name}'...")
            try:
                dataset_data = lig_proc.load_dataset(dataset_name)
                print(f"✓ Loaded {len(dataset_data)} entries")
            except Exception as e:
                print(f"⚠️  Failed: {e}")
    
    # Test similarity search (if RDKit available)
    print("\n" + "="*60)
    print("Testing Similarity Search")
    print("="*60)
    
    if registered_smiles and len(registered_smiles) > 1:
        query = registered_smiles[0]
        print(f"\nSearching similar to: {query[:30]}...")
        try:
            similar = lig_proc.search_similar_ligands(query, similarity=0.5)
            print(f"✓ Found {len(similar)} similar compounds")
            for smiles, score in similar[:2]:
                print(f"  - Score: {score:.3f} - {smiles[:30]}...")
        except Exception as e:
            print(f"⚠️  Search failed: {e}")
    
    # Test property calculation
    print("\n" + "="*60)
    print("Testing Property Calculation")
    print("="*60)
    
    if registered_smiles:
        smiles = registered_smiles[0]
        print(f"\nCalculating properties for: {smiles[:30]}...")
        props = lig_proc.calculate_properties(smiles)
        if props:
            print("✓ Properties calculated:")
            print(f"  MW: {props.get('mw', 'N/A')}")
            print(f"  LogP: {props.get('logp', 'N/A')}")
            print(f"  HBA: {props.get('hba', 'N/A')}")
            print(f"  HBD: {props.get('hbd', 'N/A')}")
        else:
            print("⚠️  Properties not available (RDKit may be missing)")
    
    # Test drug-likeness filter
    if registered_smiles:
        print("\nTesting drug-likeness filter...")
        drug_like = lig_proc.filter_drug_like(registered_smiles)
        print(f"Drug-like compounds: {len(drug_like)}/{len(registered_smiles)}")
    
    print("\n✅ All tests completed!")


if __name__ == "__main__":
    main()