"""
Minimal test script for ligand download and processing.

This script demonstrates:
1. Setting up the Protos environment
2. Downloading ligand data from ChEMBL
3. Registering ligands with LigandProcessor
4. Creating and loading datasets
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.processing.ligand import LigandProcessor


def main():
    # Setup data directory
    # Use platform-independent path
    datadir = Path(__file__).parent  # Get the directory where this script is located
    test_data_root = datadir / "data"
    test_data_root.mkdir(exist_ok=True, parents=True)
    
    # Set environment variable for data root
    os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())
    
    # Create ProtosPaths instance
    paths = ProtosPaths()
    
    # Initialize processors
    print("Initializing processors...")
    lig_proc = LigandProcessor(paths=paths)
    
    print(f"✓ LigandProcessor initialized")
    print(f"  Data path: {lig_proc.data_path}")
    print(f"  SDF directory: {lig_proc.sdf_dir}")
    
    # Initialize ChEMBL loader
    print("\n" + "="*60)
    print("Testing ChEMBL Download")
    print("="*60)
    
    # Import ChEMBL here to avoid circular import issues
    try:
        from protos.loaders.chembl_loader import ChEMBLDL
        chembl_loader = ChEMBLDL(data_root=str(test_data_root.absolute()))
    except ImportError as e:
        print(f"⚠️  Could not import ChEMBL loader: {e}")
        chembl_loader = None
    
    # Check if ChEMBL is available
    if not chembl_loader or not chembl_loader.chembl:
        print("\n⚠️  ChEMBL client not available!")
        print("To test ChEMBL functionality, install: pip install chembl_webresource_client")
        print("\nTesting with synthetic data instead...")
        
        # Create synthetic test data
        test_compounds = [
            {
                'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'chembl_id': 'CHEMBL25',
                'activity_type': 'IC50',
                'value': 100,
                'units': 'nM',
                'value_nm': 100,
                'protein_id': 'TEST_PROTEIN'
            },
            {
                'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                'chembl_id': 'CHEMBL521',
                'activity_type': 'IC50',
                'value': 250,
                'units': 'nM', 
                'value_nm': 250,
                'protein_id': 'TEST_PROTEIN'
            },
            {
                'smiles': 'COC1=CC=CC=C1OCCN',
                'chembl_id': 'CHEMBL1234',
                'activity_type': 'Ki',
                'value': 50,
                'units': 'nM',
                'value_nm': 50,
                'protein_id': 'TEST_PROTEIN'
            }
        ]
    else:
        # Test with real ChEMBL data
        print("✓ ChEMBL client available")
        
        # Test protein targets
        test_proteins = ["EGFR", "COX2"]
        
        print("\nTesting protein to ChEMBL mapping:")
        for protein in test_proteins:
            chembl_id = chembl_loader.map_protein_to_chembl_target(protein)
            print(f"  {protein} → {chembl_id}")
            
        # Try with UniProt IDs as well
        print("\nTrying with UniProt IDs:")
        uniprot_test = {"P00533": "EGFR", "P35354": "COX2"}
        for uniprot_id, gene_name in uniprot_test.items():
            chembl_id = chembl_loader.map_protein_to_chembl_target(uniprot_id)
            print(f"  {uniprot_id} ({gene_name}) → {chembl_id}")
        
        # Download ligands for a protein (limit to 5 for testing)
        print(f"\nDownloading ligands for EGFR (limit=5)...")
        chembl_loader.limit = 5
        
        try:
            compounds_dict = chembl_loader.download_protein_ligands(
                "EGFR",
                min_pchembl=6.0,  # Good potency
                save_sdf=True
            )
            test_compounds = compounds_dict.get("EGFR", [])
            print(f"✓ Downloaded {len(test_compounds)} compounds")
            
            # If no compounds found, try with UniProt ID
            if not test_compounds:
                print("\nTrying with UniProt ID P00533...")
                compounds_dict = chembl_loader.download_protein_ligands(
                    "P00533",
                    min_pchembl=5.0,  # Lower threshold
                    save_sdf=True
                )
                test_compounds = compounds_dict.get("P00533", [])
                print(f"✓ Downloaded {len(test_compounds)} compounds")
                
        except Exception as e:
            print(f"⚠️  Download failed: {e}")
            print("Using synthetic data instead...")
            test_compounds = []
            
        # If still no compounds, use synthetic data
        if not test_compounds:
            print("\nNo compounds found from ChEMBL. Using synthetic test data...")
            test_compounds = [
                {
                    'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                    'chembl_id': 'CHEMBL25',
                    'activity_type': 'IC50',
                    'value': 100,
                    'units': 'nM',
                    'value_nm': 100,
                    'protein_id': 'EGFR'
                },
                {
                    'smiles': 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
                    'chembl_id': 'CHEMBL521',
                    'activity_type': 'IC50',
                    'value': 250,
                    'units': 'nM', 
                    'value_nm': 250,
                    'protein_id': 'EGFR'
                }
            ]
    
    # Register compounds with LigandProcessor
    print("\n" + "="*60)
    print("Testing LigandProcessor Entity Management")
    print("="*60)
    
    registered_smiles = []  # Initialize here to avoid UnboundLocalError
    
    if test_compounds:
        print(f"\nRegistering {len(test_compounds)} compounds...")
        
        for i, compound in enumerate(test_compounds):
            try:
                smiles = compound['smiles']
                
                # Save entity
                lig_proc.save_entity(smiles, {
                    'smiles': smiles,
                    'chembl_id': compound.get('chembl_id', f'TEST_{i}')
                })
                
                registered_smiles.append(smiles)
                print(f"  ✓ Registered: {compound.get('chembl_id', f'TEST_{i}')}")
                
                # Calculate and display properties
                props = lig_proc.calculate_properties(smiles)
                if props:
                    print(f"    MW: {props.get('mw', 'N/A'):.1f}, LogP: {props.get('logp', 'N/A'):.2f}")
                
            except Exception as e:
                print(f"  ⚠️  Failed to register compound {i}: {e}")
        
        print(f"\n✓ Registered {len(registered_smiles)} compounds")
    
    # Test entity operations
    print("\n" + "="*60)
    print("Testing Entity Operations")
    print("="*60)
    
    # List all entities
    all_entities = lig_proc.list_entities()
    print(f"\nTotal ligand entities: {len(all_entities)}")
    
    # Show first few
    if all_entities:
        print("\nFirst 3 entities:")
        for entity in all_entities[:3]:
            print(f"  - {entity[:50]}..." if len(entity) > 50 else f"  - {entity}")
    
    # Test loading an entity
    if registered_smiles:
        test_smiles = registered_smiles[0]
        print(f"\nLoading entity: {test_smiles[:30]}...")
        
        entity_data = lig_proc.load_entity(test_smiles)
        if entity_data:
            print("✓ Entity loaded successfully")
            print(f"  ChEMBL ID: {entity_data.get('chembl_id', 'N/A')}")
            if 'properties' in entity_data:
                print(f"  MW: {entity_data['properties'].get('mw', 'N/A')}")
    
    # Create a dataset
    print("\n" + "="*60)
    print("Testing Dataset Management")
    print("="*60)
    
    if registered_smiles:
        dataset_name = "test_ligands"
        print(f"\nCreating dataset '{dataset_name}' with {len(registered_smiles)} compounds...")
        
        try:
            lig_proc.create_dataset(
                dataset_name,
                registered_smiles[:3],  # Use first 3 compounds
                metadata={
                    "description": "Test dataset for ligand download",
                    "source": "ChEMBL" if chembl_loader.chembl else "Synthetic"
                }
            )
            print("✓ Dataset created successfully")
        except Exception as e:
            print(f"⚠️  Failed to create dataset: {e}")
    
    # List datasets
    print("\nListing all datasets:")
    datasets = lig_proc.list_datasets()
    for ds in datasets:
        print(f"  - {ds}")
    
    # Load dataset
    if datasets and dataset_name in datasets:
        print(f"\nLoading dataset '{dataset_name}'...")
        try:
            dataset_contents = lig_proc.load_dataset(dataset_name)
            print(f"✓ Dataset loaded with {len(dataset_contents)} entries")
            
            # Show first entry
            if dataset_contents:
                first_key = list(dataset_contents.keys())[0]
                first_entry = dataset_contents[first_key]
                print(f"\nFirst entry:")
                print(f"  SMILES: {first_key[:30]}...")
                print(f"  ChEMBL ID: {first_entry.get('chembl_id', 'N/A')}")
        except Exception as e:
            print(f"⚠️  Failed to load dataset: {e}")
    
    # Test similarity search
    print("\n" + "="*60)
    print("Testing Similarity Search")
    print("="*60)
    
    if registered_smiles and len(registered_smiles) > 1:
        query_smiles = registered_smiles[0]
        print(f"\nSearching for compounds similar to: {query_smiles[:30]}...")
        
        try:
            similar = lig_proc.search_similar_ligands(query_smiles, similarity=0.3)
            print(f"✓ Found {len(similar)} similar compounds")
            
            for smiles, score in similar[:3]:
                print(f"  - Similarity: {score:.3f} - {smiles[:30]}...")
        except Exception as e:
            print(f"⚠️  Similarity search failed: {e}")
    
    # Test drug-likeness filter
    print("\n" + "="*60)
    print("Testing Drug-Likeness Filter")
    print("="*60)
    
    if registered_smiles:
        print(f"\nChecking drug-likeness for {len(registered_smiles)} compounds...")
        
        drug_like = lig_proc.filter_drug_like(registered_smiles)
        print(f"✓ Drug-like compounds: {len(drug_like)}/{len(registered_smiles)}")
        
        drug_like_strict = lig_proc.filter_drug_like(registered_smiles, strict=True)
        print(f"✓ Drug-like (strict): {len(drug_like_strict)}/{len(registered_smiles)}")
    
    print("\n" + "="*60)
    print("✓ All tests completed!")
    print("="*60)


if __name__ == "__main__":
    main()