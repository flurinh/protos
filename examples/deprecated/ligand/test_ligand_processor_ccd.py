#!/usr/bin/env python3
"""
Test ligand processor integration with CCD loader.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.processing.molecule import MoleculeProcessor

def test_ligand_processor_ccd():
    """Test ligand processor CCD integration."""
    print("=== Testing Ligand Processor CCD Integration ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize paths
    paths = ProtosPaths(data_root=str(data_dir.absolute()))
    
    # Create ligand processor
    processor = MoleculeProcessor(name="test_ccd", paths=paths)
    
    # 1. Test getting CCD ligand
    print("1. Testing get_ccd_ligand()...")
    test_ligands = ['ATP', 'NAD', 'HEM', 'FAD', 'COA']
    
    for ligand_id in test_ligands:
        print(f"\n   Getting {ligand_id}...")
        ligand_data = processor.get_ccd_ligand(ligand_id, download_if_missing=False)
        
        if ligand_data:
            print(f"   ✓ Found {ligand_id}")
            print(f"     Name: {ligand_data.get('name', 'N/A')}")
            print(f"     Formula: {ligand_data.get('formula', 'N/A')}")
            print(f"     Type: {ligand_data.get('type', 'N/A')}")
            print(f"     SMILES: {ligand_data['smiles'][:50]}..." if 'smiles' in ligand_data else "     No SMILES")
            
            # Test molecular properties calculation
            if 'properties' in ligand_data and ligand_data['properties']:
                print(f"     MW: {ligand_data['properties'].get('mw', 'N/A')}")
                print(f"     LogP: {ligand_data['properties'].get('logp', 'N/A')}")
        else:
            print(f"   ✗ Failed to get {ligand_id}")
    
    # 2. Test creating a CCD dataset
    print("\n\n2. Testing create_ccd_dataset()...")
    dataset_name = "common_cofactors"
    cofactor_ids = ["ATP", "ADP", "NAD", "FAD", "COA", "SAM"]
    
    print(f"   Creating dataset '{dataset_name}' with: {', '.join(cofactor_ids)}")
    successful = processor.create_ccd_dataset(dataset_name, cofactor_ids, download_if_missing=False)
    
    print(f"   Successfully added {len(successful)} ligands")
    
    # 3. Test loading the created dataset
    print("\n3. Testing load_dataset()...")
    try:
        dataset = processor.load_dataset(dataset_name)
        print(f"   ✓ Loaded dataset with {len(dataset)} entries")
        
        # Show first entry
        if dataset:
            first_key = list(dataset.keys())[0]
            first_entry = dataset[first_key]
            print(f"   First entry: {first_key[:30]}...")
            print(f"     CCD ID: {first_entry.get('ccd_id', 'N/A')}")
            print(f"     Name: {first_entry.get('name', 'N/A')}")
    except Exception as e:
        print(f"   ✗ Failed to load dataset: {e}")
    
    # 4. Test converting CCD ligand to structure format
    print("\n4. Testing convert_to_cif_dataframe()...")
    try:
        # Get ATP data
        atp_data = processor.get_ccd_ligand('ATP', download_if_missing=False)
        if atp_data and 'smiles' in atp_data:
            # Convert to CIF DataFrame
            df = processor.convert_to_cif_dataframe(atp_data['smiles'], res_name='ATP')
            if df is not None:
                print(f"   ✓ Converted ATP to CIF DataFrame")
                print(f"     Atoms: {len(df)}")
                print(f"     Elements: {', '.join(df['element'].unique())}")
                print(f"     Residue: {df['res_name'].iloc[0]}")
            else:
                print("   ✗ Failed to convert ATP to CIF format")
        else:
            print("   ✗ No ATP data available for conversion")
    except Exception as e:
        print(f"   ✗ Error during conversion: {e}")
    
    # 5. Test database statistics
    print("\n5. Testing get_database_statistics()...")
    stats = processor.get_database_statistics()
    
    for db_name, db_stats in stats.items():
        print(f"\n   {db_name}:")
        print(f"     Downloaded: {db_stats['downloaded']}")
        if db_stats['downloaded']:
            if 'component_count' in db_stats:
                print(f"     Components: {db_stats['component_count']}")
            if 'molecule_count' in db_stats:
                print(f"     Molecules: {db_stats['molecule_count']}")
        print(f"     Description: {db_stats['description']}")
    
    print("\n✅ All tests completed!")
    return True


if __name__ == "__main__":
    test_ligand_processor_ccd()