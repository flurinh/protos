#!/usr/bin/env python3
"""
Quick test of CCD integration only.
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.io.paths.path_config import ProtosPaths
from protos.processing.molecule import MoleculeProcessor

def test_ccd_integration():
    """Test CCD integration."""
    print("=== Testing CCD Integration ===\n")
    
    # Set up paths
    data_dir = Path("/mnt/c/Users/hidbe/PycharmProjects/protos/data")
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    
    # Initialize processor
    processor = MoleculeProcessor(name="test_ccd")
    
    # Test getting specific components
    test_components = ["ATP", "NAD", "FAD", "HEM", "ADP"]
    
    for ccd_id in test_components:
        print(f"\nGetting CCD component {ccd_id}...")
        try:
            ligand = processor.get_ccd_ligand(ccd_id, download_if_missing=False)
            if ligand:
                print(f"  ✓ Successfully loaded {ccd_id}")
                print(f"    SMILES: {ligand.get('smiles', 'N/A')[:60]}...")
                print(f"    Name: {ligand.get('name', 'N/A')}")
                if 'properties' in ligand:
                    props = ligand['properties']
                    print(f"    MW: {props.get('mw', 'N/A'):.2f}")
                    print(f"    Formula: {ligand.get('formula', 'N/A')}")
            else:
                print(f"  ✗ Failed to load {ccd_id}")
        except Exception as e:
            print(f"  ✗ Error loading {ccd_id}: {e}")
    
    # Test creating a dataset from CCD components
    print("\n\nCreating dataset from CCD components...")
    try:
        success = processor.create_ccd_dataset(
            "nucleotides", 
            ["ATP", "GTP", "CTP", "UTP"],
            download_if_missing=False
        )
        if success:
            print("  ✓ Successfully created nucleotides dataset")
            
            # Load and check the dataset
            dataset = processor.load_dataset("nucleotides")
            print(f"  Dataset contains {len(dataset)} ligands")
            for name, data in dataset.items():
                print(f"    - {name}: {data.get('name', 'Unknown')}")
        else:
            print("  ✗ Failed to create dataset")
    except Exception as e:
        print(f"  ✗ Error creating dataset: {e}")
    
    print("\n✅ CCD integration test completed!")
    return True


if __name__ == "__main__":
    success = test_ccd_integration()
    sys.exit(0 if success else 1)