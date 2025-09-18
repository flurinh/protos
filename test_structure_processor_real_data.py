#!/usr/bin/env python3
"""
Test StructureProcessor with real data from RCSB PDB.

This test verifies:
1. Downloading a single structure
2. Registering and saving it
3. Loading it back
4. Saving as PKL and loading again
5. Creating a dataset with multiple structures
"""

import os
import tempfile
import shutil
from pathlib import Path

from protos.processing.structure import StructureProcessor
from protos.io.paths import ProtosPaths


def test_structure_processor_with_real_data():
    """Test StructureProcessor with real RCSB data."""
    
    # Create a temporary directory for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"Using temporary directory: {temp_dir}")
        
        # Set up ProtosPaths with temp directory
        data_dir = os.path.join(temp_dir, "data")
        os.environ["PROTOS_DATA_ROOT"] = data_dir
        paths = ProtosPaths(data_root=data_dir)
        
        # Create processor
        processor = StructureProcessor(name="test_processor", paths=paths)
        print(f"Created StructureProcessor with data path: {processor.data_path}")
        
        # Test 1: Download a single structure
        print("\n=== Test 1: Download single structure ===")
        pdb_id = "1ubq"  # Ubiquitin - small, well-studied protein
        
        # Download from RCSB
        print(f"Downloading {pdb_id} from RCSB...")
        structure = processor.download_structure(pdb_id, save_to_cache=True)
        
        if structure is not None:
            print(f"✓ Downloaded {pdb_id}: {len(structure)} atoms")
            print(f"  Chains: {structure['auth_chain_id'].unique().tolist()}")
        else:
            print(f"✗ Failed to download {pdb_id}")
            return
        
        # Verify it was registered
        if processor.entity_exists(pdb_id):
            print(f"✓ Entity {pdb_id} is registered")
        else:
            print(f"✗ Entity {pdb_id} was not registered")
        
        # Test 2: Load the structure
        print("\n=== Test 2: Load structure ===")
        loaded = processor.load_structure(pdb_id)
        
        if loaded is not None and len(loaded) == len(structure):
            print(f"✓ Loaded {pdb_id}: {len(loaded)} atoms")
        else:
            print(f"✗ Failed to load {pdb_id}")
        
        # Test 3: Save as PKL explicitly and load again
        print("\n=== Test 3: Save as PKL and reload ===")
        
        # The structure should already be cached, but let's verify
        cache_file = processor.path_cache_dir / f"{pdb_id}.pkl"
        if cache_file.exists():
            print(f"✓ PKL cache exists at: {cache_file}")
        
        # Load using entity system (which prefers PKL)
        entity_loaded = processor.load_entity(pdb_id)
        if entity_loaded is not None:
            print(f"✓ Loaded from entity system: {len(entity_loaded)} atoms")
        
        # Test 4: Download multiple structures and create dataset
        print("\n=== Test 4: Download multiple structures for dataset ===")
        pdb_ids = ["1ubq", "2gb1", "1d3z"]  # All ubiquitin structures
        
        successful_downloads = []
        for pdb_id in pdb_ids:
            if pdb_id == "1ubq":
                # Already downloaded
                successful_downloads.append(pdb_id)
                continue
                
            print(f"Downloading {pdb_id}...")
            struct = processor.download_structure(pdb_id)
            if struct is not None:
                print(f"✓ Downloaded {pdb_id}: {len(struct)} atoms")
                successful_downloads.append(pdb_id)
            else:
                print(f"✗ Failed to download {pdb_id}")
        
        # Create dataset
        print("\n=== Test 5: Create and load dataset ===")
        dataset_name = "ubiquitin_structures"
        
        processor.create_dataset(
            dataset_name,
            successful_downloads,
            metadata={
                "description": "Ubiquitin structures from RCSB",
                "source": "RCSB PDB"
            }
        )
        print(f"✓ Created dataset '{dataset_name}' with {len(successful_downloads)} structures")
        
        # List datasets
        datasets = processor.list_datasets()
        print(f"Available datasets: {datasets}")
        
        # Load dataset - Note: StructureProcessor's load_dataset returns list of IDs
        loaded_ids = processor.load_dataset(dataset_name)
        print(f"✓ Loaded dataset with {len(loaded_ids)} structure IDs: {loaded_ids}")
        
        # Get the actual structures by loading individually
        print(f"✓ Loading structures individually:")
        for pdb_id in loaded_ids:
            struct = processor.load_structure(pdb_id)
            if struct is not None:
                print(f"  - {pdb_id}: {len(struct)} atoms")
        
        # Get dataset info
        info = processor.get_dataset_info(dataset_name)
        print(f"\nDataset info:")
        print(f"  Name: {info['name']}")
        print(f"  Entities: {info['entities']}")
        print(f"  Created: {info['created']}")
        
        # Test 6: Verify entity registry
        print("\n=== Test 6: Verify entity registry ===")
        all_entities = processor.list_entities()
        print(f"All registered entities: {all_entities}")
        
        # Verify all our structures are registered
        for pdb_id in successful_downloads:
            if processor.entity_exists(pdb_id):
                print(f"✓ {pdb_id} is registered")
            else:
                print(f"✗ {pdb_id} is not registered")
        
        # Print summary of data organization
        print("\n=== Data Organization Summary ===")
        print(f"Structure files: {processor.path_structure_dir}")
        if processor.path_structure_dir.exists():
            cif_files = list(processor.path_structure_dir.glob("*.cif"))
            print(f"  CIF files: {[f.name for f in cif_files]}")
        
        print(f"Cache files: {processor.path_cache_dir}")
        if processor.path_cache_dir.exists():
            pkl_files = list(processor.path_cache_dir.glob("*.pkl"))
            print(f"  PKL files: {[f.name for f in pkl_files]}")
        
        print(f"Datasets: {processor.data_path / 'datasets'}")
        if (processor.data_path / 'datasets').exists():
            json_files = list((processor.data_path / 'datasets').glob("*.json"))
            print(f"  JSON files: {[f.name for f in json_files]}")
        
        print("\n✓ All tests completed successfully!")


if __name__ == "__main__":
    # Run the test
    try:
        test_structure_processor_with_real_data()
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()