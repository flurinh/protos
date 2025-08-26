"""
Example: Data Management Workflow in Protos

This example demonstrates:
1. Downloading structures with auto-registration
2. Manually adding and registering files
3. Scanning for unregistered files
4. Creating datasets from registered entities
"""

import os
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from protos.commands.register_data import (
    register_structure_file,
    scan_for_unregistered_files,
    create_dataset_from_entities
)
from protos.commands.download_with_registration import (
    download_and_register_structure,
    bulk_download_structures
)
from protos.io.data_access import GlobalRegistry
from protos.processing.structure.struct_base_processor import CifBaseProcessor


def example_download_workflow():
    """Example 1: Download structures with automatic registration."""
    print("\n=== Example 1: Download with Auto-Registration ===")
    
    # Download a single structure
    pdb_id = "1ubq"
    print(f"\nDownloading {pdb_id}...")
    entity_id = download_and_register_structure(pdb_id)
    
    if entity_id:
        print(f"✓ Downloaded and registered: {pdb_id} -> {entity_id}")
        
        # Now we can immediately use it
        proc = CifBaseProcessor(name="example")
        structure = proc.load_structure(pdb_id)
        print(f"✓ Loaded structure with {len(structure)} atoms")
    else:
        print(f"✗ Failed to download {pdb_id}")
    
    # Bulk download
    print("\nBulk downloading structures...")
    pdb_ids = ["1a2y", "1a2x", "1crn"]
    results = bulk_download_structures(pdb_ids, max_workers=3)
    
    successful = sum(1 for v in results.values() if v is not None)
    print(f"✓ Downloaded {successful}/{len(pdb_ids)} structures")


def example_manual_registration():
    """Example 2: Manually add and register files."""
    print("\n=== Example 2: Manual File Registration ===")
    
    # Simulate manually adding a file
    # In real use, users would copy their files to data/structure/mmcif/
    print("\nFor manual registration:")
    print("1. Copy your CIF file to: data/structure/mmcif/")
    print("2. Register it:")
    print("   protos register data/structure/mmcif/my_structure.cif --type structure")
    
    # Or use Python API
    # file_path = Path("data/structure/mmcif/my_structure.cif")
    # if file_path.exists():
    #     entity_id = register_structure_file(file_path)
    #     print(f"✓ Registered: {file_path.name} -> {entity_id}")


def example_scan_and_register():
    """Example 3: Scan for unregistered files."""
    print("\n=== Example 3: Scan for Unregistered Files ===")
    
    from protos.io.paths import ProtosPaths
    paths = ProtosPaths()
    
    # Scan data directories
    print("\nScanning for unregistered files...")
    results = scan_for_unregistered_files(
        paths.data_root,
        file_type="all"
    )
    
    print(f"\nScan Results:")
    print(f"  Registered files: {len(results['registered'])}")
    print(f"  Unregistered files: {len(results['unregistered'])}")
    print(f"  Orphaned entries: {len(results['orphaned_entries'])}")
    
    if results['unregistered']:
        print("\nUnregistered files found:")
        for file_path in results['unregistered'][:5]:  # Show first 5
            print(f"  - {file_path}")
        if len(results['unregistered']) > 5:
            print(f"  ... and {len(results['unregistered']) - 5} more")


def example_create_dataset():
    """Example 4: Create a dataset from registered entities."""
    print("\n=== Example 4: Create Dataset ===")
    
    # Get some registered entities
    registry = GlobalRegistry()
    
    # Find all registered structure entities
    all_entities = []
    for entity_id, entity_data in registry.entity_registry.entities.items():
        if 'structure' in entity_data.get('formats', {}):
            all_entities.append(entity_id)
    
    if len(all_entities) >= 3:
        # Create a dataset with first 3 entities
        dataset_entities = all_entities[:3]
        
        dataset_path = create_dataset_from_entities(
            entity_ids=dataset_entities,
            dataset_name="example_kinases",
            dataset_type="structure",
            metadata={
                "description": "Example dataset of kinase structures",
                "created_by": "data_management_example.py",
                "purpose": "demonstration"
            }
        )
        
        print(f"✓ Created dataset: {dataset_path}")
        print(f"  Contains {len(dataset_entities)} structures")
        
        # Now load the dataset
        proc = CifBaseProcessor(name="dataset_example")
        proc.load_dataset("example_kinases")
        print(f"✓ Loaded dataset with {len(proc.pdb_ids)} structures")
    else:
        print("✗ Not enough registered structures to create dataset")
        print("  Run the download example first!")


def main():
    """Run all examples."""
    print("=== Protos Data Management Examples ===")
    print("\nThis example demonstrates the complete data management workflow.")
    
    # Example 1: Download with auto-registration
    example_download_workflow()
    
    # Example 2: Manual registration
    example_manual_registration()
    
    # Example 3: Scan for unregistered files
    example_scan_and_register()
    
    # Example 4: Create dataset
    example_create_dataset()
    
    print("\n=== Summary ===")
    print("Key concepts demonstrated:")
    print("1. Downloads automatically register files")
    print("2. Manual additions require explicit registration")
    print("3. Scan helps identify unregistered files")
    print("4. Datasets organize entities for analysis")
    print("\nSee docs/DATA_MANAGEMENT.md for more details!")


if __name__ == "__main__":
    main()