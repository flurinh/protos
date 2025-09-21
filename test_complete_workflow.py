#!/usr/bin/env python3
"""End-to-end demonstration of the Protos structure workflow."""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC_DIR = PROJECT_ROOT / "src"
if SRC_DIR.exists():
    sys.path.insert(0, str(SRC_DIR))

import protos


def configure_data_root() -> Path:
    """Point Protos at the home-directory data root."""

    data_root = PROJECT_ROOT / "data"
    data_root.mkdir(parents=True, exist_ok=True)
    protos.set_data_path(str(data_root))
    return data_root


def main():
    print("=== Complete Protos Workflow Test ===\n")

    print("Step 1: Initialize Protos")
    print("-" * 50)

    data_root = configure_data_root()

    from protos.io.paths import get_protos_paths

    paths = get_protos_paths()
    print("✓ Protos initialized")
    print(f"  Data directory: {paths.data_root}")
    print()
    
    # Step 2: Download and register structures
    print("Step 2: Download and Register Data")
    print("-" * 50)
    
    from protos.io.ingest.structure_loader import StructureLoader
    loader = StructureLoader()
    print("✓ StructureLoader created")
    
    # Download from RCSB PDB
    print("\nDownloading from RCSB PDB...")
    pdb_ids = ["1ubq", "7vvl"]
    
    for pdb_id in pdb_ids:
        print(f"\n  Downloading {pdb_id}...")
        try:
            result = loader.download_and_register(pdb_id)
            if result:
                print(f"  ✓ Successfully downloaded and registered '{result}'")
            else:
                print(f"  ✗ Failed to download {pdb_id}")
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    # Download from AlphaFold
    print("\nDownloading from AlphaFold...")
    uniprot_id = "P00533"
    
    try:
        result = loader.download_and_register_alphafold(
            uniprot_id, 
            name="EGFR_HUMAN"
        )
        if result:
            print(f"  ✓ Successfully downloaded and registered '{result}'")
        else:
            print(f"  ✗ Failed to download AlphaFold structure")
    except Exception as e:
        print(f"  ✗ Error: {e}")
    
    # Bulk download with dataset creation
    print("\n\nBulk download with dataset...")
    batch_ids = ["3sn6", "5d5a", "6b73"]
    
    try:
        success, failed = loader.download_batch(
            batch_ids,
            dataset_name="gpcr_structures"
        )
        print(f"  ✓ Batch complete: {len(success)} succeeded, {len(failed)} failed")
        if success:
            print(f"    Successful: {success}")
        if failed:
            print(f"    Failed: {failed}")
    except Exception as e:
        print(f"  ✗ Batch error: {e}")
    
    print()
    
    # Step 3: Load and process with StructureProcessor
    print("Step 3: Load and Process Data")
    print("-" * 50)
    
    from protos.processing.structure import StructureProcessor
    processor = StructureProcessor()
    print("✓ StructureProcessor created")
    
    # List available data
    entities = processor.list_entities()
    datasets = processor.list_datasets()
    
    print(f"\nAvailable structures: {len(entities)}")
    for entity in entities[:5]:  # Show first 5
        print(f"  - {entity}")
    if len(entities) > 5:
        print(f"  ... and {len(entities) - 5} more")
    
    print(f"\nAvailable datasets: {len(datasets)}")
    for dataset in datasets:
        print(f"  - {dataset}")
    
    # Load a single structure
    print("\n\nLoading single structure...")
    if "1ubq" in entities:
        try:
            structure = processor.load_entity("1ubq")
            if structure is not None:
                print(f"✓ Loaded '1ubq':")
                print(f"  Shape: {structure.shape}")
                print(f"  Atoms: {len(structure)}")
                print(f"  Columns: {list(structure.columns)[:5]}...")
                
                # Basic analysis
                chains = structure.index.get_level_values('structure_id').unique()
                print(f"  Structures in DataFrame: {list(chains)}")
            else:
                print("✗ Failed to load structure data")
        except Exception as e:
            print(f"✗ Error loading structure: {e}")
    
    # Load a dataset
    print("\n\nLoading dataset...")
    if "gpcr_structures" in datasets:
        try:
            dataset_dict = processor.load_dataset("gpcr_structures", return_format='dict')
            print(f"✓ Loaded dataset 'gpcr_structures' (dict view):")
            print(f"  Structures in dataset: {len(dataset_dict)}")

            for name, df in dataset_dict.items():
                if df is not None:
                    print(f"  - {name}: {len(df)} atoms")
                else:
                    print(f"  - {name}: Failed to load")

            stacked_df = processor.load_dataset("gpcr_structures", return_format='stacked')
            if stacked_df is not None:
                print("\n  Stacked DataFrame view:")
                print(f"    Shape: {stacked_df.shape}")
                print(f"    Index names: {stacked_df.index.names}")
                print(f"    Columns sample: {list(stacked_df.columns)[:5]}...")
            else:
                print("  ⚠ No stacked DataFrame available (no structures loaded)")
        except Exception as e:
            print(f"✗ Error loading dataset: {e}")
    
    # Summary
    print("\n\n" + "="*50)
    print("WORKFLOW COMPLETE")
    print("="*50)
    print(f"\nData location: {paths.data_root}")
    print(f"Structures registered: {len(entities)}")
    print(f"Datasets created: {len(datasets)}")
    
    print("\nNext steps:")
    print("1. Use processor methods for analysis")
    print("2. Export results with processor.export_entity()")
    print("3. Create visualizations with protos.visualization")


if __name__ == "__main__":
    main()
