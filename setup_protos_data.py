#!/usr/bin/env python3
"""
Setup script for Protos data directory structure.

This script:
1. Creates the proper data directory structure
2. Copies reference data from the repository
3. Downloads structure datasets if requested
"""

import os
import shutil
import argparse
from pathlib import Path
import json

def create_directory_structure(data_root):
    """Create the complete Protos data directory structure."""
    
    directories = [
        # Structure directories
        "structure/mmcif",
        "structure/alignments",
        "structure/structure_dataset",
        
        # GRN directories
        "grn/ref",
        "grn/tables",
        "grn/datasets",
        "grn/configs",
        
        # Sequence directories
        "sequence/fasta",
        "sequence/fasta/raw",
        "sequence/fasta/processed",
        "sequence/alignments",
        
        # Embedding directories
        "embedding/models",
        "embedding/embeddings",
        
        # Property directories
        "property/annotations",
        "property/calculated",
        
        # Ligand directories
        "ligand/binding_sites",
        "ligand/structures",
        
        # Graph directories
        "graph/networks",
        "graph/analysis",
        
        # Output directories
        "output",
        "temp",
        "cache"
    ]
    
    print(f"Creating directory structure at: {data_root}")
    
    for dir_path in directories:
        full_path = Path(data_root) / dir_path
        full_path.mkdir(parents=True, exist_ok=True)
        print(f"  ✓ Created: {dir_path}")
    
    return True


def copy_reference_data(repo_root, data_root):
    """Copy reference data from repository to data directory."""
    
    # Define what to copy
    copy_mappings = [
        # GRN reference tables
        ("src/protos/reference_data/grn/ref", "grn/ref"),
        ("tests/test-data/grn/ref", "grn/ref"),  # Also copy test reference
        
        # GRN configs
        ("src/protos/reference_data/grn/config.json", "grn/configs/config.json"),
        ("tests/test-data/grn/configs", "grn/configs"),
        
        # Test sequences
        ("tests/test-data/sequence/fasta", "sequence/fasta/processed"),
        
        # Example structure datasets
        ("src/protos/reference_data/structure/structure_dataset", "structure/structure_dataset"),
    ]
    
    print("\nCopying reference data...")
    
    for src_rel, dst_rel in copy_mappings:
        src_path = Path(repo_root) / src_rel
        dst_path = Path(data_root) / dst_rel
        
        if not src_path.exists():
            print(f"  ⚠ Source not found: {src_rel}")
            continue
            
        try:
            if src_path.is_file():
                dst_path.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src_path, dst_path)
                print(f"  ✓ Copied file: {src_rel} -> {dst_rel}")
            else:
                # Copy directory contents
                if dst_path.exists():
                    # Merge with existing
                    for item in src_path.iterdir():
                        if item.is_file():
                            shutil.copy2(item, dst_path / item.name)
                else:
                    shutil.copytree(src_path, dst_path, dirs_exist_ok=True)
                print(f"  ✓ Copied directory: {src_rel} -> {dst_rel}")
        except Exception as e:
            print(f"  ✗ Error copying {src_rel}: {e}")
    
    return True


def create_example_fasta(data_root):
    """Create example FASTA files for testing."""
    
    # Create a test dataset FASTA
    test_fasta = """
>BR_test Bacteriorhodopsin test sequence
MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVKGMGVSDPDAKKFYAITT
LVPAIAFTMYLSMLLGYGLTMVPFGGEQNPIYWARYADWLFTTPLLLLDLALLVDADQGT
ILALVGADGIMIGTGLVGALTKVYSYRFVWWAISTAAMLYILYVLFFGFTSKAESMRPEV
ASTFKVLRNVTVVLWSAYPVVWLIGSEGAGIVPLNIETLLFMVLDVSAKVGFGLILLRSR
AIFGEAEAPEPSAGDGAAATSD

>HR_test Halorhodopsin test sequence  
MLFPTMTDTLSVQAARMKDVNITLMIGLFVAYFLISGLFYWIKKGFGPADPRNIPIFGVM
GAIMAAAFISWWSAVAVMDKQGTILALVGADGVIIGGALGGALIETTDKRIWWAITTAAM
LYILYFFAGTSKAESIRPDVAKSLAVHVIWIIWSAYPVIWMIGGEGAGTVPADMISILLD
VAAKVGFGLILLRSRSILQSAEAVEPSAGDAAAEAD

>ChR2_test Channelrhodopsin-2 test sequence
MDYGGALSAVGRELLFVTNPVVVNGSIVIPEMKTQVEEFYSSLQSSGMSVISLVTFSLFS
SLGGWIGMTLHWCAAWGIIQVLQALLVWQVLILERQYWILVKNIPVLTQVVNVLDLVPLL
LNDLALLVDAAQGTILLGVGADIIMGVGLVGALTQVFSYRFVWWAISTAALLYILYVLFF
GFTSKADMRPEVASTFQVQRNVVVVLWCAFPVAVWLIGSEGATVAPDMKTIWDVLAKVGF
GLILLRSRSILGAGGEAPEPAAGGDALAD
""".strip()
    
    # Write test FASTA file
    fasta_path = Path(data_root) / "sequence" / "fasta" / "processed" / "mo_test.fasta"
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(fasta_path, 'w') as f:
        f.write(test_fasta)
    
    print(f"\n✓ Created example FASTA: {fasta_path}")
    
    return True


def setup_environment(data_root):
    """Set up environment variables."""
    
    env_file = Path.home() / ".protos_env"
    
    env_content = f"""# Protos environment configuration
export PROTOS_DATA_ROOT="{data_root}"
export PROTOS_REF_DATA_ROOT="{data_root}"

# Add to your .bashrc or .zshrc:
# source ~/.protos_env
"""
    
    with open(env_file, 'w') as f:
        f.write(env_content)
    
    print(f"\n✓ Created environment file: {env_file}")
    print("  Add this line to your .bashrc or .zshrc:")
    print(f"  source {env_file}")
    
    # Also set for current session
    os.environ["PROTOS_DATA_ROOT"] = str(data_root)
    os.environ["PROTOS_REF_DATA_ROOT"] = str(data_root)
    
    return True


def download_structure_dataset(data_root, dataset_name="microbial_opsins"):
    """Download a structure dataset (placeholder for actual download)."""
    
    print(f"\n📥 Downloading structure dataset: {dataset_name}")
    print("  Note: Actual download implementation needed")
    print("  For now, creating example dataset configuration...")
    
    # Create example dataset configuration
    dataset_config = {
        "name": dataset_name,
        "description": "Example microbial opsin structures",
        "pdb_ids": ["1UAZ", "3DDL", "4PXK", "6CMO", "7ZOU"],
        "filters": {
            "resolution": 3.0,
            "method": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
        }
    }
    
    dataset_path = Path(data_root) / "structure" / "structure_dataset" / f"{dataset_name}.json"
    
    with open(dataset_path, 'w') as f:
        json.dump(dataset_config, f, indent=2)
    
    print(f"  ✓ Created dataset config: {dataset_path}")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Setup Protos data directory structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Setup in default location (~/protos_data)
  python setup_protos_data.py
  
  # Setup in custom location
  python setup_protos_data.py --data-root /path/to/data
  
  # Setup and download datasets
  python setup_protos_data.py --download-datasets
  
  # Just copy reference data to existing directory
  python setup_protos_data.py --skip-create
"""
    )
    
    parser.add_argument(
        '--data-root',
        type=str,
        default=str(Path.home() / "protos_data"),
        help='Root directory for Protos data (default: ~/protos_data)'
    )
    
    parser.add_argument(
        '--repo-root',
        type=str,
        default=str(Path(__file__).parent),
        help='Root directory of Protos repository'
    )
    
    parser.add_argument(
        '--skip-create',
        action='store_true',
        help='Skip directory creation (only copy reference data)'
    )
    
    parser.add_argument(
        '--download-datasets',
        action='store_true',
        help='Download example structure datasets'
    )
    
    parser.add_argument(
        '--create-examples',
        action='store_true',
        default=True,
        help='Create example FASTA files (default: True)'
    )
    
    args = parser.parse_args()
    
    # Convert to absolute paths
    data_root = Path(args.data_root).absolute()
    repo_root = Path(args.repo_root).absolute()
    
    print(f"=== Protos Data Setup ===")
    print(f"Repository root: {repo_root}")
    print(f"Data root: {data_root}")
    
    # Step 1: Create directory structure
    if not args.skip_create:
        create_directory_structure(data_root)
    
    # Step 2: Copy reference data
    copy_reference_data(repo_root, data_root)
    
    # Step 3: Create example files
    if args.create_examples:
        create_example_fasta(data_root)
    
    # Step 4: Set up environment
    setup_environment(data_root)
    
    # Step 5: Download datasets if requested
    if args.download_datasets:
        download_structure_dataset(data_root)
    
    print("\n✅ Setup complete!")
    print(f"\nTo use Protos CLI tools:")
    print(f"  export PROTOS_DATA_ROOT={data_root}")
    print(f"  python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_test")
    
    return 0


if __name__ == "__main__":
    exit(main())