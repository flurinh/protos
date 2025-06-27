#!/usr/bin/env python3
"""
Enhanced setup script to copy real reference data to test-data directory.

This ensures tests have access to real data for rigorous testing.
"""

import os
import shutil
import json
from pathlib import Path


def copy_reference_data():
    """Copy essential reference data files to test-data directory."""
    
    print("Copying reference data to test-data...")
    
    # Define paths
    ref_data_dir = Path("src/protos/reference_data")
    test_data_dir = Path("tests/test-data")
    
    if not ref_data_dir.exists():
        print(f"Reference data directory not found: {ref_data_dir}")
        return
    
    # Ensure test-data exists
    test_data_dir.mkdir(parents=True, exist_ok=True)
    
    # Copy essential files with directory structure
    files_to_copy = [
        # GRN reference tables
        ("grn/ref/mo_ref.csv", "grn/ref/mo_ref.csv"),
        ("grn/ref/gpcrdb_ref.csv", "grn/ref/gpcrdb_ref.csv"),
        
        # GRN configs
        ("grn/configs/binding_domain.json", "grn/configs/binding_domain.json"),
        ("grn/configs/binding_domain2.json", "grn/configs/binding_domain2.json"),
        ("grn/configs/config.json", "grn/configs/config.json"),
        ("grn/configs/motif.json", "grn/configs/motif.json"),
        
        # Structure files
        ("structure/mmcif/1uaz.cif", "structure/mmcif/1uaz.cif"),
        ("structure/mmcif/3ddl.cif", "structure/mmcif/3ddl.cif"),
        ("structure/mmcif/example.cif", "structure/mmcif/example.cif"),
        
        # Dataset definitions
        ("structure/structure_dataset/microbial_opsins.json", "structure/structure_dataset/microbial_opsins.json"),
        ("structure/structure_dataset/test_mo.json", "structure/structure_dataset/test_mo.json"),
        
        # Sequence data
        ("sequence/fasta/test_mo.fasta", "sequence/fasta/test_mo.fasta"),
        ("sequence/fasta/example.fasta", "sequence/fasta/example.fasta"),
    ]
    
    copied_count = 0
    for src_rel, dst_rel in files_to_copy:
        src_path = ref_data_dir / src_rel
        dst_path = test_data_dir / dst_rel
        
        if src_path.exists():
            # Create parent directories
            dst_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Copy file
            shutil.copy2(src_path, dst_path)
            print(f"  ✓ Copied {src_rel}")
            copied_count += 1
        else:
            print(f"  ⚠ Source not found: {src_rel}")
    
    print(f"\nCopied {copied_count} files from reference data to test-data")
    
    # Update registries to reflect new data
    update_test_registries(test_data_dir)


def update_test_registries(test_data_dir: Path):
    """Update registry files to include copied reference data."""
    
    print("\nUpdating registry files...")
    
    # Update GRN registry
    grn_registry = test_data_dir / "grn" / "registry.json"
    grn_data = {}
    if grn_registry.exists():
        with open(grn_registry) as f:
            grn_data = json.load(f)
    
    grn_data.update({
        "mo_ref": {
            "path": "grn/ref/mo_ref.csv",
            "format": "csv",
            "description": "Microbial opsin reference GRN table",
            "metadata": {
                "processor_type": "grn",
                "family": "microbial_opsins"
            }
        },
        "gpcrdb_ref": {
            "path": "grn/ref/gpcrdb_ref.csv",
            "format": "csv",
            "description": "GPCR reference GRN table",
            "metadata": {
                "processor_type": "grn",
                "family": "gpcr"
            }
        }
    })
    
    with open(grn_registry, 'w') as f:
        json.dump(grn_data, f, indent=2)
    print("  ✓ Updated grn/registry.json")
    
    # Update structure registry
    struct_registry = test_data_dir / "structure" / "registry.json"
    struct_data = {}
    if struct_registry.exists():
        with open(struct_registry) as f:
            struct_data = json.load(f)
    
    struct_data.update({
        "1uaz": {
            "path": "structure/mmcif/1uaz.cif",
            "format": "cif",
            "description": "Test structure 1UAZ",
            "metadata": {
                "processor_type": "structure",
                "pdb_id": "1uaz"
            }
        },
        "3ddl": {
            "path": "structure/mmcif/3ddl.cif",
            "format": "cif",
            "description": "Test structure 3DDL",
            "metadata": {
                "processor_type": "structure",
                "pdb_id": "3ddl"
            }
        }
    })
    
    with open(struct_registry, 'w') as f:
        json.dump(struct_data, f, indent=2)
    print("  ✓ Updated structure/registry.json")


if __name__ == "__main__":
    copy_reference_data()