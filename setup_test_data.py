#!/usr/bin/env python3
"""
Setup script to create test data structure using the simplified path system.

This script creates the test data directory structure without requiring
all the protos dependencies.
"""

import os
import sys
import json
from pathlib import Path

# Add the src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

# Import only the simplified path configuration
from protos.io.paths.path_config import ProtosPaths


def setup_test_data():
    """Create test data structure in tests/test-data/"""
    
    print("Setting up test data structure...")
    
    # Define test data directory
    test_data_dir = os.path.join(os.path.dirname(__file__), "tests", "test-data")
    
    # Initialize ProtosPaths with test directory
    print(f"Creating data structure in: {test_data_dir}")
    paths = ProtosPaths(data_root=test_data_dir, create_dirs=True)
    
    # Verify directories were created
    print("\nDirectories created:")
    for item in ['structure', 'grn', 'sequence', 'graph', 'property', 'embedding']:
        path = os.path.join(test_data_dir, item)
        if os.path.exists(path):
            print(f"  ✓ {item}/")
    
    # Create sample files
    print("\nCreating sample files:")
    
    # 1. Global registry
    global_registry = paths.get_global_registry_path()
    registry_data = {
        "test_structure_1": {
            "path": "structure/mmcif/1abc.cif",
            "metadata": {
                "processor_type": "structure",
                "dataset_type": "cif",
                "description": "Test structure 1"
            },
            "timestamp": "2024-01-01T00:00:00"
        },
        "test_grn_ref": {
            "path": "grn/ref/test_ref.csv",
            "metadata": {
                "processor_type": "grn",
                "dataset_type": "csv",
                "description": "Test GRN reference table"
            },
            "timestamp": "2024-01-01T00:00:00"
        }
    }
    with open(global_registry, 'w') as f:
        json.dump(registry_data, f, indent=2)
    print(f"  ✓ Created global_registry.json")
    
    # 2. GRN config
    grn_config_dir = paths.get_grn_subdir_path('configs_dir')
    config_file = os.path.join(grn_config_dir, 'config.json')
    with open(config_file, 'w') as f:
        json.dump({
            "test_family": {
                "standard": {
                    "tm1": ["1x28", "1x64"],
                    "tm2": ["2x31", "2x71"]
                },
                "strict": {
                    "tm1": ["1x49", "1x59"],
                    "tm2": ["2x37", "2x50"]
                }
            }
        }, f, indent=2)
    print(f"  ✓ Created grn/configs/config.json")
    
    # 3. Test FASTA
    fasta_dir = paths.get_sequence_subdir_path('fasta_dir')
    fasta_file = os.path.join(fasta_dir, 'test_proteins.fasta')
    with open(fasta_file, 'w') as f:
        f.write(">test_protein_1\n")
        f.write("MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASED\n")
        f.write(">test_protein_2\n")
        f.write("MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKD\n")
    print(f"  ✓ Created sequence/fasta/test_proteins.fasta")
    
    # 4. GRN reference table
    grn_ref_dir = os.path.join(paths.get_processor_path('grn'), 'ref')
    os.makedirs(grn_ref_dir, exist_ok=True)
    ref_table = os.path.join(grn_ref_dir, 'test_ref.csv')
    with open(ref_table, 'w') as f:
        f.write(",1x50,1x51,1x52,2x50,2x51,2x52,n.1,n.2,c.1,c.2\n")
        f.write("ref_protein1,A,L,V,G,S,T,M,K,R,K\n")
        f.write("ref_protein2,A,M,V,G,T,S,-,K,K,R\n")
    print(f"  ✓ Created grn/ref/test_ref.csv")
    
    # 5. Structure dataset
    struct_dataset_dir = paths.get_structure_subdir_path('dataset_dir')
    dataset_file = os.path.join(struct_dataset_dir, 'test_dataset.json')
    with open(dataset_file, 'w') as f:
        json.dump({
            "dataset_id": "test_dataset",
            "name": "Test Structure Dataset",
            "description": "Test dataset for simplified path system",
            "content": ["1abc", "2xyz", "3def"],
            "metadata": {
                "created_by": "setup_test_data.py",
                "structure_count": 3
            }
        }, f, indent=2)
    print(f"  ✓ Created structure/structure_dataset/test_dataset.json")
    
    # 6. Processor registries
    for processor in ['structure', 'grn', 'sequence']:
        registry_path = paths.get_registry_path(processor)
        with open(registry_path, 'w') as f:
            json.dump({
                f"test_{processor}_dataset": {
                    "path": f"{processor}/test_data.dat",
                    "metadata": {
                        "processor_type": processor,
                        "test": True
                    }
                }
            }, f, indent=2)
        print(f"  ✓ Created {processor}/registry.json")
    
    # 7. README
    readme_path = os.path.join(test_data_dir, "README.md")
    with open(readme_path, 'w') as f:
        f.write("# Test Data Directory\n\n")
        f.write("This directory contains test data for the Protos framework using the simplified path system.\n\n")
        f.write("## Structure\n\n")
        f.write("- All data is contained in a single directory hierarchy\n")
        f.write("- No separation between 'reference' and 'user' data\n")
        f.write("- Automatic directory creation when processors are initialized\n\n")
        f.write("## Usage\n\n")
        f.write("```python\n")
        f.write("from protos.io.paths.path_config_simplified import ProtosPaths\n\n")
        f.write("# Initialize with this test directory\n")
        f.write("paths = ProtosPaths(data_root='tests/test-data')\n")
        f.write("```\n")
    print(f"  ✓ Created README.md")
    
    print(f"\nTest data structure created successfully in: {test_data_dir}")
    print("\nDirectory structure:")
    print("tests/test-data/")
    print("├── structure/")
    print("│   ├── mmcif/")
    print("│   ├── alignments/")
    print("│   ├── structure_dataset/")
    print("│   └── registry.json")
    print("├── grn/")
    print("│   ├── tables/")
    print("│   └── registry.json")
    print("├── sequence/")
    print("│   ├── fasta/")
    print("│   ├── alignments/")
    print("│   └── registry.json")
    print("├── graph/")
    print("├── property/")
    print("├── embedding/")
    print("├── global_registry.json")
    print("└── README.md")


if __name__ == "__main__":
    setup_test_data()