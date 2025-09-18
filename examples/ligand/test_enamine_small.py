"""
Quick test of Enamine small dataset functionality.
Uses the diversity_1k dataset (only 1,000 compounds).
"""

import os
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

from protos.processing.ligand import LigandProcessor

# Set up environment
test_data_root = Path(__file__).parent / "data"
os.environ["PROTOS_DATA_ROOT"] = str(test_data_root.absolute())

# Initialize processor
print("Initializing LigandProcessor...")
lig_proc = LigandProcessor()

# List available Enamine datasets
print("\n=== Available Enamine Datasets ===")
datasets = lig_proc.list_enamine_datasets()
print("Small datasets (good for testing):")
for name, info in datasets.items():
    if '1k' in name.lower() or '1K' in info['size']:
        print(f"  - {name}: {info['description']} ({info['size']})")

# Check which is the default
print(f"\nDefault dataset for similarity search: diversity_1k")

# Get info about the small dataset
print("\n=== Dataset Information ===")
info = lig_proc.get_enamine_dataset_info('diversity_1k')
if info:
    print(f"Name: {info['name']}")
    print(f"Description: {info['description']}")
    print(f"Size: {info['size']}")
    print(f"Downloaded: {info.get('downloaded', False)}")

# Example usage (won't actually download without valid credentials)
print("\n=== Example Usage ===")
print("To search for similar compounds:")
print("""
similar = lig_proc.search_enamine_by_similarity(
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    dataset='diversity_1k',    # Small test dataset (default)
    similarity=0.7
)
""")

print("\nNote: Actual download requires valid Enamine credentials in .env file")
print("See .env.example for format")