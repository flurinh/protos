#!/usr/bin/env python3
"""
Quick test to check the status of local database downloads.
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

print("=== Database Status Check ===")
print(f"Data root: {test_data_root}")

try:
    # Initialize processor
    lig_proc = LigandProcessor()
    
    # Check database status
    stats = lig_proc.get_database_statistics()
    
    print("\n=== Database Status ===")
    for db_name, info in stats.items():
        print(f"\n{db_name}:")
        print(f"  Description: {info['description']}")
        
        if db_name == 'Enamine':
            print(f"  Credentials set: {info.get('credentials_set', False)}")
            print(f"  Available datasets: {info.get('available_datasets', 0)}")
            print(f"  Downloaded datasets: {len(info.get('downloaded_datasets', []))}")
        else:
            print(f"  Downloaded: {info.get('downloaded', False)}")
            if info.get('downloaded'):
                print(f"  Path: {info.get('path', 'N/A')}")
    
    # Check ProtosPaths is managing everything
    print(f"\nDatabase directory: {lig_proc.paths.get_processor_path('ligand')}/databases/")
    
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()

print("\n✅ Status check complete!")