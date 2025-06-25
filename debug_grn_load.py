#!/usr/bin/env python3
"""Debug GRN loading issue"""

import os
os.environ["PROTOS_DATA_ROOT"] = "data"

from protos.processing.grn.grn_base_processor import GRNBaseProcessor

# Create processor without preload
processor = GRNBaseProcessor(
    name="debug",
    data_root="data",
    processor_data_dir="grn"
)

print(f"Data path: {processor.data_path}")

# Try to load manually
try:
    print("\nCalling load_data...")
    result = processor.load_data("ref/mo_ref", file_format="csv")
    print(f"Result type: {type(result)}")
    print(f"Result: {result}")
    
    if hasattr(result, 'shape'):
        print(f"Shape: {result.shape}")
    if hasattr(result, 'columns'):
        print(f"Columns: {result.columns[:5].tolist()}")
    
except Exception as e:
    import traceback
    print(f"Error: {e}")
    traceback.print_exc()

# Check the actual file
import pandas as pd
file_path = "data/grn/ref/mo_ref.csv"
print(f"\nDirect load from {file_path}:")
try:
    df = pd.read_csv(file_path, index_col=0)
    print(f"Loaded shape: {df.shape}")
    print(f"First 5 columns: {df.columns[:5].tolist()}")
    print(f"First 5 rows: {df.index[:5].tolist()}")
except Exception as e:
    print(f"Error: {e}")