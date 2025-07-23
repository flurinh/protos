#!/usr/bin/env python3
"""Test the updated get_grn_interval function."""

from protos.processing.grn.grn_utils import get_grn_interval, GRNConfigManager

# Test 1: Auto-generate GRNs
print("=== Test 1: Auto-generate GRNs ===")
grns = get_grn_interval("1.28", "1.64", grns_str=None)
print(f"Generated {len(grns)} GRNs from 1.28 to 1.64:")
print(f"First 5: {grns[:5]}")
print(f"Last 5: {grns[-5:]}")

# Test 2: Filter from provided list
print("\n=== Test 2: Filter from provided list ===")
# Create a list with some GRNs including the special ones
test_list = ["1.25", "1.28", "1.30", "1.50", "1.551", "1.60", "1.64", "1.70", "2.50"]
filtered = get_grn_interval("1.28", "1.64", grns_str=test_list)
print(f"Filtered {len(filtered)} GRNs from provided list:")
print(f"Result: {filtered}")

# Test 3: Test with actual config
print("\n=== Test 3: Test with GRN config ===")
try:
    config = GRNConfigManager()
    grn_config = config.get_config(protein_family='gpcr_a', strict=False)
    if grn_config and 'tm1' in grn_config:
        tm1_start, tm1_end = grn_config['tm1']
        tm1_grns = get_grn_interval(tm1_start, tm1_end, grns_str=None)
        print(f"TM1 interval: {tm1_start} to {tm1_end}")
        print(f"Generated {len(tm1_grns)} GRNs")
        print(f"Sample: {tm1_grns[::5]}")  # Every 5th GRN
    else:
        print("No config found for gpcr_a")
except Exception as e:
    print(f"Error loading config: {e}")

# Test 4: Generate all standard GRNs for a protein family
print("\n=== Test 4: Generate all standard GRNs ===")
try:
    config = GRNConfigManager()
    grn_config = config.get_config(protein_family='gpcr_a', strict=False)
    if grn_config:
        all_grns = []
        for region, (start, end) in grn_config.items():
            region_grns = get_grn_interval(start, end, grns_str=None)
            all_grns.extend(region_grns)
            print(f"{region}: {len(region_grns)} GRNs ({start} to {end})")
        
        # Remove duplicates
        all_grns = list(set(all_grns))
        print(f"\nTotal unique GRNs: {len(all_grns)}")
    else:
        print("No config found")
except Exception as e:
    print(f"Error: {e}")