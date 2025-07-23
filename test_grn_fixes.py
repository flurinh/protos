#!/usr/bin/env python3
"""Test the GRN annotation fixes."""

from protos.processing.grn.grn_utils import generate_all_grns_from_config, GRNConfigManager
from protos.processing.grn.grn_table_utils import expand_annotation
import pandas as pd

# Test 1: Check if generate_all_grns_from_config works
print("=== Test 1: Generate all GRNs from config ===")
config = GRNConfigManager()
grn_config = config.get_config(protein_family='gpcr_a', strict=False)
print(f"Config keys: {list(grn_config.keys())}")

all_grns = generate_all_grns_from_config(grn_config)
print(f"Total GRNs generated: {len(all_grns)}")
print(f"First 10 GRNs: {all_grns[:10]}")
print(f"Last 10 GRNs: {all_grns[-10:]}")

# Check for TM1 GRNs
tm1_grns = [g for g in all_grns if g.startswith('1.')]
print(f"\nTM1 GRNs ({len(tm1_grns)}): {tm1_grns[:5]}...{tm1_grns[-5:]}")

# Test 2: Test expand_annotation with a simple case
print("\n=== Test 2: Test expand_annotation ===")

# Create a simple test case
test_row = pd.Series({
    '1.50': 'N24',
    '2.50': 'D52',
    '3.50': 'R105',
    '4.50': 'W132',
    '5.50': 'P178',
    '6.48': 'W307',
    '7.50': 'P348'
})

query_seq = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRTPLNYILLNLAVADLFMVFGGFTTTLYTSLHGYFVFGPTGCNLEGFFATLGGEIALWSLVVLAIERYVVVCKPMSNFRFGENHAIMGVAFTWVMALACAAPPLVGWSRYIPEGMQCSCGIDYYTPHEETNNESFVIYMFVVHFIIPLIVIFFCYGQLVFTVKEAAAQQQESATTQKAEKEVTRMVIIMVIAFLICWLPYAGVAFYIFTHQGSDFGPIFMTIPAFFAKTSAVYNPVIYIMMNKQFRNCMVTTLCCGKNPLGDDEASTTVSKTETSQVAPA"

# Create a mock alignment (query aligned to itself for simplicity)
alignment = [query_seq, '|' * len(query_seq), query_seq]

# Run expand_annotation
grn_list, rn_list, missing = expand_annotation(
    test_row, 
    query_seq, 
    alignment,
    protein_family='gpcr_a',
    verbose=2
)

print(f"\nExpansion results:")
print(f"Total annotated: {len(grn_list)}")
print(f"Missing positions: {len(missing)}")

# Check for invalid GRN formats
invalid_grns = [g for g in grn_list if '.' in g and len(g.split('.')[1]) > 3]
if invalid_grns:
    print(f"\nWARNING: Found invalid GRN formats: {invalid_grns}")
else:
    print("\nAll GRN formats are valid!")

# Check for loop annotations
loop_grns = [g for g in grn_list if '.' in g and len(g.split('.')[0]) == 2]
print(f"\nLoop GRNs found: {len(loop_grns)}")
if loop_grns:
    print(f"Examples: {loop_grns[:5]}")