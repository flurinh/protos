# CLI Scripts Update Summary

## Date: 2025-06-25

This document summarizes the updates made to the GRN CLI scripts to use the new `GRNBaseProcessor` instead of the legacy `GRNProcessor`.

## Overview

All GRN CLI scripts have been updated to:
1. Use `GRNBaseProcessor` for data loading/saving
2. Support configurable data roots via `--data-root` argument
3. Handle both dot notation (3.50) and x notation (3x50)
4. Provide better documentation with examples

## Updated Scripts

### 1. assign_grns.py

**Key Changes:**
- Now uses `GRNBaseProcessor` with proper path configuration
- Added `--data-root` argument for custom data directory
- Handles both dot and x notation for Schiff base checking
- Uses dataset IDs with subdirectory support (e.g., `ref/mo_ref`)

**New Usage:**
```bash
# Assign GRNs using default data root
python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_dataset

# Use custom data root
python -m protos.cli.grn.assign_grns -p gpcr_a -d gpcr_dataset --data-root /path/to/data

# Use more cores
python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_dataset -n 16
```

### 2. clean_grn_table.py

**Key Changes:**
- Added support for dataset IDs vs file paths
- New `--use-files` flag to work with file paths
- Uses `GRNBaseProcessor` by default for dataset operations

**New Usage:**
```bash
# Clean using dataset IDs
python -m protos.cli.grn.clean_grn_table -i ref/mo_ref -o ref/mo_ref_cleaned

# Clean using file paths
python -m protos.cli.grn.clean_grn_table -i input.csv -o output.csv --use-files

# With custom data root
python -m protos.cli.grn.clean_grn_table -i my_table -o my_table_clean --data-root /path/to/data
```

### 3. diagnose_grn.py

**Key Changes:**
- Uses `GRNBaseProcessor` for loading tables
- Handles both dot and x notation for Schiff base checking
- Added `--use-file` flag for direct file operations
- Improved diagnostic reporting

**New Usage:**
```bash
# Diagnose using dataset ID
python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t ref/mo_ref

# Diagnose with file path
python -m protos.cli.grn.diagnose_grn -p gpcr_a -t grn_table.csv --use-file

# Save diagnostic report
python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t mo_table -o report.json

# Skip certain checks
python -m protos.cli.grn.diagnose_grn -p gpcr_a -t ref_table --no-loops --no-schiff
```

### 4. expand_annotation.py

**Note:** This file is a function implementation, not a CLI script. No changes were needed.

### 5. visualize_grn.py

**Key Changes:**
- Uses `GRNBaseProcessor` for loading tables
- Handles both dot and x notation for helix coloring
- Added `--use-file` flag for direct file operations
- Fixed sorting to work with both notation formats

**New Usage:**
```bash
# Visualize using dataset ID
python -m protos.cli.grn.visualize_grn -t ref/mo_ref -o grn_heatmap.png

# Visualize with file path
python -m protos.cli.grn.visualize_grn -t grn_table.csv -o heatmap.png --use-file

# Show distribution
python -m protos.cli.grn.visualize_grn -t ref_table -o distribution.png -v distribution

# Filter specific GRNs and proteins
python -m protos.cli.grn.visualize_grn -t my_table -o filtered.png -g 3.50,7.50 -p BR,1UAZ
```

## Common Features Across All Scripts

1. **Path Configuration:**
   - All scripts now respect `PROTOS_DATA_ROOT` environment variable
   - `--data-root` argument overrides the environment variable
   - Proper absolute path handling

2. **Notation Compatibility:**
   - Scripts handle both dot notation (3.50) and x notation (3x50)
   - Automatic detection of notation format in data
   - Conversion between formats where needed

3. **Dataset vs File Operations:**
   - Default behavior uses dataset IDs with `GRNBaseProcessor`
   - Optional flags (`--use-file`, `--use-files`) for direct file operations
   - Consistent interface across all scripts

4. **Improved Documentation:**
   - All scripts now have comprehensive help text
   - Examples in epilog sections
   - Clear parameter descriptions

## Migration Guide

For users upgrading from the old CLI scripts:

1. **Data Directory:** Ensure `PROTOS_DATA_ROOT` is set or use `--data-root`
2. **Reference Tables:** Now use subdirectory paths like `ref/mo_ref` instead of direct paths
3. **Output Paths:** The scripts create directories as needed for output files
4. **Notation:** Scripts automatically handle both notations, no changes needed

## Testing the Updates

To test the updated scripts:

```bash
# Set up environment
export PROTOS_DATA_ROOT=/path/to/your/data

# Test assignment
python -m protos.cli.grn.assign_grns -p microbial_opsins -d test_dataset

# Test cleaning
python -m protos.cli.grn.clean_grn_table -i ref/mo_ref -o ref/mo_ref_clean

# Test diagnosis
python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t ref/mo_ref

# Test visualization
python -m protos.cli.grn.visualize_grn -t ref/mo_ref -o test_viz.png
```

## Known Issues

1. The GRN parser still internally expects x notation - this needs to be updated in a future release
2. Some terminal GRN positions (like '-001', '101') may show parsing warnings
3. The expand_annotation.py file is a function implementation, not a CLI script

## Future Work

1. Update GRN parser to natively handle dot notation
2. Add support for batch processing in CLI scripts
3. Improve error messages for common issues
4. Add progress bars for long-running operations