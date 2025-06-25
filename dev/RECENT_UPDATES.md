# Recent Updates to Protos Library

## Date: 2025-06-25

### Summary of Changes

This document summarizes the significant changes made to the Protos library, focusing on GRN-Structure integration, test improvements, and bug fixes.

## 1. GRN Notation Update

### Key Change: Dot Notation is Now Standard
- The GRN system uses **dot notation** (e.g., `1.50`, `7.50`) as the standard format
- The older 'x' notation (e.g., `1x50`, `7x50`) is deprecated and will be removed
- All new code should use dot notation exclusively

## 2. Structure-GRN Integration

### New Methods Added to CifBaseProcessor

#### `get_seq_dict()` 
Extracts protein sequences from structure data.
```python
sequences = processor.get_seq_dict()
# Returns: {'1UAZ_A': 'MLELLPTAVEGVSQ...', '1UAZ_B': '...'}
```

#### `get_grn_dict()`
Extracts GRN annotations from structure data (if present).
```python
grn_dict = processor.get_grn_dict()
# Returns: {
#   '1UAZ': {
#     'A': {
#       '1.50': 'R82',
#       '7.50': 'K216'
#     }
#   }
# }
```

#### `assign_grns()`
Automatically assigns GRN annotations based on sequence similarity.
```python
grn_assignments = processor.assign_grns(
    protein_family='microbial_opsins',
    similarity_threshold=0.2,
    use_mmseqs=True,
    save_results=True
)
```

### Bug Fixes in Integration
1. **Import Error Fixed**: Removed non-existent `GRNAssignment` class import
2. **Path Resolution Fixed**: GRN reference tables now correctly load from `ref/` subdirectory
3. **Column Compatibility**: Added support for both `auth_comp_id` and `res_name3l` columns
4. **Sequence Extraction**: Fixed to properly get unique residues using `drop_duplicates`
5. **Attribute Error**: Fixed `data_dirs` attribute error by using correct `data_path`

## 3. Testing Improvements

### New Test Files Created

#### `test_grn_processor_with_real_data.py`
- Comprehensive tests for GRN processor using real reference data
- 8 tests covering all major functionality
- Uses actual microbial opsin reference data from `test-data/grn/ref/mo_grn.csv`

#### `test_cifbase_grn_integration.py`
- Integration tests for GRN-Structure functionality
- Tests sequence extraction, GRN assignment, and mapping

#### `examples/test_grn_struct_step_by_step.py`
- Step-by-step demonstration of the full workflow
- Downloads real structures, extracts sequences, assigns GRNs, and maps back

### Test Data Organization

Three separate data sources are now properly distinguished:
1. **reference_data**: In `src/protos/reference_data/` - Reference datasets
2. **test-data**: In `tests/test-data/` - Test fixtures mimicking real architecture
3. **data/**: Production data directory for CLI scripts and real use cases

## 4. Core Library Fixes

### BaseProcessor Path Handling
Fixed PosixPath vs string comparison issues in `register_dataset` method:
```python
# Before
if file_path.startswith(self.data_path):

# After  
if file_path.startswith(str(self.data_path)):
```

### GRN Processor Updates
1. Fixed `save_grn_table` method signature issues
2. Clarified `filter_data_by_occurances` expects number of proteins, not fraction
3. Updated `get_grn_dict` documentation to clarify it returns list of positions

## 5. Documentation Updates

### New Documentation
- `docs/GRN_STRUCTURE_INTEGRATION.md` - Comprehensive guide for integration features
- `docs/RECENT_UPDATES.md` - This file, tracking all recent changes

### Updated Documentation
- Updated existing docs to reflect dot notation as standard
- Added troubleshooting sections for common issues
- Included real-world examples with actual protein data

## 6. CLI Scripts Updates

### ✅ All CLI Scripts Updated (2025-06-25)
All GRN CLI scripts have been successfully migrated to use `GRNBaseProcessor`:

#### Updated Scripts:
1. **assign_grns.py**
   - Now uses GRNBaseProcessor with proper path configuration
   - Added `--data-root` argument for custom data directories
   - Handles both dot and x notation for Schiff base checking
   - Uses dataset IDs with subdirectory support (e.g., `ref/mo_ref`)

2. **clean_grn_table.py**
   - Added support for dataset IDs vs file paths
   - New `--use-files` flag to work with file paths
   - Uses GRNBaseProcessor by default for dataset operations

3. **diagnose_grn.py**
   - Uses GRNBaseProcessor for loading tables
   - Handles both dot and x notation in diagnostics
   - Added `--use-file` flag for direct file operations
   - Improved diagnostic reporting

4. **expand_annotation.py**
   - Note: This is a function implementation file, not a CLI script
   - No changes were needed

5. **visualize_grn.py**
   - Uses GRNBaseProcessor for loading tables
   - Handles both dot and x notation for helix coloring
   - Added `--use-file` flag for direct file operations
   - Fixed sorting to work with both notation formats

#### Common Features Added:
- All scripts respect `PROTOS_DATA_ROOT` environment variable
- `--data-root` argument overrides environment variable
- Consistent interface across all scripts
- Comprehensive help text with examples
- Support for both dataset IDs and file paths

See [CLI_UPDATES.md](CLI_UPDATES.md) for detailed usage examples.

## 7. CLI Testing Results (2025-06-25)

### Testing Setup
- Created proper data directory structure in `data/`
- Set up reference GRN tables and test FASTA files
- All CLI tools now use local `data/` directory

### Test Results
1. **diagnose_grn.py**: ✅ Working perfectly
   - Successfully diagnosed mo_ref table with 124 proteins
   - Detected all GRNs in dot notation format

2. **assign_grns.py**: ⚠️ Functional but needs better test data
   - MMseqs2 alignment working
   - Test sequences too divergent - no strict GRNs found
   - Need sequences closer to reference opsins

### Issues Discovered
1. **BaseProcessor Bug**: `load_data()` returns file path for subdirectory datasets
   - Added workaround in GRNBaseProcessor
   - Needs proper fix in BaseProcessor

2. **Path Issues**: Fixed FASTA path from `fasta/processed/` to `sequence/fasta/processed/`

See [CLI_TESTING_STATUS.md](CLI_TESTING_STATUS.md) for detailed testing log.

## 8. Known Issues and Future Work

### To Be Fixed
1. GRN parser still expects 'x' notation internally - needs update to handle dot notation natively
2. CLI scripts need migration to new processor architecture
3. Some terminal GRN positions (like '-001', '101') show parsing warnings

### Successfully Tested Workflows
1. **Structure → GRN**: Extract sequences and assign GRNs ✓
2. **GRN → Structure**: Map GRN annotations to structure residues ✓
3. **Real Data Integration**: Works with actual PDB structures ✓

## 8. Example Usage

### Complete Working Example
```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
import os
from pathlib import Path

# Set up paths
data_dir = Path("/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# Initialize processor
processor = CifBaseProcessor(
    name="mo_analysis",
    data_root=str(data_dir.absolute())
)

# Load structures
processor.load_structure("1uaz")  # Bacteriorhodopsin

# Extract sequences
sequences = processor.get_seq_dict()
print(f"Extracted {len(sequences)} sequences")

# Assign GRNs
grn_assignments = processor.assign_grns(
    protein_family='microbial_opsins',
    similarity_threshold=0.2,
    use_mmseqs=True
)

# Get GRN annotations
grn_dict = processor.get_grn_dict()

# Find Schiff base lysine
k750 = processor.data[processor.data['grn'] == '7.50']
print(f"Found {len(k750)} atoms at position 7.50")
```

## Conclusion

The Protos library now has robust GRN-Structure integration with comprehensive testing using real data. The main remaining work is updating the CLI scripts to use the new architecture and completing the transition to dot notation throughout the codebase.