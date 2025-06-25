# CLI Testing Status and Current Intent

## Date: 2025-06-25

## Executive Summary

We've successfully updated all GRN CLI scripts to use the new `GRNBaseProcessor` architecture and have begun testing them with real data. The CLI tools are now functional but revealed some issues that need addressing to make them truly flexible and useful for researchers.

## Current Intent

Our goal is to make the Protos CLI tools **flexible and powerful** for real-world use:

1. **Quick FASTA Processing**: Researchers should be able to quickly load FASTA files and assign GRNs
2. **Protein Family Detection**: Automatically detect protein family from sequence similarity when not specified
3. **Robust Workflow**: Handle edge cases gracefully and provide meaningful feedback
4. **Real Data Focus**: Always test with real biological data, not just synthetic test cases

## Testing Results

### ✅ What's Working

1. **Directory Structure**: Successfully set up proper data directory structure
   ```
   data/
   ├── grn/
   │   ├── ref/          # Reference GRN tables
   │   ├── tables/       # Calculated GRN tables
   │   ├── datasets/     # Output datasets
   │   └── configs/      # Configuration files
   ├── sequence/
   │   └── fasta/
   │       └── processed/  # Input FASTA files
   └── structure/
       └── structure_dataset/
   ```

2. **CLI Tools Status**:
   - **diagnose_grn.py**: ✅ Working - Successfully diagnoses GRN tables
   - **assign_grns.py**: ⚠️ Runs but needs better test data
   - **clean_grn_table.py**: ✅ Updated but not tested
   - **visualize_grn.py**: ✅ Updated but not tested

### ⚠️ Issues Found

1. **BaseProcessor Bug**: 
   - `load_data()` returns file path string instead of loaded data for subdirectory datasets
   - Workaround added to GRNBaseProcessor to handle this case

2. **Test Data Quality**:
   - Current test sequences (BR_test, HR_test, ChR2_test) are too divergent from reference
   - No strict GRNs found in alignment - need better test sequences

3. **Path Assumptions**:
   - Fixed: FASTA path was `fasta/processed/` but should be `sequence/fasta/processed/`

## Test Execution Log

### 1. Setup Data Directory
```bash
python setup_protos_data.py --data-root data --download-datasets
```
✅ Successfully created directory structure and copied reference data

### 2. Test Diagnose CLI
```bash
export PROTOS_DATA_ROOT=data
python -m protos.cli.grn.diagnose_grn -p microbial_opsins -t ref/mo_ref
```
✅ Output:
- 124 proteins, 119 GRN positions
- All positions in dot notation (1.44, 1.45, etc.)
- No invalid GRNs or issues found

### 3. Test Assign GRNs CLI
```bash
python -m protos.cli.grn.assign_grns -p microbial_opsins -d mo_test -n 2
```
⚠️ Result:
- MMseqs2 alignment worked
- Found matches but "No strict GRNs found in alignment"
- Need better test sequences that are closer to known microbial opsins

## Next Steps

### Immediate Actions
1. **Create Better Test Data**: 
   - Extract actual sequences from known microbial opsins (1UAZ, 3DDL, etc.)
   - These will have higher similarity to reference sequences

2. **Test Remaining CLI Tools**:
   - Test clean_grn_table.py with a table that needs cleaning
   - Test visualize_grn.py to generate plots

3. **Fix BaseProcessor Issue**:
   - The current workaround works but the underlying issue should be fixed
   - Need to investigate why load_data returns path for subdirectory datasets

### Future Enhancements (As Requested)
1. **Flexible FASTA Loading**:
   - Add command to process any FASTA file without pre-registration
   - Support drag-and-drop or file path input

2. **Auto-detect Protein Family**:
   - Run sequence against all reference families
   - Report best match and confidence score
   - Use detected family for GRN assignment

3. **Batch Processing**:
   - Process multiple FASTA files in one command
   - Generate summary report of all assignments

4. **Direct Workflow**:
   ```bash
   # Dream CLI usage:
   protos grn assign my_sequences.fasta --auto-detect-family --output my_grns.csv
   ```

## Code Changes Made

1. **GRNBaseProcessor** (line 238-242):
   ```python
   # Handle case where load_data returns a file path instead of data
   if isinstance(df, str):
       self.logger.debug(f"load_data returned path: {df}, loading directly")
       import pandas as pd
       df = pd.read_csv(df, index_col=0, low_memory=low_memory, **kwargs)
   ```

2. **assign_grns.py** (line 174):
   ```python
   # Fixed path from 'fasta/processed' to 'sequence/fasta/processed'
   fasta_path = data_root / 'sequence' / 'fasta' / 'processed' / f'{dataset}.fasta'
   ```

## Conclusion

The CLI tools are now functional with the new architecture, but we need:
1. Better test data that matches reference sequences
2. More flexible commands for real-world use
3. Fix for the BaseProcessor subdirectory issue

The foundation is solid - we just need to build on it to make the tools truly useful for researchers who want quick, flexible GRN assignment capabilities.