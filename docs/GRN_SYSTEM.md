# GRN System Overview

## What is GRN?

Generic Residue Numbering (GRN) is a standardized system for numbering residues across protein families, enabling consistent structural comparison and analysis. In Protos, GRN is primarily used for GPCRs and microbial opsins.

## GRN Format

### Current Standard: Dot Notation
As of the latest version, Protos uses **dot notation** as the standard format:
- **Helix positions**: `TM.BW` (e.g., "3.50" = Transmembrane helix 3, Ballesteros-Weinstein position 50)
- **N-terminus**: `n.1`, `n.2`, etc. (sequential from N-terminal)
- **C-terminus**: `c.1`, `c.2`, etc. (sequential from C-terminal)
- **Loop regions**: `<closer_helix><further_helix>.<distance>` (e.g., "12.003" = loop between TM1 and TM2, 3 residues from closest helix)

### Deprecated: X Notation
- **X notation**: `TMxBW` (e.g., "3x50") - This format is deprecated and will be removed in future versions
- The system still accepts x notation for backward compatibility but converts it internally

## Core Components

### 1. **GRNProcessor** (`grn_processor.py`)
- Loads and manages GRN reference tables
- Provides sequence extraction and validation
- Handles GRN table I/O operations

### 2. **GRNBaseProcessor** (`grn_base_processor.py`)
- Extends BaseProcessor with GRN-specific functionality
- Manages protein family configurations
- Integrates with Protos path system

### 3. **Assignment Module** (`grn_assignment.py`, `grn_table_utils.py`)
- Core assignment logic for mapping GRNs to new sequences
- Handles N/C-terminal annotations
- Fills gaps and assigns loop regions

### 4. **Configuration System**
- **config.json**: Defines secondary structure boundaries per protein family
- **Standard mode**: Extended boundaries including variable regions
- **Strict mode**: Core conserved positions only

## GRN Assignment Workflow

```
1. Input Sequence
      ↓
2. MMseqs2 Search → Find best reference match
      ↓
3. BLOSUM62 Alignment → Detailed pairwise alignment
      ↓
4. Initial Transfer → Map GRNs from reference based on alignment
      ↓
5. Expand Annotation → Fill gaps, annotate terminals and loops
      ↓
6. Quality Control → Validate sequential numbering
      ↓
7. Output GRN Table → Wide format with protein IDs as rows, GRNs as columns
```

## Reference Tables

Located in `data/grn/ref/`:
- **curated_grn.csv**: Manually curated reference
- **gpcrdb_ref.csv**: GPCR reference from GPCRdb
- **mo_ref.csv**: Microbial opsins reference
- **ref.csv**: General reference table

## Table Format

Wide format with:
- **Rows**: Protein/structure IDs
- **Columns**: GRN positions
- **Values**: Residue+position strings (e.g., "K270")

Example:
```
        1.50   2.50   3.50   7.53
7BMH    M62    I90    R115   K270
1ABC    L22    I50    R82    K225
```

## Key Functions

### Assignment Functions
- `assign_grns()`: Main entry point for batch processing
- `annotate_gpcr()`: Single sequence annotation
- `expand_annotation()`: Core expansion logic
- `init_row_from_alignment()`: Initialize GRN row from alignment

### Utility Functions
- `parse_grn_str2float()`: Convert GRN string to sortable float
- `sort_grns()`: Sort GRN positions correctly
- `validate_grn_string()`: Validate GRN format

## Current Implementation Notes

1. **Mock implementations**: Some functions in `grn_assignment.py` have simplified implementations that should be replaced with schema version logic
2. **Path dependencies**: Complex path resolution can cause issues - always use explicit paths
3. **Protein families**: Currently supports `gpcr_a`, `microbial_opsins`, and `iLBP`
4. **Performance**: Uses caching for repeated alignments to improve speed

## Usage Example

```python
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_table_utils import annotate_gpcr

# Initialize processor
grnp = GRNProcessor(name="mo_processor", protein_family="microbial_opsins")
grnp.load_grn_table("mo_ref")

# Annotate a single sequence
grn_row = annotate_gpcr(
    grnp=grnp,
    query_name="my_protein",
    query_seq="MGLDAVCKVLVDDDD...",
    protein_family="microbial_opsins"
)
```