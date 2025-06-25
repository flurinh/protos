# GRN Assignment Workflow Documentation

## Overview
The GRN (Generic Residue Numbering) assignment system is a critical component of the Protos library that enables consistent residue numbering across protein families. This system assigns standardized numbers to residues based on structural alignment to reference sequences.

## File Architecture

### Directory Structure
```
data/grn/
├── ref/                    # Reference GRN tables (protein family specific)
│   ├── curated_grn.csv    # Curated reference GRN table
│   ├── gpcrdb_ref.csv     # GPCR reference table
│   ├── mo_ref.csv         # Microbial opsins reference table
│   └── ref.csv            # General reference table
├── tables/                # Output GRN tables from assignment process
├── assignments/           # Intermediate assignment results
└── configs/              # Configuration files

src/protos/reference_data/grn/
├── config.json           # Secondary structure boundaries (strict/standard)
├── ref/                  # Reference tables (fallback location)
└── registry.json         # Registry of available GRN tables
```

### Configuration Files

#### config.json Structure
Contains secondary structure boundary definitions for different protein families:
- **gpcr_a**: Class A GPCRs with TM1-TM7 + H8
- **microbial_opsins**: Microbial opsins with TM1-TM7
- **iLBP**: Intracellular lipid-binding proteins

Each family has:
- **standard**: Broader boundaries for secondary structure elements
- **strict**: Core conserved regions only

Example:
```json
{
  "gpcr_a": {
    "standard": {
      "tm1": ["1x28", "1x64"],
      "tm2": ["2x31", "2x71"],
      ...
    },
    "strict": {
      "tm1": ["1x49", "1x59"],
      "tm2": ["2x37", "2x50"],
      ...
    }
  }
}
```

## GRN Table Format

GRN tables use a wide format with:
- **Row index**: Protein/structure IDs (e.g., "7BMH", "1ABC")
- **Columns**: GRN positions in dot notation (e.g., "3.50", "7.53")
- **Values**: Residue+position strings (e.g., "K270", "M62")

Example:
```
        1.50   2.50   3.50   7.53
7BMH    M62    I90    R115   K270
1ABC    L22    I50    R82    K225
4HYJ    -      V50    R80    Y240
```

## Assignment Process

### 1. Input Processing
- **Input**: FASTA file(s) with query sequences
- **Reference**: Protein family-specific GRN reference table
- **Config**: Secondary structure boundaries for the family

The reference tables contain pre-assigned GRN numbers for well-characterized proteins:
- **GPCR**: Uses GPCRdb reference numbering
- **Microbial Opsins**: Uses microbial opsin-specific numbering
- **Format**: Wide table with protein IDs as rows, GRN positions as columns

### 2. Initial Sequence Search (MMseqs2)
```python
# Find best matching reference sequence using fast sequence search
output = mmseqs2_align2(query_seqs=query_dict, ref_seqs=ref_dict)
best_matches = output.loc[output.groupby('query_id')['e_value'].idxmin()]
```

**Key Parameters**:
- E-value threshold for matches
- Coverage requirements
- Identity percentage cutoffs

### 3. Pairwise Alignment (BLOSUM62)
```python
# Perform detailed pairwise alignment with best reference match
aligner = init_aligner(open_gap_score=-22)
alignment = align_blosum62(query_seq, ref_seq, aligner)
formatted_alignment = format_alignment(alignment)
```

**Alignment Format**:
```
Query:    MLELLPTAVEGVSQAQITGRPEWIWLALGTALMGLGTLYFLVK
          |||||||| ||||||||||||||||||||||||||||||||||
Reference: MLELLPTA-EGVSQAQITGRPEWIWLALGTALMGLGTLYFLVK
```

### 4. Initial GRN Transfer Based on Strict Boundaries
```python
# Create position mapping from reference
ref_row = grnp.data.loc[best_match, :]
ref_dict = {grn: res for grn, res in ref_row.to_dict().items() if res != '-'}
seq_pos2grn = dict([(i + 1, grn) for i, grn in enumerate(list(ref_dict.keys()))])

# Transfer GRNs from reference based on alignment
new_row = init_row_from_alignment(alignment, seq_pos2grn)

# Filter by strict GRN positions only
grn_config_strict = config.get_config(strict=True)
grns_str_strict = init_grn_intervals(grn_config_strict)
new_row_strict = new_row.loc[grns_str_strict]
```

**Strict vs Standard Boundaries**:
- **Strict**: Core conserved positions only (e.g., 1x49-1x59 for TM1)
- **Standard**: Extended boundaries including variable regions (e.g., 1x28-1x64 for TM1)

### 5. Annotation Expansion Process
The `expand_annotation` function performs comprehensive GRN assignment in multiple stages:

#### 5.1 Correctly Aligned GRNs
```python
# Get GRNs that align well with the query
aligned_grns = get_correctly_aligned_grns(
    all_query_gene_numbers, aligned_grns, alignment, max_alignment_gap=1
)
```

#### 5.2 N-terminal Annotation
```python
# Assign negative numbers to N-terminal residues
n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(
    aligned_grns, missing_gene_numbers, grns_float_std
)
# Results in: n.1, n.2, n.3, etc.
```

#### 5.3 C-terminal Annotation
```python
# Assign high numbers to C-terminal residues
c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(
    aligned_grns, missing_gene_numbers, query_gene_len, grns_float_std
)
# Results in: c.1, c.2, c.3, etc.
```

#### 5.4 Missing Standard GRNs
```python
# Fill in missing conserved positions within standard boundaries
missing_std_grns = [grn for grn in grns_str_std if grn not in present_grns]
grns, missing = assign_missing_std_grns(
    missing_std_grns, present_seq_nr_grn_list, query_seq, missing, grns_str_std
)
```

#### 5.5 Loop and Gap Annotation
```python
# Annotate loops between TM helices and gaps
nloop, gaps, cloop = annotate_gaps_and_loops(
    expanded_grn_list, missing, query_seq, grn_config_std, grns_str_std
)
# Results in decimal positions: 12.001, 12.002, 23.001, etc.
```

**Loop Naming Convention**:
- `12.xxx`: Loop between TM1 and TM2
- `23.xxx`: Loop between TM2 and TM3
- Sequential numbering within each loop

### 6. Quality Control and Validation

#### 6.1 Sequential Numbering Check
```python
def is_sequential(result):
    """Check if residue numbers are sequential."""
    current = 1
    for value in result.values():
        if int(value[1:]) != current:
            return False
        current += 1
    return True
```

#### 6.2 Conserved Residue Validation
```python
# Check family-specific conserved residues
if protein_family == 'microbial_opsins':
    # Schiff base lysine at 7x50
    df = df[df['7x50'].str.contains('K')]
elif protein_family == 'gpcr_a':
    # DRY motif arginine
    df = df[df['3x50'].str.contains('R')]
```

#### 6.3 Alignment Quality Metrics
- **Coverage**: Percentage of reference positions aligned
- **Identity**: Percentage of identical residues in alignment
- **Gap Penalty**: Number and size of gaps introduced

### 7. Output Generation and Format
```python
# Convert results to wide format GRN table
df = pd.DataFrame.from_dict(final_results, orient='index')
# Sort columns by GRN position
df = df.loc[:, sort_grns_str(df.columns.tolist())]
# Save with protein IDs as index
df.to_csv(f"data/grn/tables/{dataset}.csv", index=True)
```

**Output Format Example**:
```
protein_id,1.50,2.50,3.50,...,7.50
QUERY001,M22,I90,R115,...,K296
QUERY002,L22,V90,R115,...,K296
```

### 8. Error Detection and Handling

#### 8.1 Validation by Re-assignment
A robust validation approach:
1. Extract strict positions from reference sequence
2. Align full query to strict-only reference
3. Run GRN assignment
4. Compare with expected results

```python
# Extract strict positions only
strict_table = extract_strict_residues(reference_table, strict_config)
strict_seq = get_sequence_from_grn_row(strict_table.iloc[0])

# Align and assign
alignment = align_blosum62(query_seq, strict_seq, aligner)
new_row = init_row_from_alignment(alignment, seq_pos2grn)

# Validate assignments match expected
for grn, assigned_val in new_row.items():
    expected_val = reference_table.loc[ref_id, grn]
    assert assigned_val[0] == expected_val[0]  # Check amino acid matches
```

#### 8.2 Common Failure Modes
1. **No suitable reference found**: E-value too high
2. **Poor alignment quality**: Too many gaps or mismatches
3. **Missing conserved residues**: Indicates wrong family or mutant
4. **Non-sequential numbering**: Gaps in critical regions
5. **Boundary violations**: Assignments outside expected ranges

#### 8.3 Error Reporting
```python
# Log detailed error information
if len(missing) > 0:
    logger.warning(f"Missing residues in {query_id}: {missing}")
if not is_sequential(grn_dict):
    logger.error(f"Non-sequential numbering in {query_id}")
```

## Key Functions

### Core Assignment Functions
- `assign_grns()`: Main entry point for batch assignment
- `expand_annotation()`: Core GRN expansion algorithm
- `init_row_from_alignment()`: Transfer GRNs based on alignment
- `annotate_gpcr()`: Single sequence annotation

### Utility Functions
- `calculate_missing_gene_numbers()`: Find gaps in numbering
- `assign_missing_std_grns()`: Fill standard GRN positions
- `annotate_gaps_and_loops()`: Handle non-conserved regions
- `sort_grns_str()`: Sort GRN positions correctly

## Usage Examples

### Command Line
```bash
# Assign GRNs to microbial opsins dataset
python -m protos.cli.grn.assign_grns \
    -p microbial_opsins \
    -d Bacteriorhodopsins \
    -n 8
```

### Python API
```python
from protos.cli.grn.assign_grns import assign_grns

# Assign GRNs to a dataset
result_df = assign_grns(
    protein_family='gpcr_a',
    dataset='my_gpcrs',
    num_cores=8
)
```

### Single Sequence
```python
from protos.processing.grn.grn_processor import GRNProcessor
from protos.processing.grn.grn_table_utils import annotate_gpcr

# Load reference GRN table
grnp = GRNProcessor('ref', 'data/grn/ref/')

# Annotate single sequence
annotate_gpcr(
    grnp=grnp,
    query_name='my_protein',
    query_seq='MDWLVGYGFGG...',
    protein_family='gpcr_a'
)
```

## Algorithm Details

### GRN Number Format
- **Standard positions**: "TM.BW" (e.g., "3.50" = TM3, position 50)
- **X notation**: "TMxBW" (e.g., "3x50" = same as above)
- **N-terminus**: "n.1", "n.2", etc. (parsed as -1.0, -2.0)
- **C-terminus**: "c.1", "c.2", etc. (parsed as 101.0, 102.0)
- **Loops**: Decimal positions (e.g., "45.500", "45.510")

### Boundary Constraints
- **Strict**: Core conserved positions only
- **Standard**: Extended boundaries including variable regions
- Used to determine which positions to transfer during alignment

### Gap Handling
- **Max alignment gap**: Default 1 (configurable)
- Gaps in conserved regions are annotated with decimal GRNs
- Large gaps may indicate alignment issues

## Error Handling

### Common Issues
1. **No alignment found**: Query sequence too divergent
2. **Missing conserved residues**: May indicate wrong family
3. **Sequential numbering failure**: Gaps in critical regions
4. **Empty results**: Check input format and reference availability

### Validation Steps
1. Check E-value threshold for initial alignment
2. Verify conserved residue presence (family-specific)
3. Ensure sequential numbering in output
4. Validate against known structures if available

## Performance Considerations

### Parallel Processing
- Default: 8 cores for batch processing
- Adjustable via `num_cores` parameter
- Each sequence processed independently

### Caching
- Reference sequences loaded once
- Alignment results can be cached
- GRN configurations loaded on demand

### Memory Usage
- Scales with number of sequences
- Large alignments may require more memory
- Batch processing for very large datasets

## Future Improvements

1. **Multi-family support**: Handle chimeric proteins
2. **Custom boundaries**: User-defined structure elements  
3. **Alignment refinement**: Iterative improvement
4. **Confidence scores**: Per-position assignment confidence
5. **Web service**: REST API for GRN assignment