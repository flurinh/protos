# GRN System: Generic Residue Numbering

This document provides a comprehensive overview of the Generic Residue Numbering (GRN) system implemented in the Protos framework, including utilities, processors, and workflows for annotating protein sequences.

## Table of Contents

1. [Introduction to GRNs](#introduction-to-grns)
2. [GRN Notation Formats](#grn-notation-formats)
3. [GRN Utility Functions](#grn-utility-functions)
4. [GRN Assignment Utilities](#grn-assignment-utilities)
5. [GRN Processor Classes](#grn-processor-classes)
6. [Complete Workflow](#complete-workflow)
7. [Implementation Details](#implementation-details)
8. [Usage Examples](#usage-examples)
9. [Troubleshooting](#troubleshooting)

## Introduction to GRNs

Generic Residue Numbering (GRN) is a standardized positional reference system for protein structures, particularly useful for membrane proteins like GPCRs and microbial rhodopsins. The system allows for consistent referencing of structurally equivalent positions across different proteins even when their sequences vary significantly.

### Benefits of GRN System

- **Structural Consistency**: Identifies structurally equivalent positions regardless of sequence variations
- **Cross-Protein Comparison**: Enables direct comparison between different proteins
- **Evolutionary Analysis**: Facilitates tracking of conserved positions across evolutionary distances
- **Standardized Communication**: Provides a common language for describing protein positions in literature

### Historical Context

The GRN system has evolved from the Ballesteros-Weinstein numbering system for GPCRs, which has been extended to other membrane protein families. The implementation in Protos refines this system with more precise notation for loop regions and comprehensive coverage of entire proteins.

## GRN Notation Formats

The GRN system uses distinct notation formats for different regions of proteins:

### Transmembrane Region Format: `NxYY`

- **Format**: `<helix_number>x<position>`
- **Example**: `3x50` (helix 3, position 50)
- **Range**: Typically 1-8 for helix number, 01-99 for position
- **Float Representation**: 3.50 (helix + position/100)

### N-Terminal Format: `n.XX`

- **Format**: `n.<distance>`
- **Example**: `n.10` (10 positions before first transmembrane helix)
- **Float Representation**: -0.10 (negative values for N-terminal)

### C-Terminal Format: `c.XX`

- **Format**: `c.<distance>`
- **Example**: `c.5` (5 positions after last transmembrane helix)
- **Float Representation**: 100.05 (values over 100 for C-terminal)

### Loop Region Format: `AB.CCC`

- **Format**: `<closer_helix><further_helix>.<distance>`
- **Example**: `12.003` (loop between helix 1-2, closer to helix 1, distance 3)
- **Components**:
  - First digit: Closer helix number (1-8)
  - Second digit: Further helix number (1-8)
  - Three-digit decimal: Distance from closer helix (001-999)
- **Float Representation**: 12.003 (helix pair + distance/1000)

### Legacy Format Variations

The codebase historically used various formats for loop regions:
- `12x05` (loop using x notation)
- `12.5` (loop without zero padding)
- Various other inconsistent formats

The current implementation standardizes these to the `AB.CCC` format while maintaining backward compatibility.

## GRN Utility Functions

Core utility functions are defined in `protos.processing.schema.grn_utils_updated`, providing the fundamental operations for working with GRNs:

### `parse_grn_str2float(grn: str) -> float`

Converts a GRN string to its float representation:

```python
# Examples
parse_grn_str2float('1x50') == 1.50   # Standard TM region
parse_grn_str2float('n.10') == -0.10  # N-terminal
parse_grn_str2float('c.5') == 100.05  # C-terminal
parse_grn_str2float('12.003') == 12.003  # Loop region
```

Implementation details:
- N-terminal: Converts to negative values (`n.10` → -0.10)
- C-terminal: Adds to 100 (`c.5` → 100.05)
- TM regions: Decimal conversion (`1x50` → 1.50)
- Loop regions: Preserves format with precise decimal (`12.003` → 12.003)

### `parse_grn_float2str(grn_float: float) -> str`

Converts a float representation back to GRN string format:

```python
# Examples
parse_grn_float2str(1.50) == '1x50'    # Standard TM region
parse_grn_float2str(-0.10) == 'n.10'   # N-terminal
parse_grn_float2str(100.05) == 'c.5'   # C-terminal
parse_grn_float2str(12.003) == '12.003'  # Loop region
```

Implementation details:
- Handles rounding to avoid floating-point precision issues
- Formats with appropriate leading zeros where needed
- Preserves the three-digit precision for loop regions

### `normalize_grn_format(grn: str) -> str`

Converts legacy GRN formats to the standardized notation:

```python
# Examples
normalize_grn_format('12x05') == '12.005'  # Legacy loop format
normalize_grn_format('12.5') == '12.005'   # Non-zero-padded loop
normalize_grn_format('1.50') == '1x50'     # Dot instead of x in TM region
```

Implementation details:
- Uses regular expressions to identify format patterns
- Applies appropriate conversions based on recognized patterns
- Preserves already-standardized formats

### `validate_grn_string(grn: str) -> Tuple[bool, str]`

Validates a GRN string against standard patterns:

```python
# Examples
validate_grn_string('1x50') == (True, "Valid standard GRN format")
validate_grn_string('9x50') == (False, "Invalid helix number: 9 (expected 1-8)")
validate_grn_string('12.003') == (True, "Valid loop GRN format")
```

Implementation details:
- Checks against regular expression patterns from schema definitions
- Applies format-specific validation rules
- Provides detailed error messages for invalid formats
- Attempts normalization for legacy formats

### `sort_grns(grn_list: List[Union[str, float]]) -> List[Union[str, float]]`

Sorts a list of GRNs in standard order:

```python
# Examples
sort_grns(['3x50', '1x50', 'n.10', 'c.5', '12.003']) 
# Results in: ['n.10', '1x50', '3x50', '12.003', 'c.5']
```

Implementation details:
- Converts all GRNs to float representation for consistent sorting
- Orders N-terminal (negative values) first
- Followed by TM regions in numerical order
- Then loop regions
- Finally C-terminal regions (values > 100)

## GRN Assignment Utilities

Assignment utilities in `protos.processing.schema.grn_assignment_utils` provide specialized functions for mapping GRNs to protein sequences:

### `assign_gene_nr(seq: str) -> List[str]`

Assigns position numbers to amino acids in a sequence:

```python
# Example
assign_gene_nr("MEAKL") == ["M1", "E2", "A3", "K4", "L5"]
```

### `get_closest_present_seqnr(missing_seqnr: int, present_seq_nr_grn_list: List[Tuple[str, str]], loop_side: str) -> Tuple[Tuple[str, str], int]`

Finds the closest annotated residue to a missing position:

```python
# present_pairs = [("A10", "1x50"), ("L20", "2x50"), ("K30", "3x50")]
get_closest_present_seqnr(25, present_pairs, loop_side='n')
# Returns: (("L20", "2x50"), -5)  # Closest on N-terminal side, distance -5
```

### `annotate_loop_region(interval: List[int], present_seq_nr_grn_list: List[Tuple[str, str]], query_seq: str) -> List[Tuple[str, str]]`

Annotates residues in a loop region between two helices:

```python
# interval = [20, 21, 22]  # Residue positions to annotate
# present_pairs = residues with known GRNs
annotate_loop_region(interval, present_pairs, query_seq)
# Returns: [("A20", "12.003"), ("V21", "12.004"), ("G22", "12.005")]
```

Implementation details:
- Identifies the helices bordering the loop
- Determines which helix is closer to each loop residue
- Calculates appropriate distance values
- Formats GRNs using the standardized AB.CCC format

### `calculate_missing_ntail_grns(aligned_grns: Dict[str, str], missing_gene_numbers: List[str], grns_float: List[float]) -> Tuple[List[Tuple[str, str]], int]`

Assigns GRNs to N-terminal residues:

```python
# Missing positions before first aligned residue
calculate_missing_ntail_grns(aligned_grns, missing_gene_numbers, grns_float)
# Returns: N-terminal GRN pairs and first gene number
```

Implementation details:
- Calculates distance from first aligned position
- Assigns negative values (n.XX format) to N-terminal residues
- Maintains proper ordering relative to first transmembrane position

### `calculate_missing_ctail_grns(aligned_grns: Dict[str, str], missing_gene_numbers: List[str], query_gene_len: int, grns_float: List[float]) -> Tuple[List[Tuple[str, str]], Optional[int]]`

Assigns GRNs to C-terminal residues:

```python
# Missing positions after last aligned residue
calculate_missing_ctail_grns(aligned_grns, missing_gene_numbers, query_gene_len, grns_float)
# Returns: C-terminal GRN pairs and last gene number
```

Implementation details:
- Calculates distance from last aligned position
- Assigns values > 100 (c.XX format) to C-terminal residues
- Maintains proper ordering relative to last transmembrane position

### `valid_jump(prev_ref_grn: Optional[str], curr_ref_grn: str, prev_query_key: Optional[str], curr_query_key: str, max_alignment_gap: int) -> bool`

Checks if an alignment jump is valid based on GRN and sequence continuity:

```python
valid_jump("1x50", "1x51", "A10", "A11", max_alignment_gap=1) == True
valid_jump("1x50", "1x52", "A10", "A12", max_alignment_gap=1) == False
```

### `get_correctly_aligned_grns(all_query_gene_numbers: List[str], reference_grn_dict: Dict[str, str], alignment: Tuple[str, str, str], max_alignment_gap: int = 1) -> Dict[str, str]`

Extracts valid GRN mappings from a sequence alignment:

```python
# Example
get_correctly_aligned_grns(all_query_gene_numbers, reference_grn_dict, alignment)
# Returns: Dictionary mapping query positions to GRNs based on alignment
```

Implementation details:
- Parses the alignment string to match query and reference positions
- Applies valid_jump to ensure alignment continuity
- Handles gaps in alignment correctly
- Returns a clean mapping of query positions to reference GRNs

## GRN Processor Classes

The GRN processor classes provide high-level interfaces for working with GRN tables:

### GRNBaseProcessor Class

Located in `protos.processing.grn.grn_base_processor`, this is the primary class for GRN table operations:

#### Initialization

```python
processor = GRNBaseProcessor(
    name="my_processor",
    dataset="ref",  # Reference dataset to load
    data_root="/path/to/data",
    processor_data_dir="grn"
)
```

#### Key Methods

##### Loading and Saving

- `load_grn_table(dataset_id, normalize_formats=True)`: Loads a GRN table from a dataset
  - Handles various index formats
  - Optionally normalizes GRN formats
  - Creates bidirectional mappings between notation styles
  
- `save_grn_table(dataset_id, normalize_formats=True)`: Saves the current GRN table
  - Optionally normalizes GRN formats before saving
  - Handles index formatting correctly

##### Table Manipulation

- `filter_by_ids(ids_to_keep)`: Filters the table to include only specified IDs
- `apply_interval(grn_interval)`: Limits the table to specific GRN positions
- `filter_data_by_occurances(threshold)`: Keeps only GRNs that occur in at least threshold proteins
- `sort_columns()`: Sorts columns by GRN position in standard order

##### GRN Format Handling

- Built-in handling of both dot notation (3.50) and x notation (3x50)
- Automatic detection of preferred notation in datasets
- Bidirectional conversion between notation styles

##### Data Access

- `get_seq_dict()`: Gets sequences from the GRN table
- `get_grn_dict(notation=None)`: Gets a dictionary mapping proteins to their GRNs

### Legacy GRNProcessor Class

Located in `protos.processing.grn.grn_processor`, this class is maintained for backward compatibility:

#### Key Differences from GRNBaseProcessor

- Does not integrate with BaseProcessor for path handling
- Uses direct file paths rather than dataset IDs
- Lacks some of the format normalization capabilities
- Less robust error handling

## Complete Workflow

The full process of GRN assignment follows these steps:

### 1. Input Preparation

- **FASTA Sequence File**: Contains protein sequences to be annotated
- **Reference GRN Table**: Contains known GRNs for reference proteins
- **Configuration**: Protein family-specific settings

### 2. Alignment to Reference Sequences

```python
# 1. Load query sequences
query_dict = read_fasta("data/fasta/processed/my_dataset.fasta")

# 2. Load reference sequences with known GRNs
grnp = GRNBaseProcessor(dataset="ref")
ref_dict = grnp.get_seq_dict()

# 3. Perform sequence alignment
output = mmseqs2_align2(query_seqs=query_dict, ref_seqs=ref_dict)
best_matches = output.loc[output.groupby('query_id')['e_value'].idxmin()][
    ['query_id', 'target_id']].values.tolist()

# 4. Get pairwise alignments
alignments = get_pairwise_alignment(query_dict, ref_dict, best_matches)
```

### 3. Initial GRN Transfer

```python
# 5. Transfer GRNs from reference to query based on alignment
new_rows = get_aligned_grns(grnp, alignments, best_matches, grns_str_strict)
```

### 4. Expansion to Full Sequence

```python
# 6. For each query sequence, expand GRNs to all positions
def annotate_sequence(query_id, query_seq, new_row, protein_family):
    # Format the alignment
    alignment = align_blosum62(query_seq, new_row_seq)
    alignment = format_alignment(alignment)
    
    # Expand annotations to full sequence
    grns, rns, missing = expand_annotation(
        new_row, query_seq, alignment, protein_family)
    
    return query_id, dict(zip(grns, rns))
```

The `expand_annotation` function performs several crucial steps:

```python
def expand_annotation(new_row, query_seq, alignment, protein_family):
    # 1. Get initial GRNs from alignment
    aligned_grns = get_correctly_aligned_grns(
        all_query_gene_numbers, reference_grn_dict, alignment)
    
    # 2. Identify missing positions
    missing_gene_numbers = calculate_missing_gene_numbers(
        all_gene_numbers, aligned_grns)
    
    # 3. Annotate N-terminal region
    n_tail_list, first_gene_number_int = calculate_missing_ntail_grns(
        aligned_grns, missing_gene_numbers, grns_float)
    
    # 4. Annotate C-terminal region
    c_tail_list, last_gene_number_int = calculate_missing_ctail_grns(
        aligned_grns, missing_gene_numbers, query_gene_len, grns_float)
    
    # 5. Identify and fill in missing standard GRNs
    std_grns, missing = assign_missing_std_grns(
        missing_std_grns, present_seq_nr_grn_list, query_seq, missing, grns_str)
    
    # 6. Annotate loops and gaps
    nloop, gaps, cloop = annotate_gaps_and_loops(
        present_seq_nr_grn_list, missing, query_seq, grn_config, grns_str)
    
    # 7. Combine all regions and sort
    all_grns = n_tail_list + list(aligned_grns.items()) + std_grns + \
               gaps + nloop + cloop + c_tail_list
    
    # 8. Extract results
    grns = [g[1] for g in all_grns]
    rns = [g[0] for g in all_grns]
    
    return grns, rns, missing
```

### 5. Parallel Processing

```python
# 7. Process sequences in parallel
results = main_parallel_execution(
    query_dict, new_rows, protein_family=protein_family, num_cores=num_cores)
```

### 6. Output Generation

```python
# 8. Create and save the final GRN table
df = pd.DataFrame.from_dict(results, orient='index')
df = df.loc[:, sort_grns_str(df.columns.tolist())]
df.to_csv(f"data/grn/datasets/{dataset}.csv", index=True)
```

## Implementation Details

### Float Representation Logic

The internal float representation of GRNs is a key design decision that enables consistent sorting and comparison:

- **N-terminal**: Negative values (`n.10` → -0.10)
- **TM regions**: Helix + position/100 (`1x50` → 1.50)
- **Loop regions**: 10*helix1 + helix2 + distance/1000 (`12.003` → 12.003)
- **C-terminal**: 100 + position/100 (`c.5` → 100.05)

This ensures proper ordering: N-terminal → TM helices → loops → C-terminal.

### Loop Format Evolution

The loop format has evolved from:

1. **Legacy Notation**: Using x with varied precision (`12x5`, `23x01`) 
2. **Intermediate Notation**: Using unpadded decimals (`12.5`, `23.1`)
3. **Current Standard**: Using three-digit padded decimals (`12.005`, `23.001`)

The current standard clarifies which helix is closer and precisely defines the distance.

### Backward Compatibility

Backward compatibility is maintained through several mechanisms:

1. **Format Normalization**: Legacy formats are automatically converted during load/save
2. **Dual Notation Support**: Both x and dot notation are supported
3. **Deprecation Warnings**: Legacy functions show warnings but continue to function
4. **Fallback Methods**: New implementations fall back to legacy methods on error

### Performance Considerations

- **Parallel Processing**: The workflow uses ProcessPoolExecutor for parallel sequence annotation
- **Normalization Toggle**: Format normalization can be disabled for performance-critical operations
- **Caching**: Bidirectional mappings between notation styles are cached for repeated use

## Usage Examples

### Basic GRN Table Loading and Manipulation

```python
from protos.processing.grn.grn_base_processor import GRNBaseProcessor

# Initialize processor with reference dataset
processor = GRNBaseProcessor(dataset="ref")

# Get all available datasets
available_datasets = processor.list_available_datasets()
print(f"Available datasets: {available_datasets}")

# Load another dataset
processor.load_grn_table(dataset_id="my_dataset")

# Filter to specific proteins
processor.filter_by_ids(["protein1", "protein2", "protein3"])

# Filter to conserved positions
processor.filter_data_by_occurances(threshold=5)  # Keep GRNs present in at least 5 proteins

# Extract sequences
sequences = processor.get_seq_dict()

# Save the modified table
processor.save_grn_table(dataset_id="my_filtered_dataset")
```

### Running GRN Assignment for New Sequences

```python
from protos.processing.grn.run_grn_assignment import main_parallel_execution
from protos.processing.sequence.fasta_utils import read_fasta

# Load sequences
query_dict = read_fasta("my_sequences.fasta")

# Set up workflow
# (See Complete Workflow section for detailed steps)

# Process sequences
results = main_parallel_execution(
    query_dict, new_rows, protein_family="gpcr_a", num_cores=8)

# Create output table
df = pd.DataFrame.from_dict(results, orient='index')
df.to_csv("my_grn_table.csv")
```

### Working with Different GRN Formats

```python
from protos.processing.schema.grn_utils_updated import normalize_grn_format

# Normalize legacy formats
loop_grn = normalize_grn_format("12x05")  # Returns "12.005"
tm_grn = normalize_grn_format("1.50")     # Returns "1x50"

# Convert between representations
from protos.processing.schema.grn_utils_updated import parse_grn_str2float, parse_grn_float2str

float_value = parse_grn_str2float("12.003")  # Returns 12.003
grn_string = parse_grn_float2str(12.003)     # Returns "12.003"
```

## Troubleshooting

### Common Issues

#### Invalid GRN Formats

**Symptom**: Warning messages about invalid GRN formats

**Solution**: 
- Use `normalize_grn_format()` to convert legacy formats
- Ensure loop regions use the correct AB.CCC format
- Verify that loop helices are in range 1-8
- Check that loop distances use three digits

#### Missing GRN Annotations

**Symptom**: Many positions remain unannotated after expansion

**Solution**:
- Check alignment quality between query and reference
- Verify that the reference protein has comprehensive GRN annotations
- Ensure the correct protein family configuration is used
- Try different reference proteins for better coverage

#### Notation Confusion

**Symptom**: Unexpected GRN formats in output

**Solution**:
- The processor auto-detects notation style (dot vs x)
- Use the `notation` parameter in `get_grn_dict()` to force a specific notation
- Check the `using_dot_notation` attribute to see the detected style

### Debugging Tips

1. **Inspect Alignments**: Poor alignments lead to poor GRN assignment
2. **Check Reference GRNs**: Ensure reference proteins have complete annotations
3. **Validate Format**: Use `validate_grn_string()` to check format validity
4. **Log Format Conversions**: Enable logging to track format normalizations
5. **Examine Missing Positions**: The `expand_annotation()` function returns missing positions

---

*This documentation is maintained as part of the Protos framework. For updates and issues, please refer to the project repository.*