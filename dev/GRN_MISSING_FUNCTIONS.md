# Missing GRN Assignment Functions Analysis

## Overview
After reviewing `grn_assignment_deprecated.py`, I've identified numerous critical functions missing from the current `grn_assignment.py` implementation. These functions are essential for a complete GRN assignment workflow.

## Missing Core Functions

### 1. **Alignment Processing Functions**

#### `get_correctly_aligned_grns()`
- **Purpose**: Extract correctly aligned GRN positions from a sequence alignment
- **Critical for**: Transferring GRN numbers from reference to query through alignment
- **Details**: 
  - Validates alignment jumps to avoid incorrect assignments
  - Handles gaps in both query and reference sequences
  - Uses match quality indicators ('|' for exact match, '.' for similar)

#### `valid_jump()`
- **Purpose**: Validate if a GRN assignment jump is valid
- **Used by**: `get_correctly_aligned_grns()`
- **Validates**:
  - Transmembrane helix transitions
  - Maximum allowed gaps between assignments
  - Sequence continuity

### 2. **N-terminal and C-terminal Extension Functions**

#### `calculate_missing_ntail_grns()`
- **Purpose**: Calculate GRN assignments for N-terminal residues
- **Handles**: Residues before the first aligned GRN position
- **Uses**: Linear extrapolation from first assigned position

#### `calculate_missing_ctail_grns()`
- **Purpose**: Calculate GRN assignments for C-terminal residues
- **Handles**: Residues after the last aligned GRN position
- **Uses**: Linear extrapolation from last assigned position

### 3. **Gap and Loop Annotation Functions**

#### `_get_seq_nr_intervals()`
- **Purpose**: Group consecutive missing sequence numbers into intervals
- **Used for**: Identifying continuous gaps or loops

#### `_is_valid_gap()`
- **Purpose**: Determine if missing positions represent a valid gap
- **Returns**: Information about gap validity and position assignment

#### `_check_interval_is_gap()`
- **Purpose**: Check if an interval represents a gap vs. a loop
- **Uses**: Secondary structure boundaries to make determination

#### `_annotate_missing_rns()`
- **Purpose**: Annotate missing residue numbers as gaps or loops
- **Handles**: Both internal gaps and inter-helix loops

### 4. **Utility Functions**

#### `_remove_loop_grns()`
- **Purpose**: Remove loop GRNs from a list
- **Filters**: GRNs containing 'c' or 'n' (loop indicators)

#### `_get_closest_present_seqnr()`
- **Purpose**: Find closest assigned sequence number
- **Used for**: Gap/loop position assignment

#### `_get_closest_present_grn()`
- **Purpose**: Find closest assigned GRN position
- **Used for**: Interpolating missing positions

#### `_check_for_duplicate_grns()`
- **Purpose**: Validate no duplicate GRN assignments
- **Returns**: Validation status and list of duplicates

#### `sort_grn_rn_pairs()`
- **Purpose**: Sort GRN-residue number pairs
- **Ensures**: Proper ordering of assignments

## Complete GRN Assignment Workflow

### Phase 1: Initial Alignment
1. **Load reference GRN table** with known assignments
2. **Perform sequence alignment** (query vs reference)
3. **Extract aligned positions** using `get_correctly_aligned_grns()`
   - Validate jumps with `valid_jump()`
   - Handle alignment gaps

### Phase 2: Terminal Extensions
1. **N-terminal extension**
   - Use `calculate_missing_ntail_grns()`
   - Extrapolate from first aligned position
   - Respect TM1 boundary constraints

2. **C-terminal extension**
   - Use `calculate_missing_ctail_grns()`
   - Extrapolate from last aligned position
   - Respect TM7/H8 boundary constraints

### Phase 3: Missing Standard GRNs
1. **Identify missing standard positions**
   - Compare assigned vs expected GRNs
   - Use `assign_missing_std_grns()`
   
2. **Pivot-based assignment**
   - Use `_is_valid_gap()` to check availability
   - Assign based on sequence proximity
   - Validate against boundaries

### Phase 4: Gap and Loop Annotation
1. **Group missing positions**
   - Use `_get_seq_nr_intervals()`
   - Classify as gaps or loops

2. **Annotate each interval**
   - Use `_annotate_missing_rns()`
   - Apply appropriate numbering scheme:
     - Gaps: increment by 0.001
     - N-loops: use (TM)(TM+1)x(distance)
     - C-loops: use (TM)(TM-1)x(distance)

### Phase 5: Validation and Cleanup
1. **Check for duplicates** with `_check_for_duplicate_grns()`
2. **Sort assignments** with `sort_grn_rn_pairs()`
3. **Apply boundary constraints** from config
4. **Remove invalid assignments**

## Critical Missing Components

### 1. **Alignment Quality Control**
- Match quality assessment
- Gap penalty handling
- Alignment score thresholds

### 2. **Boundary Constraint Enforcement**
- Strict vs standard boundaries
- Family-specific constraints
- Helix length validation

### 3. **Loop Detection Algorithm**
- Distinguish gaps from loops
- Apply correct numbering scheme
- Handle edge cases

### 4. **Interpolation Logic**
- Linear interpolation for gaps
- Distance-based loop numbering
- Pivot selection algorithm

## Implementation Priority

1. **High Priority** (Core functionality):
   - `get_correctly_aligned_grns()`
   - `valid_jump()`
   - `calculate_missing_ntail_grns()`
   - `calculate_missing_ctail_grns()`

2. **Medium Priority** (Gap/Loop handling):
   - `_annotate_missing_rns()`
   - `_is_valid_gap()`
   - `_get_seq_nr_intervals()`

3. **Low Priority** (Utilities):
   - Sorting functions
   - Validation utilities
   - Helper functions

## Testing Requirements

Each function needs tests for:
- Normal cases
- Edge cases (empty sequences, no alignment)
- Boundary conditions
- Error handling
- Integration with full workflow