# Utility Function Refactoring Plan

This document outlines the plan for refactoring key utility functions to use the standard schemas and interfaces defined in the `schema` module.

## Phase 2: Refactoring Priority Utility Functions

### Structure Utilities

1. **`load_structure` in struct_utils.py**

   **Current Issues:**
   - Uses hardcoded column mappings and transformations
   - Performs validation and conversion implicitly
   - Has limited error handling
   
   **Refactoring Plan:**
   - Update to use `STRUCTURE_CORE_COLUMNS` and `STRUCTURE_EXTENDED_COLUMNS` from schema_definitions
   - Add explicit validation using `validate_structure_df`
   - Implement as method for `StructureInterface` class
   - Use `cif_to_structure_df` from conversion_utilities
   - Improve error handling and logging

2. **`normalize_structures` in struct_utils.py**

   **Current Issues:**
   - Uses implicit type conversion
   - Has mixed responsibilities
   - Lacks validation
   
   **Refactoring Plan:**
   - Split into smaller, focused functions
   - Add schema validation at input and output
   - Use type annotations for clarity
   - Ensure consistency with column naming standards

3. **`get_ca_ret_coords` in struct_utils.py**

   **Current Issues:**
   - Uses indexing assumptions for DataFrames
   - Has limited error handling
   - Uses direct column access without validation
   
   **Refactoring Plan:**
   - Add parameter validation
   - Use standard column naming
   - Add schema validation for input and output
   - Improve error handling for missing data

### GRN Utilities

1. **`parse_grn_str2float` and `parse_grn_float2str` in grn_utils_updated.py**

   **Current Issues:**
   - Some redundancy between the two implementations
   - Error handling could be improved
   - Limited validation of inputs
   - **CRITICAL**: Loop region format implementation is incorrect
   
   **Refactoring Plan:**
   - Use `GRN_PATTERNS` from schema_definitions for validation
   - Improve error handling for edge cases
   - Standardize return values
   - Add comprehensive docstrings
   - **Fix loop region formatting** (see detailed plan below)

2. **`check_str_grn_valid` in grn_utils_updated.py**

   **Current Issues:**
   - Uses hardcoded regex patterns
   - Complex logic with limited documentation
   - Multiple overlapping validations
   - **CRITICAL**: Loop pattern validation is inconsistent
   
   **Refactoring Plan:**
   - Use `GRN_PATTERNS` from schema_definitions
   - Simplify validation logic
   - Add clearer error messages
   - Return validation details for debugging
   - **Fix loop pattern validation** (see detailed plan below)

3. **`get_grn_interval` in grn_utils_updated.py**

   **Current Issues:**
   - Complex function with multiple responsibilities
   - Limited validation
   - Uses direct DataFrame access without schema validation
   
   **Refactoring Plan:**
   - Split into smaller focused functions
   - Add schema validation using `validate_grn_table`
   - Add parameter validation
   - Implement as method for `GRNInterface` class
   - Improve error messages

### Sequence Utilities

1. **`init_aligner` in seq_alignment.py**

   **Current Issues:**
   - Hardcoded parameters
   - Limited configuration options
   - Direct matrix manipulation without abstractions
   
   **Refactoring Plan:**
   - Make parameters configurable
   - Use abstraction for matrix manipulation
   - Add validation for configuration
   - Improve error handling

2. **`contains_only_valid_residues` in seq_alignment.py**

   **Current Issues:**
   - Uses hardcoded set of valid residues
   - Limited error reporting
   - Limited extensibility
   
   **Refactoring Plan:**
   - Use schema-defined validation rules
   - Add better error reporting
   - Implement as method for `SequenceInterface` class
   - Make extensible for different residue sets

3. **`msa_blosum62` in seq_alignment.py**

   **Current Issues:**
   - Complex function with multiple responsibilities
   - Limited validation
   - Uses multiple internal functions
   
   **Refactoring Plan:**
   - Split into smaller focused functions
   - Add schema validation for input and output
   - Standardize alignment result format
   - Return results in schema-compatible format
   - Implement as method for `SequenceInterface` class

## GRN Loop Region Format Refactoring

### 1. Current Issues with GRN Loop Region Formatting

The current implementation of GRN utilities has a critical inconsistency in how loop regions are represented:

- **Current Format**: Varies between `[1-8][1-8]x[0-9][0-9]*` (e.g., "12x05") and imprecise decimal notations
- **Required Format**: `[1-8][1-8].[0-9]{3}` where:
  - First digit: Closer helix number (1-8)
  - Second digit: Further away helix number (1-8)
  - Three-digit decimal: Distance from closer helix (001-999)
  - Examples: "12.003" (loop between helix 1-2, closer to 1, distance 3), "65.011" (loop between helix 5-6, closer to 6, distance 11)

This inconsistency affects parsing, validation, and assignment of GRNs for loop regions.

### 2. Revised Schema Definition

Update the `GRN_PATTERNS` in `schema_definitions.py`:

```python
# Current definition
GRN_PATTERNS = {
    'standard': r'^(\d+)x(\d+)$',  # e.g., 1x50
    'n_term': r'^n\.(\d+)$',       # e.g., n.10
    'c_term': r'^c\.(\d+)$',       # e.g., c.5
    'loop': r'^(\d+)\.(\d+)$'      # e.g., 2.45
}

# Updated definition
GRN_PATTERNS = {
    'standard': r'^(\d+)x(\d+)$',  # e.g., 1x50
    'n_term': r'^n\.(\d+)$',       # e.g., n.10
    'c_term': r'^c\.(\d+)$',       # e.g., c.5
    'loop': r'^([1-8])([1-8])\.(\d{3})$'  # e.g., 12.003, 65.011
}
```

Add descriptive documentation to clarify the loop format:

```python
# Documentation for GRN formats
GRN_FORMAT_DOCS = {
    'standard': "Standard GRN format: <helix>x<position> (e.g., 1x50)",
    'n_term': "N-terminal format: n.<position> (e.g., n.10)",
    'c_term': "C-terminal format: c.<position> (e.g., c.5)",
    'loop': """Loop region format: <closer helix><further helix>.<distance> where:
            - First digit: Closer helix (1-8)
            - Second digit: Further helix (1-8)
            - Three-digit decimal: Distance from closer helix (001-999)
            Examples: 12.003 (between helix 1-2, closer to 1, distance 3)
                     65.011 (between helix 5-6, closer to 6, distance 11)"""
}
```

### 3. Revised Parsing Functions

Refactor `parse_grn_str2float` and `parse_grn_float2str` in both `grn_utils.py` and `grn_utils_updated.py`:

#### 3.1 Parsing GRN Strings to Floats

The float representation of loop GRNs should encode:
- Which helices are involved (first two digits of integer part)
- Which helix is closer (encoded in the floating point value)
- Distance from closer helix (encoded in the decimal part)

```python
# Current issue:
# For loop GRNs like "12.003", the function interprets this as helix 12, position 3/1000
# This is incorrect - it should be loop between helix 1-2, closer to 1, distance 3

# Corrected Implementation:
def parse_grn_str2float(grn: str) -> float:
    """
    Convert a GRN string to its float representation.
    
    Handles:
    - Standard: '1x50' -> 1.50
    - N-terminal: 'n.10' -> -0.10
    - C-terminal: 'c.5' -> 100.05
    - Loop: '12.003' -> 12.003 (loop between helix 1-2, closer to 1, distance 3)
    """
    try:
        # N-terminal region
        if 'n.' in grn:
            # Parse position number
            position = int(grn.split('n.')[1])
            return -0.01 * position
            
        # C-terminal region
        elif 'c.' in grn:
            # Parse position number
            position = int(grn.split('c.')[1])
            return 100.0 + 0.01 * position
            
        # Loop region with format AB.CCC
        elif '.' in grn and len(grn.split('.')[1]) == 3:
            # Parse helix pair and distance
            helix_pair = grn.split('.')[0]
            distance = int(grn.split('.')[1]) / 1000.0
            
            if len(helix_pair) == 2:
                # Extract closer and further helix
                closer_helix = int(helix_pair[0])
                further_helix = int(helix_pair[1])
                
                # Create float representation:
                # Integer part: 10*smaller_helix + larger_helix
                # Decimal part: distance value / 1000
                return float(f"{min(closer_helix, further_helix)}{max(closer_helix, further_helix)}.{int(distance*1000):03d}")
            else:
                raise ValueError(f"Invalid loop format: {grn}, expected format: AB.CCC")
                
        # Standard GRN format (TM regions)
        elif 'x' in grn:
            # Parse helix and position
            helix_str, position_str = grn.split('x')
            helix = int(helix_str)
            position = int(position_str)
            return helix + position / 100.0
            
        # Invalid or unrecognized format
        else:
            raise ValueError(f"Unrecognized GRN format: {grn}")
            
    except (ValueError, IndexError) as e:
        # Log the error with the GRN string for debugging
        logging.error(f"Error parsing GRN string '{grn}': {e}")
        return 0.0
```

#### 3.2 Converting Floats to GRN Strings

```python
def parse_grn_float2str(grn_float: float) -> str:
    """
    Convert a GRN float representation to its string format.
    
    Handles:
    - Standard: 1.50 -> '1x50'
    - N-terminal: -0.10 -> 'n.10'
    - C-terminal: 100.05 -> 'c.5'
    - Loop: 12.003 -> '12.003' (loop between helix 1-2, closer to 1, distance 3)
    """
    # Round to 3 decimal places to avoid floating point issues
    grn_float = round(grn_float, 3)
    
    # N-terminal region (negative values)
    if grn_float < 0:
        # Convert to n.XX format
        position = int(abs(grn_float) * 100)
        return f"n.{position}"
    
    # C-terminal region (100+)
    elif grn_float >= 100:
        # For values like 100.05, convert to c.5 format
        position = int(round((grn_float - 100) * 100))
        return f"c.{position}"
    
    # Loop region (values between 10 and 100)
    elif grn_float >= 10:
        # Extract the parts
        int_part = int(grn_float)
        decimal_part = round((grn_float - int_part) * 1000)
        
        # Get helix numbers
        helix1 = int(int_part / 10)  # First digit
        helix2 = int(int_part % 10)  # Second digit
        
        # Format with proper zero padding for the distance
        return f"{helix1}{helix2}.{decimal_part:03d}"
    
    # Standard transmembrane region
    else:
        # Split into helix and position parts
        helix = int(grn_float)
        position = int(round((grn_float - helix) * 100))
        
        # Format with proper zero padding
        return f"{helix}x{position:02d}"
```

### 4. Validation Functions

Refactor the `validate_grn_string` function to properly validate loop GRNs:

```python
def validate_grn_string(grn: str) -> Tuple[bool, str]:
    """Validate a GRN string against standard patterns."""
    # Empty or None check
    if not grn:
        return False, "Empty or None GRN string"
    
    # Loop format validation
    loop_pattern = re.compile(GRN_PATTERNS['loop'])
    if loop_pattern.match(grn):
        # Parse the components
        helix_pair = grn.split('.')[0]
        distance_str = grn.split('.')[1]
        
        # Validate helices
        closer_helix = int(helix_pair[0])
        further_helix = int(helix_pair[1])
        
        # Validate helix values
        if not (1 <= closer_helix <= 8):
            return False, f"Invalid closer helix: {closer_helix} (expected 1-8)"
        if not (1 <= further_helix <= 8):
            return False, f"Invalid further helix: {further_helix} (expected 1-8)"
            
        # Validate distance format (3 digits)
        if len(distance_str) != 3:
            return False, f"Invalid distance format: {distance_str} (expected 3 digits)"
            
        # Validate distance value
        distance = int(distance_str)
        if not (0 <= distance <= 999):
            return False, f"Invalid distance: {distance} (expected 0-999)"
            
        return True, "Valid loop GRN format"
    
    # Other GRN format validation...
```

### 5. GRN Assignment Functions

Update the assignment functions to use the correct loop format:

```python
def _annotate_missing_rns(interval, present_seq_nr_grn_list, query_seq, grn_config, grns_str):
    # Loop handling section:
    for seqnr in n_interval:
        closest, min_dist = _get_closest_present_seqnr(seqnr, present_seq_nr_grn_list, loop_side='n')
        region = closest[1][0]
        # Form the loop GRN with the closest helix and next helix
        next_region = str(int(region) + 1)
        # Distance from closest helix (format as three digits)
        distance = f"{abs(min_dist):03d}"
        # Create loop GRN in format: <closer helix><further helix>.<distance>
        grn = f"{region}{next_region}.{distance}"
        
        if grn not in known_grns:
            seq_id = query_seq[seqnr - 1] + str(seqnr)
            nloop.append((seq_id, grn))
```

### 6. Backward Compatibility

To ensure backward compatibility with existing data:

```python
def normalize_grn_format(grn: str) -> str:
    """
    Normalize a GRN string to the standardized format.
    
    Converts legacy formats to the new standard:
    - '12x05' -> '12.005' (loop with x notation)
    - '12.5' -> '12.005' (loop without zero padding)
    - '1.2' -> '1x20' (standard GRN with dot instead of x)
    
    Args:
        grn: GRN string to normalize
        
    Returns:
        Normalized GRN string
    """
    # Already in standard format
    if re.match(GRN_PATTERNS['standard'], grn) or \
       re.match(GRN_PATTERNS['n_term'], grn) or \
       re.match(GRN_PATTERNS['c_term'], grn) or \
       re.match(GRN_PATTERNS['loop'], grn):
        return grn
    
    # Legacy loop format with x
    if re.match(r'^([1-8])([1-8])x(\d+)$', grn):
        helix_pair = grn.split('x')[0]
        distance_str = grn.split('x')[1]
        # Convert to new format with three-digit distance
        return f"{helix_pair}.{int(distance_str):03d}"
    
    # Legacy loop format without zero padding
    if re.match(r'^([1-8])([1-8])\.(\d+)$', grn) and len(grn.split('.')[1]) < 3:
        helix_pair = grn.split('.')[0]
        distance_str = grn.split('.')[1]
        # Convert to new format with three-digit distance
        return f"{helix_pair}.{int(distance_str):03d}"
    
    # Standard GRN with dot instead of x
    if re.match(r'^([1-8])\.(\d+)$', grn) and len(grn.split('.')[0]) == 1:
        helix = grn.split('.')[0]
        position = grn.split('.')[1]
        # Convert to standard x format
        return f"{helix}x{int(position):02d}"
    
    # Return as is if we can't normalize
    return grn
```

### 7. Implementation Priority Order

1. **Update Schema Definitions**
   - First update `GRN_PATTERNS` in `schema_definitions.py`
   - Add documentation for the loop format
   - Add validation helper functions

2. **Implement Normalization Functions**
   - Create `normalize_grn_format()` for backward compatibility
   - Test with various legacy formats

3. **Update Parsing Functions**
   - Update `parse_grn_str2float()` to handle the loop format correctly
   - Update `parse_grn_float2str()` to generate correct loop GRNs
   - Ensure consistency across all utilities

4. **Update Assignment Functions**
   - Modify `_annotate_missing_rns()` and related functions
   - Ensure correct loop GRN generation during assignment

5. **Update Tests**
   - Create tests for the new loop format
   - Test backward compatibility
   - Compare results with expected behavior

### 8. Testing Strategy

1. **Unit Tests**
   - Test each parsing function with:
     - Standard format GRNs
     - N-terminal/C-terminal GRNs
     - Loop GRNs in the new format
     - Legacy format loop GRNs

2. **Integration Tests**
   - Test full assignment workflow with sample sequences
   - Compare results with expected GRN assignments
   - Test migration of existing GRN tables

3. **Edge Cases**
   - Test with all possible helix combinations (1-8)
   - Test with minimum and maximum distances
   - Test with invalid inputs and ensure proper error handling

### 9. Documentation Updates

1. **Add Format Documentation**
   - Clear explanation of the loop format
   - Examples of valid and invalid GRNs
   - Migration guide for legacy formats

2. **Update Docstrings**
   - Add Examples section to all parsing functions
   - Document float representation logic

3. **Add Deprecation Warnings**
   - Warn when encountering legacy formats
   - Log format conversions

## Implementation Steps

1. **For Each Function:**
   - Create new version in appropriate module
   - Add parameter validation with type hints
   - Add schema validation for input/output
   - Add docstrings following standard format
   - Implement tests for the new version
   - Add deprecation warning to old version

2. **For Integration:**
   - Update processors to use new functions
   - Add validation calls for data
   - Ensure backward compatibility
   - Document usage patterns

## Testing Strategy

1. **Unit Tests for New Functions:**
   - Test with valid inputs
   - Test with invalid inputs
   - Test with edge cases
   - Compare results with old functions

2. **Integration Tests:**
   - Use existing pipelines with new functions
   - Verify results match expected outputs
   - Check performance impact

## First Refactoring Target

We'll start with refactoring the GRN parsing functions since they're foundational to many other operations:

1. `parse_grn_str2float` and `parse_grn_float2str` 
2. `check_str_grn_valid`
3. `normalize_grn_format` (new function for backward compatibility)

These functions are relatively self-contained and will provide immediate benefits when standardized, especially with the corrected loop format handling.