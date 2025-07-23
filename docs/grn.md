# Generic Residue Numbering (GRN) System

## Overview

The Generic Residue Numbering (GRN) system provides a standardized way to refer to equivalent positions across protein families, particularly useful for GPCRs and other membrane proteins. This system allows researchers to compare structures and sequences even when the actual residue numbers differ.

## GRN Format Specifications

### Standard Transmembrane (TM) Positions

**Format**: `<helix>.<position>`

- `<helix>`: Single digit (1-8) representing the helix number
- `<position>`: Two-digit position within the helix (01-99)

**Examples**:
- `1.50`: Position 50 in helix 1 (most conserved position)
- `3.32`: Position 32 in helix 3 (DRY motif in GPCRs)
- `7.50`: Position 50 in helix 7 (NPxxY motif in GPCRs)

### N-terminal Positions

**Format**: `n.<position>`

- `<position>`: Distance from the start of TM1

**Examples**:
- `n.1`: First residue before TM1
- `n.10`: 10th residue before TM1

### C-terminal Positions

**Format**: `c.<position>`

- `<position>`: Distance from the end of the last helix (TM7 or H8)

**Examples**:
- `c.1`: First residue after the last helix
- `c.25`: 25th residue after the last helix

### Loop Positions

**Format**: `<nearest_helix><distant_helix>.<distance>`

- `<nearest_helix>`: The helix the residue is closer to
- `<distant_helix>`: The other helix bounding the loop
- `<distance>`: Three-digit distance from the nearest helix (001-999)

**Examples**:
- `12.003`: 3rd residue from helix 1 in the loop between helices 1 and 2
- `34.001`: 1st residue from helix 3 in the loop between helices 3 and 4
- `43.002`: 2nd residue from helix 4 in the loop between helices 3 and 4

### Loop Annotation Pattern

For a loop with 5 residues between helices 3 and 4:
```
34.001  (1st residue from helix 3)
34.002  (2nd residue from helix 3)
34.003  (3rd residue from helix 3)
43.002  (2nd residue from helix 4)
43.001  (1st residue from helix 4)
```

## Special GRN Formats in Reference Data

### Intermediate Positions

In reference alignment tables, positions between standard GRNs may use additional decimal places:

- `2.551`: Position between 2.55 and 2.56
- `5.461`: Position between 5.46 and 5.47

These are valid in pre-aligned reference data but should not be generated for new annotations.

## Protein Family Configurations

### GPCR Class A (gpcr_a)

Standard GRN ranges:
- **TM1**: 1.28 - 1.64
- **TM2**: 2.31 - 2.71  
- **TM3**: 3.20 - 3.60
- **TM4**: 4.34 - 4.69
- **TM5**: 5.36 - 5.69
- **TM6**: 6.24 - 6.61
- **TM7**: 7.30 - 7.56
- **H8**: 8.37 - 8.72

Strict (conserved) positions:
- **TM1**: 1.49 - 1.59
- **TM2**: 2.37 - 2.50
- **TM3**: 3.42 - 3.51
- **TM4**: 4.47 - 4.60
- **TM5**: 5.47 - 5.50
- **TM6**: 6.43 - 6.50
- **TM7**: 7.45 - 7.53
- **H8**: 8.47 - 8.58

### Microbial Opsins (mo)

Standard GRN ranges:
- **TM1**: 1.36 - 1.61
- **TM2**: 2.39 - 2.67
- **TM3**: 3.41 - 3.62
- **TM4**: 4.35 - 4.61
- **TM5**: 5.39 - 5.68
- **TM6**: 6.33 - 6.59
- **TM7**: 7.33 - 7.62

Note: Microbial opsins typically don't have an H8 helix.

## GRN Assignment Algorithm

The GRN assignment process follows these steps:

1. **Initial Alignment**: Align query sequence to reference sequences with known GRN annotations
2. **Strict Position Transfer**: Transfer GRN numbers for conserved positions
3. **N-terminal Assignment**: Assign n.1, n.2, etc. for residues before TM1
4. **C-terminal Assignment**: Assign c.1, c.2, etc. for residues after the last helix
5. **Standard GRN Assignment**: Fill in missing standard GRN positions within TM helices
6. **Loop/Gap Annotation**: Assign loop positions for residues between TM helices

## Usage Examples

### Getting GRN Interval
```python
from protos.processing.grn.grn_utils import get_grn_interval

# Auto-generate all positions in TM1
tm1_grns = get_grn_interval("1.28", "1.64")
# Returns: ['1.28', '1.29', '1.30', ..., '1.64']

# Filter from a provided list
ref_grns = ['1.25', '1.28', '1.50', '1.551', '1.64', '1.70']
tm1_filtered = get_grn_interval("1.28", "1.64", grns_str=ref_grns)
# Returns: ['1.28', '1.50', '1.64']
```

### Annotating a Sequence
```python
from protos.processing.grn import GRNProcessor

# Load reference GRN table
grnp = GRNProcessor(name="gpcr_grn")
grnp.load_grn_table("gpcrdb_ref")

# Annotate a new sequence
result = grnp.annotate_sequence(
    name="my_receptor",
    sequence="MNASW...VIFI",
    protein_family="gpcr_a"
)
```

## Best Practices

1. **Use Standard Notation**: Always use dot notation (not the deprecated 'x' notation)
2. **Validate GRN Strings**: Use `validate_grn_string()` to check format validity
3. **Sort GRN Lists**: Use `sort_grns_str()` for proper helix-loop-helix ordering
4. **Handle Missing Data**: Use '-' for gaps in GRN tables
5. **Preserve Reference Formats**: Don't modify GRN formats in reference data (e.g., keep '2.551')

## Common Issues and Solutions

### Issue: Invalid GRN format generated
**Solution**: Ensure using the updated functions that generate proper dot notation

### Issue: Loops not annotated
**Solution**: Check that the gap between helices is properly identified and use the correct loop format

### Issue: Missing standard GRNs
**Solution**: Verify the protein family configuration is loaded and contains the expected ranges

## See Also

- [GRN Processor Documentation](processors/grn_processor.md)
- [Annotation Processing](processing/annotations.md)
- [Sequence Alignment](processing/alignments.md)