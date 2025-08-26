# GRN Processor

The GRN Processor manages Generic Residue Numbering (GRN) data, providing standardized residue numbering schemes for protein families. This enables consistent position referencing across homologous proteins despite sequence variations.

## Overview

The GRN system provides:
- Standardized position numbering for protein families (GPCRs, microbial opsins, etc.)
- Complete residue annotation including TM regions, loops, N/C-terminals
- Cross-species residue mapping and conservation analysis
- Integration with structural and sequence data
- Support for multiple numbering schemes

For detailed GRN format specifications, see the [GRN Documentation](../grn.md).

## Basic Usage

### Initialization

```python
from protos.processing.grn import GRNProcessor

# Create processor
grnp = GRNProcessor(name="grn_analysis")

# Load reference GRN table
grnp.load_grn_table("gpcrdb_ref")  # For GPCRs
# or
grnp.load_grn_table("mo_ref")      # For microbial opsins
```

### GRN Table Format

GRN tables are pandas DataFrames with:
- **Index**: Protein/sequence identifiers
- **Columns**: GRN positions (e.g., "1.50", "3.50", "n.10", "c.5", "34.001")
- **Values**: Amino acid + position (e.g., "K296", "D83", "-" for gaps)

Example structure:
```python
        1.50   2.50   3.50   ...  7.50   n.1    c.1
RHO     N55    D83    R135   ...  N302   M1     A348
OPSD    N55    D83    R135   ...  N302   M1     V348
ADRB2   N51    D79    R131   ...  Y316   M1     L412
```

## Core Operations

### Annotating New Sequences

```python
# Simple annotation
result = grnp.annotate_sequence(
    name="my_receptor",
    sequence="MNGTEGPNFYVPFSNKTGVVRSPFEYPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRT...",
    protein_family="gpcr_a"
)
```

### Detailed Annotation with Verbosity

```python
from protos.processing.grn.grn_table_utils import annotate_grnp

# Annotate with detailed progress
grn_row = annotate_grnp(
    grnp=grnp,
    query_name="NEW_RECEPTOR",
    query_seq=sequence,
    protein_family="gpcr_a",
    verbose=2  # 0=silent, 1=major steps, 2=detailed
)

# The annotation process:
# 1. Finds best reference match using MMseqs2
# 2. Transfers strict GRN positions via alignment
# 3. Assigns N-terminal positions (n.1, n.2, ...)
# 4. Assigns C-terminal positions (c.1, c.2, ...)
# 5. Fills in missing standard GRN positions
# 6. Annotates loops between helices (e.g., 34.001)
```

### Working with GRN Data

```python
# Get sequence from GRN table
sequence = grnp.get_sequence("5HT2B_HUMAN")

# Get specific positions
dry_motif = grnp.get_grn_positions("5HT2B_HUMAN", ["3.49", "3.50", "3.51"])

# Filter by GRN range
tm3_residues = grnp.filter_by_grn_range("3.20", "3.60")

# Get all sequences as dictionary
seq_dict = grnp.get_seq_dict()

# Save results
grnp.save_grn_table("my_annotated_table")
```

## Advanced Features

### GRN Configuration Management

```python
from protos.processing.grn.grn_utils import GRNConfigManager

# Load configuration for specific protein family
config = GRNConfigManager()
grn_config = config.get_config(protein_family="gpcr_a", strict=False)

# Get GRN interval for a region
from protos.processing.grn.grn_utils import get_grn_interval

# Generate all positions in TM3
tm3_grns = get_grn_interval("3.20", "3.60")  # Auto-generates with 0.01 step

# Or filter from a reference list
ref_grns = ["3.20", "3.32", "3.50", "3.551", "3.60"]
tm3_filtered = get_grn_interval("3.20", "3.60", grns_str=ref_grns)
```

### Conservation Analysis

```python
# Find conserved positions
conserved = grnp.find_conserved_positions(threshold=0.9)

# Analyze specific motifs
motifs = {
    "DRY": ["3.49", "3.50", "3.51"],
    "NPxxY": ["7.49", "7.50", "7.51", "7.52", "7.53"],
    "CWxP": ["6.47", "6.48", "6.50"]
}

for motif_name, positions in motifs.items():
    residues = grnp.get_grn_positions("RHO_HUMAN", positions)
    print(f"{motif_name}: {residues}")
```

### Structure Integration

```python
from protos.processing.structure import StructureProcessor

# Load structure with GRN assignments
sp = StructureProcessor(name="structure_grn")
sp.load_dataset("gpcr_structures")

# Get GRN assignments for structure
grn_dict = sp.get_grn_dict("5dsg_A")

# Map GRN positions to structure residues
position_3_50 = grn_dict.get("3.50")  # Returns residue number in structure
```

## Protein Family Configurations

### GPCR Class A
- **Standard range**: TM1 (1.28-1.64), TM2 (2.31-2.71), ..., H8 (8.37-8.72)
- **Strict positions**: Highly conserved positions only
- **Key motifs**: DRY (3.49-3.51), NPxxY (7.49-7.53), CWxP (6.47-6.50)

### Microbial Opsins
- **Standard range**: TM1 (1.36-1.61), TM2 (2.39-2.67), ..., TM7 (7.33-7.62)
- **No H8 helix**
- **Key positions**: Schiff base (7.43), Proton acceptor (3.28), Proton donor (3.32)

## Common Workflows

### Workflow 1: Annotate Multiple Sequences

```python
# Load sequences
sequences = {
    "receptor1": "MNGTEGPNF...",
    "receptor2": "MEGNPNYFT...",
    "receptor3": "MPPGWNNTA..."
}

# Annotate all sequences
annotated_table = pd.DataFrame()

for name, seq in sequences.items():
    row = grnp.annotate_sequence(
        name=name,
        sequence=seq,
        protein_family="gpcr_a"
    )
    annotated_table = pd.concat([annotated_table, row.to_frame().T])

# Save complete table
grnp.data = annotated_table
grnp.save_grn_table("my_family_grn")
```

### Workflow 2: Opsin Classification

```python
# Classify opsins based on key GRN positions
def classify_opsin(protein_name):
    row = grnp.data.loc[protein_name]
    
    # Check key functional positions
    if row["3.28"] == "D" and row["7.43"] == "K":
        if row["3.32"] == "D":
            return "proton_pump"
        elif row["3.32"] == "T":
            return "chloride_pump"
    elif row["3.28"] == "E":
        return "channel"
    
    return "unknown"

# Apply to all opsins
opsin_types = {name: classify_opsin(name) for name in grnp.data.index}
```

### Workflow 3: Export for Analysis

```python
# Export specific regions
tm_regions = []
for i in range(1, 8):
    start = f"{i}.{grn_config[f'tm{i}'][0].split('.')[1]}"
    end = f"{i}.{grn_config[f'tm{i}'][1].split('.')[1]}"
    tm_cols = grnp.filter_columns_by_range(start, end)
    tm_regions.append(grnp.data[tm_cols])

# Export conserved positions only
conserved_pos = grnp.find_conserved_positions(0.95)
conserved_table = grnp.data[conserved_pos]
conserved_table.to_csv("conserved_positions.csv")
```

## Best Practices

1. **Always specify protein family**: Different families have different GRN configurations
   ```python
   # Good
   result = grnp.annotate_sequence(name="seq1", sequence=seq, protein_family="gpcr_a")
   
   # Bad - will use default or fail
   result = grnp.annotate_sequence(name="seq1", sequence=seq)
   ```

2. **Handle missing positions**: Not all positions can be assigned
   ```python
   grn_list, rn_list, missing = expand_annotation(...)
   if missing:
       print(f"Warning: {len(missing)} positions could not be annotated")
   ```

3. **Use verbosity for debugging**: Set verbose=2 to understand annotation decisions
   ```python
   result = annotate_grnp(grnp, "test", seq, verbose=2)
   ```

4. **Validate GRN formats**: Use proper dot notation
   ```python
   from protos.processing.grn.grn_utils import validate_grn_string
   
   is_valid, message = validate_grn_string("3.50")  # Valid
   is_valid, message = validate_grn_string("3x50")  # Invalid (old format)
   ```

## Troubleshooting

### Issue: No standard GRNs found
**Solution**: Check that the protein family configuration is loaded correctly
```python
config = GRNConfigManager()
available_families = config.configs.keys()
print(f"Available families: {available_families}")
```

### Issue: Poor annotation coverage
**Solution**: Try different reference sequences or protein families
```python
# Try microbial opsin family if GPCR annotation fails
result_mo = grnp.annotate_sequence(name="seq1", sequence=seq, protein_family="mo")
```

### Issue: Invalid GRN formats in output
**Solution**: Update to latest version that uses dot notation exclusively
```python
# Check for old x notation
if 'x' in grn_string and '.' not in grn_string:
    print("Warning: Old GRN format detected")
```

## Summary

The GRN Processor provides:
- Complete residue annotation for protein families
- Flexible configuration for different numbering schemes
- Integration with structure and sequence analysis
- Conservation and motif analysis capabilities
- Standardized data format for cross-study comparison

For more details on GRN formats and specifications, see the [GRN Documentation](../grn.md).