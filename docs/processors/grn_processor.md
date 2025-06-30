# GRNBaseProcessor

The GRNBaseProcessor manages Generic Residue Numbering (GRN) data, providing standardized residue numbering schemes for protein families. This enables consistent position referencing across homologous proteins despite sequence variations.

## Overview

GRN (Generic Residue Numbering) provides:
- Standardized position numbering for protein families
- Cross-species residue mapping
- Conservation analysis across homologs
- Integration with structural and sequence data
- Support for GPCR and other numbering schemes

## GRN Concept

GRN uses a hierarchical numbering system:
- **Helix/Region** (1-7 for GPCRs, customizable for other families)
- **Position** (typically 50 as reference point per helix)
- **Format**: `helix.position` (e.g., "3.50", "7.50")

Example:
```
Position 3.50 = Position 50 in helix/region 3
- Highly conserved DRY motif region in GPCRs
- K296 in rhodopsin
- K231 in β2-adrenergic receptor
```

## Basic Usage

### Initialization

```python
from protos.processing.grn.grn_base_processor import GRNBaseProcessor

# Create processor
gp = GRNBaseProcessor(name="grn_analysis")

# With options
gp = GRNBaseProcessor(
    name="gpcr_grn",
    preload=True,          # Load reference tables on init
    verbose=True           # Detailed logging
)
```

### Loading GRN Tables

```python
# Load existing GRN table
gp.load_grn_table("rhodopsin_family")

# Load with validation
gp.load_grn_table("microbial_opsins", validate=True)

# Check available tables
tables = gp.list_grn_tables()
```

## GRN Table Format

GRN tables are structured with:
- **Rows**: Sequence/protein identifiers
- **Columns**: GRN positions (helix.position format)
- **Values**: Residue and position (e.g., "K296", "D83")

```python
# Example GRN table structure
grn_df = pd.DataFrame({
    '1.50': ['N55', 'N55', 'T62'],      # Position 50 in TM1
    '2.50': ['D83', 'D83', 'V90'],      # Position 50 in TM2  
    '3.50': ['R135', 'R135', 'L129'],   # DRY motif region
    '4.50': ['W161', 'W161', 'W171'],   # Conserved tryptophan
    '5.50': ['P215', 'Y215', 'T205'],   # Position 50 in TM5
    '6.50': ['F261', 'F261', 'P238'],   # Position 50 in TM6
    '7.50': ['N302', 'K296', 'K257']    # NPxxY motif/Schiff base
}, index=['RHO_HUMAN', 'OPSD_HUMAN', 'BACR_HALSA'])
```

## Core Operations

### Creating GRN Tables

```python
# Create from alignment
alignment = {
    'PROT1': 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTK',
    'PROT2': 'MVLSPADKTNVKGAWGKVGGHAGEYGAEALERMFLSFPTTK',
    'PROT3': 'MVLSPADKTNVKAGWGKVGAHAGEYGAEALERMFLSFPTTK'
}

# Define key positions
key_positions = {
    '1.50': 10,  # V in all sequences
    '2.50': 20,  # W/A in sequences
    '3.50': 28,  # Y in all sequences
}

# Create GRN table
grn_table = gp.create_grn_from_alignment(alignment, key_positions)
```

### Saving GRN Tables

```python
# Save GRN table
gp.save_grn_table("my_family_grn", grn_data)

# Save with metadata
gp.save_grn_table(
    "gpcr_class_a",
    grn_data,
    metadata={
        "family": "Class A GPCRs",
        "reference": "RHO_HUMAN",
        "created": "2024-01-15",
        "positions": 285  # Number of annotated positions
    }
)
```

### Sequence Operations

```python
# Extract sequences from GRN
sequences = gp.get_seq_dict()
# {'RHO_HUMAN': 'NDRRMFLSFP...', 'OPSD_HUMAN': 'NDRYMFLSFP...'}

# Get specific position across all sequences
position_3_50 = gp.get_position_column('3.50')
# ['R135', 'R135', 'L129']

# Get conserved positions
conserved = gp.find_conserved_positions(threshold=0.9)
# ['1.50', '4.50', '6.50']  # Positions conserved in >90% sequences
```

### Analysis Functions

```python
# Conservation analysis
conservation = gp.calculate_conservation()
# DataFrame with conservation scores per position

# Position frequency
freq_3_50 = gp.get_position_frequency('3.50')
# {'R': 0.67, 'L': 0.33}

# Find variable positions
variable = gp.find_variable_positions(min_variants=3)

# Clustering by GRN patterns
clusters = gp.cluster_by_grn_patterns()
```

## Advanced Features

### GRN Annotation

```python
# Annotate new sequence with GRN
new_sequence = "MNGTEGPNFYVPFSNKTGVVRSPFEAPQYYLAEPWQFSMLAAYMFLLIMLGFPINFLTLYVTVQHKKLRT"

# Align to reference and transfer GRN
grn_annotation = gp.annotate_sequence(
    sequence=new_sequence,
    sequence_id="NEW_OPSIN",
    reference="BACR_HALSA"
)

# Add to existing table
gp.add_sequence_to_grn(
    sequence_id="NEW_OPSIN",
    grn_annotation=grn_annotation
)
```

### Motif Analysis

```python
# Define functional motifs
motifs = {
    "DRY": ["3.49", "3.50", "3.51"],
    "NPxxY": ["7.49", "7.50", "7.51", "7.52", "7.53"],
    "CWxP": ["6.47", "6.48", "6.50"]
}

# Analyze motif conservation
motif_analysis = gp.analyze_motifs(motifs)

# Find sequences with specific motif patterns
dry_intact = gp.find_sequences_with_motif(
    positions=motifs["DRY"],
    pattern=["D", "R", "Y"]
)
```

### Integration with Structures

```python
# Map GRN to structure
from protos.processing.structure.struct_base_processor import CifBaseProcessor
cp = CifBaseProcessor()

# Load structure
structure = cp.load_structure("5uig")  # Rhodopsin structure

# Map GRN positions to structure residues
grn_mapping = gp.map_grn_to_structure(
    grn_sequence_id="RHO_HUMAN",
    structure_pdb="5uig",
    chain="A"
)

# Extract GRN positions from structure
position_coords = gp.extract_grn_positions_from_structure(
    structure,
    grn_positions=["3.50", "6.50", "7.50"]
)
```

### Family-Specific Schemes

```python
# Load predefined numbering schemes
gp.load_numbering_scheme("ballesteros-weinstein")  # GPCRs
gp.load_numbering_scheme("kabat")                  # Antibodies
gp.load_numbering_scheme("imgt")                   # Immunoglobulins

# Create custom scheme
custom_scheme = {
    "helices": 7,
    "reference_positions": {
        1: 50, 2: 50, 3: 50, 4: 50, 
        5: 50, 6: 50, 7: 50
    },
    "conserved_positions": {
        "1.50": "N",
        "2.50": "D", 
        "3.50": "R",
        "7.50": "Y"
    }
}
gp.define_custom_scheme("my_family", custom_scheme)
```

## Working with Microbial Opsins

### Opsin-Specific Analysis

```python
# Load microbial opsin GRN
gp.load_grn_table("microbial_opsins")

# Analyze key functional positions
key_positions = {
    "proton_acceptor": "3.28",      # D85 in BR
    "proton_donor": "3.32",         # D96 in BR  
    "schiff_base": "7.43",          # K216 in BR
    "proton_release": "3.16"        # E194/E204 in BR
}

# Check conservation of functional residues
for name, pos in key_positions.items():
    conservation = gp.get_position_frequency(pos)
    print(f"{name} ({pos}): {conservation}")
```

### Opsin Classification

```python
# Classify opsins by GRN patterns
def classify_opsin(sequence_id):
    """Classify opsin type based on GRN patterns."""
    grn_row = gp.get_grn_for_sequence(sequence_id)
    
    # Check key positions
    if grn_row["3.28"] == "D" and grn_row["7.43"] == "K":
        if grn_row["3.32"] == "D":
            return "proton_pump"
        elif grn_row["3.32"] == "T":
            return "chloride_pump"
    elif grn_row["3.28"] == "E":
        return "channel"
    else:
        return "unknown"

# Apply classification
opsin_types = {}
for seq_id in gp.data.index:
    opsin_types[seq_id] = classify_opsin(seq_id)
```

## Best Practices

### 1. Validate GRN Tables

```python
# Always validate after creation/loading
issues = gp.validate_grn_table()
if issues:
    for issue in issues:
        print(f"Warning: {issue}")

# Check format consistency
gp.standardize_grn_format()  # Convert all to residue+number format
```

### 2. Handle Missing Positions

```python
# Check for gaps in annotation
missing = gp.find_missing_positions()

# Fill gaps with placeholder
gp.fill_missing_positions(placeholder="-")

# Or interpolate from neighbors
gp.interpolate_missing_positions(method="nearest")
```

### 3. Reference Selection

```python
# Choose appropriate reference sequence
# Usually the best-characterized member
reference_options = {
    "rhodopsins": "RHO_HUMAN",      # Bovine/human rhodopsin
    "bacteriorhodopsins": "BACR_HALSA",  # Well-studied
    "channelrhodopsins": "CHR2_CHLRE"    # ChR2
}

# Set reference for alignment
gp.set_reference_sequence(reference_options["rhodopsins"])
```

### 4. Conservation-Based Filtering

```python
# Focus on conserved positions for analysis
conserved_positions = gp.find_conserved_positions(threshold=0.8)

# Create reduced GRN table
reduced_grn = gp.data[conserved_positions]

# Use for clustering or classification
clusters = gp.cluster_sequences(
    positions=conserved_positions,
    method="hierarchical"
)
```

## Common Workflows

### Workflow 1: Family Analysis

```python
def analyze_protein_family(sequences, reference_id):
    """Complete GRN analysis for protein family."""
    
    # Create GRN from sequences
    grn_table = gp.create_grn_from_sequences(
        sequences,
        reference=reference_id
    )
    
    # Save GRN table
    gp.data = grn_table
    gp.save_grn_table("family_analysis")
    
    # Analyze conservation
    conservation = gp.calculate_conservation()
    
    # Find motifs
    motifs = gp.find_conserved_motifs(min_length=3)
    
    # Cluster subfamilies
    clusters = gp.cluster_by_grn_patterns()
    
    return {
        "grn_table": grn_table,
        "conservation": conservation,
        "motifs": motifs,
        "clusters": clusters
    }
```

### Workflow 2: Structure-Function Mapping

```python
def map_function_to_structure(grn_table, functional_data, structure_pdb):
    """Map functional data to structural positions via GRN."""
    
    # Load GRN
    gp.data = grn_table
    
    # Map functional positions
    functional_positions = {}
    for seq_id, function in functional_data.items():
        if seq_id in gp.data.index:
            # Get GRN positions for functional residues
            positions = gp.get_functional_positions(seq_id, function)
            functional_positions[seq_id] = positions
    
    # Map to structure
    cp = CifBaseProcessor()
    structure = cp.load_structure(structure_pdb)
    
    # Highlight functional positions
    for grn_pos in functional_positions.values():
        struct_residues = gp.grn_to_structure_residues(grn_pos, structure_pdb)
        # Annotate structure...
    
    return functional_positions
```

### Workflow 3: Novel Sequence Annotation

```python
def annotate_novel_sequences(novel_sequences, reference_grn_table):
    """Annotate new sequences with established GRN."""
    
    # Load reference GRN
    gp.load_grn_table(reference_grn_table)
    
    annotated = {}
    for seq_id, sequence in novel_sequences.items():
        # Find best reference
        best_ref = gp.find_closest_reference(sequence)
        
        # Transfer annotation
        annotation = gp.transfer_grn_annotation(
            sequence,
            reference=best_ref,
            method="profile_alignment"
        )
        
        annotated[seq_id] = annotation
        
        # Add to GRN table
        gp.add_sequence_to_grn(seq_id, annotation)
    
    # Save updated table
    gp.save_grn_table(f"{reference_grn_table}_extended")
    
    return annotated
```

## Troubleshooting

### Common Issues

1. **Invalid GRN format**
```python
# Check and fix format
invalid = gp.find_invalid_positions()
if invalid:
    gp.fix_grn_format(invalid)
```

2. **Alignment issues**
```python
# Use more sensitive alignment
gp.annotate_sequence(
    sequence,
    method="hmmer",  # More sensitive than pairwise
    e_value=0.001
)
```

3. **Missing reference positions**
```python
# Handle sequences missing key positions
if pd.isna(grn_row["7.50"]):  # Missing critical position
    # Try alternative references
    alt_annotation = gp.annotate_with_multiple_references(sequence)
```

## Summary

GRNBaseProcessor provides:
- Standardized residue numbering for protein families
- Conservation and motif analysis
- Integration with structure and sequence data  
- Family-specific numbering schemes
- Annotation transfer to novel sequences
- Clustering and classification capabilities

The GRN system enables consistent analysis across diverse homologous sequences, facilitating structure-function studies and evolutionary analysis.