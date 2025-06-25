# GRN-Structure Integration Guide

## Important Updates

### Dot Notation is Now Standard
The GRN system now uses **dot notation** (e.g., `1.50`, `7.50`) as the standard format. The older 'x' notation (e.g., `1x50`, `7x50`) is deprecated and will be removed in future versions.

### Recent Changes
- Fixed GRN assignment to use correct reference table paths (`ref/` subdirectory)
- Added column name compatibility for different CIF parsers (`auth_comp_id` vs `res_name3l`)
- Fixed sequence extraction to properly handle unique residues
- Improved path resolution for cross-platform compatibility

## Overview

This guide explains how to use the GRN (Generic Residue Numbering) integration features in CifBaseProcessor to:
1. Extract sequences from protein structures
2. Assign GRN numbers to sequences
3. Map GRN annotations back to structure data

## Understanding the Data Format

### Structure Data Format

The CifBaseProcessor works with structure data in a DataFrame format with the following key columns:

```python
# Essential columns for sequence extraction
'pdb_id'         # PDB identifier (e.g., '1UAZ')
'auth_chain_id'  # Chain identifier (e.g., 'A', 'B')
'auth_seq_id'    # Residue sequence number (1, 2, 3...)
'auth_comp_id'   # 3-letter amino acid code (e.g., 'MET', 'ALA')

# Additional structure columns
'atom_name'      # Atom name (e.g., 'CA', 'N', 'C')
'x', 'y', 'z'    # 3D coordinates
'atom_id'        # Unique atom identifier

# GRN column (added after assignment)
'grn'           # GRN position (e.g., '7.50', '3.50')
```

### GRN Table Format

GRN reference tables have:
- **Rows**: Protein identifiers (e.g., 'BR', '1UAZ')
- **Columns**: GRN positions (e.g., '1.50', '2.50', '3.50', '7.50')
- **Values**: Residue+position (e.g., 'K296', 'D85')

Example:
```
       1.50  2.50  3.50  ...  7.50
BR     R82   D85   D96   ...  K216
1UAZ   R82   D85   D96   ...  K216
```

## Step-by-Step Workflow

### 1. Load Structure Data

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Initialize processor
processor = CifBaseProcessor(
    name="my_analysis",
    data_root="/path/to/data",
    processor_data_dir="structure"
)

# Load structures
processor.load_dataset("microbial_opsins")
# OR load individual structure
processor.load_structure("1UAZ")
```

### 2. Extract Sequences

```python
# Extract sequences from structure data
sequences = processor.get_seq_dict()
# Returns: {'1UAZ_A': 'MLELLPTAVEGVSQ...', '1UAZ_B': '...'}

# Save sequences (optional)
from protos.io.fasta_utils import write_fasta
write_fasta(sequences, "extracted_sequences.fasta")
```

### 3. Assign GRNs

```python
# Assign GRNs based on sequence similarity
grn_assignments = processor.assign_grns(
    protein_family='microbial_opsins',  # or 'gpcr_a'
    similarity_threshold=0.2,            # 20% identity minimum
    use_mmseqs=True,                     # Use fast MMseqs2 search
    save_results=True                    # Save GRN table
)

# This automatically:
# 1. Finds similar sequences in reference GRN table
# 2. Aligns sequences and transfers GRN numbers
# 3. Adds 'grn' column to structure data
# 4. Saves GRN assignments to file
```

### 4. Extract GRN Annotations

```python
# Get GRN annotations from structure
grn_dict = processor.get_grn_dict()

# Returns nested dictionary:
# {
#   '1UAZ': {
#     'A': {
#       '1.50': 'R82',
#       '3.50': 'D96',
#       '7.50': 'K216'
#     }
#   }
# }

# Access specific GRN position
k216_grn = grn_dict['1UAZ']['A']['7.50']  # Returns 'K216'
```

### 5. Query Structure by GRN

```python
# Find all residues at position 7.50 (Schiff base lysine)
schiff_base = processor.data[processor.data['grn'] == '7.50']

# Get coordinates of specific GRN position
k750_atoms = processor.data[
    (processor.data['grn'] == '7.50') & 
    (processor.data['atom_name'] == 'CA')
]

# Calculate distances from retinal binding lysine
for idx, atom in k750_atoms.iterrows():
    coords = np.array([atom['x'], atom['y'], atom['z']])
    print(f"{atom['pdb_id']}:{atom['auth_chain_id']} K{atom['auth_seq_id']} at {coords}")
```

## Key GRN Positions

For microbial opsins:
- **1.50**: Conserved position in TM1
- **2.50**: Conserved position in TM2
- **3.50**: D85 in bacteriorhodopsin (proton acceptor)
- **6.48**: Conserved tryptophan
- **7.50**: K216 in bacteriorhodopsin (Schiff base lysine)

## Complete Example

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
import numpy as np

# Initialize
processor = CifBaseProcessor(name="mo_analysis")

# Load bacteriorhodopsin
processor.load_structure("1UAZ")

# Extract sequences
sequences = processor.get_seq_dict()
print(f"Found {len(sequences)} chains")

# Assign GRNs
grn_assignments = processor.assign_grns(
    protein_family='microbial_opsins',
    similarity_threshold=0.2
)

# Get GRN annotations
grn_dict = processor.get_grn_dict()

# Find key functional residues
key_residues = processor.data[
    processor.data['grn'].isin(['3.50', '7.50'])
]

# Analyze
for _, res in key_residues.iterrows():
    print(f"{res['grn']}: {res['auth_comp_id']}{res['auth_seq_id']}")
```

## Troubleshooting

### No sequences extracted
- Check that `auth_comp_id` contains valid 3-letter amino acid codes
- Verify `auth_seq_id` is numeric and sequential
- Ensure data has both `pdb_id` and `auth_chain_id`

### No GRN assignments
- Lower similarity threshold (try 0.1 instead of 0.2)
- Check that reference GRN table exists for protein family
- Try with `use_mmseqs=False` to use BioPython instead
- Verify sequences are from the expected protein family

### GRN mapping issues
- Ensure residue numbers in GRN table match structure
- Check that 3-letter codes match between structure and GRN assignment
- Verify chain IDs are consistent

## Data Format Examples

### Minimal structure data for testing:
```python
import pandas as pd

test_data = pd.DataFrame({
    'pdb_id': ['TEST1'] * 5,
    'auth_chain_id': ['A'] * 5,
    'auth_seq_id': [1, 2, 3, 4, 5],
    'auth_comp_id': ['MET', 'ALA', 'GLY', 'VAL', 'LEU']
})
```

### After GRN assignment:
```python
# Same data with added 'grn' column
test_data['grn'] = [None, '1.50', None, None, '2.50']
```

### GRN dictionary format:
```python
{
    'TEST1': {
        'A': {
            '1.50': 'A2',
            '2.50': 'L5'
        }
    }
}
```