# GRN Annotation Approaches in Protos

This document explains three different approaches to annotate sequences with Generic Residue Numbers (GRNs) in the Protos framework.

## Overview

GRN (Generic Residue Numbering) provides a standardized way to number residues across protein families, enabling consistent comparison of homologous positions. In GPCRs, for example, position 7.50 is the most conserved residue in helix 7 across the family.

## Approach 1: CLI Command

The CLI provides a user-friendly command-line interface for GRN annotation:

### Command Structure
```bash
protos grn assign -p <protein_family> [options]
```

### Examples
```bash
# Assign GRNs to all registered microbial opsin sequences
protos grn assign -p microbial_opsins

# Assign GRNs to specific sequences  
protos grn assign -p gpcr_a -s P12345 P67890

# Use custom reference table and output name
protos grn assign -p gpcr_a --reference gpcrdb_ref -o my_gpcr_grns

# Use more cores for parallel processing
protos grn assign -p gpcr_a -n 16
```

### How it Works
1. **Load sequences**: Either specific sequence IDs or all registered sequences
2. **Load reference GRN table**: Based on protein family (e.g., `gpcrdb_ref` for GPCRs)
3. **Find best matches**: Uses MMseqs2 for fast similarity search
4. **Perform alignments**: Aligns query sequences to best-matching references
5. **Transfer GRNs**: Maps GRN positions through the alignment
6. **Save results**: Stores annotated GRN table with entity tracking

### Key Features
- Parallel processing support
- Entity registry integration
- Automatic reference table selection
- MMseqs2 for fast similarity search

## Approach 2: Direct Code (Using GRN Processors)

For programmatic access, you can use the GRN processors directly:

### Basic Example
```python
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
from protos.processing.sequence.seq_processor import SeqProcessor
from protos.processing.grn.grn_assignment import assign_grns_to_sequences
from protos.processing.sequence.seq_alignment import mmseqs2_align2

# Initialize processors
grn_proc = GRNBaseProcessor(name="my_grn_analysis")
seq_proc = SeqProcessor(name="my_sequences")

# Load sequences
sequences = seq_proc.load_sequences("my_sequences.fasta")

# Load GRN reference table  
grn_proc.load_grn_table("ref/gpcrdb_ref")  # For GPCRs

# Get reference sequences
ref_sequences = grn_proc.get_seq_dict()

# Find similar sequences using MMseqs2
hits = mmseqs2_align2(sequences, ref_sequences)

# Filter by similarity and assign GRNs
# (Implementation would follow similar logic to CLI)
```

### Key Components
- **GRNBaseProcessor**: Manages GRN tables and reference data
- **grn_assignment.py**: Contains core GRN assignment functions
- **grn_table_utils.py**: Utilities for GRN table manipulation
- **seq_alignment.py**: Sequence alignment tools (MMseqs2, BioPython)

### Workflow
1. Load query sequences and reference GRN table
2. Perform similarity search (MMseqs2)
3. Align sequences to best matches
4. Transfer GRN annotations through alignment
5. Apply boundary constraints (helix/loop regions)
6. Save annotated table

## Approach 3: Structure-Based GRN Assignment (CifBaseProcessor)

The CifBaseProcessor includes built-in GRN assignment for structure data:

### Method: `assign_grns()`
```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Initialize structure processor
struct_proc = CifBaseProcessor(name="gpcr_structures")

# Load structures
struct_proc.load_dataset("gpcr_dataset")

# Assign GRNs to all chains
grn_assignments = struct_proc.assign_grns(
    protein_family='gpcr_a',           # Protein family
    similarity_threshold=0.2,           # 20% identity threshold
    grn_table_name='gpcrdb_ref',       # Reference table (auto-detected if None)
    use_mmseqs=True,                   # Use fast similarity search
    save_results=True                  # Save assignments
)

# Result: Dictionary of {pdb_chain: grn_series}
# Example: {"1U19_A": pd.Series({"7.50": "K296", "3.50": "D85", ...})}
```

### How it Works
1. **Extract sequences** from loaded structures
2. **Load reference GRN table** (auto-detects based on family)
3. **Find similar sequences** using MMseqs2
4. **Perform pairwise alignments** for similar chains
5. **Transfer GRN annotations** through alignments
6. **Update structure data** with GRN positions
7. **Save results** as GRN annotations

### Additional Methods

#### `get_grn_dict()`
Extracts existing GRN annotations from structure data:
```python
grn_dict = struct_proc.get_grn_dict()
# Returns: {pdb_id: {chain_id: {grn_position: "ResName1L+SeqPos"}}}
# Example: {"1U19": {"A": {"7.50": "K296", "3.50": "D85"}}}
```

#### `apply_grn_interval()`
Filters structure data to specific GRN positions:
```python
# Keep only binding pocket residues
binding_pocket_grns = ["3.32", "3.33", "5.46", "6.48", "7.39", "7.43", "7.50"]
struct_proc.apply_grn_interval(binding_pocket_grns)
```

## Comparison of Approaches

| Feature | CLI Command | Direct Code | Structure-Based |
|---------|------------|-------------|-----------------|
| **Use Case** | Batch processing, automation | Custom pipelines | Structure analysis |
| **Input** | FASTA files or sequence IDs | Any sequence dict | Loaded structures |
| **Flexibility** | Limited to CLI options | Full control | Integrated with structure data |
| **Parallel Processing** | Built-in | Manual setup | Not implemented |
| **Entity Integration** | Automatic | Manual | Automatic |
| **Best For** | Quick analysis, standard workflows | Complex pipelines | Structure-sequence integration |

## Key Limitations Identified

1. **GRN Reference Format**: The `gpcrdb_ref.csv` format may not match expected column names (uses "12.50" instead of "12.500")
2. **Missing Methods**: Some expected methods like `struct_proc.add_grn_annotations()` are not implemented
3. **Ligand Analysis**: No built-in ligand processor for analyzing protein-ligand interactions by GRN position

## Recommendations

1. **For simple batch processing**: Use the CLI command
2. **For custom analysis pipelines**: Use direct code approach
3. **For structure-based analysis**: Use CifBaseProcessor's `assign_grns()` method
4. **For GPCRs specifically**: Ensure the reference table format matches (may need GPCRdb API integration)