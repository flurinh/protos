# Alignments

Protos provides comprehensive alignment capabilities for both sequences and structures, enabling comparative analysis across homologous proteins. The framework integrates multiple alignment algorithms and provides unified interfaces for different alignment types.

## Overview

Alignment types supported:
- **Sequence alignments**: Pairwise and multiple sequence alignments
- **Structure alignments**: 3D coordinate superposition
- **Structure-based sequence alignments**: Using structural information
- **Profile alignments**: HMM-based alignments

## Sequence Alignments

### Pairwise Sequence Alignment

```python
from protos.processing.sequence.seq_processor import SeqProcessor
sp = SeqProcessor()

# Global alignment (Needleman-Wunsch)
alignment = sp.align_sequences(
    seq1="MKTAYIAKQRQISFVKSHFSRQLE",
    seq2="MKTAYIAKQKQISFVKSHFSRQLE",
    method="needle",
    matrix="BLOSUM62",
    gap_open=-10,
    gap_extend=-1
)

# Local alignment (Smith-Waterman)
local_alignment = sp.align_sequences(
    seq1, seq2,
    method="water",
    matrix="PAM250"
)

# Access alignment details
print(f"Identity: {alignment['identity']}%")
print(f"Similarity: {alignment['similarity']}%")
print(f"Gaps: {alignment['gaps']}")
print(f"Score: {alignment['score']}")
```

### Multiple Sequence Alignment (MSA)

```python
# Prepare sequences
sequences = {
    "PROT1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
    "PROT2": "MVLSPADKTNVKGAWGKVGGHAGEYGAEAL",
    "PROT3": "MVLSPADKTNVKAGWGKVGAHAGEYGAEAL",
    "PROT4": "MVLSPADKTNVKASWGKVGGHAGGYGAEAL"
}

# MUSCLE alignment
msa = sp.multiple_sequence_alignment(
    sequences,
    method="muscle",
    params={
        "maxiters": 16,
        "diags": True  # Faster for similar sequences
    }
)

# Clustal Omega alignment
msa_clustal = sp.multiple_sequence_alignment(
    sequences,
    method="clustalo",
    params={
        "iterations": 5,
        "threads": 4
    }
)

# MAFFT alignment (for large datasets)
msa_mafft = sp.multiple_sequence_alignment(
    sequences,
    method="mafft",
    algorithm="auto"  # Automatically choose best algorithm
)
```

### Alignment Analysis

```python
# Calculate conservation
conservation = sp.calculate_conservation(msa)

# Find conserved regions
conserved_regions = sp.find_conserved_regions(
    msa,
    min_length=5,
    min_conservation=0.9
)

# Extract consensus sequence
consensus = sp.get_consensus_sequence(
    msa,
    threshold=0.7,  # 70% agreement required
    ambiguous='X'   # Character for ambiguous positions
)

# Identify variable regions
variable_regions = sp.find_variable_regions(
    msa,
    window_size=10,
    variability_threshold=0.5
)
```

### Profile Alignments

```python
# Build HMM profile from MSA
profile = sp.build_hmm_profile(
    msa,
    name="my_family"
)

# Align new sequence to profile
profile_alignment = sp.align_to_profile(
    sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFL",
    profile=profile,
    e_value=0.001
)

# Search database with profile
hits = sp.search_with_profile(
    profile=profile,
    database=sp.list_entities(),
    e_value_threshold=1e-5
)
```

## Structure Alignments

### Basic Structure Alignment

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
cp = CifBaseProcessor()

# Align two structures
alignment = cp.align_structures(
    structure1="1ubq",
    structure2="2gb1",
    method="kabsch",  # Kabsch algorithm
    atoms="CA"        # Alpha carbons only
)

# Get transformation matrix
rotation = alignment['rotation']
translation = alignment['translation']
rmsd = alignment['rmsd']

# Apply transformation
aligned_coords = cp.apply_transformation(
    structure="2gb1",
    rotation=rotation,
    translation=translation
)
```

### Multiple Structure Alignment

```python
# Align multiple structures
structures = ["1ubq", "2gb1", "1crn", "3nir"]

# Using FoldMason
multi_alignment = cp.align_multiple_structures(
    structures,
    method="foldmason",
    reference=0,  # Use first structure as reference
    output_name="aligned_structures"
)

# Using iterative alignment
iterative_alignment = cp.iterative_structure_alignment(
    structures,
    method="progressive",
    guide_tree="upgma"
)

# Extract core regions
core_residues = cp.find_structural_core(
    multi_alignment,
    coverage=0.8  # Present in 80% of structures
)
```

### Advanced Structure Alignment

```python
# Flexible alignment (allows conformational differences)
flexible_alignment = cp.flexible_align(
    structure1="open_form",
    structure2="closed_form",
    n_fragments=5,  # Number of rigid fragments
    flexibility_threshold=5.0  # Å
)

# Domain-based alignment
domain_alignment = cp.align_domains(
    structure1="1abc",
    structure2="2def",
    domain_definitions={
        "domain1": (1, 150),
        "domain2": (151, 300)
    }
)

# Align with sequence constraints
constrained_alignment = cp.align_with_constraints(
    structure1="1ubq",
    structure2="2gb1",
    sequence_alignment=msa,  # Pre-computed sequence alignment
    gap_penalty=2.0
)
```

## Structure-Based Sequence Alignment

### Extracting Sequence Alignments from Structures

```python
# Get structure-based sequence alignment
struct_seq_alignment = cp.get_structure_sequence_alignment(
    structures=["1ubq", "2gb1", "1crn"],
    method="3dcoffee"
)

# Improve MSA with structural information
improved_msa = sp.refine_alignment_with_structure(
    msa=msa,
    structures={
        "PROT1": "1ubq",
        "PROT2": "2gb1"
    }
)

# Identify structurally equivalent positions
equivalent_positions = cp.find_equivalent_positions(
    structure1="1ubq",
    structure2="2gb1",
    distance_threshold=3.0  # Å
)
```

## Integration with GRN

### GRN-Guided Alignments

```python
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
gp = GRNBaseProcessor()

# Align using GRN positions as anchors
grn_alignment = gp.align_with_grn_anchors(
    sequences=sequences,
    grn_positions=["1.50", "2.50", "3.50", "7.50"],
    anchor_weight=2.0  # Increase weight for GRN positions
)

# Transfer GRN annotation through alignment
grn_transfer = gp.transfer_grn_via_alignment(
    source_sequence="BACR_HALSA",
    target_sequence="NEW_OPSIN",
    source_grn=grn_table
)

# Validate GRN consistency in alignment
validation = gp.validate_grn_alignment(
    alignment=msa,
    grn_table=grn_table,
    required_positions=["3.50", "7.50"]  # Critical positions
)
```

## Alignment Utilities

### Format Conversion

```python
# Convert between alignment formats
# Supported: clustal, fasta, stockholm, phylip, nexus

# Read alignment
alignment = sp.read_alignment("alignment.aln", format="clustal")

# Convert to different format
sp.write_alignment(
    alignment,
    "alignment.sto",
    format="stockholm"
)

# Convert to DataFrame
alignment_df = sp.alignment_to_dataframe(alignment)
```

### Alignment Editing

```python
# Trim alignment
trimmed = sp.trim_alignment(
    alignment,
    start=10,
    end=250,
    remove_gaps=True
)

# Remove gappy columns
clean_alignment = sp.remove_gappy_columns(
    alignment,
    gap_threshold=0.5  # Remove columns with >50% gaps
)

# Remove similar sequences
nr_alignment = sp.remove_redundant_sequences(
    alignment,
    identity_threshold=0.95
)

# Extract sub-alignment
sub_alignment = sp.extract_sub_alignment(
    alignment,
    sequences=["PROT1", "PROT3", "PROT5"]
)
```

### Alignment Visualization

```python
# Generate alignment visualization
sp.visualize_alignment(
    alignment,
    output="alignment.png",
    color_scheme="clustal",  # or "hydrophobicity", "charge"
    show_consensus=True,
    show_conservation=True
)

# Create alignment logo
sp.create_sequence_logo(
    alignment,
    output="logo.png",
    positions=range(50, 100)  # Focus on specific region
)

# Export for visualization tools
sp.export_for_jalview(alignment, "alignment.jvp")
sp.export_for_pymol(alignment, structures, "alignment.pse")
```

## Best Practices

### 1. Algorithm Selection

Choose appropriate algorithms based on sequence similarity:

```python
def select_alignment_method(sequences):
    """Select best alignment method based on sequence properties."""
    
    # Calculate average pairwise identity
    identities = []
    seq_list = list(sequences.values())
    
    for i in range(len(seq_list)):
        for j in range(i+1, len(seq_list)):
            alignment = sp.align_sequences(seq_list[i], seq_list[j])
            identities.append(alignment['identity'])
    
    avg_identity = np.mean(identities)
    
    if avg_identity > 70:
        return "muscle"  # Fast for similar sequences
    elif avg_identity > 40:
        return "clustalo"  # Good balance
    else:
        return "mafft"  # Better for divergent sequences
```

### 2. Parameter Optimization

```python
# Test different parameters
gap_penalties = [(-10, -1), (-12, -2), (-8, -1)]
matrices = ["BLOSUM62", "BLOSUM45", "PAM250"]

best_score = -float('inf')
best_params = None

for gap_open, gap_extend in gap_penalties:
    for matrix in matrices:
        alignment = sp.align_sequences(
            seq1, seq2,
            matrix=matrix,
            gap_open=gap_open,
            gap_extend=gap_extend
        )
        
        if alignment['score'] > best_score:
            best_score = alignment['score']
            best_params = (matrix, gap_open, gap_extend)
```

### 3. Quality Assessment

```python
def assess_alignment_quality(alignment):
    """Evaluate alignment quality metrics."""
    
    quality_metrics = {
        'conservation': sp.calculate_conservation(alignment),
        'gap_fraction': sp.calculate_gap_fraction(alignment),
        'column_occupancy': sp.calculate_column_occupancy(alignment),
        'sequence_identity': sp.calculate_average_identity(alignment)
    }
    
    # Check for problematic features
    if quality_metrics['gap_fraction'] > 0.3:
        print("Warning: High gap content")
    
    if quality_metrics['sequence_identity'] < 0.2:
        print("Warning: Low sequence identity")
    
    return quality_metrics
```

### 4. Large-Scale Alignments

```python
# For very large alignments
def align_large_dataset(sequences, chunk_size=100):
    """Align large sequence sets efficiently."""
    
    # First, cluster sequences
    clusters = sp.cluster_sequences(
        sequences,
        identity_threshold=0.9
    )
    
    # Align representatives
    representatives = {}
    for cluster_id, members in clusters.items():
        representatives[f"cluster_{cluster_id}"] = sequences[members[0]]
    
    repr_alignment = sp.multiple_sequence_alignment(
        representatives,
        method="mafft",
        algorithm="linsi"  # Accurate mode
    )
    
    # Add remaining sequences
    for cluster_id, members in clusters.items():
        for member in members[1:]:
            repr_alignment = sp.add_to_alignment(
                repr_alignment,
                sequence_name=member,
                sequence=sequences[member]
            )
    
    return repr_alignment
```

## Common Workflows

### Workflow 1: Homology Modeling Preparation

```python
def prepare_for_homology_modeling(target_seq, template_pdb):
    """Prepare alignment for homology modeling."""
    
    # Load template structure
    template_structure = cp.load_structure(template_pdb)
    template_seq = cp.extract_sequence(template_structure)
    
    # Perform careful alignment
    alignment = sp.align_sequences(
        target_seq,
        template_seq,
        method="needle",
        matrix="BLOSUM62",
        gap_open=-10,
        gap_extend=-1
    )
    
    # Identify gaps in alignment
    gaps = sp.find_alignment_gaps(alignment)
    
    # Check structural context of gaps
    gap_regions = []
    for gap in gaps:
        region_structure = cp.extract_region(
            template_structure,
            start=gap['start']-5,
            end=gap['end']+5
        )
        
        secondary_structure = cp.assign_secondary_structure(region_structure)
        gap_regions.append({
            'gap': gap,
            'secondary_structure': secondary_structure,
            'loop_region': 'C' in secondary_structure  # Coil/loop
        })
    
    return {
        'alignment': alignment,
        'gaps': gap_regions,
        'identity': alignment['identity'],
        'coverage': alignment['coverage']
    }
```

### Workflow 2: Family Conservation Analysis

```python
def analyze_family_conservation(family_sequences, structures=None):
    """Comprehensive conservation analysis for protein family."""
    
    # Create MSA
    msa = sp.multiple_sequence_alignment(
        family_sequences,
        method="mafft"
    )
    
    # Calculate conservation
    conservation = sp.calculate_conservation(msa)
    
    # Identify conserved motifs
    motifs = sp.find_conserved_motifs(
        msa,
        min_length=3,
        min_conservation=0.8
    )
    
    # If structures available, add structural context
    if structures:
        structural_conservation = {}
        
        for position in range(len(conservation)):
            struct_contexts = []
            
            for seq_name, structure in structures.items():
                if seq_name in msa:
                    context = cp.get_structural_context(
                        structure,
                        residue_position=position,
                        radius=8.0
                    )
                    struct_contexts.append(context)
            
            structural_conservation[position] = {
                'sequence_conservation': conservation[position],
                'avg_burial': np.mean([c['burial'] for c in struct_contexts]),
                'secondary_structure': most_common([c['ss'] for c in struct_contexts])
            }
    
    return {
        'alignment': msa,
        'conservation': conservation,
        'motifs': motifs,
        'structural_conservation': structural_conservation if structures else None
    }
```

### Workflow 3: Structure Superposition Pipeline

```python
def structure_superposition_pipeline(structure_list):
    """Complete structure alignment and analysis pipeline."""
    
    # Initial pairwise alignments
    pairwise_rmsds = {}
    for i, struct1 in enumerate(structure_list):
        for j, struct2 in enumerate(structure_list[i+1:], i+1):
            alignment = cp.align_structures(struct1, struct2)
            pairwise_rmsds[(struct1, struct2)] = alignment['rmsd']
    
    # Select reference (most similar to all others)
    avg_rmsds = {}
    for struct in structure_list:
        rmsds = [rmsd for (s1, s2), rmsd in pairwise_rmsds.items() 
                 if struct in (s1, s2)]
        avg_rmsds[struct] = np.mean(rmsds)
    
    reference = min(avg_rmsds, key=avg_rmsds.get)
    
    # Align all to reference
    aligned_structures = {}
    transformations = {}
    
    for struct in structure_list:
        if struct != reference:
            alignment = cp.align_structures(
                reference,
                struct,
                method="kabsch",
                atoms="CA"
            )
            
            aligned_structures[struct] = cp.apply_transformation(
                struct,
                alignment['rotation'],
                alignment['translation']
            )
            transformations[struct] = alignment
        else:
            aligned_structures[struct] = cp.load_structure(struct)
    
    # Calculate structural divergence
    divergence = cp.calculate_structural_divergence(aligned_structures)
    
    # Identify flexible regions
    flexible_regions = cp.find_flexible_regions(
        aligned_structures,
        rmsd_threshold=3.0
    )
    
    return {
        'reference': reference,
        'aligned_structures': aligned_structures,
        'transformations': transformations,
        'divergence': divergence,
        'flexible_regions': flexible_regions
    }
```

## Troubleshooting

### Common Issues

1. **Poor alignment quality**
```python
# Try different scoring matrices
matrices = ["BLOSUM45", "BLOSUM62", "BLOSUM80", "PAM250"]
best_alignment = None
best_score = -float('inf')

for matrix in matrices:
    alignment = sp.align_sequences(seq1, seq2, matrix=matrix)
    if alignment['score'] > best_score:
        best_alignment = alignment
        best_score = alignment['score']
```

2. **Memory issues with large alignments**
```python
# Use iterative approach
sp.align_sequences_iterative(
    sequences,
    chunk_size=50,
    method="mafft",
    memory_efficient=True
)
```

3. **Structural alignment failures**
```python
# Try different atom selections
atom_selections = ["CA", "backbone", "CB"]
for atoms in atom_selections:
    try:
        alignment = cp.align_structures(s1, s2, atoms=atoms)
        break
    except:
        continue
```

## Summary

Protos alignment capabilities include:
- Multiple sequence alignment algorithms
- Structure superposition methods
- Structure-based sequence alignment
- Profile and HMM alignments
- Integration with GRN system
- Comprehensive analysis tools

The alignment functionality enables comparative analysis across sequences and structures, supporting a wide range of structural biology applications.