# SeqProcessor

The SeqProcessor manages protein and nucleotide sequence data, providing comprehensive functionality for sequence loading, manipulation, alignment, and analysis.

## Overview

SeqProcessor handles:
- FASTA file parsing and writing
- Sequence validation and cleaning
- Multiple sequence alignments
- Sequence search and clustering
- Integration with sequence databases
- Cross-format operations with other processors

## Basic Usage

### Initialization

```python
from protos.processing.sequence.seq_processor import SeqProcessor

# Create processor instance
sp = SeqProcessor(name="sequence_analysis")

# With options
sp = SeqProcessor(
    name="my_sequences",
    verbose=True,           # Detailed logging
    validate_on_load=True   # Automatic sequence validation
)
```

### Loading Sequences

```python
# Load single sequence by name
sequence = sp.load_sequence("P12345")

# Load from FASTA file
sequences = sp.load_fasta("proteins.fasta")

# Load multiple sequences
seq_dict = sp.load_sequences(["P12345", "Q9Y6K1", "O14786"])

# Check if sequence exists
if sp.has_entity("MY_PROTEIN"):
    sequence = sp.load_sequence("MY_PROTEIN")
```

### Saving Sequences

```python
# Save single sequence
sp.save_sequence("NEW_PROTEIN", "MKTAYIAKQRQISFVKSHFSRQLE...")

# Save multiple sequences
sequences = {
    "PROT1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
    "PROT2": "MVLSPADKTNVKGAWGKVGGHAGEYGAEAL",
    "PROT3": "MVLSPADKTNVKAGWGKVGAHAGEYGAEAL"
}
sp.save_sequences(sequences, "my_proteins.fasta")

# Save with metadata in header
sp.save_sequence(
    "EGFR_HUMAN",
    sequence="MRPSGTAGAALLALLAAL...",
    description="Epidermal growth factor receptor OS=Homo sapiens GN=EGFR"
)
```

## Sequence Operations

### Basic Manipulations

```python
# Get sequence length
length = sp.get_sequence_length("P12345")

# Extract subsequence
subseq = sp.extract_subsequence("P12345", start=100, end=200)

# Reverse sequence
reversed_seq = sp.reverse_sequence("P12345")

# Complement (for nucleotides)
complement = sp.complement_sequence("NT_SEQUENCE", type="dna")

# Translate (DNA/RNA to protein)
protein = sp.translate_sequence("mRNA_SEQ")
```

### Sequence Validation

```python
# Validate protein sequence
is_valid = sp.validate_sequence(
    sequence="MKTAYIAKQRQISFVKSHFSRQLE",
    seq_type="protein"
)

# Clean sequence (remove invalid characters)
clean_seq = sp.clean_sequence(
    "MKT-AYI*AKQ RQIS",
    remove_gaps=True,
    remove_stops=True
)

# Check for ambiguous residues
has_ambiguous = sp.has_ambiguous_residues(sequence)

# Find stop codons (in protein sequences)
stop_positions = sp.find_stop_codons(sequence)
```

### Sequence Properties

```python
# Calculate composition
composition = sp.calculate_composition("P12345")
# {'A': 45, 'C': 12, 'D': 23, ...}

# Calculate molecular weight
mw = sp.calculate_molecular_weight("P12345")
# 28453.2

# Calculate isoelectric point
pi = sp.calculate_isoelectric_point("P12345")
# 6.82

# Get hydrophobicity profile
hydro = sp.calculate_hydrophobicity("P12345", window=7)

# Calculate sequence complexity
complexity = sp.calculate_complexity("P12345")
```

## Alignment Operations

### Pairwise Alignment

```python
# Align two sequences
alignment = sp.align_sequences(
    seq1="MKTAYIAKQRQISFVKSHFSRQLE",
    seq2="MKTAYIAKQKQISFVKSHFSRQLE",
    method="needle"  # Global alignment
)

# Local alignment
local_align = sp.align_sequences(seq1, seq2, method="water")

# Get alignment scores
print(f"Identity: {alignment['identity']}%")
print(f"Similarity: {alignment['similarity']}%")
print(f"Score: {alignment['score']}")
```

### Multiple Sequence Alignment

```python
# Align multiple sequences
sequences = {
    "SEQ1": "MVLSPADKTNVKAAWGKVGAHAGEYGAEAL",
    "SEQ2": "MVLSPADKTNVKGAWGKVGGHAGEYGAEAL",
    "SEQ3": "MVLSPADKTNVKAGWGKVGAHAGEYGAEAL"
}

# Using MUSCLE
msa = sp.multiple_sequence_alignment(
    sequences,
    method="muscle",
    output="alignment.aln"
)

# Using Clustal Omega
msa = sp.multiple_sequence_alignment(
    sequences,
    method="clustalo",
    params={"iterations": 3}
)

# Extract consensus sequence
consensus = sp.get_consensus_sequence(msa, threshold=0.7)
```

### Alignment Analysis

```python
# Calculate conservation scores
conservation = sp.calculate_conservation(msa)

# Find conserved regions
conserved_regions = sp.find_conserved_regions(
    msa,
    min_length=5,
    min_conservation=0.9
)

# Get gap statistics
gap_stats = sp.analyze_gaps(msa)

# Trim alignment
trimmed = sp.trim_alignment(
    msa,
    gap_threshold=0.5,  # Remove columns with >50% gaps
    terminal_gaps=True   # Remove terminal gaps
)
```

## Database Integration

### UniProt Integration

```python
# Download from UniProt
sp.download_from_uniprot(["P12345", "Q9Y6K1", "O14786"])

# Search UniProt
results = sp.search_uniprot(
    query="kinase AND organism:human",
    limit=100
)

# Get UniProt metadata
metadata = sp.get_uniprot_metadata("P12345")
```

### NCBI Integration

```python
# Download from NCBI
sp.download_from_ncbi(["NP_001005484", "NP_005219"])

# BLAST search
blast_results = sp.blast_sequence(
    sequence="MKTAYIAKQRQISFVKSHFSRQLE",
    database="nr",
    e_value=0.001,
    max_hits=100
)
```

## Sequence Analysis

### Motif and Domain Analysis

```python
# Find sequence motifs
motifs = sp.find_motifs("P12345", database="prosite")

# Find domains
domains = sp.find_domains("P12345", database="pfam")

# Custom pattern search
pattern = "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H"  # Zinc finger
matches = sp.find_pattern(sequence, pattern)

# Find repeats
repeats = sp.find_repeats(sequence, min_length=3, min_copies=2)
```

### Clustering and Similarity

```python
# Calculate sequence similarity matrix
sequences = sp.load_sequences(["PROT1", "PROT2", "PROT3"])
similarity_matrix = sp.calculate_similarity_matrix(sequences)

# Cluster sequences
clusters = sp.cluster_sequences(
    sequences,
    method="hierarchical",
    threshold=0.7
)

# Find similar sequences in database
similar = sp.find_similar_sequences(
    query_sequence,
    database=sp.list_entities(),
    identity_threshold=0.3
)

# Build sequence profile
profile = sp.build_sequence_profile(aligned_sequences)
```

### Evolutionary Analysis

```python
# Calculate evolutionary distances
distances = sp.calculate_evolutionary_distances(
    sequences,
    model="JTT"  # Jones-Taylor-Thornton
)

# Build phylogenetic tree
tree = sp.build_phylogenetic_tree(
    sequences,
    method="neighbor_joining"
)

# Calculate dN/dS ratio (for coding sequences)
dnds = sp.calculate_dnds(seq1, seq2)

# Identify orthologs
orthologs = sp.find_orthologs(
    query="EGFR_HUMAN",
    species=["mouse", "rat", "zebrafish"]
)
```

## Integration Features

### Structure Integration

```python
# Extract sequences from structures
from protos.processing.structure.struct_base_processor import CifBaseProcessor
cp = CifBaseProcessor()

# Get sequences from all loaded structures
structure_sequences = cp.get_seq_dict()

# Save to sequence processor
sp.save_sequences(structure_sequences, "structure_derived.fasta")

# Compare structure vs UniProt sequence
struct_seq = cp.extract_sequence("1atp")
uniprot_seq = sp.load_sequence("P00568")
alignment = sp.align_sequences(struct_seq, uniprot_seq)
```

### GRN Integration

```python
# Create sequences from GRN table
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
gp = GRNBaseProcessor()

# Extract sequences from GRN
grn_sequences = gp.get_seq_dict()
sp.save_sequences(grn_sequences, "grn_sequences.fasta")

# Align sequences for GRN annotation
aligned = sp.multiple_sequence_alignment(grn_sequences)
```

### Embedding Preparation

```python
# Prepare sequences for embedding
sequences = sp.load_sequences(["PROT1", "PROT2", "PROT3"])

# Filter by length for embedding models
filtered = sp.filter_sequences(
    sequences,
    min_length=50,
    max_length=1000
)

# Remove redundancy
nr_sequences = sp.remove_redundancy(
    filtered,
    identity_threshold=0.9
)

# Pass to embedding processor
from protos.processing.embedding.embedding_processor import EmbeddingProcessor
ep = EmbeddingProcessor()
embeddings = ep.embed_sequences(nr_sequences)
```

## Advanced Features

### Batch Processing

```python
# Process sequences in batches
def batch_analysis(sequence_list, batch_size=100):
    results = []
    
    for i in range(0, len(sequence_list), batch_size):
        batch = sequence_list[i:i+batch_size]
        
        # Load batch
        sequences = sp.load_sequences(batch)
        
        # Analyze
        for seq_id, sequence in sequences.items():
            result = {
                "id": seq_id,
                "length": len(sequence),
                "mw": sp.calculate_molecular_weight(sequence),
                "pi": sp.calculate_isoelectric_point(sequence)
            }
            results.append(result)
    
    return pd.DataFrame(results)
```

### Custom File Formats

```python
# Read custom format
def read_custom_format(file_path):
    sequences = {}
    with open(file_path) as f:
        # Custom parsing logic
        pass
    
    # Save in standard format
    sp.save_sequences(sequences, "converted.fasta")
    return sequences

# Write custom format
def write_custom_format(sequences, output_path):
    with open(output_path, 'w') as f:
        for seq_id, sequence in sequences.items():
            # Custom writing logic
            f.write(f">{seq_id}\n{sequence}\n")
```

### Sequence Modification

```python
# Add tags
tagged_seq = sp.add_tag(
    sequence="MKTAYIAKQRQISFVKSHFSRQLE",
    tag="HHHHHH",  # His-tag
    position="N"   # N-terminus
)

# Remove signal peptides
mature_seq = sp.remove_signal_peptide("P12345")

# Add linkers
linked = sp.add_linker(
    seq1="MKTAYIAK",
    seq2="QRQISFVK",
    linker="GGGGS"
)

# Mutate sequence
mutated = sp.mutate_sequence(
    sequence="MKTAYIAKQRQISFVK",
    mutations={"K3A": 3, "Q9E": 9}  # K3A, Q9E mutations
)
```

## Best Practices

### 1. Sequence Validation

Always validate sequences before analysis:

```python
# Validate on load
sequences = sp.load_sequences(seq_list, validate=True)

# Manual validation
for seq_id, sequence in sequences.items():
    if not sp.validate_sequence(sequence):
        sp.logger.warning(f"Invalid sequence: {seq_id}")
        # Clean or skip
        sequences[seq_id] = sp.clean_sequence(sequence)
```

### 2. Memory Management

For large sequence datasets:

```python
# Use generators for large files
def process_large_fasta(file_path):
    for seq_id, sequence in sp.read_fasta_iterator(file_path):
        # Process one sequence at a time
        result = analyze_sequence(sequence)
        yield seq_id, result

# Chunk processing
chunk_size = 1000
for chunk in sp.read_fasta_chunks(file_path, chunk_size):
    process_chunk(chunk)
```

### 3. Alignment Parameters

Choose appropriate alignment parameters:

```python
# For similar sequences
alignment = sp.align_sequences(
    seq1, seq2,
    gap_open=-10,
    gap_extend=-1
)

# For divergent sequences
alignment = sp.align_sequences(
    seq1, seq2,
    matrix="BLOSUM45",  # More permissive
    gap_open=-5,
    gap_extend=-0.5
)
```

### 4. Error Handling

```python
# Handle missing sequences gracefully
missing_sequences = []
for seq_id in sequence_list:
    try:
        seq = sp.load_sequence(seq_id)
        if seq:
            process_sequence(seq)
    except Exception as e:
        sp.logger.error(f"Failed to load {seq_id}: {e}")
        missing_sequences.append(seq_id)
```

## Common Workflows

### Workflow 1: Homolog Analysis

```python
def analyze_homologs(query_id, identity_threshold=0.3):
    """Find and analyze sequence homologs."""
    
    # Load query
    query_seq = sp.load_sequence(query_id)
    
    # Find homologs
    homologs = sp.find_similar_sequences(
        query_seq,
        identity_threshold=identity_threshold
    )
    
    # Align all homologs
    sequences = {query_id: query_seq}
    for homolog in homologs:
        sequences[homolog['id']] = sp.load_sequence(homolog['id'])
    
    msa = sp.multiple_sequence_alignment(sequences)
    
    # Analyze conservation
    conservation = sp.calculate_conservation(msa)
    
    # Find conserved motifs
    motifs = sp.find_conserved_regions(msa)
    
    return {
        "homologs": homologs,
        "alignment": msa,
        "conservation": conservation,
        "motifs": motifs
    }
```

### Workflow 2: Mutation Analysis

```python
def analyze_mutations(wild_type_id, mutant_sequences):
    """Analyze effects of mutations."""
    
    # Load wild-type
    wt_seq = sp.load_sequence(wild_type_id)
    
    results = []
    for mutant_id, mut_seq in mutant_sequences.items():
        # Align to identify mutations
        alignment = sp.align_sequences(wt_seq, mut_seq)
        
        # Find mutation positions
        mutations = sp.identify_mutations(alignment)
        
        # Analyze changes
        for mut in mutations:
            result = {
                "mutant": mutant_id,
                "position": mut['position'],
                "wild_type": mut['from'],
                "mutant": mut['to'],
                "property_change": sp.analyze_property_change(
                    mut['from'], mut['to']
                )
            }
            results.append(result)
    
    return pd.DataFrame(results)
```

### Workflow 3: Domain Architecture

```python
def analyze_domain_architecture(sequence_ids):
    """Analyze domain composition of proteins."""
    
    architectures = {}
    
    for seq_id in sequence_ids:
        sequence = sp.load_sequence(seq_id)
        
        # Find domains
        domains = sp.find_domains(sequence)
        
        # Sort by position
        domains.sort(key=lambda x: x['start'])
        
        # Create architecture string
        architecture = "-".join([d['name'] for d in domains])
        
        architectures[seq_id] = {
            "domains": domains,
            "architecture": architecture,
            "domain_count": len(domains)
        }
    
    # Group by architecture
    architecture_groups = {}
    for seq_id, info in architectures.items():
        arch = info['architecture']
        if arch not in architecture_groups:
            architecture_groups[arch] = []
        architecture_groups[arch].append(seq_id)
    
    return architectures, architecture_groups
```

## Troubleshooting

### Common Issues

1. **Invalid sequence characters**
```python
# Clean sequences automatically
sequence = sp.clean_sequence(raw_sequence, strict=True)
```

2. **Memory issues with large alignments**
```python
# Use iterative alignment for large sets
sp.align_sequences_iterative(
    sequences,
    chunk_size=50,
    output="large_alignment.aln"
)
```

3. **Slow database searches**
```python
# Use local database copies
sp.setup_local_blast_db("nr", update=True)
results = sp.blast_sequence(query, database="local_nr")
```

## Summary

SeqProcessor provides:
- Comprehensive sequence data management
- Multiple alignment algorithms
- Database integration
- Sequence analysis tools
- Cross-format operations
- Batch processing capabilities

The processor handles all aspects of sequence data from basic storage to complex evolutionary analyses, integrating seamlessly with other Protos components.