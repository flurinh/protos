# CifBaseProcessor

The CifBaseProcessor manages 3D protein structure data, providing comprehensive functionality for loading, analyzing, and manipulating structural information from PDB and mmCIF files.

## Overview

CifBaseProcessor handles:
- Loading structures from PDB/mmCIF files
- Parsing atomic coordinates and metadata
- Chain and residue management
- Structure-based sequence extraction
- Integration with structure alignment tools
- Cross-format operations with other processors

## Basic Usage

### Initialization

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Create processor instance
cp = CifBaseProcessor(name="structure_analysis")

# With options
cp = CifBaseProcessor(
    name="my_structures",
    verbose=True,              # Enable detailed logging
    apply_standard_names=True  # Standardize atom/residue names
)
```

### Loading Structures

```python
# Load by PDB ID
structure = cp.load_structure("1ubq")

# Load with type enforcement
structure = cp.load_structure("1ubq", apply_dtypes=True)

# Load specific chains
structure_chain_a = cp.load_structure("1ubq", chain="A")

# Check if structure exists
if cp.has_entity("1ubq"):
    structure = cp.load_structure("1ubq")
```

## Data Format

Structures are loaded as pandas DataFrames with standardized columns:

```python
# Structure DataFrame schema
structure_df = pd.DataFrame({
    'record_type': str,      # ATOM, HETATM
    'atom_serial': int,      # Atom serial number
    'atom_name': str,        # Atom name (CA, CB, etc.)
    'alt_loc': str,          # Alternate location indicator
    'residue_name': str,     # Residue name (ALA, GLY, etc.)
    'chain_id': str,         # Chain identifier
    'residue_number': int,   # Residue sequence number
    'insertion_code': str,   # Insertion code
    'x': float,              # X coordinate
    'y': float,              # Y coordinate  
    'z': float,              # Z coordinate
    'occupancy': float,      # Occupancy factor
    'b_factor': float,       # Temperature factor
    'element': str,          # Element symbol
    'charge': str,           # Formal charge
    
    # Additional mmCIF fields
    'auth_chain_id': str,    # Author chain ID
    'auth_seq_id': int,      # Author sequence number
    'auth_comp_id': str,     # Author compound ID
    'pdb_id': str            # PDB identifier
})
```

## Core Operations

### Structure Analysis

```python
# Get structure information
chains = cp.get_chains("1ubq")
# ['A']

residue_count = cp.get_residue_count("1ubq")
# 76

# Get specific atoms
ca_atoms = structure[structure['atom_name'] == 'CA']

# Get specific residues
active_site = structure[structure['residue_number'].isin([35, 36, 37])]

# Calculate properties
center_of_mass = cp.calculate_center_of_mass(structure)
radius_of_gyration = cp.calculate_radius_of_gyration(structure)
```

### Sequence Extraction

```python
# Extract sequence from structure
sequence = cp.extract_sequence("1ubq")
# Returns: "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

# Extract from specific chain
sequence_a = cp.extract_sequence("1ubq", chain="A")

# Get all sequences as dictionary
sequences = cp.get_seq_dict()
# {'1ubq_A': 'MQIFVKTLTG...'}

# Save extracted sequences
cp.save_sequences(sequences, "extracted_sequences.fasta")
```

### Dataset Management

```python
# Create dataset
cp.create_dataset(
    dataset_id="kinases",
    name="Protein kinase structures",
    description="Human kinase structures for drug discovery",
    content=["1atp", "2src", "3erk", "4bcr"]
)

# Load dataset
cp.load_dataset("kinases")

# Access loaded data
for pdb_id in cp.pdb_ids:
    structure = cp.data[cp.data['pdb_id'] == pdb_id]
    # Analyze structure...

# List available datasets
datasets = cp.list_datasets()
```

## Advanced Features

### Structure Filtering

```python
# Filter by resolution
high_res_structures = cp.filter_by_resolution(max_resolution=2.0)

# Filter by experimental method
xray_structures = cp.filter_by_method("X-RAY DIFFRACTION")

# Filter by organism
human_structures = cp.filter_by_organism("Homo sapiens")

# Complex filtering
filtered = cp.filter_structures(
    resolution_range=(1.0, 2.5),
    methods=["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"],
    min_length=100,
    max_length=500
)
```

### Coordinate Operations

```python
# Extract coordinate array
coords = cp.get_coordinates("1ubq", atoms="CA")
# Returns numpy array of shape (n_residues, 3)

# Calculate distances
distances = cp.calculate_distance_matrix("1ubq", atoms="CA")

# Align structures
aligned = cp.align_structures("1ubq", "2gb1", method="kabsch")

# Calculate RMSD
rmsd = cp.calculate_rmsd("1ubq", "2gb1", atoms="backbone")
```

### Chain Operations

```python
# Extract single chain
chain_a = cp.extract_chain("1ubq", chain="A")

# Rename chains
cp.rename_chain("1ubq", old_chain="A", new_chain="B")

# Merge structures
merged = cp.merge_structures(["1ubq", "2gb1"], new_name="merged_complex")

# Split by chains
chains = cp.split_by_chains("1abc")
# Returns: {"1abc_A": chain_a_data, "1abc_B": chain_b_data}
```

### Residue Operations

```python
# Get residue information
residues = cp.get_residues("1ubq", chain="A")

# Extract residue range
segment = cp.extract_residue_range("1ubq", start=10, end=20, chain="A")

# Mutate residue (in silico)
mutated = cp.mutate_residue("1ubq", chain="A", residue_num=35, new_residue="ALA")

# Find missing residues
missing = cp.find_missing_residues("1ubq")
```

## Integration with Other Tools

### Structure Alignment with FoldMason

```python
# Align multiple structures
alignment_result = cp.align_multiple_structures(
    ["1ubq", "2gb1", "1crn"],
    method="foldmason",
    output_name="aligned_structures"
)

# Access alignment
aligned_coords = alignment_result['coordinates']
alignment_scores = alignment_result['scores']
```

### Structure Download

```python
# Download from PDB
cp.download_structure("1ubq", source="pdb")

# Download from AlphaFold
cp.download_structure("P12345", source="alphafold")

# Batch download
structures = ["1ubq", "2gb1", "1crn", "3nir"]
cp.download_structures(structures, source="pdb")
```

### Format Conversion

```python
# Convert PDB to CIF
cp.convert_format("1ubq", from_format="pdb", to_format="cif")

# Save in different format
structure = cp.load_structure("1ubq")
cp.save_structure_as_pdb("1ubq", "output.pdb")
```

## Working with Metadata

### Structure Metadata

```python
# Get metadata
metadata = cp.get_structure_metadata("1ubq")
# {
#     "resolution": 1.8,
#     "method": "X-RAY DIFFRACTION",
#     "space_group": "P 21 21 21",
#     "cell_dimensions": {...},
#     "deposition_date": "1995-07-10",
#     "authors": ["Vijay-Kumar, S.", "Bugg, C.E."],
#     "organism": "Homo sapiens"
# }

# Add custom metadata
cp.add_metadata("1ubq", {
    "project": "Ubiquitin study",
    "analyzed_by": "J. Smith",
    "analysis_date": "2024-01-15"
})
```

### Quality Metrics

```python
# Calculate structure quality metrics
quality = cp.assess_structure_quality("1ubq")
# {
#     "resolution": 1.8,
#     "r_factor": 0.176,
#     "r_free": 0.221,
#     "clashscore": 4.2,
#     "ramachandran_outliers": 0.0,
#     "rotamer_outliers": 1.3
# }
```

## Visualization Support

```python
# Generate PyMOL script
cp.generate_pymol_script("1ubq", 
    highlight_residues=[35, 36, 37],
    color_by="chain",
    output="visualize_1ubq.pml"
)

# Export for visualization
cp.export_for_visualization("1ubq",
    format="pdb",
    include_heterogens=False,
    output="1ubq_clean.pdb"
)

# Create structure summary report
cp.create_structure_report("1ubq", output="1ubq_report.html")
```

## Best Practices

### 1. Type Enforcement

Always use `apply_dtypes=True` for numerical operations:

```python
# Ensure correct data types
structure = cp.load_structure("1ubq", apply_dtypes=True)

# Critical for coordinate operations
coords = structure[['x', 'y', 'z']].values  # Now guaranteed numeric
```

### 2. Chain Handling

Be explicit about chains:

```python
# Check available chains first
chains = cp.get_chains("1abc")

# Process each chain
for chain in chains:
    chain_data = cp.extract_chain("1abc", chain=chain)
    # Process chain...
```

### 3. Memory Management

For large structures or datasets:

```python
# Load structures one at a time
for pdb_id in large_dataset:
    structure = cp.load_structure(pdb_id)
    # Process and save results
    result = analyze_structure(structure)
    save_result(result)
    # Clear from memory
    del structure
```

### 4. Error Handling

```python
# Handle missing structures gracefully
structures_analyzed = []
for pdb_id in pdb_list:
    try:
        structure = cp.load_structure(pdb_id)
        if structure is not None and not structure.empty:
            result = analyze_structure(structure)
            structures_analyzed.append(pdb_id)
    except Exception as e:
        cp.logger.warning(f"Failed to process {pdb_id}: {e}")
        continue
```

## Common Workflows

### Workflow 1: Structure Quality Analysis

```python
def analyze_structure_quality(pdb_ids):
    """Comprehensive structure quality analysis."""
    results = []
    
    for pdb_id in pdb_ids:
        # Load structure
        structure = cp.load_structure(pdb_id, apply_dtypes=True)
        if structure is None:
            continue
            
        # Get metadata
        metadata = cp.get_structure_metadata(pdb_id)
        
        # Calculate metrics
        result = {
            "pdb_id": pdb_id,
            "resolution": metadata.get("resolution"),
            "method": metadata.get("method"),
            "chain_count": len(cp.get_chains(pdb_id)),
            "residue_count": cp.get_residue_count(pdb_id),
            "missing_residues": len(cp.find_missing_residues(pdb_id))
        }
        
        results.append(result)
    
    return pd.DataFrame(results)
```

### Workflow 2: Homolog Analysis

```python
def analyze_homologs(reference_pdb, homolog_pdbs):
    """Compare reference structure to homologs."""
    # Load reference
    ref_structure = cp.load_structure(reference_pdb)
    ref_sequence = cp.extract_sequence(reference_pdb)
    
    results = []
    for homolog in homolog_pdbs:
        # Load and align
        homolog_structure = cp.load_structure(homolog)
        alignment = cp.align_structures(reference_pdb, homolog)
        
        # Calculate similarity metrics
        rmsd = alignment["rmsd"]
        seq_identity = calculate_sequence_identity(
            ref_sequence,
            cp.extract_sequence(homolog)
        )
        
        results.append({
            "homolog": homolog,
            "rmsd": rmsd,
            "sequence_identity": seq_identity,
            "aligned_residues": alignment["n_aligned"]
        })
    
    return pd.DataFrame(results)
```

### Workflow 3: Active Site Analysis

```python
def analyze_active_site(pdb_id, active_site_residues):
    """Detailed analysis of active site region."""
    structure = cp.load_structure(pdb_id, apply_dtypes=True)
    
    # Extract active site
    active_site = structure[
        structure['residue_number'].isin(active_site_residues)
    ]
    
    # Find nearby residues (within 5Å)
    nearby = cp.find_nearby_residues(
        pdb_id,
        active_site_residues,
        distance=5.0
    )
    
    # Calculate properties
    analysis = {
        "active_site_atoms": len(active_site),
        "nearby_residues": len(nearby),
        "solvent_accessible": cp.calculate_sasa(active_site),
        "pocket_volume": cp.calculate_pocket_volume(active_site)
    }
    
    return analysis
```

## Troubleshooting

### Common Issues

1. **Missing coordinates**
```python
# Check for missing coordinates
structure = cp.load_structure("1abc")
missing = structure[structure[['x', 'y', 'z']].isna().any(axis=1)]
if not missing.empty:
    print(f"Warning: {len(missing)} atoms have missing coordinates")
```

2. **Chain ID mismatches**
```python
# Use auth_chain_id for author-defined chains
structure = cp.load_structure("1abc")
# Check both chain_id and auth_chain_id
print(structure[['chain_id', 'auth_chain_id']].value_counts())
```

3. **Memory issues with large structures**
```python
# Process in chunks
chunk_size = 10000
for chunk in pd.read_csv(large_structure_file, chunksize=chunk_size):
    process_chunk(chunk)
```

## Summary

CifBaseProcessor provides:
- Comprehensive structure data management
- Rich analysis capabilities
- Integration with structural biology tools
- Cross-format operations
- Metadata and quality assessment
- Dataset organization

The processor abstracts the complexity of structure file formats while providing powerful tools for structural analysis and manipulation.