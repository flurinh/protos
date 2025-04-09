# Processor Classes Documentation

This document provides detailed information about the key processor classes that facilitate data interactions between different datatypes in the Protos framework.

## Overview

The framework uses four main processor classes that serve as the backbone for data manipulation:

1. **GRNProcessor (GRNP)**: Manages Generic Residue Numbering (GRN) for protein sequences
2. **EMBProcessor (EMBP)**: Handles protein sequence embeddings 
3. **CifProcessor (CP)**: Processes protein structures from CIF/PDB files
4. **PropertyProcessor (PP)**: Manages protein properties and metadata

These classes are designed to work together seamlessly, providing a cohesive data flow for protein analysis:

```
┌───────────────┐       ┌───────────────┐       ┌───────────────┐
│               │       │               │       │               │
│  CifProcessor ├───────┤  GRNProcessor ├───────┤ EMBProcessor  │
│     (CP)      │       │    (GRNP)     │       │    (EMBP)     │
│               │       │               │       │               │
└───────┬───────┘       └───────┬───────┘       └───────┬───────┘
        │                       │                       │
        │                       │                       │
        │                       ▼                       │
        │               ┌───────────────┐               │
        └───────────────┤PropertyProcessor│──────────────┘
                        │      (PP)      │
                        └───────────────┘
```

## Path Management in Processors

All processor classes inherit from `BaseProcessor` which implements the core path management system:

1. **Automatic Processor Type Detection**:
   - When a processor is instantiated, it automatically determines its type from the class name
   - Example: `CifProcessor` becomes "structure", `GRNProcessor` becomes "grn"

2. **Path Resolver Creation**:
   - A `ProtosPaths` instance is created to handle path resolution
   - Resolves reference data (package resources) and user data (runtime created)

3. **Default Path Construction**:
   - Data paths follow a convention: `<data_root>/<processor_type>/<subdir>`
   - Example: `data/structure/mmcif` for structure files

4. **Cross-Platform Compatibility**:
   - Paths use `pathlib.Path` for consistent handling across Windows/Linux
   - Handles path separators (\ vs. /) automatically

5. **Reference vs. User Data**:
   - File operations check both reference and user data sources
   - First tries user data, then falls back to reference data
   - Example: `load_structure()` tries `data/structure/mmcif/<file>` then `protos/reference_data/structure/mmcif/<file>`

## GRNProcessor (GRNP)

### Purpose
The GRNProcessor class provides a framework for working with Generic Residue Numbering (GRN) systems, which allow for consistent residue numbering across different protein sequences. This is especially important for GPCRs and opsins where structural alignment is critical.

### Key Attributes
- `data`: DataFrame containing GRN mapping information for multiple proteins
- `ids`: List of protein identifiers
- `grns`: List of GRN positions
- `maps`: Dictionary mapping between GRN and sequence positions
- `features`: DataFrame for associated protein features
- `path_resolver`: ProtosPaths instance for path management

### Primary Methods
1. **Loading and Saving**
   - `load_grn_table(dataset)`: Load a GRN table from a CSV file
   - `save_grn_table(dataset)`: Save the current GRN table
   - `load_and_merge_grn_tables(datasets)`: Merge multiple GRN tables

2. **Data Manipulation**
   - `get_seq_dict()`: Convert GRN table to a dictionary of sequences
   - `filter_data_by_occurances(threshold)`: Filter GRN positions by occurrence frequency
   - `filter_by_ids(ids)`: Filter data to include only specified proteins
   - `apply_interval(grn_interval)`: Limit analysis to a specific GRN interval

3. **GRN Conversion**
   - `get_grn_dict()`: Get a dictionary of GRN positions for each protein

### Usage Example
```python
# Initialize GRN processor (paths automatically set up)
grnp = GRNProcessor()

# Paths are automatically set:
# - User data: data/grn/tables/
# - Reference data: protos/reference_data/grn/tables/

# Load a GRN table (checks both user and reference data)
grnp.load_grn_table('gpcrs')

# Get sequences from GRN table
seq_dict = grnp.get_seq_dict()

# Filter by specific proteins
grnp.filter_by_ids(['protein1', 'protein2'])

# Limit to specific GRN positions
grnp.apply_interval(['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50'])
```

## EMBProcessor (EMBP)

### Purpose
The EMBProcessor class handles protein sequence embeddings, providing a framework to generate, store, and manipulate embedding vectors for protein sequences. These embeddings capture the semantic information of protein sequences and can be used for various machine learning tasks.

### Key Attributes
- `emb_dict`: Dictionary mapping protein IDs to embedding matrices
- `model_manager`: Manages the embedding model (ESM, Ankh, etc.)
- `ids`: List of protein identifiers
- `device`: Device for computation (CPU, CUDA)
- `path_resolver`: ProtosPaths instance for path management

### Primary Methods
1. **Embedding Generation**
   - `emb_grnp(grnp)`: Generate embeddings for sequences in a GRNProcessor
   - `emb_seq_dict(seq_dict)`: Generate embeddings for a dictionary of sequences

2. **Data Management**
   - `save_dataset(dataset_name)`: Save embeddings to a pickle file
   - `load_dataset(dataset_name)`: Load embeddings from a pickle file
   - `filter_by_ids(ids)`: Filter embeddings to include only specified proteins

3. **Embedding Transformations**
   - `get_grn_embeddings(grn_list, grnp)`: Extract embeddings for specific GRN positions
   - `map_embeddings_to_grns(grn_list, grnp)`: Map embeddings to GRN positions
   - `aggr_embeddings(embeddings, operation)`: Aggregate embeddings using various operations

### Usage Example
```python
# Initialize embedding processor with a model (paths automatically set up)
model_manager = EMBModelManager(model_name='esm2', device='cuda')
embp = EMBProcessor(model_manager=model_manager)

# Paths are automatically set:
# - User data: data/embedding/datasets/
# - Reference data: protos/reference_data/embedding/datasets/

# Generate embeddings for sequences in GRN processor
embp.emb_grnp(grnp)

# Save the embeddings (to user data path)
embp.save_dataset('gpcr_embeddings')

# Extract embeddings for binding pocket residues
binding_pocket_grns = ['3.32', '3.33', '5.42', '5.46', '6.48', '6.51', '7.39']
binding_pocket_embs = embp.get_grn_embeddings(binding_pocket_grns, grnp)
```

## CifProcessor (CP)

### Purpose
The CifProcessor class handles protein structure data, processing PDB/CIF files to extract coordinates, chain information, and other structural data. It facilitates structure analysis and visualization.

### Key Attributes
- `data`: DataFrame containing atom coordinates and structural information
- `pdb_ids`: List of PDB identifiers
- `structure_filenames`: List of filenames for structure files
- `structure_info`: Dictionary with structure metadata
- `path_resolver`: ProtosPaths instance for path management

### Path Handling
CifProcessor uses the refactored path system to manage file locations:

1. **Automatic Path Configuration**:
   - Determines processor type (`structure`) from class name
   - Sets up a ProtosPaths resolver for path management
   - Creates standard directories if they don't exist

2. **Cross-Platform Compatibility**:
   - Uses `pathlib.Path` for consistent path handling across Windows/Linux
   - Properly handles file separators (\ vs. /)
   - Works with relative and absolute paths

3. **Reference vs. User Data**:
   - Automatically checks both reference and user data locations
   - Tries user data location first, then falls back to reference data

### Primary Methods
1. **Structure Loading**
   - `load_structure(pdb_id)`: Loads a structure from PDB/CIF file, checking both user and reference data
   - `get_available_pdb_files()`: Lists available PDB files from both user and reference directories
   - `load_data(version)`: Loads preprocessed structural data

2. **Structure Filtering and Manipulation**
   - `filter_by_ids(pdb_ids)`: Filters structures to include only specified PDBs
   - `reset_index()`: Resets the index of the structure DataFrame
   - `get_chains(pdb_id)`: Gets available chains for a structure

3. **Structure Analysis**
   - `get_chain_residues(pdb_id, chain_id)`: Extracts residues from a specific chain
   - `get_ca_coordinates(pdb_id, chain_id)`: Gets alpha carbon coordinates
   - `get_backbone_coordinates(pdb_id, chain_id)`: Gets backbone atom coordinates
   - `extract_binding_pocket(pdb_id, ligand, distance)`: Extracts atoms in binding pocket around a ligand
   
4. **Structure Export and Temporary Files**
   - `to_cif(pdb_id)`: Converts structure data to CIF format text
   - `save_cif(pdb_id, output_path, versioned, force_overwrite)`: Saves structure to CIF file
   - `save_temp_cif(pdb_id, suffix)`: Creates temporary CIF file with unique name
   - `cleanup_temp_files(older_than_hours)`: Removes temporary CIF files
   - `export_structures(pdb_ids, output_dir, validate)`: Exports multiple structures to CIF files

### Usage Example
```python
# Initialize structure processor (path handling is automatic)
cp = CifProcessor()

# Paths are automatically set:
# - User data: data/structure/mmcif/
# - Reference data: protos/reference_data/structure/mmcif/

# Load a specific structure (checks both user and reference data)
cp.load_structure('6CMO')

# Get chains for the structure
chains = cp.get_chains('6CMO')

# Extract CA coordinates for chain A
ca_coords = cp.get_ca_coordinates('6CMO', 'A')

# Create temporary files for alignment (saved to data/structure/temp_cif/)
temp_file1 = cp.save_temp_cif('6CMO', suffix='for_alignment')
temp_file2 = cp.save_temp_cif('7ZOU', suffix='for_alignment')
```

## PropertyProcessor (PP)

### Purpose
The PropertyProcessor class manages protein properties and metadata, providing a framework to associate various properties with protein identifiers. It allows for filtering, merging, and analyzing protein properties.

### Key Attributes
- `identity`: DataFrame containing protein identifiers and basic information
- `properties`: DataFrame containing protein properties
- `metadata`: Dictionary with dataset metadata
- `path_resolver`: ProtosPaths instance for path management

### Primary Methods
1. **Data Management**
   - `set_dataframe(dataframe_type, dataframe)`: Set identity or properties DataFrame
   - `load_dataset(dataset)`: Load a dataset from file
   - `save_dataset(dataset)`: Save current data to file

2. **Filtering and Transformation**
   - `filter(df_type, filters, map_to_other)`: Apply filters to dataframes
   - `add_new_column(new_data, data_type)`: Add new property data
   - `synchronize_dataframes(common_ids)`: Ensure consistent protein IDs across dataframes

3. **Analysis and Visualization**
   - `plot_property_distribution(property_name)`: Visualize property distributions
   - `correlate_properties(prop1, prop2)`: Calculate correlation between properties
   - `group_by_property(property_name)`: Group proteins by a specific property

### Usage Example
```python
# Initialize property processor (path handling is automatic)
pp = PropertyProcessor()

# Paths are automatically set:
# - User data: data/property/datasets/
# - Reference data: protos/reference_data/property/datasets/

# Load a dataset (checks both user and reference data)
pp.load_dataset('gpcr_properties')

# Add identity information
pp.set_dataframe('identity', identity_df)

# Add property information
pp.set_dataframe('properties', properties_df)

# Filter by protein family
filtered_df = pp.filter('properties', {'family__eq': 'GPCR_A'}, inplace=True)

# Add new experimental data
pp.add_new_column(new_experiment_data, data_type='properties')
```

## Data Flow Between Processors

### GRNP to EMBP
The GRNProcessor provides sequence information that the EMBProcessor uses to generate embeddings:

```python
# Generate embeddings for sequences in GRN processor
grnp = GRNProcessor(dataset='gpcrs')
embp = EMBProcessor(model_manager=model_manager)
embp.emb_grnp(grnp)
```

### CP to GRNP
The CifProcessor provides structural information that can be mapped to GRN positions:

```python
# Extract sequence from structure
cp = CifProcessor()
structure = cp.load_structure('6CMO')
seq = cp.get_sequence('6CMO', 'A')

# Map to GRN positions
grnp = GRNProcessor()
grnp.assign_grns_to_sequence(seq, pdb_id='6CMO')
```

### GRNP/EMBP/CP to PP
The PropertyProcessor can integrate data from all other processors:

```python
# Create property processor
pp = PropertyProcessor()

# Add GRN information as properties
grnp = GRNProcessor(dataset='gpcrs')
grn_data = {'protein_id': grnp.ids, 'has_grn': [True] * len(grnp.ids)}
pp.add_new_column(pd.DataFrame(grn_data), data_type='properties')

# Add embedding information
embp = EMBProcessor()
embp.load_dataset('gpcr_embeddings')
emb_data = {'protein_id': embp.ids, 'has_embedding': [True] * len(embp.ids)}
pp.add_new_column(pd.DataFrame(emb_data), data_type='properties')

# Add structure information
cp = CifProcessor()
struct_ids = cp.get_available_pdb_files()
struct_data = {'protein_id': struct_ids, 'has_structure': [True] * len(struct_ids)}
pp.add_new_column(pd.DataFrame(struct_data), data_type='properties')
```

## Best Practices

1. **Consistent Identifiers**: Use consistent protein identifiers across all processors
2. **Data Validation**: Validate input data before processing
3. **Memory Management**: Use filtering methods to reduce memory usage
4. **Error Handling**: Implement proper error handling for missing data
5. **Documentation**: Document custom processing pipelines
6. **Path Handling**: Use PathLib for cross-platform compatibility
7. **Checking Both Sources**: Always check both reference and user data when loading files

## Common Workflows

### Structure-Based Embedding Analysis
```python
# 1. Load structure
cp = CifProcessor()
cp.load_structure('6CMO')

# 2. Extract sequence and assign GRNs
seq = cp.get_sequence('6CMO', 'A')
grnp = GRNProcessor()
grnp.assign_grns_to_sequence(seq, pdb_id='6CMO')

# 3. Generate embeddings
embp = EMBProcessor(model_manager=model_manager)
embp.emb_grnp(grnp)

# 4. Map embeddings to structure
binding_pocket_grns = ['3.32', '3.33', '5.42', '5.46', '6.48', '6.51', '7.39']
binding_pocket_embs = embp.get_grn_embeddings(binding_pocket_grns, grnp)
```

### Structure Alignment with External Tools
```python
# 1. Load structures
cp = CifProcessor()
cp.load_structure('6CMO')
cp.load_structure('7ZOU')

# 2. Create temporary files for alignment
temp_file1 = cp.save_temp_cif('6CMO', suffix='for_alignment')
temp_file2 = cp.save_temp_cif('7ZOU', suffix='for_alignment')

# 3. Align with FoldMason
from protos.processing.structure.foldmason import FoldMason
fm = FoldMason(use_wsl=True)  # Set to False if running natively on Linux/Mac
alignment = fm.easy_msa(
    input_files=[temp_file1, temp_file2],
    output_prefix="my_alignment",
    tmp_folder="tmp"
)

# 4. Align with GTalign
from protos.processing.structure.gtalign import GTalign
gta = GTalign()
result = gta.align(
    query_structures=temp_file1,
    reference_structures=temp_file2,
    tm_score_threshold=0.0
)

# 5. Clean up temporary files
cp.cleanup_temp_files(older_than_hours=1)
```

### Property-Based Filtering and Analysis
```python
# 1. Load property data
pp = PropertyProcessor(dataset='gpcr_properties')

# 2. Filter by specific property
pp.filter('properties', {'ligand_type__eq': 'peptide'}, inplace=True)

# 3. Load GRN data for filtered proteins
grnp = GRNProcessor()
grnp.load_grn_table('gpcrs')
grnp.filter_by_ids(pp.identity['protein_id'].tolist())

# 4. Generate embeddings for filtered set
embp = EMBProcessor(model_manager=model_manager)
embp.emb_grnp(grnp)
```

## Extending the Processors

The processor classes are designed to be extensible. Common extension patterns include:

1. **Subclassing**: Create specialized subclasses for specific protein families
2. **Pipeline Creation**: Create pipeline classes that combine multiple processors
3. **New Processors**: Develop new processors for specific data types
4. **Custom Methods**: Add custom methods to existing processors

Example of a custom pipeline:
```python
class ProteinAnalysisPipeline:
    def __init__(self, protein_ids):
        self.protein_ids = protein_ids
        self.grnp = GRNProcessor()
        self.embp = EMBProcessor()
        self.cp = CifProcessor()
        self.pp = PropertyProcessor()
        
    def load_data(self, grn_dataset, emb_dataset, prop_dataset):
        self.grnp.load_grn_table(grn_dataset)
        self.grnp.filter_by_ids(self.protein_ids)
        
        self.embp.load_dataset(emb_dataset)
        self.embp.filter_by_ids(self.protein_ids)
        
        self.pp.load_dataset(prop_dataset)
        self.pp.filter('identity', {'protein_id__is_in': self.protein_ids}, inplace=True)
        
    def analyze_binding_pocket(self, binding_pocket_grns):
        # Extract binding pocket embeddings
        binding_pocket_embs = self.embp.get_grn_embeddings(binding_pocket_grns, self.grnp)
        
        # Get properties for these proteins
        binding_props = self.pp.properties[self.pp.properties['protein_id'].isin(self.protein_ids)]
        
        # Return combined analysis
        return binding_pocket_embs, binding_props
```