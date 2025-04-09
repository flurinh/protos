# Protos File Format Specifications

This document defines the standard file formats used within the Protos framework. All processors and utilities should adhere to these specifications for consistency and interoperability.

## Table of Contents

- [Data Registry](#data-registry)
- [Sequence Formats](#sequence-formats)
  - [FASTA](#fasta)
- [Structure Formats](#structure-formats)
  - [PDB](#pdb)
  - [mmCIF](#mmcif)
- [Table Formats](#table-formats)
  - [GRN Tables](#grn-tables)
  - [Property Tables](#property-tables)
- [Embedding Formats](#embedding-formats)
  - [Embedding Pickle](#embedding-pickle)
- [Graph Formats](#graph-formats)
  - [Graph JSON](#graph-json)
- [Temporary Files](#temporary-files)
- [External Tool Integration](#external-tool-integration)
- [Conversion Guidelines](#conversion-guidelines)

## Data Registry

The data registry (`registry.json`) maps dataset identifiers to file paths and metadata:

```json
{
  "dataset_id": {
    "path": "tables/dataset.csv",  // Relative path within processor directory
    "metadata": {
      "processor_type": "grn",
      "dataset_type": "table",
      "source": "reference",  // reference or user
      "rows": 250,
      "columns": ["3.50", "6.48", "7.49"],
      "created_at": "2023-08-15T12:34:56",
      "description": "GRN reference table for opsins"
    },
    "timestamp": "2023-08-15T12:34:56"
  }
}
```

Processor-specific registries are stored at `data/<processor_type>/registry.json` while the global registry is at `data/global_registry.json`. Both reference and user data have their own registries, but the path resolution system automatically checks both locations when loading files.

## Sequence Formats

### FASTA

FASTA files store biological sequences:

```
>sequence_id [optional description]
AGCTAGCTAGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCTAGCTAGCT
>another_sequence_id
AGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

**Format specifications:**
- Each sequence typically starts with `>` followed by an identifier
- Sequence can span multiple lines
- No blank lines within a sequence
- Standard alphabet: A-Z for proteins, ACGTU for nucleic acids
- For robust parsing, files without `>` prefix may be treated as a single unnamed sequence
- Multiple sequences can be stored in a single file
- Description text after the identifier is optional and may contain metadata

## Structure Formats

### PDB

PDB (Protein Data Bank) files store protein structures:

**Format specifications:**
- Fixed-width columnar format that requires precise positioning
- HEADER, TITLE, AUTHOR records for metadata
- ATOM records for atomic coordinates with:
  - Atom serial number (columns 7-11)
  - Atom name (columns 13-16)
  - Residue name (columns 18-20)
  - Chain identifier (column 22)
  - Residue sequence number (columns 23-26)
  - X, Y, Z coordinates (columns 31-54)
  - Occupancy and temperature factor (columns 55-66)
  - Element symbol (columns 77-78)
- HETATM records for non-standard residues (same format as ATOM)
- TER records for chain terminators
- END record at file end
- Important: Some tools (like GTalign) require strict adherence to column positions

### mmCIF

mmCIF (macromolecular Crystallographic Information File) files store detailed structural data:

**Format specifications:**
- Data organized in key-value pairs or tables
- Categories prefixed with underscore (e.g., `_atom_site.group_PDB`)
- Multiple entries formatted as loops with `loop_` prefix
- Block structure delimited by `data_` entries
- More comprehensive than PDB format with explicit data typing
- Important: Some tools may require specific mmCIF tags and formatting conventions

## Table Formats

### GRN Tables

GRN (Generic Residue Numbering) tables map sequence positions to standardized numbering schemes:

```csv
protein_id,3.50,6.48,7.49
opsin1,M125,Y261,P296
opsin2,M123,Y259,P294
```

**Format specifications:**
- CSV format with protein_id as the first column
- Each additional column represents a GRN position
- Values are formatted as `<amino acid><sequence position>` (e.g., "M125" for Methionine at position 125)
- Empty cells indicate absent positions
- Amino acid is provided to facilitate validation and cross-referencing

### Property Tables

Property tables store metadata and properties for proteins:

**Identity Table (`dataset_identity.csv`):**
```csv
protein_id,uniprot_id,name,organism
opsin1,O12345,Rhodopsin,Homo sapiens
opsin2,O23456,Blue opsin,Bos taurus
```

**Properties Table (`dataset_properties.csv`):**
```csv
protein_id,activation,spectral_sensitivity,ligand_binding
opsin1,light,500,all-trans-retinal
opsin2,light,460,11-cis-retinal
```

**Metadata File (`dataset_metadata.pkl`):**
```python
{
    "description": "Dataset of opsin proteins with spectral properties",
    "version": "1.0",
    "created": "2023-08-15T12:34:56",
    "source": "UniProt database",
    "columns": {
        "activation": {"type": "categorical", "values": ["light", "chemical"]},
        "spectral_sensitivity": {"type": "numeric", "unit": "nm"},
        "ligand_binding": {"type": "categorical"}
    }
}
```

**Format specifications:**
- Paired CSV files sharing protein_id as the key
- Identity table contains core identifiers
- Properties table contains measured/functional data
- Property tables can be extended with custom columns
- Optional metadata file (pickle format) containing dataset information, column descriptions, and other metadata

## Embedding Formats

### Embedding Pickle

Embedding files store vector representations of sequences or structures:

```python
{
  "protein_id1": np.array([...]),  # Shape: [sequence_length, embedding_dim]
  "protein_id2": np.array([...])
}
```

**Format specifications:**
- Pickle (.pkl) format containing a dictionary
- Keys are protein identifiers
- Values are NumPy arrays of embeddings
- Embedding dimension should be consistent

## Graph Formats

### Graph JSON

Graph files store protein structure graphs:

```json
{
  "protein_id": {
    "nodes": [
      {"id": 0, "residue": "ALA", "position": 1, "grn": "1.50", "coords": [1.0, 2.0, 3.0]},
      {"id": 1, "residue": "GLY", "position": 2, "grn": null, "coords": [4.0, 5.0, 6.0]}
    ],
    "edges": [
      {"source": 0, "target": 1, "weight": 0.95, "type": "covalent"},
      {"source": 0, "target": 2, "weight": 0.82, "type": "spatial"}
    ]
  }
}
```

**Format specifications:**
- JSON format containing graph definitions
- Nodes represent amino acid residues
- Edges represent connections between residues
- Node attributes include residue type, position, and spatial coordinates
- Edge attributes include connection type and weight

## Temporary Files

Temporary files are used for intermediate processing and external tool integration.

**Format specifications:**
- Use the `TempFileManager` for creating and tracking temporary files
- Prefer context managers for automatic cleanup
- Use standardized prefixes (`protos_`) for easy identification
- Group temporary files in dedicated subdirectories when possible

**Example usage with context manager:**
```python
# Create a temporary file that is automatically cleaned up
with temp_manager.temp_file(suffix=".fasta") as temp_path:
    # Write data to the file
    with open(temp_path, 'w') as f:
        f.write(">seq1\nAGCT")
    
    # Use the file
    process_fasta(temp_path)
    
# File is automatically deleted when the context exits
```

**Example usage with explicit cleanup:**
```python
# Create a temporary file without automatic cleanup
temp_path = temp_manager.create_temp_file(suffix=".fasta")

try:
    # Write data to the file
    with open(temp_path, 'w') as f:
        f.write(">seq1\nAGCT")
    
    # Use the file
    process_fasta(temp_path)
    
finally:
    # Explicitly clean up the file
    temp_manager.cleanup(temp_path)
```

## External Tool Integration

External bioinformatics tools often require specific file formats and temporary workspaces.

**Format specifications:**
- Use the `ExternalToolHelper` for preparing tool inputs and outputs
- Create temporary workspaces for tools with multiple intermediate files
- Convert between internal data representations and tool-specific formats
- Handle tool-specific format requirements and extensions

**Example usage with FoldMason:**
```python
# Create a tool helper for FoldMason
foldmason_helper = ExternalToolHelper("foldmason")

# Prepare input sequences
with foldmason_helper.tool_input(sequences, format_type="fasta") as input_file:
    # Create a workspace
    with foldmason_helper.tool_workspace() as workspace:
        # Prepare output file
        with foldmason_helper.tool_output(suffix=".db", expected_format="sqlite") as (output_file, read_output):
            # Run FoldMason
            subprocess.run(["foldmason", "easy_msa", 
                           "--input", input_file,
                           "--output", output_file,
                           "--tmp", workspace])
            
            # Read results
            results = read_output()
```

**Example usage with GTalign:**
```python
# Create a tool helper for GTalign
gtalign_helper = ExternalToolHelper("gtalign")

# Create a temporary copy of a structure file
temp_structure = gtalign_helper.create_temp_copy(structure_file)

# Run GTalign with the temporary file
results = subprocess.run(["gtalign", "align", 
                         "--query", temp_structure,
                         "--reference", reference_file])

# Clean up
temp_manager.cleanup(temp_structure)
```

## Conversion Guidelines

When converting between formats, follow these principles:

1. **Lossless conversion**: Preserve all original data when possible
2. **Identifiers**: Maintain consistent identifiers across formats
3. **Metadata**: Include appropriate metadata in the new format
4. **Validation**: Validate converted data against format specifications
5. **Documentation**: Note any conversion limitations or special handling

## Path Resolution

The Protos path system automatically handles file location for all formats described in this document:

1. **Processor-Type Detection**: 
   - Each processor type (structure, grn, sequence, etc.) has its own directory
   - Processor classes automatically determine their type from class name

2. **Reference vs. User Data**:
   - Reference data: Read-only data distributed with package (`protos/reference_data/`)
   - User data: Read-write data created at runtime (`data/`)

3. **Automatic Path Resolution**:
   - When loading files, system checks both reference and user locations
   - Example: `cp.load_structure("1abc")` checks:
     - First: `data/structure/mmcif/1abc.cif` (user data)
     - Then: `protos/reference_data/structure/mmcif/1abc.cif` (reference data)

4. **Cross-Platform Compatibility**:
   - Paths use `pathlib.Path` for consistent handling across Windows/Linux
   - File separators (\ vs. /) are handled automatically

---

Note: This document will be updated as new formats are introduced or existing formats are refined. All data handling code should refer to these specifications to ensure consistency across the Protos framework.