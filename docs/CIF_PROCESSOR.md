# CIF Processor Documentation

## Overview

The CIF (Crystallographic Information File) processor system in Protos provides a comprehensive framework for loading, manipulating, and storing protein structures in mmCIF format. This document outlines the complete architecture of the CIF processing system, from base classes to specific utilities.

## Core Components

The CIF processor system is built on several interconnected components:

### 1. BaseProcessor Framework

All processors in Protos, including the CIF processor, inherit from the `BaseProcessor` abstract base class. This design provides:

- Consistent API across different processor types
- Standardized path resolution and data management
- Integration with the dataset and registry systems
- Configurable data flow between processors

### 2. CIF-specific Processors

#### StructBaseProcessor / CifBaseProcessor

This processor extends the base framework with structure-specific functionality:

- Loading and processing mmCIF structure files
- Managing collections of protein structures
- Providing structure-specific interfaces for data extraction
- Facilitating operations on protein chains, residues, and atoms
- Integration with structure alignment tools

#### StructProcessor / CifProcessor

Built on top of the base processor, this implementation adds:

- Advanced filtering of structure data
- Structure-specific data manipulations
- Domain-specific functions for structural biology
- Calculation of structural properties

### 3. CIF Handling (io/cif_handler.py)

The `CifHandler` class handles the actual file operations:

- Reading mmCIF files into internal data structures
- Writing processed data back to mmCIF format
- Validation of mmCIF data for consistency
- Converting between different representation formats
- Merging structures from different sources
- Extracting biological assemblies

### 4. CIF Utilities (io/cif_utils.py)

The utilities module provides lower-level functions for:

- Parsing mmCIF format into pandas DataFrames
- Converting DataFrames back to mmCIF format
- Validating CIF data structure and content
- Performing common operations on CIF data
- Efficient handling of large structure files

### 5. Path Management (io/paths/)

The path system handles:

- Resolution of file paths for structures
- Management of directory structure
- Distinction between reference and user data
- Platform-independent path handling
- Environment variable configuration

### 6. Registry and Dataset System

The registry and dataset systems provide:

- Tracking of available structure datasets
- Standard dataset definitions
- User-created dataset management
- Metadata for datasets and structures
- Cross-referencing between different data types

## Data Flow

The typical flow of data through the CIF processor system is:

1. **Input**: mmCIF files are read using the CifHandler
2. **Parsing**: CIF data is converted to pandas DataFrames
3. **Processing**: StructProcessor applies transformations and analyses
4. **Organization**: Results are organized into datasets
5. **Output**: Processed data is written back to mmCIF or other formats

## Key Methods and Attributes

### StructBaseProcessor

**Key Methods:**
- `load_structure(pdb_id)`: Load a structure by PDB ID
- `load_structures(pdb_ids)`: Load multiple structures
- `load_dataset(dataset_id)`: Load a predefined dataset
- `get_chains(pdb_id)`: Get chain IDs for a structure
- `get_sequence(pdb_id, chain_id)`: Get sequence for a specific chain
- `filter_by_ids(pdb_ids)`: Filter loaded structures by PDB ID
- `save_data(filename)`: Save processor data
- `create_dataset(dataset_id, name, description, content)`: Create a new dataset

**Key Attributes:**
- `data`: DataFrame containing structure data
- `chain_dict`: Dictionary mapping chain IDs to sequences
- `pdb_ids`: List of loaded PDB IDs
- `structure_filenames`: Paths to structure files

### CifHandler

**Key Methods:**
- `read(file_path)`: Read a mmCIF file
- `write(data, file_path)`: Write data to mmCIF format
- `validate_data(data)`: Ensure data meets mmCIF format requirements
- `merge_structures(structure1, structure2)`: Combine structures
- `extract_bioassembly(structure, assembly_id)`: Extract biological assembly

## Directory Structure

The CIF processor uses a standardized directory structure:

```
<data_root>/
└── structure/
    ├── mmcif/                # Raw mmCIF files
    ├── structure_dataset/    # Dataset definitions
    │   └── datasets.json     # Registry of available datasets
    ├── alignments/           # Structure alignment results
    ├── temp_cif/             # Temporary files
    └── registry.json         # Registry of all structure data
```

## Working with Datasets

The dataset system allows working with collections of structures:

```python
# Load a standard dataset
processor = StructBaseProcessor()
processor.load_dataset("gpcr_structures")

# Create a custom dataset
processor.create_dataset(
    dataset_id="my_structures",
    name="My Structure Collection",
    description="A custom set of structures",
    content=["1abc", "2xyz", "3pqr"]
)

# Save and share dataset
processor.save_dataset("my_structures")
```

## Path Configuration

The path system can be configured in several ways:

1. **Environment Variables**:
   ```
   PROTOS_DATA=/path/to/data
   PROTOS_REFERENCE_DATA=/path/to/reference/data
   ```

2. **Programmatic Configuration**:
   ```python
   from protos.io.paths.path_config import configure_paths
   
   configure_paths(
       data_dir="/path/to/data",
       reference_data_dir="/path/to/reference/data"
   )
   ```

## Integration with Other Components

The CIF processor system integrates with other Protos components:

- **Sequence Processing**: Extract and process protein sequences
- **GRN System**: Annotate structures with Generic Residue Numbers
- **Property Processing**: Associate properties with structures
- **Embedding**: Generate embeddings from structures
- **Visualization**: Render structures for analysis

## Example Workflows

### Basic Structure Loading and Analysis

```python
from protos.processing.structure.struct_processor import StructProcessor

# Initialize processor
processor = StructProcessor()

# Load structures
processor.load_structures(["1abc", "2xyz"])

# Extract information
for pdb_id in processor.pdb_ids:
    chains = processor.get_chains(pdb_id)
    for chain in chains:
        sequence = processor.get_sequence(pdb_id, chain)
        print(f"{pdb_id} {chain}: {sequence[:10]}...")
```

### Structure Filtering and Modification

```python
# Filter to keep only specified chains
processor.filter_chains({"1abc": ["A", "B"]})

# Extract binding pocket
binding_site = processor.extract_binding_site("1abc", chain_id="A", center_residue="50", radius=10.0)

# Save modified structure
processor.save_structure("1abc", "/path/to/output.cif")
```

### Working with Datasets

```python
# Create dataset from structures with resolution < 2.5Å
high_res = processor.filter_by_resolution(max_resolution=2.5)
processor.create_dataset("high_resolution", "High Resolution Structures", 
                         "Structures with resolution < 2.5Å", high_res)

# Load and use dataset
processor.load_dataset("high_resolution")
```

## Best Practices

1. **Initialization**: Always initialize processors with a name for better tracking
2. **Path Management**: Use the built-in path resolution rather than hardcoded paths
3. **Data Validation**: Validate structure data before saving
4. **Dataset Usage**: Use datasets for organizing collections of structures
5. **Chain Handling**: Be explicit about chain IDs when working with multi-chain proteins
6. **Memory Management**: For large datasets, use batch processing or filtering

## Troubleshooting

Common issues and solutions:

1. **File Not Found**: Ensure proper path configuration
2. **Invalid CIF Format**: Validate data before writing
3. **Memory Issues**: Use filtering to reduce data size
4. **Path Resolution**: Check environment variables and path configuration
5. **Registry Errors**: Ensure registry files are not corrupted

## Temporary File Handling

The CIF processor system includes robust support for temporary file operations, which is essential for external tool integration, structure alignments, and other operations that require intermediate file storage.

### TempFileManager

The core temporary file functionality is provided by the `TempFileManager` class in `protos.io.formats`:

```python
from protos.io.formats import temp_manager

# Create a temporary file
with temp_manager.temp_file(suffix=".cif") as tmp_file:
    # Write structure data to temporary file
    processor.save_structure("1abc", tmp_file)
    
    # Use the temporary file (e.g., for alignment)
    run_alignment_tool(tmp_file, other_structure)
    
    # File is automatically deleted when the context is exited
```

### Dedicated Temporary CIF Directory

The path system includes a standard location for temporary CIF files:

```
<data_root>/structure/temp_cif/
```

This directory is used for temporary operations that need to persist between function calls.

### Save Temporary CIF Files

While there isn't a built-in `save_temp_cif` method in the core `StructBaseProcessor`, the recommended pattern for saving temporary CIF files is:

```python
def save_temp_cif(processor, pdb_id: str, suffix: str = None, output_dir: Path = None) -> Path:
    """
    Save a temporary CIF file for a given PDB ID.
    
    Args:
        processor: The CifBaseProcessor instance
        pdb_id (str): The PDB ID to generate the CIF file for
        suffix (Optional[str]): An optional suffix for the filename
        output_dir (Optional[Path]): Directory to save the CIF file; defaults to temp_dir if None

    Returns:
        Path: The path to the saved CIF file
    """
    from protos.io.cif_handler import CifHandler
    import time
    from pathlib import Path
    
    # Check if PDB ID exists in loaded data
    if pdb_id not in processor.pdb_ids:
        raise ValueError(f"PDB ID {pdb_id} not found in the data")
    
    # Resolve output directory
    if output_dir is None:
        if not hasattr(processor, "temp_dir"):
            temp_dir = Path(processor.data_root) / "structure" / "temp_cif"
            processor.temp_dir = temp_dir
        output_dir = processor.temp_dir
    
    # Ensure output directory exists
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate filename with timestamp for uniqueness
    timestamp = int(time.time())
    suffix_str = f"_{suffix}" if suffix else ""
    filename = f"{pdb_id}{suffix_str}_{timestamp}.cif"
    temp_path = output_dir / filename
    
    # Get structure data for this PDB ID
    structure_data = processor.data[processor.data['pdb_id'] == pdb_id].copy()
    
    # Use CifHandler to write the structure data
    cif_handler = CifHandler()
    cif_handler.write(structure_data, str(temp_path))
    
    return temp_path

# Usage example:
# Save a structure to a temporary file
temp_path = save_temp_cif(processor, "1abc")

# Use the temporary file
run_external_tool(temp_path)
```

This function can be used for external tools like structure alignment programs that need temporary files as input. In the example below, we show how to use this function with GTalign for structure alignment:

```python
# Initialize GTalign
from protos.processing.structure.gtalign import GTalign
gtalign = GTalign()

# Save reference structure to a temporary file
ref_temp_path = save_temp_cif(processor, reference_pdb_id, suffix="for_alignment")

# Save target structure to a temporary file
target_temp_path = save_temp_cif(processor, target_pdb_id, suffix="for_alignment")

# Run alignment
result = gtalign.align(
    query_structures=str(ref_temp_path),
    reference_structures=[str(target_temp_path)],
    output_dir=alignment_output_dir,
    tm_score_threshold=0.5,
    speed=9,
    depth=2,
    verbose=True
)
```

You can add this utility function to your own processors or extend the `StructBaseProcessor` class to include it as a method. This pattern is used in various alignment utilities throughout the codebase.

### ExternalToolHelper

For integration with external structural biology tools, the `ExternalToolHelper` class provides dedicated methods:

```python
from protos.io.formats import ExternalToolHelper

# Create a helper for a specific tool
helper = ExternalToolHelper("alignment_tool")

# Setup input and output files
input_file = helper.tool_input("input.cif", processor.get_structure("1abc"))
output_file = helper.tool_output("result.txt")

# Run the tool
helper.run_command(f"alignment_tool -i {input_file} -o {output_file}")

# Access the results
with open(output_file, 'r') as f:
    results = f.read()
    
# Cleanup is handled automatically when helper goes out of scope
```

### Best Practices for Temporary Files

1. **Use Context Managers**: Whenever possible, use context managers to ensure cleanup
2. **Include Timestamps**: When creating persistent temporary files, include timestamps in filenames
3. **Cleanup Regularly**: Implement periodic cleanup of older temporary files
4. **Track Created Files**: Keep a record of created temporary files for later cleanup
5. **Use Standard Directories**: Store temporary files in the standard `temp_cif` directory
6. **Validate Before Use**: Always validate temporary files before using them with external tools

## Advanced Topics

- **Custom Data Processors**: Extending the StructProcessor
- **Cross-Processor Integration**: Combining structure data with other processors
- **Performance Optimization**: Strategies for handling large datasets
- **Validation Pipelines**: Ensuring data quality across processing steps