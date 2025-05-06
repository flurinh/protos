# Protos Overview

<p align="center">
  <img src="logo.png" alt="Protos Logo">
</p>


## What is Protos?

Protos is a Python library designed to **standardize and execute complex computational workflows** essential for structural biology research. It provides integrated capabilities for handling diverse biological data types – including sequences, 3D structures, alignments, and associated properties – through a unified framework.

The core function of Protos is to manage multi-step analyses by breaking them down into defined tasks handled by modular components.

## Core Architecture: Processors & Interoperability

Protos utilizes a modular architecture built upon distinct Python components called **'Processors'**. Each Processor is specialized for a specific domain, such as:

*   **`CifBaseProcessor`**: Manages 3D structure data
*   **`SeqProcessor`**: Handles sequence data and alignments
*   **`GRNBaseProcessor`**: Implements Generic Residue Numbering systems
*   **`LigProcessor`**: Deals with ligand information and interactions
*   **`EMBProcessor`**: Manages protein embeddings
*   **`PropertyProcessor`**: Integrates metadata and calculated properties

A key feature is the **interoperability** between these Processors. Outputs from one (e.g., selected residues from `CifBaseProcessor`) can directly serve as inputs for another (e.g., for GRN mapping by `GRNBaseProcessor`, followed by sequence analysis by `SeqProcessor`), enabling flexible construction of sophisticated analysis pipelines.

The relationships and primary data flow between these core processors are visualized below:

![Protos Processor Overview Diagram](resources/protos_overview.png)
*(Diagram showing connections between CP, SP, GRNP, LP, EMBP, and PP)*

---

## Working with Protos

Protos was developed to solve real-world computational challenges in structural biology by providing a robust framework for building complex workflows. The following examples demonstrate how Protos is currently used in production workflows.

### Important Note on Path Management

Protos uses a specialized path management system that requires proper initialization to avoid common issues. Production code must explicitly set data paths rather than rely on defaults:

```python
# Explicit path configuration is essential to avoid path resolution bugs
import os
from pathlib import Path
from protos.io.paths.path_config import ProtosPaths

# Set environment variables for consistent path resolution
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())

# Initialize paths with explicit directories
paths = ProtosPaths(
    user_data_root=str(data_dir.absolute()),
    ref_data_root=str(data_dir.absolute()),
    create_dirs=True
)

print(f"Path resolver user_data_root: {paths.user_data_root}")
print(f"Path resolver ref_data_root: {paths.ref_data_root}")
```

This approach ensures consistent behavior across different environments and prevents the path resolution bug documented in our technical documentation.

## Current Production Use Cases

### Example 1: Structure Data Processing and Filtering

This example demonstrates how to load and filter structure data from standardized datasets, which is a common first step in structural biology workflows.

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.io.paths.path_config import ProtosPaths
import os
from pathlib import Path

# Define paths explicitly
data_dir = Path("/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())

# Initialize processor with explicit paths
processor = CifBaseProcessor(
    name="structure_processor",
    data_root=str(data_dir.absolute()),
    processor_data_dir="structure"
)

# Ensure critical directories exist
mmcif_dir = data_dir / "structure" / "mmcif"
dataset_dir = data_dir / "structure" / "structure_dataset"
os.makedirs(mmcif_dir, exist_ok=True)
os.makedirs(dataset_dir, exist_ok=True)

# Load a dataset with proper data type enforcement
processor.load_dataset("my_dataset", apply_dtypes=True)
print(f"Loaded dataset with {len(processor.pdb_ids)} structures")

# Filter structures by chain ID
filtered_structures = {}
for pdb_id in processor.pdb_ids:
    df = processor.data[processor.data['pdb_id'] == pdb_id].copy()
    df_chain = df[df['auth_chain_id'] == 'A'].copy()
    
    if not df_chain.empty:
        # Ensure coordinates are numeric (critical for calculations)
        for coord in ['x', 'y', 'z']:
            df_chain[coord] = pd.to_numeric(df_chain[coord], errors='coerce')
            
        filtered_structures[pdb_id] = {'df': df_chain}

print(f"Filtered to {len(filtered_structures)} structures with chain A")
```

This example handles real-world considerations like ensuring proper data types for coordinates, which is critical for numerical operations. It also demonstrates the explicit path configuration required for reliable operation.

### Example 2: Structure Comparison with GRN Annotation

This example shows how to align structures, assign Generic Residue Numbers (GRNs), and perform comparative analysis - a common workflow in analyzing protein families.

```python
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.grn.grn_base_processor import GRNBaseProcessor
import os
import pandas as pd
import numpy as np
from pathlib import Path

# Set up paths properly
data_dir = Path("/path/to/data")
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())

# Create processors with explicit paths
cp = CifBaseProcessor(
    name="structure_processor",
    data_root=str(data_dir.absolute()),
    processor_data_dir="structure"
)

grnp = GRNBaseProcessor(
    name="grn_processor",
    data_root=str(data_dir.absolute()),
    processor_data_dir="grn"
)

# Override key paths to ensure consistency
cp.path_structure_dir = os.path.join(str(data_dir.absolute()), "structure", "mmcif")
cp.path_dataset_dir = os.path.join(str(data_dir.absolute()), "structure", "structure_dataset")
cp.path_alignment_dir = os.path.join(str(data_dir.absolute()), "structure", "alignments")

# Load reference GRN table
grnp.load_grn_table("reference_grn_table")

# Load structures
pdb_ids = ['6CMO', '7ZOU', '5XEZ']
cp.load_structures(pdb_ids, apply_dtypes=True)
print(f"Loaded {len(cp.pdb_ids)} structures")

# Align structures
alignment_info = cp.perform_multiple_alignment(pdb_ids, reference_id='6CMO')
print("Structures aligned to reference")

# Process each structure
results = {}
for pdb_id in pdb_ids:
    # Extract sequence from structure
    chain_a_df = cp.data[(cp.data['pdb_id'] == pdb_id) & (cp.data['auth_chain_id'] == 'A')]
    sequence = cp.extract_sequence_from_dataframe(chain_a_df)
    
    # Assign GRNs to the sequence
    grn_assignment = grnp.assign_grns(
        sequence=sequence,
        protein_id=f"{pdb_id}_A",
        reference_id="reference_protein",
        use_cached=True
    )
    
    # Identify ligand binding pocket
    pocket_residues = cp.extract_binding_pocket(
        pdb_id=pdb_id,
        chain_id='A',
        ligand_name='RET',
        distance_cutoff=4.0
    )
    
    # Map pocket residues to GRNs
    pocket_grns = grnp.get_grns_from_residues(
        residue_ids=pocket_residues,
        protein_id=f"{pdb_id}_A"
    )
    
    # Calculate distances using aligned coordinates
    reference_position = "3.50"
    ref_position_coord = cp.get_ca_coordinate_by_grn(
        pdb_id=pdb_id,
        chain_id='A',
        grn=reference_position,
        alignment_info=alignment_info
    )
    
    pocket_coords = cp.get_ca_coordinates_by_grns(
        pdb_id=pdb_id,
        chain_id='A',
        grns=pocket_grns,
        alignment_info=alignment_info
    )
    
    # Calculate average distance
    distances = []
    for grn, coord in pocket_coords.items():
        dist = np.linalg.norm(ref_position_coord - coord)
        distances.append(dist)
    
    avg_distance = np.mean(distances)
    results[pdb_id] = {
        'avg_distance': avg_distance,
        'pocket_size': len(pocket_grns),
        'pocket_grns': pocket_grns
    }

# Rank and output results
ranked_results = sorted(results.items(), key=lambda x: x[1]['avg_distance'])
for pdb_id, result in ranked_results:
    print(f"{pdb_id}: Avg distance = {result['avg_distance']:.2f}Å, Pocket size = {result['pocket_size']}")
```

This example demonstrates important production practices:
1. Explicit path configuration for reliability
2. Type validation for coordinates to prevent calculation errors
3. Proper use of GRN assignment with caching for performance
4. Structure alignment with reference selection
5. Complex multi-step analysis workflow

### Example 3: Opsin Analysis Workflow

This example shows a real production workflow used for analyzing opsin protein structures, demonstrating how Protos is used in scientific research:

```python
"""
Opsin Analysis Workflow Example - Based on actual production code
This workflow demonstrates:
1. Loading and processing structure data from multiple sources
2. Applying proper path configuration for reliability
3. Implementing multi-stage caching for performance
4. Structure filtering and property annotation
5. Integration with advanced analysis functions
"""
import os
import numpy as np
import pandas as pd
from pathlib import Path
import pickle

# Import Protos components
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.io.paths.path_config import ProtosPaths

# Set up paths with proper configuration
data_dir = Path("/path/to/data")
output_dir = Path("/path/to/results")
cache_dir = output_dir / "cache"
os.makedirs(cache_dir, exist_ok=True)

# Configure environment variables - CRITICAL for reliable path resolution
os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())

# Define cache file paths
raw_cache_file = cache_dir / "raw_structures.pkl"
processed_cache_file = cache_dir / f"processed_structures_A.pkl"

# STEP 1: Load structure data (with caching)
if os.path.exists(raw_cache_file):
    # Load from cache if available
    try:
        print(f"Loading raw structure data from cache: {raw_cache_file}")
        with open(raw_cache_file, 'rb') as f:
            raw_data = pickle.load(f)
            cp_mo_exp = raw_data.get('cp_mo_exp')
            cp_mo_pred = raw_data.get('cp_mo_pred')
            datasets = raw_data.get('datasets', {})
    except Exception as e:
        print(f"Cache loading failed: {str(e)}")
        raw_data = None
else:
    raw_data = None

# Load from scratch if cache failed or doesn't exist
if raw_data is None:
    print("Creating new processor instances and loading datasets")
    
    # Initialize paths handler with explicit paths
    paths = ProtosPaths(
        user_data_root=str(data_dir.absolute()),
        ref_data_root=str(data_dir.absolute()),
        create_dirs=True
    )
    
    # Initialize processors with explicit paths
    cp_mo_exp = CifBaseProcessor(
        name="mo_exp_processor",
        data_root=str(data_dir.absolute()),
        processor_data_dir="structure"
    )
    
    cp_mo_pred = CifBaseProcessor(
        name="mo_pred_processor",
        data_root=str(data_dir.absolute()),
        processor_data_dir="structure"
    )
    
    # Ensure critical directories exist
    mmcif_dir = data_dir / "structure" / "mmcif"
    dataset_dir = data_dir / "structure" / "structure_dataset"
    os.makedirs(mmcif_dir, exist_ok=True)
    os.makedirs(dataset_dir, exist_ok=True)
    
    # Override key paths to ensure consistency
    for processor in [cp_mo_exp, cp_mo_pred]:
        processor.path_structure_dir = os.path.join(str(data_dir.absolute()), "structure", "mmcif")
        processor.path_dataset_dir = os.path.join(str(data_dir.absolute()), "structure", "structure_dataset")
        processor.path_alignment_dir = os.path.join(str(data_dir.absolute()), "structure", "alignments")
    
    # Load datasets
    try:
        cp_mo_exp.load_dataset("mo_exp", apply_dtypes=True)
        print(f"Loaded experimental dataset with {len(cp_mo_exp.pdb_ids)} structures")
    except Exception as e:
        print(f"Failed to load experimental dataset: {str(e)}")
    
    try:
        cp_mo_pred.load_dataset("mo_pred", apply_dtypes=True)
        print(f"Loaded predicted dataset with {len(cp_mo_pred.pdb_ids)} structures")
    except Exception as e:
        print(f"Failed to load predicted dataset: {str(e)}")
    
    # Ensure retinal ligand has consistent naming ('RET' vs 'LIG')
    if hasattr(cp_mo_pred, 'data') and cp_mo_pred.data is not None:
        cp_mo_pred.data.loc[cp_mo_pred.data['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
        print("Renamed 'LIG' to 'RET' in predicted structures")
    
    # Cache raw data
    raw_data_to_save = {
        'cp_mo_exp': cp_mo_exp,
        'cp_mo_pred': cp_mo_pred,
        'datasets': {
            "mo_exp": {"id": "mo_exp", "pdb_ids": cp_mo_exp.pdb_ids if hasattr(cp_mo_exp, 'pdb_ids') else []},
            "mo_pred": {"id": "mo_pred", "pdb_ids": cp_mo_pred.pdb_ids if hasattr(cp_mo_pred, 'pdb_ids') else []}
        }
    }
    
    with open(raw_cache_file, 'wb') as f:
        pickle.dump(raw_data_to_save, f)
    print(f"Saved raw data to cache: {raw_cache_file}")

# STEP 2: Filter structures by chain and process
processed_structures = {}
chain_id = 'A'

# Define helper function to filter by chain and retinal
def filter_structures_by_chain_and_retinal(processor, chain='A', retinal_name='RET', cutoff=6.0):
    filtered_structures = {}
    
    for pdb_id in processor.pdb_ids:
        df_pdb = processor.data[processor.data['pdb_id'] == pdb_id]
        df_chain = df_pdb[df_pdb['auth_chain_id'] == chain].copy()
        
        if df_chain.empty:
            continue
        
        # Select retinal atoms
        df_ret = df_pdb[(df_pdb['res_name3l'] == retinal_name) & 
                        (df_pdb['auth_chain_id'] == chain)].copy()
        
        # Ensure coordinates are numeric
        for df in [df_chain, df_ret]:
            for coord in ['x', 'y', 'z']:
                if coord in df.columns:
                    df[coord] = pd.to_numeric(df[coord], errors='coerce')
        
        filtered_structures[pdb_id] = {
            'df': df_chain,
            'df_ret': df_ret,
            'chain_id': chain
        }
        
        # Extract CA atoms
        df_ca = df_chain[df_chain['atom_name'] == 'CA'].copy()
        filtered_structures[pdb_id]['df_ca'] = df_ca
    
    return filtered_structures

# Process both experimental and predicted structures
exp_structures = filter_structures_by_chain_and_retinal(cp_mo_exp, chain=chain_id)
pred_structures = filter_structures_by_chain_and_retinal(cp_mo_pred, chain=chain_id)

# Combine all structures
processed_structures.update(exp_structures)
processed_structures.update(pred_structures)

print(f"Processed {len(exp_structures)} experimental and {len(pred_structures)} predicted structures")

# STEP 3: Create structure mapping (experimental to predicted)
structure_mapping = {}
for exp_id in exp_structures.keys():
    # Create base name without suffixes
    base_name = exp_id.split('_')[0] if '_' in exp_id else exp_id
    
    # Look for matching predicted structure
    for pred_id in pred_structures.keys():
        if base_name in pred_id:
            structure_mapping[exp_id] = pred_id
            break

print(f"Created {len(structure_mapping)} experimental-predicted structure pairs")

# STEP 4: Save processed results
processed_data = {
    'processed_structures': processed_structures,
    'structure_mapping': structure_mapping
}

with open(processed_cache_file, 'wb') as f:
    pickle.dump(processed_data, f)
print(f"Saved {len(processed_structures)} processed structures to cache")

# Example usage of processed data:
# Calculate structure statistics
stats = {}
for pdb_id, data in processed_structures.items():
    if 'df_ca' in data and not data['df_ca'].empty:
        df_ca = data['df_ca']
        
        # Calculate center of mass
        com = df_ca[['x', 'y', 'z']].mean().values
        
        # Calculate radius of gyration
        distances = np.sqrt(np.sum((df_ca[['x', 'y', 'z']].values - com)**2, axis=1))
        rg = np.mean(distances)
        
        stats[pdb_id] = {
            'residue_count': len(df_ca),
            'center_of_mass': com,
            'radius_of_gyration': rg
        }

print(f"\nAnalyzed {len(stats)} structures")
for pdb_id, stat in list(stats.items())[:3]:  # Show top 3 examples
    print(f"{pdb_id}: {stat['residue_count']} residues, Rg = {stat['radius_of_gyration']:.2f}Å")

# Final results
results = {
    'processed_structures': processed_structures,
    'structure_mapping': structure_mapping,
    'statistics': stats
}
```

This example demonstrates a real-world scientific workflow using Protos:
1. Multi-stage caching strategy for efficient data handling
2. Proper path configuration with environment variables and explicit paths
3. Structure filtering with chain and ligand selection
4. Data transformation and type enforcement for calculations
5. Creating and managing mappings between related structures
6. Calculating structural properties from processed data

## Future Direction: Protos-MCP

While the current Protos library provides a powerful toolkit for structural biology, it requires programming expertise to use effectively. We envision a future where Protos can be accessed through natural language via integration with Large Language Models (LLMs).

The **Protos-MCP** project aims to expose Protos capabilities as discrete tools that LLMs can utilize through the Model Context Protocol (MCP). This would allow researchers to interact with complex bioinformatics workflows using natural language, making the power of Protos accessible to a broader audience.

![Protos-MCP Interaction Flow](resources/MCP_integration.png)
*(Diagram showing User -> LLM -> MCP Server -> Protos Library -> MCP Server -> LLM -> User)*

### Required Improvements for MCP Integration

Before Protos can be effectively integrated with MCP, several key improvements are needed:

#### 1. Path Management Refinement
The current path handling system has a critical bug related to `DataSource.AUTO` that must be fixed:
```python
# Current buggy implementation
def _resolve_data_root(self, source: DataSource) -> str:
    if source == DataSource.USER:
        return self.user_data_root
    elif source == DataSource.REFERENCE:
        return self.ref_data_root
    else:  # AUTO
        # Default to reference data for reading, but this should be refined
        return self.ref_data_root  # <-- BUG: Should use user_data_root in many cases
```

The path management needs to be simplified for MCP use with:
- Automatic initialization from a single project directory
- Consistent default behavior that prioritizes user data
- Hidden implementation details that don't leak into API

#### 2. Creating a Stateless API Layer

MCP tools require stateless functions that don't rely on maintaining processor state:

```python
# Current approach (stateful)
cp = CifBaseProcessor()
cp.load_structure("1abc")
sequence = cp.get_sequence("1abc", "A")

# Required approach for MCP (stateless)
@mcp_tool
def get_structure_sequence(pdb_id: str, chain_id: str) -> str:
    """Get the amino acid sequence from a structure."""
    cp = CifBaseProcessor()
    cp.load_structure(pdb_id) 
    return cp.get_sequence(pdb_id, chain_id)
```

#### 3. Standardizing Error Handling and Input Validation

MCP requires predictable error handling and validation:
- Consistent exception types with informative messages
- Input parameter validation for all external functions
- Graceful failure modes instead of crashing
- Error messages designed for LLM comprehension

#### 4. Streamlining Common Workflows

The current multi-stage caching pattern should be automated:
```python
# Common pattern in current code (too verbose for MCP)
if cache_file.exists():
    try:
        with open(cache_file, 'rb') as f:
            data = pickle.load(f)
    except Exception as e:
        data = None

# Desired pattern for MCP
@mcp_tool
def align_structures(pdb_ids: List[str], use_cache: bool = True) -> Dict:
    """Align multiple protein structures with automatic caching."""
    return alignment_manager.perform_alignment(pdb_ids, use_cache=use_cache)
```

These improvements would make Protos suitable for MCP integration by providing a cleaner, more reliable interface that LLMs can confidently use through natural language interaction.

## Installation

### Prerequisites

- **Python**: Protos requires Python 3.9 or later.
- **Package Manager**: pip (or optionally conda) for installing Python packages.

### Standard Installation (Development)

Since Protos is currently under active development, we recommend installing from the source repository:

```bash
git clone https://github.com/your-organization/protos.git
cd protos
pip install -e .
```

### Environment Variables

Protos relies on environment variables for path configuration:

```bash
# Essential configuration for reliable path resolution
export PROTOS_DATA_ROOT=/path/to/your/data
export PROTOS_REF_DATA_ROOT=/path/to/your/data
```

### External Dependencies

Protos integrates with several external bioinformatics tools:

- **MMseqs2**: For sequence searching and clustering
- **GTalign**: For GPU-accelerated protein structure alignment

These tools need to be installed separately and made available in your system PATH.
