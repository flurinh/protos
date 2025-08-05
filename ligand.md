# Protos Ligand Processor: Capabilities and Usage Guide

## Overview

The Protos Ligand Processor provides comprehensive functionality for handling small molecule data, including structure-based drug design, bioactivity analysis, and protein-ligand interaction profiling. It integrates seamlessly with the Structure Processor to enable sophisticated analyses of protein-ligand complexes.

## Core Capabilities

### 1. Ligand Data Management

The LigandProcessor follows Protos' entity-registry system with human-readable naming:

- **Entity Storage**: Ligands are stored as SDF files with SMILES-based filenames
- **Multiple Identifiers**: Supports SMILES, ChEMBL IDs, InChI Keys as aliases
- **Automatic Registration**: Entities are tracked in the central registry
- **Dataset Management**: Create collections of ligands for batch analysis

### 2. Molecular Property Calculation

When RDKit is available, the processor calculates:

- **Basic Properties**: MW, LogP, HBA, HBD, TPSA, rotatable bonds
- **Lipinski Rule of Five**: Drug-likeness assessment
- **Aromatic rings, heavy atoms, formal charge**
- **Custom property filters** for lead optimization

### 3. Structure-Ligand Analysis (Key Feature)

The most powerful capabilities come from integration with the Structure Processor:

#### Ligand Extraction from Structures
```python
from protos.analysis.structure_ligand_analysis import extract_all_ligands

# Extract all non-water/ion ligands from a structure
ligands = extract_all_ligands(struct_processor, "1ATP", exclude_common=True)
# Returns: List of ligand dictionaries with coordinates and metadata
```

#### Binding Site Analysis
```python
from protos.analysis.structure_ligand_analysis import get_binding_site

# Identify residues within 5Å of ligand
binding_site = get_binding_site(struct_processor, "1ATP", ligand_atoms, cutoff=5.0)
# Returns: DataFrame of binding residues sorted by distance
```

#### Detailed Interaction Profiling
```python
from protos.analysis.structure_ligand_analysis import calculate_ligand_interactions

# Get comprehensive interaction analysis
interactions = calculate_ligand_interactions(struct_processor, "1ATP", ligand_atoms, detailed=True)
# Returns: 
# - Hydrogen bonds (donor/acceptor pairs with distances)
# - Hydrophobic contacts
# - Water-mediated interactions
# - π-π stacking and salt bridges
# - Summary statistics
```

#### Comparative Analysis
```python
from protos.analysis.structure_ligand_analysis import compare_ligand_binding_sites

# Compare binding sites across multiple structures
similarity_matrix = compare_ligand_binding_sites(
    struct_processor,
    [("1ATP", "ATP", "A"), ("2ATP", "ADP", "A")],
    cutoff=5.0
)
# Returns: Jaccard similarity matrix of binding sites
```

### 4. ChEMBL Integration

With `chembl_webresource_client` installed:

- **Download bioactivity data** for target proteins
- **Map between identifier systems** (UniProt, PDB, gene names)
- **Built-in gene name aliases** for common drug targets (EGFR, HER2, COX2, etc.)
- **Smart target selection** - prioritizes SINGLE PROTEIN over chimeric proteins
- **Filter by activity type** (IC50, Ki, Kd) and potency
- **Automatic SDF file generation** from ChEMBL compounds
- **Caching system** for protein mappings and activity data

### 5. Similarity Search

With RDKit available:

- **Molecular fingerprint generation**
- **Tanimoto similarity calculation**
- **Search within datasets or globally**
- **Configurable similarity thresholds**

## Utilities Included

### ligand_utils.py Functions

1. **SMILES Handling**
   - `sanitize_smiles_filename()`: Convert SMILES to safe filenames
   - `validate_smiles()`: Check validity and get canonical form
   - `smiles_to_inchi()`: Generate InChI identifiers

2. **Property Calculation**
   - `calculate_molecular_properties()`: Full property profile
   - `is_drug_like()`: Lipinski and custom filters
   - `parse_activity_value()`: Convert activity units to nM

3. **Protein Mapping**
   - `extract_protein_mapping()`: Identify protein ID types
   - Supports: UniProt, PDB, gene names, ChEMBL targets
   - Pre-configured aliases for common drug targets:
     * Kinases: EGFR, HER2, ABL, SRC, VEGFR2, MET, ALK, BTK, JAK2, BRAF
     * Cell cycle: CDK4, CDK6
     * PI3K/mTOR pathway: PI3K, PIK3CA, MTOR
     * Epigenetics: HDAC
     * Immuno-oncology: PD1, PDL1, CTLA4
     * Inflammation: COX2 (PTGS2)

4. **File Generation**
   - `create_sdf_from_smiles()`: Generate SDF with properties

### structure_ligand_analysis.py Functions

1. **Ligand Extraction**
   - `extract_all_ligands()`: Get all ligands from structure
   - `get_ligand_by_id()`: Extract specific ligand
   - Automatic exclusion of water, ions, crystallization agents

2. **Binding Site Analysis**
   - `get_binding_site()`: Find interacting residues
   - `estimate_binding_site_volume()`: Calculate pocket volume
   - `identify_key_residues()`: Find critical interactions

3. **Interaction Detection**
   - Uses `LigandInteractionAnalyzer` for detailed analysis
   - Hydrogen bonds with geometry validation
   - Hydrophobic contact identification
   - Water bridge detection
   - Distance-based interaction profiling

4. **Comparative Analysis**
   - `compare_ligand_binding_sites()`: Cross-structure comparison
   - `find_conserved_interactions()`: Identify key interactions
   - Jaccard similarity for binding site comparison

5. **Export Functions**
   - `export_ligand_sdf()`: Save ligand with 3D coordinates
   - `create_ligand_interaction_report()`: Comprehensive analysis report

## Recent Updates (August 2025)

### Improvements to ChEMBL Integration
1. **Fixed ChEMBL client initialization** - Now uses proper lazy loading
2. **Enhanced gene name mapping**:
   - Built-in aliases for 20+ common drug targets
   - Fallback to partial matching (icontains) when exact match fails
   - Smart target selection prioritizes SINGLE PROTEIN over chimeric proteins
3. **Fixed entity registry integration**:
   - Uses public API methods only
   - Proper EntityInfo attribute handling
   - Efficient format-filtered entity listing

### Bug Fixes
1. **Path handling** - Fixed TypeError in _load_ligand_data
2. **Similarity search** - Fixed EntityInfo.name → EntityInfo.original_id
3. **List entities** - Now uses proper EntityRegistry API with format filtering
4. **SDF file handling** - Fixed "Bad input file" errors by:
   - Checking if SDF files exist before attempting to read
   - Gracefully handling missing RDKit by falling back to metadata
   - Preventing SDF creation errors when RDKit is not available

## Current Limitations & Missing Functionality

### 1. Coordinate-to-Structure Conversion
- **Issue**: Cannot generate proper bond orders from PDB coordinates
- **Impact**: Ligands extracted from structures lack chemical connectivity
- **Workaround**: Use template matching with known ligand databases

### 2. PDB Chemical Component Dictionary Integration
- **Missing**: Direct access to PDB CCD for ligand templates
- **Impact**: Limited ability to handle novel ligands
- **Needed**: API integration or local CCD database

### 3. Conformational Analysis
- **Missing**: 3D conformer generation and analysis
- **Impact**: Cannot assess ligand flexibility or binding poses
- **Needed**: RDKit conformer generation integration

### 4. Docking Integration
- **Missing**: Interface to docking programs (AutoDock, Vina, etc.)
- **Impact**: Cannot predict binding poses for new ligands
- **Potential**: Add docking wrapper classes

### 5. Advanced Interaction Types
- **Partially Implemented**: Basic H-bonds and hydrophobic contacts
- **Missing**: 
  - Cation-π interactions
  - Halogen bonds
  - Metal coordination
  - Detailed π-π stacking geometry

### 6. QSAR/ML Features
- **Missing**: Interaction fingerprints for machine learning
- **Impact**: Cannot build predictive models from interaction data
- **Needed**: Standardized feature extraction

## Recommended Workflow

### 1. Structure-Based Ligand Analysis
```python
# 1. Load structure with ligands
struct_proc.load_structures(["1ATP"])

# 2. Extract and analyze all ligands
ligands = extract_all_ligands(struct_proc, "1ATP")
for ligand in ligands:
    # Get binding site
    binding = get_binding_site(struct_proc, "1ATP", ligand['atoms'])
    
    # Calculate interactions
    interactions = calculate_ligand_interactions(
        struct_proc, "1ATP", ligand['atoms']
    )
    
    # Store properties
    prop_proc.assign_property(
        f"1ATP_{ligand['res_name3l']}",
        'binding_residues',
        len(binding['residues'])
    )
```

### 2. ChEMBL Bioactivity Integration
```python
# 1. Download active compounds for target
compounds = lig_proc.get_protein_ligands(
    "EGFR", 
    min_pchembl=6.0,
    activity_types=["IC50", "Ki"]
)

# 2. Register compounds
for compound in compounds:
    lig_proc.save_entity(compound['smiles'], compound)

# 3. Create dataset for analysis
lig_proc.create_dataset("egfr_inhibitors", 
                       [c['smiles'] for c in compounds])
```

### 3. Comparative Binding Site Analysis
```python
# Compare ATP binding sites across kinases
ligand_list = [
    ("1ATP", "ATP", "A"),
    ("1BYQ", "ATP", "A"),
    ("1PKG", "ATP", "A")
]

similarity = compare_ligand_binding_sites(
    struct_proc, ligand_list, cutoff=5.0
)
```

## Future Development Priorities

1. **PDB CCD Integration**: Enable automatic ligand template retrieval
2. **Bond Order Assignment**: Implement heuristics for connectivity inference
3. **Conformer Analysis**: Add 3D conformer generation and clustering
4. **ML Features**: Standardized interaction fingerprints
5. **Extended Interactions**: Complete interaction type coverage
6. **Pharmacophore Analysis**: 3D pharmacophore generation
7. **Fragment Analysis**: Decompose ligands into chemical fragments

## Dependencies

### Required
- Protos core (BaseProcessor, EntityRegistry, DatasetManager)
- NumPy, Pandas, SciPy

### Optional but Recommended
- **RDKit**: Molecular property calculation, SMILES handling, similarity search
- **chembl_webresource_client**: ChEMBL database access
- **BioPython**: Enhanced structure manipulation

### Installation
```bash
# Basic installation
pip install protos

# With RDKit support
conda install -c conda-forge rdkit

# With ChEMBL support
pip install chembl_webresource_client
```

## Key Integration Points

The Ligand Processor is designed to work seamlessly with:

1. **Structure Processor**: Extract and analyze ligands from protein structures
2. **Property Processor**: Store bioactivity and calculated properties
3. **GRN Processor**: Map binding sites to generic residue numbers
4. **Sequence Processor**: Link ligand binding to sequence features

This integration enables sophisticated workflows like structure-based drug design, chemogenomics analysis, and binding site comparison across protein families.