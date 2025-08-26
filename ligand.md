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

### 5. Local Database Access (NEW)

The LigandProcessor now provides seamless access to three major ligand databases that are downloaded once and cached locally. All database paths are managed by ProtosPaths - no hardcoded paths!

#### PDB Chemical Component Dictionary (CCD)
Access all small molecules found in PDB structures:

```python
# Get specific ligand from CCD
atp = lig_proc.get_ccd_ligand('ATP')  # Auto-downloads CCD on first use
print(f"ATP SMILES: {atp['smiles']}")
print(f"ATP formula: {atp['formula']}")

# Create dataset from CCD components
cofactors = lig_proc.create_ccd_dataset(
    "cofactors",
    ["ATP", "NAD", "FAD", "COA", "HEM", "PLP"]
)
```

#### QM9 Quantum Chemistry Dataset
Search 134k molecules by quantum properties:

```python
# Find molecules with specific HOMO-LUMO gap
molecules = lig_proc.search_qm9_by_properties({
    'gap': (0.1, 0.3),      # eV
    'dipole': (0, 2.0),     # Debye
    'homo': (-7.0, -5.0)    # eV
}, limit=100)

# Get specific molecule with all quantum properties
mol = lig_proc.get_qm9_molecule(12345)
print(f"Quantum properties: {mol['quantum_properties']}")
```

#### Enamine REAL Database
Commercial compound library (requires Enamine subscription):

```python
# Search by similarity in Enamine subsets
similar = lig_proc.search_enamine_by_similarity(
    "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin
    dataset='diversity_1k',    # Small test dataset (default)
    similarity=0.7
)

# Available datasets:
# - diversity_1k: 1,000 diverse compounds (test)
# - hit2lead_1k: 1,000 lead-like compounds (test)
# - diversity_10k: 10,000 diverse compounds
# - fragments_5k: 5,000 fragments (MW < 300)
# - kinase_focused: 50,000 kinase-targeted compounds
# - gpcr_focused: 50,000 GPCR-targeted compounds
```

**Note**: Enamine is a commercial database. The loader includes placeholder URLs
that need to be replaced with actual download links from your Enamine subscription.
Contact Enamine (https://enamine.net) for access to their REAL database.

#### Database Statistics
Check which databases are available:

```python
stats = lig_proc.get_database_statistics()
for db_name, info in stats.items():
    print(f"{db_name}: {'Downloaded' if info['downloaded'] else 'Not downloaded'}")
    if info['downloaded']:
        print(f"  Path: {info['path']}")
```

### 6. Similarity Search

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
5. **PDB to ChEMBL mapping** - Fixed PDB code mapping by:
   - Adding PDB_ALIASES dictionary for common PDB structures
   - Updating API calls to avoid deprecated paths
   - Example: 1M17 now correctly maps to CHEMBL203 (EGFR)

### New Features (August 5, 2025)

1. **SDF File Operations** - Full support for Structure Data Format files:
   ```python
   # Load SDF file
   molecules = lig_proc.load_sdf_file('drug_library', as_entities=True)
   
   # Save molecules to SDF
   sdf_path = lig_proc.save_sdf_file('egfr_inhibitors', molecules)
   
   # Convert DataFrame to SDF
   lig_proc.save_sdf_file('compounds', df)
   ```

2. **Structure Format Conversion** - Convert ligands to 3D structures:
   ```python
   # Convert to CIF format (DEFAULT) - fully compatible with StructureProcessor
   cif_path = lig_proc.convert_to_structure_format(smiles)  # CIF is default
   
   # Specify chain and residue name for CIF
   cif_path = lig_proc.convert_to_structure_format(smiles, 
                                                   chain_id='L', 
                                                   res_name='ATP')
   
   # Also supports PDB format (if specifically needed)
   pdb_path = lig_proc.convert_to_structure_format(smiles, output_format='pdb')
   
   # And MOL2 format
   mol2_path = lig_proc.convert_to_structure_format(smiles, output_format='mol2')
   ```

3. **CIF DataFrame Integration** - Direct integration with StructureProcessor:
   ```python
   # Convert ligand to CIF DataFrame format
   ligand_df = lig_proc.convert_to_cif_dataframe(smiles, 
                                                 chain_id='L', 
                                                 res_name='LIG')
   
   # Merge with protein structure
   complex_df = pd.concat([protein_df, ligand_df], ignore_index=True)
   
   # Save using CIF handler
   from protos.io.cif_handler import CifHandler
   handler = CifHandler()
   handler.write('complex.cif', complex_df)
   ```

4. **Format Handler Integration** - SDF files are now part of Protos format registry:
   - Automatic format detection
   - Consistent read/write interface
   - Validation support
   - Works with both .sdf and .mol extensions

5. **SDF Utilities** (in protos.io.sdf_utils):
   - `read_sdf_file()` - Read with optional sanitization
   - `write_sdf_file()` - Write with property selection
   - `sdf_to_dataframe()` - Convert to pandas DataFrame
   - `dataframe_to_sdf()` - Convert from DataFrame
   - `merge_sdf_files()` - Combine multiple SDFs
   - `validate_sdf_file()` - Check file validity
   - `filter_sdf_by_property()` - Filter by property values
   - `extract_unique_properties()` - List all properties

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

### 4. Creating Protein-Ligand Complexes (NEW)
```python
# Download ligand from ChEMBL
compounds = lig_proc.get_protein_ligands("EGFR", limit=1)
smiles = compounds[0]['smiles']

# Convert to CIF format (default)
ligand_df = lig_proc.convert_to_cif_dataframe(smiles, 
                                              chain_id='L',
                                              res_name='INH')

# Load protein structure
struct_proc.load_structure('1M17')
protein_df = struct_proc.data

# Create complex
complex_df = pd.concat([protein_df, ligand_df], ignore_index=True)

# Save complex as CIF
from protos.io.cif_handler import CifHandler
handler = CifHandler()
handler.write('1M17_with_inhibitor.cif', complex_df)
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