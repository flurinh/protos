# Structure Utils Refactoring Plan

**Status: COMPLETED** ✅

**Completion Date**: 2025-01-19

## Overview

The current structure utilities are scattered across multiple files with mixed concerns. This document outlines a comprehensive plan to refactor and reorganize these utilities following clean separation of concerns.

## Current State Analysis

### 1. `struct_utils.py` (Main offender - mixed concerns)

**Format/IO Functions:**
- `load_structure()` - Loads CIF files and converts to DataFrame
- `get_protein_structure_df()` - Extracts coordinates from gemmi model

**Analysis Functions:**
- `get_distance_matrix()` - Calculates pairwise distances
- `calculate_membrane_normal()` - PCA-based membrane normal calculation
- `orient_structure_in_membrane()` - Iterative membrane orientation
- `flip_protein()` - 180° rotation around x-axis
- `apply_rotation_matrix()` - General rotation application
- `calc_rotation_matrix()` - Rodrigues' rotation formula
- `align_proteins_on_retinal()` - Retinal-based alignment
- `annotate_helix_numbers()` - Helix annotation
- `orient_structures_n_terminus_up()` - N-terminus orientation
- `normalize_structures()` - Structure normalization

**Constants/Schema:**
- `CIF_COLUMNS` - Column names for CIF parsing
- `STRUCT_COLUMNS` - Internal DataFrame column names
- `STRUCT_COLUMN_DTYPE` - Data types for columns

### 2. `structure_utils.py` (Better organized)

Already follows better practices:
- Imports schemas from `schema_definitions.py`
- Uses standardized interfaces
- Better separation of concerns

**Functions:**
- `three_to_one_letter_code()` - Amino acid conversion
- `load_structure()` - Refactored version with schemas
- `cif_block_to_df()` - CIF parsing with schema validation
- `extract_torsion_angles()` - Torsion angle calculation
- `normalize_structures()` - Structure normalization
- `get_ca_ret_coords()` - CA/retinal coordinate extraction

### 3. `struct_alignment.py`

Pure analysis functions:
- `align_structures()` - Structure alignment algorithm
- `get_structure_alignment()` - Alignment with transformation
- `structure_comparison_ava()` - All-vs-all comparison
- `structure_comparison_1va()` - One-vs-all comparison

### 4. `cif_handler.py` & `cif_utils.py`

Well-organized format handling:
- `CifHandler` - Format handler class
- `cif_to_df()` - CIF parsing
- `df_to_cif()` - CIF writing
- Column mappings and validation

## Refactoring Plan

### ✅ COMPLETED REFACTORING

The refactoring has been successfully completed with the following changes:

### Phase 1: Create New Structure

```
protos/
├── analysis/
│   └── structure/
│       ├── __init__.py
│       ├── geometry.py          # Geometric calculations
│       ├── alignment.py         # Structure alignment
│       ├── membrane.py          # Membrane-specific analysis
│       └── comparison.py        # Structure comparison
│
└── io/
    └── formats/
        ├── cif_handler.py       # Keep as-is
        ├── cif_utils.py         # Consolidate all CIF utilities
        └── structure_schema.py  # Structure-specific schemas
```

### Phase 2: File Contents After Refactoring

#### `protos/io/formats/structure_schema.py`
```python
"""
Structure data schemas and constants.
"""

# Column definitions
STRUCTURE_COLUMNS = [...]
STRUCTURE_COLUMN_DTYPE = {...}
CIF_COLUMNS = [...]

# Validation functions
def validate_structure_df(df: pd.DataFrame) -> bool:
    """Validate structure DataFrame against schema."""
    pass
```

#### `protos/io/formats/cif_utils.py` (Enhanced)
```python
"""
Complete CIF format utilities.
Consolidates all CIF parsing and writing functionality.
"""

# Move from struct_utils.py:
def load_structure_from_cif(path: str) -> pd.DataFrame:
    """Load CIF file and return standardized DataFrame."""
    # Refactored version of load_structure()
    pass

def extract_structure_from_model(model, chain_id: str) -> pd.DataFrame:
    """Extract structure data from gemmi model."""
    # Refactored version of get_protein_structure_df()
    pass

# Keep existing functions from cif_utils.py
```

#### `protos/analysis/structure/geometry.py`
```python
"""
Geometric analysis functions for protein structures.
"""

def calculate_distance_matrix(coordinates: Union[np.ndarray, pd.DataFrame]) -> np.ndarray:
    """Calculate pairwise distance matrix."""
    # From struct_utils.get_distance_matrix()
    pass

def calculate_rotation_matrix(vec_from: np.ndarray, vec_to: np.ndarray) -> np.ndarray:
    """Calculate rotation matrix using Rodrigues' formula."""
    # From struct_utils.calc_rotation_matrix()
    pass

def apply_rotation(df: pd.DataFrame, rotation_matrix: np.ndarray) -> pd.DataFrame:
    """Apply rotation matrix to coordinates."""
    # From struct_utils.apply_rotation_matrix()
    pass

def apply_translation(df: pd.DataFrame, translation: np.ndarray) -> pd.DataFrame:
    """Apply translation to coordinates."""
    pass

def flip_structure(df: pd.DataFrame, axis: str = 'x') -> pd.DataFrame:
    """Flip structure around specified axis."""
    # From struct_utils.flip_protein()
    pass
```

#### `protos/analysis/structure/membrane.py`
```python
"""
Membrane protein analysis functions.
"""

def calculate_membrane_normal(df_ca: pd.DataFrame) -> np.ndarray:
    """Calculate membrane normal using PCA."""
    # From struct_utils.calculate_membrane_normal()
    pass

def orient_in_membrane(df: pd.DataFrame, max_iterations: int = 12) -> pd.DataFrame:
    """Orient structure with membrane normal along z-axis."""
    # From struct_utils.orient_structure_in_membrane()
    pass

def orient_n_terminus_up(structures: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """Orient structures with N-terminus pointing up."""
    # From struct_utils.orient_structures_n_terminus_up()
    pass

def annotate_transmembrane_helices(df: pd.DataFrame, k: int = 7, w: int = 5) -> pd.DataFrame:
    """Annotate transmembrane helix numbers."""
    # From struct_utils.annotate_helix_numbers()
    pass
```

#### `protos/analysis/structure/alignment.py`
```python
"""
Structure alignment and superposition functions.
"""

# Move existing functions from struct_alignment.py
def align_structures(...):
    pass

def get_structure_alignment(...):
    pass

# Add from struct_utils.py:
def align_on_retinal(structures: Dict[str, pd.DataFrame], 
                    reference_id: str = '4PXK') -> Dict[str, pd.DataFrame]:
    """Align structures based on retinal position."""
    # From struct_utils.align_proteins_on_retinal()
    pass

def extract_ca_and_ligand_coords(df: pd.DataFrame, 
                                 chain_id: str,
                                 ligand_id: str = 'RET') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Extract CA atoms and ligand coordinates."""
    # From struct_utils.get_ca_ret_coords()
    pass
```

#### `protos/analysis/structure/comparison.py`
```python
"""
Structure comparison and analysis functions.
"""

def compare_all_vs_all(structures: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Perform all-vs-all structure comparison."""
    # From struct_alignment.structure_comparison_ava()
    pass

def compare_one_vs_all(structures: Dict[str, pd.DataFrame], 
                      reference: str) -> pd.DataFrame:
    """Compare one structure against all others."""
    # From struct_alignment.structure_comparison_1va()
    pass

def normalize_structures(structures: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    """Normalize structure coordinates and B-factors."""
    # From struct_utils.normalize_structures()
    pass
```

### Phase 3: Migration Completed ✅

#### Files Created:

1. **`protos/io/formats/structure_schema.py`** ✅
   - Column definitions with `structure_id` (replaced `pdb_id`)
   - Data type mappings
   - Validation functions

2. **`protos/io/formats/cif_utils.py`** (Enhanced) ✅
   - Added `load_structure_from_cif()` - replaces old `load_structure()`
   - Added `extract_structure_from_model()` - replaces `get_protein_structure_df()`
   - Simplified column mapping (e.g., `Cartn_x` → `x`)
   - Support for `grn` (generic residue number)

3. **`protos/analysis/structure/geometry.py`** ✅
   - `calculate_distance_matrix()`
   - `calculate_rotation_matrix()`
   - `apply_rotation()`, `apply_translation()`
   - `flip_structure()`

4. **`protos/analysis/structure/membrane.py`** ✅
   - `calculate_membrane_normal()`
   - `orient_in_membrane()`
   - `orient_n_terminus_up()`
   - `annotate_transmembrane_helices()`

5. **`protos/analysis/structure/alignment.py`** ✅
   - `align_structures()`
   - `get_structure_alignment()`
   - `align_on_retinal()`
   - `extract_ca_and_ligand_coords()`

6. **`protos/analysis/structure/comparison.py`** ✅
   - `compare_all_vs_all()`
   - `compare_one_vs_all()`
   - `normalize_structures()`

#### Files Deleted:

- **`struct_utils.py`** ✅ - All functions migrated
- **`struct_alignment.py`** ✅ - Moved to `analysis/structure/alignment.py`

### Phase 4: Import Updates ✅

#### Import Updates Completed:

1. **`structure_processor.py`** ✅
   - Now imports from `protos.io.formats.cif_utils` and `structure_schema`
   - Uses `structure_id` parameter consistently
   - Uses `load_structure_from_cif()` instead of old `load_structure_util()`
   - Properly handles `pdb_id` to `structure_id` mapping in `_ensure_canonical()`

2. **`structure_processor_deprecated.py`** ✅
   - Updated imports to new locations
   - Note: This file still exists but uses old patterns
   
3. **`visualization/structure_vis.py`** ✅
   - Updated from `struct_alignment` to `protos.analysis.structure.alignment`

4. **`to_be_updated/binding_domain_alignment.py`** ✅
   - Updated alignment imports

#### Old imports:
```python
from protos.processing.structure.struct_utils import load_structure, get_distance_matrix
from protos.processing.structure.struct_alignment import get_structure_alignment
```

#### New imports:
```python
from protos.io.formats.cif_utils import load_structure_from_cif
from protos.analysis.structure.geometry import calculate_distance_matrix
from protos.analysis.structure.alignment import get_structure_alignment
```

## Additional Improvements

### Column Standardization ✅

During the refactoring, we also standardized column naming across all structure files:

1. **Replaced `pdb_id` with `structure_id`** throughout the codebase
2. **Removed non-existent columns** like `auth_atom_name`
3. **Simplified CIF column mapping**:
   - `Cartn_x` → `x`, `Cartn_y` → `y`, `Cartn_z` → `z`
   - `label_atom_id` → `atom_name`
   - `label_comp_id` → `res_name`
   - `auth_asym_id` → `auth_chain_id`

4. **Added support for**:
   - `grn` (Generic Residue Number)
   - `phi`, `psi`, `omega` (torsion angles from gemmi)

### Column Mapping Updates ✅

The CIF column mapping in `cif_utils.py` now properly maps CIF columns to simplified internal names:
- Removed non-existent `auth_atom_name` mapping
- Kept actual CIF column names that exist in mmCIF format
- Added support for calculated columns (phi, psi, omega) and grn

## Benefits of Refactoring

1. **Clear Separation of Concerns**
   - I/O operations in `io/formats/`
   - Analysis functions in `analysis/structure/`
   - No mixing of file operations with analysis

2. **Better Organization**
   - Related functions grouped together
   - Easier to find and maintain code
   - Logical import paths

3. **Improved Testability**
   - Pure functions without I/O dependencies
   - Easier to mock and test components

4. **Extensibility**
   - Clear places to add new analysis functions
   - Format handling separate from analysis

5. **Consistency**
   - Standardized function signatures
   - Consistent error handling
   - Unified coordinate system handling

## Implementation Summary

### ✅ All Components Completed:

1. **Core functionality**
   - `structure_schema.py` - Schemas defined with `structure_id`
   - `geometry.py` - All geometric operations migrated
   - Enhanced `cif_utils.py` - I/O consolidated with simplified mappings

2. **Domain-specific**
   - `membrane.py` - All membrane protein analysis functions
   - `alignment.py` - Structure alignment with retinal support

3. **Extended features**
   - `comparison.py` - RMSD calculations and normalization
   - Old files deleted (`struct_utils.py`, `struct_alignment.py`)

## Testing Strategy

1. **Unit Tests** for each new module
2. **Integration Tests** to ensure refactored code produces same results
3. **Performance Tests** to ensure no regression
4. **Migration Tests** to verify old code still works during transition

## Deprecation Status

### ✅ Completed:

1. **Phase 1** ✅: New structure created, functions migrated
2. **Phase 2** ✅: Old files removed (no deprecation period needed as internal refactoring)
3. **Phase 3** ✅: All internal usage updated
4. **Phase 4** ✅: Deprecated code removed

### Files Removed:
- `protos/processing/structure/struct_utils.py`
- `protos/processing/structure/struct_alignment.py`

### Files Still Using Old Patterns:
- `protos/processing/structure/structure_utils.py` - Uses old schema_definitions imports but follows better practices

## Notes

- Preserve all functionality during refactoring
- Maintain backward compatibility during transition
- Document all changes in CHANGELOG
- Update user documentation with new import paths