# Structure Analysis Functions Implementation TODO

## Overview
Separate read-only analysis functions from data-modifying operations in StructureProcessor. Analysis functions will be organized into focused modules under `protos/analysis/structure/`.

## 1. Directory Structure Creation

### Create the following structure:
```
protos/
└── analysis/
    └── structure/
        ├── __init__.py
        ├── geometry.py       # Distance and coordinate calculations
        ├── alignment.py      # Structure alignment algorithms
        ├── quality.py        # Structure validation and quality checks
        ├── properties.py     # Structural property calculations
        └── sequence.py       # Sequence extraction and analysis
```

## 2. Module Implementation Details

### 2.1 `geometry.py` - Geometric Analysis Functions

**Functions to implement:**

1. **`calculate_distance_matrix(coords: np.ndarray) -> np.ndarray`**
   - Calculate pairwise distance matrix between all atoms
   - Input: Nx3 coordinate array
   - Output: NxN distance matrix
   - Use efficient numpy operations

2. **`calculate_distance(coord1: np.ndarray, coord2: np.ndarray) -> float`**
   - Calculate distance between two points
   - Handle both single points and arrays of points

3. **`find_contacts(df: pd.DataFrame, cutoff: float = 5.0, exclude_same_residue: bool = True) -> pd.DataFrame`**
   - Find all atom pairs within cutoff distance
   - Options to exclude intra-residue contacts
   - Return DataFrame with atom1, atom2, distance

4. **`calculate_center_of_mass(df: pd.DataFrame, mass_column: Optional[str] = None) -> np.ndarray`**
   - Calculate center of mass for structure
   - Use atomic masses if provided, else uniform weight

5. **`calculate_radius_of_gyration(df: pd.DataFrame) -> float`**
   - Calculate radius of gyration
   - Measure of structure compactness

### 2.2 `alignment.py` - Structure Alignment Functions

**Functions to implement:**

1. **`calculate_rmsd(df1: pd.DataFrame, df2: pd.DataFrame, atom_selection: str = 'CA') -> float`**
   - Calculate RMSD between two structures
   - Options for different atom selections (CA, backbone, all)
   - Handle missing atoms gracefully

2. **`kabsch_alignment(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]`**
   - Kabsch algorithm for optimal superposition
   - Return: rotation matrix, translation vector, RMSD
   - Used by other alignment functions

3. **`align_structures(df1: pd.DataFrame, df2: pd.DataFrame, method: str = 'kabsch', atom_selection: str = 'CA') -> Tuple[np.ndarray, np.ndarray]`**
   - General structure alignment function
   - Support different methods (kabsch, quaternion)
   - Return transformation matrices

4. **`sequence_based_alignment(df1: pd.DataFrame, df2: pd.DataFrame) -> Dict[int, int]`**
   - Align structures based on sequence alignment
   - Return residue mapping between structures

### 2.3 `quality.py` - Structure Validation Functions

**Functions to implement:**

1. **`check_missing_atoms(df: pd.DataFrame, residue_templates: Optional[Dict] = None) -> pd.DataFrame`**
   - Check for missing atoms in each residue
   - Use standard residue templates
   - Return DataFrame with residue, missing_atoms

2. **`validate_bond_lengths(df: pd.DataFrame, tolerance: float = 0.5) -> pd.DataFrame`**
   - Check if bond lengths are within expected ranges
   - Flag unusual bond lengths
   - Return DataFrame with bond, length, expected, deviation

3. **`check_clashes(df: pd.DataFrame, clash_cutoff: float = 2.0) -> pd.DataFrame`**
   - Find steric clashes (atoms too close)
   - Exclude bonded atoms
   - Return DataFrame with clashing atom pairs

4. **`calculate_b_factor_statistics(df: pd.DataFrame) -> Dict[str, float]`**
   - Calculate B-factor statistics
   - Mean, std, percentiles by chain/residue type
   - Identify high B-factor regions

5. **`validate_chirality(df: pd.DataFrame) -> pd.DataFrame`**
   - Check amino acid chirality (L vs D)
   - Flag unusual chirality
   - Important for structure quality

### 2.4 `properties.py` - Structural Property Calculations

**Functions to implement:**

1. **`calculate_secondary_structure(df: pd.DataFrame) -> pd.Series`**
   - Assign secondary structure elements
   - Use DSSP-like algorithm
   - Return Series with residue -> ss_type mapping

2. **`calculate_solvent_accessibility(df: pd.DataFrame, probe_radius: float = 1.4) -> pd.Series`**
   - Calculate solvent accessible surface area
   - Per-residue SASA values
   - Use rolling ball algorithm

3. **`calculate_hydrophobic_moment(df: pd.DataFrame, window_size: int = 11) -> pd.Series`**
   - Calculate hydrophobic moment for helices
   - Sliding window approach
   - Useful for membrane proteins

4. **`identify_binding_sites(df: pd.DataFrame, ligand_df: Optional[pd.DataFrame] = None) -> List[Dict]`**
   - Identify potential binding pockets
   - Use cavity detection algorithms
   - Return list of pocket descriptions

5. **`calculate_electrostatic_potential(df: pd.DataFrame, grid_spacing: float = 1.0) -> np.ndarray`**
   - Calculate electrostatic potential on grid
   - Simple Coulomb potential
   - Return 3D grid of potential values

### 2.5 `sequence.py` - Sequence Extraction Functions

**Functions to implement:**

1. **`extract_sequence(df: pd.DataFrame, chain_id: Optional[str] = None, one_letter: bool = True) -> str`**
   - Extract amino acid sequence from structure
   - Handle missing residues with 'X'
   - Option for one or three letter codes

2. **`extract_all_sequences(df: pd.DataFrame) -> Dict[str, str]`**
   - Extract sequences for all chains
   - Return dict of chain_id -> sequence

3. **`map_structure_to_sequence(df: pd.DataFrame, sequence: str) -> Dict[int, int]`**
   - Map structure residue numbers to sequence positions
   - Handle insertions and deletions
   - Return mapping dict

4. **`identify_missing_residues(df: pd.DataFrame, expected_sequence: str) -> List[Tuple[int, str]]`**
   - Compare structure to expected sequence
   - Identify missing residues
   - Return list of (position, residue) tuples

## 3. Data Modification Methods for StructureProcessor

### 3.1 `merge_structures(structure_ids: List[str], chain_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame`
**Purpose:** Combine multiple structures into single multi-chain structure
**Implementation:**
- Load all structures
- Optionally rename chains to avoid conflicts
- Concatenate DataFrames
- Re-canonicalize with new structure_id
- Handle coordinate conflicts if needed

### 3.2 `align_structures(structure_ids: List[str], reference_id: str, method: str = 'kabsch') -> Dict[str, pd.DataFrame]`
**Purpose:** Align multiple structures and apply transformations
**Implementation:**
- Use `alignment.align_structures()` to calculate transformations
- Apply transformations using `apply_transformation()`
- Return aligned structures
- Option to save as new entities or modify in-place

### 3.3 `orient_structure(structure_id: str, method: str = 'principal_axes') -> pd.DataFrame`
**Purpose:** Orient structure using standard methods
**Implementation:**
- Methods: 'principal_axes', 'membrane_normal', 'inertia_tensor'
- Calculate orientation matrix
- Apply transformation
- Useful for consistent visualization

### 3.4 `renumber_residues(structure_id: str, start: int = 1, by_chain: bool = True) -> pd.DataFrame`
**Purpose:** Renumber residues sequentially
**Implementation:**
- Renumber auth_seq_id sequentially
- Option to restart numbering for each chain
- Maintain mapping of old -> new numbers
- Update both auth_seq_id and label_seq_id

## 4. Integration Pattern

### StructureProcessor uses analysis functions:
```python
def align_structures(self, structure_ids, reference_id):
    from protos.analysis.structure.alignment import align_structures, calculate_rmsd
    
    # Get structures
    ref_df = self.load_entity(reference_id)
    
    results = {}
    for struct_id in structure_ids:
        df = self.load_entity(struct_id)
        
        # Calculate alignment (read-only)
        rotation, translation = align_structures(df, ref_df)
        rmsd = calculate_rmsd(df, ref_df)
        
        # Apply transformation (modifies data)
        df_aligned = self.apply_transformation(struct_id, rotation, translation)
        
        results[struct_id] = {
            'aligned_structure': df_aligned,
            'rmsd': rmsd,
            'transformation': (rotation, translation)
        }
    
    return results
```

## 5. Testing Strategy

### 5.1 Unit Tests for Analysis Functions
- Test with known structures (e.g., simple geometric shapes)
- Verify mathematical correctness
- Test edge cases (empty structures, single atom, etc.)

### 5.2 Integration Tests
- Test StructureProcessor methods that use analysis functions
- Verify data modification works correctly
- Test with real PDB structures

### 5.3 Performance Tests
- Test with large structures (>100k atoms)
- Benchmark critical functions
- Optimize bottlenecks

## 6. Documentation Requirements

### 6.1 API Documentation
- Detailed docstrings for all functions
- Include mathematical formulas where relevant
- Usage examples

### 6.2 User Guide
- How to use StructureProcessor for common tasks
- When to use analysis functions vs processor methods
- Best practices for structure manipulation

## Implementation Priority:
1. Create directory structure ✓
2. Implement geometry.py (most fundamental) ✓
3. Implement alignment.py (commonly needed) ✓
4. Implement quality.py ✓
5. Implement properties.py ✓
6. Implement sequence.py ✓
7. Update StructureProcessor with data modification methods ✓
   - merge_structures() ✓
   - align_structures() ✓
   - orient_structure() ✓
   - renumber_residues() ✓
8. Write tests (pending)
9. Update documentation (pending)

## Status Summary:
- ✅ All analysis modules implemented
- ✅ All data modification methods added to StructureProcessor
- ⏳ Tests need to be created
- ⏳ Documentation needs to be updated