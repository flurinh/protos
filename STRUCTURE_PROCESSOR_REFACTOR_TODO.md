# StructureProcessor Refactoring: PKL-Canonical Architecture

## Overview

Establish PKL as the sole canonical format for structures. CIF serves only for initial ingest and optional export. No versioning, no format-specific registry keys - just "structure" pointing to PKL files.

---

## Core Architecture

### Registry Simplification
- **One registration per structure**: `format_type="structure"` → PKL file
- **No CIF registration**: Keeps registry clean and simple
- **Metadata in PKL entry**: Source, URL, ingest time, content hash

### Storage Invariants
- **Memory**: `self.frames: Dict[str, DataFrame]` with lazy `self._data` stacking
- **MultiIndex**: Always `(structure_id, atom_id)`, sorted
- **Schema**: All ints as `Int64`, coords numeric, no `pdb_id` column
- **Paths**: All via ProtosPaths only - no hardcoded paths, no fallback cascades

---

## Implementation Tasks

### Phase 1: Core Methods

#### 1.1 Canonicalization
```python
def _ensure_canonical(df: pd.DataFrame, structure_id: str) -> pd.DataFrame:
    """Transform any DataFrame to canonical form"""
    # Tasks:
    # - Remove pdb_id column if exists (map to structure_id first)
    # - Ensure structure_id and atom_id columns
    # - Set MultiIndex (structure_id, atom_id)
    # - Apply STRUCT_COLUMN_DTYPE with Int64 for integers
    # - Ensure optional columns (grn as empty string)
    # - Validate coordinate columns are numeric
```

#### 1.2 Content Hashing
```python
def compute_content_hash(df: pd.DataFrame) -> str:
    """Hash core structural content for change detection"""
    # Tasks:
    # - Select core columns: chain, seq_id, residue, atom, x, y, z
    # - Sort by index for consistency
    # - Generate sha256 hash
```

#### 1.3 Load Entity
```python
def load_entity(self, structure_id: str, *, auto_ingest: bool = True) -> Optional[pd.DataFrame]:
    """PKL-first loading with auto-ingest fallback"""
    # Tasks:
    # 1. Try registry lookup (format_type="structure")
    # 2. Try PKL in cache dir (auto-register if found)
    # 3. If auto_ingest and CIF exists: ingest_cif()
    # 4. Update frames storage
```

#### 1.4 Save Entity
```python
def save_entity(self, structure_id: str, df: pd.DataFrame, 
                format: str = 'pkl', metadata: Optional[Dict] = None) -> None:
    """Save canonical PKL and register"""
    # Tasks:
    # - Validate format='pkl' only
    # - Canonicalize DataFrame
    # - Compute content hash
    # - Write PKL to cache dir
    # - Register/update with format_type="structure"
    # - Include metadata: source, hash, schema version
```

#### 1.5 Ingest CIF
```python
def ingest_cif(self, structure_id: str, cif_path: Optional[Path] = None) -> pd.DataFrame:
    """One-time CIF to PKL conversion"""
    # Tasks:
    # - Find CIF (provided path or standard location)
    # - Parse using existing CIF parser
    # - Canonicalize DataFrame
    # - Save as PKL via save_entity()
    # - Update storage
```

#### 1.6 Export Entity
```python
def export_entity(self, structure_id: str, format: str = 'cif', 
                  out_path: Optional[Path] = None) -> Path:
    """Export PKL to other formats (no registration)"""
    # Tasks:
    # - Load from PKL (must exist)
    # - Reset index for CIF compatibility
    # - Write CIF to structure dir
    # - Return path (no registry update)
```

### Phase 2: Path Management

#### 2.1 Path Resolution
```python
@property
def path_pkl_dir(self) -> Path:
    """Get PKL directory from ProtosPaths"""
    # Use ProtosPaths to get the cache directory for PKL files
    return Path(self.paths.get_subdir_path('structure', 'cache_dir'))

@property
def path_cif_dir(self) -> Path:
    """Get CIF directory from ProtosPaths"""
    # Use ProtosPaths to get the structure directory for CIF files
    return Path(self.paths.get_subdir_path('structure', 'structure_dir'))

@property
def path_dataset_dir(self) -> Path:
    """Get dataset directory from ProtosPaths"""
    # Use ProtosPaths to get the dataset directory
    return Path(self.paths.get_subdir_path('structure', 'dataset_dir'))

@property
def path_temp_dir(self) -> Path:
    """Get temp directory from ProtosPaths"""
    # Use ProtosPaths to get the temp directory
    return Path(self.paths.get_subdir_path('structure', 'temp_dir'))
```

### Phase 3: Storage Management

#### 3.1 Frame Management
```python
def _set_frame(self, structure_id: str, df: pd.DataFrame):
    """Update frame and mark dirty"""
    self.frames[structure_id] = df
    self._dirty = True

def _remove_frame(self, structure_id: str):
    """Remove frame and mark dirty"""
    if structure_id in self.frames:
        del self.frames[structure_id]
        self._dirty = True

@property
def data(self) -> Optional[pd.DataFrame]:
    """Lazy stacked view with rebuild on access"""
    if self._dirty and self.frames:
        self._data = pd.concat(self.frames.values()).sort_index()
        self._dirty = False
    return self._data
```

### Phase 4: Dataset Operations

#### 4.1 Load Dataset
```python
def load_dataset(self, dataset_name: str, return_format: str = 'stacked'):
    """Load all dataset members from individual PKLs"""
    # Tasks:
    # - Get member list from dataset
    # - Load each via load_entity(auto_ingest=False)
    # - Update frames in bulk
    # - Return stacked or dict format
```

#### 4.2 Save Dataset  
```python
def save_dataset(self, dataset_name: str, structure_ids: List[str], metadata: Optional[Dict] = None):
    """Ensure PKLs exist for all members"""
    # Tasks:
    # - Check each member has registered PKL
    # - Save missing ones from frames if available
    # - Create logical dataset entry
```

### Phase 5: Deprecations & Cleanup

#### 5.1 Remove Legacy
- Remove `pdb_ids` property completely
- Remove any CIF-first loading paths
- Remove versioning logic
- Clean up old column names

#### 5.2 Add Deprecation Warnings
- `load_structure()` → `load_entity()`
- `save_structure()` → `save_entity()` or `export_entity()`
- Format parameter in save (only pkl allowed)

---

## Refactoring Checklist

### Immediate Changes
- [ ] Implement `_ensure_canonical()` with full schema enforcement
- [ ] Implement `compute_content_hash()` for metadata
- [ ] Update `load_entity()` to be PKL-first with registry lookup
- [ ] Restrict `save_entity()` to PKL only
- [ ] Add `ingest_cif()` for one-time conversions
- [ ] Add `export_entity()` for CIF generation without registration

### Path Management
- [ ] Add `path_pkl_dir` property using ProtosPaths
- [ ] Add `path_cif_dir` property using ProtosPaths  
- [ ] Add `path_dataset_dir` property using ProtosPaths
- [ ] Add `path_temp_dir` property using ProtosPaths
- [ ] Ensure all paths use ProtosPaths.get_subdir_path()

### Storage Updates
- [ ] Replace list storage with dict (`self.frames`)
- [ ] Implement lazy stacking with dirty flag
- [ ] Add `_set_frame()` and `_remove_frame()` helpers
- [ ] Ensure MultiIndex on all DataFrames

### Registry Updates  
- [ ] Use generic `format_type="structure"` for all registrations
- [ ] Store CIF provenance in PKL metadata only
- [ ] Remove any CIF registrations
- [ ] Add content hash to all metadata

### Dataset Handling
- [ ] Update `load_dataset()` to load from individual PKLs
- [ ] Update `save_dataset()` to ensure PKLs exist
- [ ] Remove any dataset-level PKL creation
- [ ] Keep datasets as logical collections only

### Testing Requirements
- [ ] Test canonicalization with various input formats
- [ ] Test PKL-first loading with fallbacks
- [ ] Test CIF ingest creates correct PKL
- [ ] Test export doesn't affect registry
- [ ] Test dataset operations with missing PKLs
- [ ] Test path resolution with different configurations

---

## Key Principles

1. **PKL is Truth**: All operations work from PKL files
2. **CIF is Transient**: Parse once on ingest, generate on export
3. **No Versioning**: Always overwrite in place
4. **Simple Registry**: One entry per structure, generic type
5. **Lazy Operations**: Stack only when needed
6. **Path Consistency**: Use ProtosPaths exclusively for all path resolution

---

## Migration Strategy

1. **Existing PKLs**: Check and canonicalize if needed
2. **Existing CIFs**: Ingest to PKL on first load
3. **Registry Cleanup**: Remove CIF entries, update to generic type
4. **Path Migration**: Update to use ProtosPaths properties