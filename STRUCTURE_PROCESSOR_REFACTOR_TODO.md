# StructureProcessor Refactoring Plan: PKL-Canonical Architecture

## Executive Summary

This refactoring establishes PKL as the canonical storage format for structures, with mmCIF serving strictly for ingest and optional export. The plan leverages BaseProcessor's entity system while implementing proper versioning, content hashing, and clear format roles.

### Core Principles
1. **PKL is canonical** - All structure operations work from PKL files
2. **mmCIF for I/O only** - Read once → normalize → PKL (ingest), or PKL → CIF (export)
3. **Content-based versioning** - New versions only when atoms/coordinates change
4. **Format-specific registry** - Use "STRUCTURE_PKL" and "MMCIF" as format keys
5. **Preserve BaseProcessor architecture** - No separate IO classes, leverage existing entity system

---

## Canonical Storage Policy

### File Formats and Roles

**STRUCTURE_PKL** (Canonical)
- **Role**: Primary storage format, always present for usable entities
- **Structure**: Per-structure PKL with MultiIndex `(structure_id, atom_id)`
- **Location**: `structure/cache/{structure_id}.pkl`

**MMCIF** (I/O Only)
- **Roles**:
  - `source`: Original file preserved during ingest (optional)
  - `export`: Generated from PKL on demand
- **Location**: 
  - Source: `structure/raw/{structure_id}.cif` (if preserved)
  - Export: `structure/mmcif/{structure_id}.cif`

### Registry Schema

```python
# Format-specific registrations
entity_registry.register_entity(
    name="1ubq",
    format_type="STRUCTURE_PKL",  # Not just "structure"
    file_path="structure/cache/1ubq.pkl",
    metadata={
        "schema_version": "v2",
        "io_version": "1",
        "content_hash": "sha256:abcd1234...",  # Core columns only
        "role": "canonical",
        "source": "rcsb",
        "source_url": "https://files.rcsb.org/download/1ubq.cif",
        "ingested_at": "2025-01-20T12:00:00Z",
        "atom_count": 660,
        "chains": ["A"],
    }
)

# If source was preserved
entity_registry.register_entity(
    name="1ubq", 
    format_type="MMCIF",
    file_path="structure/raw/1ubq.cif",
    metadata={
        "role": "source",
        "source": "rcsb",
        "preserved_at": "2025-01-20T12:00:00Z"
    }
)

# Relationships
entity_registry.add_relationship(
    source_name="1ubq",
    target_name="1ubq", 
    relationship_type="derived_from",
    metadata={"from_format": "MMCIF", "to_format": "STRUCTURE_PKL"}
)
```

---

## Versioning Policy

### When to Overwrite (Same Entity)
- Schema migrations (v1 → v2)
- Annotation changes (adding GRN column)
- Data type fixes (int → Int64)
- **Keep same `content_hash`**, bump `schema_version` or `io_version`

### When to Version (New Entity)
- Atom changes (added/removed atoms)
- Coordinate changes (mutations, relaxations)
- Residue renumbering
- **New `content_hash`** triggers versioning

### Version Naming
```python
def default_version_namer(structure_id: str) -> str:
    """Default: append @YYYY-MM-DD"""
    return f"{structure_id}@{datetime.now().strftime('%Y-%m-%d')}"

# Examples:
# 1ubq → 1ubq@2025-01-20 (new version)
# Relationship: 1ubq@2025-01-20 version_of 1ubq
# Alias: 1ubq → 1ubq@2025-01-20 (latest)
```

---

## Phase 0 — Content Hashing & Canonicalization

### Content Hash Implementation
```python
def compute_content_hash(df: pd.DataFrame) -> str:
    """Compute hash of core structural content"""
    # Core columns only (exclude analysis columns like grn)
    CORE_COLUMNS = [
        'auth_chain_id', 'auth_seq_id', 'res_name3l', 
        'atom_name', 'x', 'y', 'z'
    ]
    
    # Sort by index for consistency
    df_sorted = df[CORE_COLUMNS].sort_index()
    
    # Hash the core content
    import hashlib
    content_bytes = pd.util.hash_pandas_object(df_sorted).values.tobytes()
    return f"sha256:{hashlib.sha256(content_bytes).hexdigest()}"
```

### Canonical DataFrame Structure
```python
def _ensure_canonical(self, df: pd.DataFrame, structure_id: str) -> pd.DataFrame:
    """Ensure DataFrame meets canonical requirements"""
    df = df.copy()
    
    # 1. Handle legacy pdb_id
    if 'pdb_id' in df.columns and 'structure_id' not in df.columns:
        df['structure_id'] = df['pdb_id']
    if 'pdb_id' in df.columns:
        df = df.drop(columns=['pdb_id'])
    
    # 2. Required columns
    if 'structure_id' not in df.columns:
        df['structure_id'] = structure_id
    if 'atom_id' not in df.columns:
        raise ValueError("DataFrame must have atom_id column")
    
    # 3. Set MultiIndex
    if not isinstance(df.index, pd.MultiIndex):
        df = df.set_index(['structure_id', 'atom_id'])
    df = df.sort_index()
    
    # 4. Apply dtypes
    for col, dtype in STRUCT_COLUMN_DTYPE_V2.items():
        if col in df.columns:
            if dtype == int:
                df[col] = df[col].astype('Int64')
            else:
                df[col] = df[col].astype(dtype, errors='ignore')
    
    # 5. Ensure optional columns exist
    if 'grn' not in df.columns:
        df['grn'] = ''
    
    # 6. Validate coordinates
    for coord in ['x', 'y', 'z']:
        df[coord] = pd.to_numeric(df[coord], errors='coerce')
    
    return df
```

---

## Phase 1 — Core API Implementation

### Load (PKL-First)
```python
def load_entity(self, structure_id: str, *, auto_ingest: bool = True) -> Optional[pd.DataFrame]:
    """Load structure from PKL, with optional auto-ingestion from CIF"""
    
    # 1. Try PKL from registry
    pkl_info = self.entity_registry.find_entity(structure_id, "STRUCTURE_PKL")
    if pkl_info:
        pkl_path = self._resolve_path(pkl_info.file_path)
        df = pd.read_pickle(pkl_path)
        self._set_frame(structure_id, df)
        return df
    
    # 2. Try cache directory
    cache_path = self.path_cache_dir / f"{structure_id}.pkl"
    if cache_path.exists():
        df = pd.read_pickle(cache_path)
        # Auto-register if missing
        self._auto_register_pkl(structure_id, df, cache_path)
        self._set_frame(structure_id, df)
        return df
    
    # 3. Auto-ingest from CIF if available
    if auto_ingest:
        # Check for registered MMCIF
        cif_info = self.entity_registry.find_entity(structure_id, "MMCIF")
        if cif_info:
            return self.ingest_cif(structure_id, preserve_source=False)
        
        # Check for CIF file on disk
        cif_path = self.path_structure_dir / f"{structure_id}.cif"
        if cif_path.exists():
            return self.ingest_cif(structure_id, preserve_source=False)
        
        # Optionally try downloading
        if self._should_auto_download(structure_id):
            return self.download_structure(structure_id, preserve_source=False)
    
    return None
```

### Ingest CIF → PKL
```python
def ingest_cif(self, structure_id: str, *, 
               cif_path: Optional[Path] = None,
               preserve_source: bool = False,
               source: Optional[str] = None,
               source_url: Optional[str] = None) -> pd.DataFrame:
    """Parse CIF once, create canonical PKL"""
    
    # 1. Find or use provided CIF
    if cif_path is None:
        # Try registry
        cif_info = self.entity_registry.find_entity(structure_id, "MMCIF")
        if cif_info:
            cif_path = self._resolve_path(cif_info.file_path)
            source = source or cif_info.metadata.get('source')
        else:
            # Try standard location
            cif_path = self.path_structure_dir / f"{structure_id}.cif"
    
    if not cif_path.exists():
        raise FileNotFoundError(f"No CIF file found for {structure_id}")
    
    # 2. Parse CIF
    self.logger.info(f"Ingesting CIF for {structure_id}")
    df_raw = self._parse_cif_file(cif_path)
    
    # 3. Canonicalize
    df = self._ensure_canonical(df_raw, structure_id)
    content_hash = compute_content_hash(df)
    
    # 4. Save canonical PKL
    pkl_metadata = {
        "schema_version": "v2",
        "io_version": "1",
        "content_hash": content_hash,
        "role": "canonical",
        "source": source or "unknown",
        "source_url": source_url,
        "ingested_at": datetime.utcnow().isoformat(),
        "atom_count": len(df),
        "chains": sorted(df['auth_chain_id'].unique().tolist())
    }
    
    self.save_entity(structure_id, df, format='pkl', metadata=pkl_metadata)
    
    # 5. Optionally preserve source
    if preserve_source and not self._is_source_preserved(structure_id):
        source_path = self.path_structure_dir / "raw" / f"{structure_id}.cif"
        source_path.parent.mkdir(exist_ok=True)
        
        import shutil
        shutil.copy2(cif_path, source_path)
        
        # Register source CIF
        self.entity_registry.register_entity(
            name=structure_id,
            format_type="MMCIF",
            file_path=str(source_path.relative_to(self.paths.data_root)),
            metadata={
                "role": "source",
                "source": source or "unknown",
                "preserved_at": datetime.utcnow().isoformat()
            }
        )
        
        # Add relationship
        self.entity_registry.add_relationship(
            source_name=structure_id,
            target_name=structure_id,
            relationship_type="derived_from",
            metadata={"from_format": "MMCIF", "to_format": "STRUCTURE_PKL"}
        )
    
    # 6. Update storage
    self._set_frame(structure_id, df)
    
    return df
```

### Save Canonical PKL
```python
def save_entity(self, structure_id: str, df: pd.DataFrame, *,
                format: str = 'pkl', 
                metadata: Optional[Dict] = None,
                versioning: str = 'auto',
                version_namer: Optional[Callable] = None) -> str:
    """Save canonical DataFrame with versioning support"""
    
    if format != 'pkl':
        raise ValueError("Use export_entity() for non-PKL formats")
    
    # 1. Canonicalize
    df = self._ensure_canonical(df, structure_id)
    new_hash = compute_content_hash(df)
    
    # 2. Check existing PKL
    existing = self.entity_registry.find_entity(structure_id, "STRUCTURE_PKL")
    if existing:
        old_hash = existing.metadata.get('content_hash', '')
        
        if old_hash == new_hash:
            # Content unchanged - overwrite for schema/annotation updates
            self.logger.info(f"Updating {structure_id} PKL (content unchanged)")
            resolved_id = structure_id
        else:
            # Content changed - check versioning policy
            if versioning == 'never':
                self.logger.warning(f"Overwriting {structure_id} with changed content")
                resolved_id = structure_id
            else:
                # Create new version
                if version_namer:
                    resolved_id = version_namer(structure_id)
                else:
                    resolved_id = f"{structure_id}@{datetime.now():%Y-%m-%d}"
                
                self.logger.info(f"Creating new version: {resolved_id}")
    else:
        resolved_id = structure_id
    
    # 3. Write PKL
    pkl_path = self.path_cache_dir / f"{resolved_id}.pkl"
    df.to_pickle(pkl_path)
    
    # 4. Register
    pkl_metadata = {
        "schema_version": "v2",
        "io_version": "1", 
        "content_hash": new_hash,
        "role": "canonical",
        **(metadata or {})
    }
    
    self.entity_registry.register_entity(
        name=resolved_id,
        format_type="STRUCTURE_PKL",
        file_path=str(pkl_path.relative_to(self.paths.data_root)),
        metadata=pkl_metadata
    )
    
    # 5. Add version relationship if needed
    if resolved_id != structure_id and existing:
        self.entity_registry.add_relationship(
            source_name=resolved_id,
            target_name=structure_id,
            relationship_type="version_of",
            metadata={"created_at": datetime.utcnow().isoformat()}
        )
        
        # Update alias to point to latest
        self.entity_registry.add_alias(structure_id, resolved_id)
    
    return resolved_id
```

### Export to CIF
```python
def export_entity(self, structure_id: str, *,
                  format: str = 'cif',
                  out_path: Optional[Path] = None,
                  overwrite: bool = True) -> Path:
    """Export canonical PKL to other formats"""
    
    if format != 'cif':
        raise ValueError(f"Unsupported export format: {format}")
    
    # 1. Load canonical PKL
    pkl_info = self.entity_registry.find_entity(structure_id, "STRUCTURE_PKL")
    if not pkl_info:
        raise ValueError(f"No PKL found for {structure_id}")
    
    pkl_path = self._resolve_path(pkl_info.file_path)
    df = pd.read_pickle(pkl_path)
    
    # 2. Determine output path
    if out_path is None:
        out_path = self.path_structure_dir / f"{structure_id}.cif"
    
    if out_path.exists() and not overwrite:
        raise FileExistsError(f"Output file exists: {out_path}")
    
    # 3. Export
    self.logger.info(f"Exporting {structure_id} to CIF")
    df_reset = df.reset_index()
    cif_utils.dataframe_to_cif(df_reset, str(out_path))
    
    # 4. Register exported CIF
    self.entity_registry.register_entity(
        name=structure_id,
        format_type="MMCIF", 
        file_path=str(out_path.relative_to(self.paths.data_root)),
        metadata={
            "role": "export",
            "exported_at": datetime.utcnow().isoformat(),
            "exported_from_hash": pkl_info.metadata.get('content_hash')
        }
    )
    
    # 5. Add relationship
    self.entity_registry.add_relationship(
        source_name=structure_id,
        target_name=structure_id,
        relationship_type="derived_from",
        metadata={"from_format": "STRUCTURE_PKL", "to_format": "MMCIF"}
    )
    
    return out_path
```

### Download & Ingest
```python
def download_structure(self, structure_id: str, *,
                      source: str = 'rcsb',
                      preserve_source: bool = False,
                      auto_register: bool = True) -> Optional[pd.DataFrame]:
    """Download CIF and ingest to PKL"""
    
    # 1. Download to temp
    self.logger.info(f"Downloading {structure_id} from {source}")
    
    if source == 'rcsb':
        url = f"https://files.rcsb.org/download/{structure_id}.cif"
    elif source == 'alphafold':
        url = f"https://alphafold.ebi.ac.uk/files/AF-{structure_id}-F1-model_v4.cif"
    else:
        raise ValueError(f"Unknown source: {source}")
    
    # Download
    import requests
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    
    # Save to temp
    temp_path = self.temp_dir / f"{structure_id}_{int(time.time())}.cif"
    temp_path.write_bytes(response.content)
    
    try:
        # 2. Ingest
        df = self.ingest_cif(
            structure_id,
            cif_path=temp_path,
            preserve_source=preserve_source,
            source=source,
            source_url=url
        )
        
        return df
        
    finally:
        # Clean up temp
        if temp_path.exists():
            temp_path.unlink()
```

---

## Phase 2 — Storage Management

### Storage Implementation
```python
class StructureProcessor(BaseProcessor):
    def __init__(self):
        super().__init__()
        
        # Storage
        self.frames: Dict[str, pd.DataFrame] = {}
        self._data: Optional[pd.DataFrame] = None
        self._dirty: bool = False
        
    def _set_frame(self, structure_id: str, df: pd.DataFrame):
        """Update single frame and mark dirty"""
        self.frames[structure_id] = df
        self._dirty = True
        
    def _remove_frame(self, structure_id: str):
        """Remove frame and mark dirty"""
        if structure_id in self.frames:
            del self.frames[structure_id]
            self._dirty = True
    
    @property
    def data(self) -> Optional[pd.DataFrame]:
        """Lazy stacked view"""
        if self._dirty:
            self._rebuild_data()
        return self._data
    
    def _rebuild_data(self):
        """Rebuild stacked DataFrame"""
        if not self.frames:
            self._data = None
        else:
            self._data = pd.concat(self.frames.values()).sort_index()
        self._dirty = False
    
    @property
    def structure_ids(self) -> List[str]:
        """List of loaded structures"""
        return sorted(self.frames.keys())
```

---

## Phase 3 — Dataset Operations

### Dataset Loading (No Monolithic PKLs)
```python
def load_dataset(self, dataset_name: str, 
                return_format: str = 'stacked') -> Union[pd.DataFrame, Dict[str, pd.DataFrame]]:
    """Load dataset members from individual PKLs"""
    
    # Get dataset info
    dataset = self.get_dataset(dataset_name)
    structure_ids = dataset.get('entities', [])
    
    # Load each structure's PKL
    frames = {}
    for sid in structure_ids:
        df = self.load_entity(sid, auto_ingest=False)  # PKL must exist
        if df is not None:
            frames[sid] = df
        else:
            self.logger.warning(f"Missing PKL for dataset member: {sid}")
    
    # Bulk update storage
    self.frames.update(frames)
    self._dirty = True
    
    # Return format
    if return_format == 'dict':
        return frames
    else:
        return self.data  # Triggers lazy rebuild
```

### Dataset Saving (Ensure PKLs)
```python
def save_dataset(self, dataset_name: str, structure_ids: List[str],
                metadata: Optional[Dict] = None):
    """Ensure each member has PKL, save logical dataset"""
    
    # Ensure each structure has a PKL
    missing = []
    for sid in structure_ids:
        pkl_info = self.entity_registry.find_entity(sid, "STRUCTURE_PKL")
        if not pkl_info:
            if sid in self.frames:
                # Save from memory
                self.save_entity(sid, self.frames[sid], format='pkl')
            else:
                missing.append(sid)
    
    if missing:
        raise ValueError(f"Cannot save dataset, missing structures: {missing}")
    
    # Create logical dataset (no PKL file)
    self.create_dataset(dataset_name, structure_ids, metadata)
```

---

## Migration & Compatibility

### Legacy Support
```python
# Deprecated methods with warnings
def load_structure(self, pdb_id: str, **kwargs):
    warnings.warn("load_structure is deprecated, use load_entity", 
                  DeprecationWarning, stacklevel=2)
    return self.load_entity(pdb_id, **kwargs)

def save_structure(self, pdb_id: str, df: pd.DataFrame, **kwargs):
    warnings.warn("save_structure is deprecated, use save_entity or export_entity",
                  DeprecationWarning, stacklevel=2)
    fmt = kwargs.pop('format', 'pkl')
    if fmt == 'cif':
        self.save_entity(pdb_id, df, format='pkl')
        return self.export_entity(pdb_id, format='cif')
    else:
        return self.save_entity(pdb_id, df, **kwargs)

@property
def pdb_ids(self):
    raise AttributeError(
        "pdb_ids is removed. Use structure_ids instead.\n"
        "If you need PDB codes, they are now generic structure_ids"
    )
```

### Cache Migration
```python
def migrate_cache_to_canonical(self, dry_run: bool = True):
    """One-time migration of existing cache"""
    
    for pkl_file in self.path_cache_dir.glob("*.pkl"):
        try:
            # Load existing
            df = pd.read_pickle(pkl_file)
            structure_id = pkl_file.stem
            
            # Check if needs migration
            needs_migration = (
                'pdb_id' in df.columns or
                not isinstance(df.index, pd.MultiIndex) or
                df.index.names != ['structure_id', 'atom_id']
            )
            
            if needs_migration:
                self.logger.info(f"Migrating {structure_id}")
                if not dry_run:
                    # Canonicalize and save
                    df_canonical = self._ensure_canonical(df, structure_id)
                    self.save_entity(structure_id, df_canonical, 
                                   format='pkl', versioning='never')
        
        except Exception as e:
            self.logger.error(f"Failed to migrate {pkl_file}: {e}")
```

---

## Summary of Changes

### What's Different from Original Plan

1. **Format-Specific Registry Keys**: Use "STRUCTURE_PKL" and "MMCIF" instead of generic "structure"
2. **Content Hashing**: Implement `compute_content_hash()` for versioning decisions
3. **Clear CIF Roles**: `source` (preserved original) vs `export` (generated from PKL)
4. **Versioning Policy**: Automatic versioning when content changes
5. **No CIF Loading**: After ingest, always load from PKL

### What Stays the Same

1. **BaseProcessor Architecture**: No separate IO classes
2. **Entity Registry Usage**: Leverages existing system with format-specific keys
3. **ProtosPaths**: All path management unchanged
4. **Storage Design**: Dict-based frames with lazy stacked view
5. **Dataset Semantics**: Logical collections, no monolithic PKLs

### Implementation Priority

1. **Phase 0**: Implement `_ensure_canonical()` and `compute_content_hash()` (1 day)
2. **Phase 1**: Core API - `load_entity()`, `ingest_cif()`, `save_entity()`, `export_entity()` (3 days)
3. **Phase 2**: Storage management with lazy rebuilding (1 day)
4. **Phase 3**: Dataset operations (1 day)
5. **Testing**: Versioning scenarios, format relationships (2 days)
6. **Migration**: Cache migration tool (1 day)

**Total**: ~9 days for complete implementation