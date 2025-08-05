# Ligand Processor: Protos Principle Violations

## Overview

The current ligand processor violates multiple core Protos principles as defined in DATA_MANAGEMENT_UNIFIED.md. This document details these violations and shows the correct patterns from other processors.

## Major Violations

### 1. Direct Filesystem Management

**VIOLATION**: The ligand processor creates and manages directories directly:

```python
# ❌ WRONG - ligand_processor.py lines 77-87
self.data_path = Path(self.paths.get_processor_path(self.processor_type))
self.sdf_dir = self.data_path / 'sdf'
self.tables_dir = self.data_path / 'tables'
self.cache_dir = self.data_path / 'cache'
self.models_dir = self.data_path / 'models'
self.pockets_dir = self.data_path / 'pockets'

# Ensure directories exist
for dir_path in [self.sdf_dir, self.tables_dir, self.cache_dir, 
                self.models_dir, self.pockets_dir]:
    dir_path.mkdir(parents=True, exist_ok=True)
```

**CORRECT PATTERN** (from sequence_processor.py):
```python
# ✅ CORRECT - Use properties that delegate to BaseProcessor
@property
def path_fasta_dir(self):
    """Get path to FASTA files directory."""
    return self.get_subdirectory_path('fasta_dir')
```

### 2. Explicit Path Construction

**VIOLATION**: The ligand processor constructs file paths explicitly:

```python
# ❌ WRONG - ligand_processor.py line 165
sdf_path = self.sdf_dir / f"{safe_filename}.sdf"

# ❌ WRONG - ligand_processor.py line 486
file_path = Path(self.paths.data_root) / entity_info.file_path
```

**CORRECT PATTERN**: Let the entity registry and ProtosPaths handle all path operations.

### 3. Direct File I/O Without Registry

**VIOLATION**: The ligand processor reads/writes files directly:

```python
# ❌ WRONG - ligand_processor.py lines 175-176
with open(sdf_path, 'w') as f:
    f.write(sdf_content)
```

**CORRECT PATTERN**: All file operations should go through the entity registry system.

### 4. Creating ChEMBL Loader with Explicit Paths

**VIOLATION**: The ligand processor passes explicit paths to the ChEMBL loader:

```python
# ❌ WRONG - ligand_processor.py line 93
self.chembl_loader = ChEMBLDL(data_root=self.paths.data_root)
```

**CORRECT PATTERN**: Components should accept ProtosPaths instances, not raw paths.

### 5. Cache Management

**VIOLATION**: The ligand processor implements its own caching logic:

```python
# ❌ WRONG - Managing cache files directly
self.cache_dir = self.data_path / 'cache'
```

**CORRECT PATTERN**: Caching should be handled at the entity registry level or use BaseProcessor's facilities.

## Core Principles Violated

From DATA_MANAGEMENT_UNIFIED.md:

1. **"ProtosPaths is the ONLY path management system"** - Violated by creating own path attributes
2. **"Zero configuration"** - Violated by explicit directory creation
3. **"No component ever specifies or manages paths directly"** - Violated throughout
4. **"Processors are the primary interface"** - Violated by exposing filesystem details

## Correct Implementation Pattern

The refactored ligand processor follows these principles:

1. **NO directory creation** - ProtosPaths handles this
2. **NO path attributes** - Use properties or methods that delegate to BaseProcessor
3. **NO direct file I/O** - Everything goes through entity registry
4. **Focus on data operations** - Not filesystem management

## Example: Correct Entity Save

```python
# ✅ CORRECT - From refactored version
def save_entity(self, name: str, data: Union[str, Dict], metadata: Optional[Dict] = None):
    # Validate data
    is_valid, canonical_smiles = validate_smiles(smiles)
    
    # Build relative path (not absolute!)
    safe_filename = sanitize_smiles_filename(canonical_smiles)
    relative_path = f"{self.processor_type}/{self.sdf_subdir}/{safe_filename}.sdf"
    
    # Register with entity registry - it handles everything
    self.entity_registry.register_entity(
        name=canonical_smiles,
        format_type=self.processor_type,
        file_path=relative_path,  # Relative path only!
        metadata=metadata
    )
```

## Summary

The ligand processor needs a complete refactor to:
1. Remove ALL filesystem operations
2. Remove ALL path management
3. Use entity registry for ALL data storage/retrieval
4. Follow the patterns established by sequence_processor.py

The refactored version (ligand_processor_refactored.py) demonstrates the correct implementation.