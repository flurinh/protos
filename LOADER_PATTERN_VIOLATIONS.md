# Loader Pattern Violations in Protos

## Overview

The current loaders (ChEMBL, GPCR, UniProt) violate core Protos principles by managing paths, creating directories, and maintaining their own registries. This document outlines the correct pattern for loaders.

## What Loaders Should Be

Loaders should be **pure utility modules** that:
1. Fetch data from external sources
2. Parse/transform data formats
3. Return data to processors
4. Have NO filesystem operations
5. Have NO path management
6. Have NO entity registry access

## Current Violations

### ChEMBL Loader (chembl_loader.py)

**VIOLATIONS**:
```python
# ❌ Takes data_root parameter
def __init__(self, data_root: Optional[str] = None, ...):

# ❌ Creates its own ProtosPaths
self.paths = ProtosPaths(data_root=data_root)

# ❌ Manages directories
self.ligand_dir = Path(self.paths.get_processor_path('ligand'))
self.sdf_dir = self.ligand_dir / 'sdf'
self.chembl_dir = self.ligand_dir / 'chembl'

# ❌ Creates directories
for dir_path in [...]:
    dir_path.mkdir(parents=True, exist_ok=True)

# ❌ Creates its own EntityRegistry
self.entity_registry = EntityRegistry(self.paths)

# ❌ Saves files directly
with open(sdf_path, 'w') as f:
    f.write(sdf_content)

# ❌ Registers entities
self.entity_registry.register_entity(...)
```

### GPCR Loader (gpcrdb_loader.py)

**VIOLATIONS**:
```python
# ❌ Takes explicit paths
def __init__(self, path='data/', fileformat='pdb'):
    self.path_pdb = path + 'pdb/'

# ❌ Checks directory existence
if not Path(self.path_pdb).exists():
    raise OSError(f"Directory {self.path_pdb} not found")

# ❌ Downloads files directly
download_pdb(url, folder=self.path_pdb, fileformat=self.fileformat)
```

### UniProt Loader (uniprot_loader.py)

**VIOLATIONS**:
```python
# ❌ Creates directories
os.makedirs(self.metadata_dir, exist_ok=True)
os.makedirs(self.fasta_dir, exist_ok=True)

# ❌ Manages paths
self.fasta_dir = self.paths.get_subdir_path('sequence', 'fasta_dir')
self.metadata_dir = self.paths.get_subdir_path('sequence', 'metadata_dir')
```

## Correct Pattern

### Loader as Pure Utility (uniprot_utils.py)

**CORRECT PATTERN**:
```python
# ✅ Pure functions with no side effects
def submit_id_mapping(from_db, to_db, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    request.raise_for_status()
    return request.json()["jobId"]

# ✅ Returns data, doesn't save it
def get_id_mapping_results(job_id):
    # ... fetch and parse data ...
    return parsed_results
```

### Refactored ChEMBL Loader Example

**CORRECT PATTERN**:
```python
# ✅ Pure utility functions
def map_protein_to_chembl_target(protein_id: str) -> Optional[str]:
    """Map protein ID to ChEMBL target. No side effects."""
    # ... mapping logic ...
    return chembl_target_id

def query_protein_ligands(protein_id: str, ...) -> List[Dict]:
    """Query ChEMBL for ligands. Returns data only."""
    # ... query logic ...
    return ligand_data  # Just return data!

# ✅ NO class needed - just functions
# ✅ NO path operations
# ✅ NO file saving
# ✅ NO entity registration
```

## Usage Pattern

### How Processors Should Use Loaders

```python
class LigandProcessor(BaseProcessor):
    def get_protein_ligands(self, protein_id: str) -> List[Dict]:
        # 1. Use loader to fetch data
        ligand_data = chembl_loader.query_protein_ligands(protein_id)
        
        # 2. Process and save through processor's entity system
        for ligand in ligand_data:
            self.save_entity(ligand['smiles'], ligand)
        
        return ligand_data
```

## Key Principles

1. **Separation of Concerns**
   - Loaders: Fetch and parse external data
   - Processors: Manage storage and entities
   - ProtosPaths: Handle all path operations

2. **No Side Effects**
   - Loaders should be pure functions
   - No filesystem operations
   - No state management

3. **Data Flow**
   ```
   External API → Loader (fetch/parse) → Processor (store/manage) → Entity Registry
   ```

4. **Testing**
   - Loaders can be tested without filesystem
   - Just test data transformation logic
   - Mock external API calls

## Migration Strategy

1. **Phase 1**: Create pure utility versions (e.g., `chembl_loader_utils.py`)
2. **Phase 2**: Update processors to use utilities
3. **Phase 3**: Deprecate class-based loaders
4. **Phase 4**: Remove old loaders

## Benefits

1. **Cleaner Architecture**: Clear separation of responsibilities
2. **Easier Testing**: Pure functions are easier to test
3. **Flexibility**: Processors control how data is stored
4. **Consistency**: All filesystem operations go through ProtosPaths
5. **Reusability**: Utilities can be used by multiple processors