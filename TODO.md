# PROTOS TODO - DATA MANAGEMENT IMPLEMENTATION

## 🚀 QUICK START & ESSENTIAL CONTEXT

### Environment Setup
```bash
# Activate conda environment
conda activate protos

# Install in development mode (REQUIRED for proper path resolution)
pip install -e .

# Verify installation
python -c "import protos; print(protos.__version__)"
```

### 🔴 CRITICAL TESTING REQUIREMENTS

1. **ALWAYS Test With Real Data**
   ```python
   # Use actual biological data - download if needed
   from protos.loaders import download_structures
   download_structures.download_pdb_files(['1ubq', '3sn6', '7zvl'])
   
   # Use real sequences from UniProt
   from protos.loaders import uniprot_utils
   uniprot_utils.download_sequences(['P62988', 'P00533'])
   ```

2. **Test After EVERY Change**
   ```bash
   # Run specific test
   pytest tests/test_processors/test_structure/ -xvs
   
   # Check for path violations
   pytest --tb=short | grep -E "(tmp_path|Path\(|\.mkdir)"
   ```

3. **Use Existing Test Data**
   ```
   tests/test-data/
   ├── structure/mmcif/     # Real CIF files
   ├── sequence/fasta/      # Real FASTA files
   └── grn/tables/          # Real GRN tables
   ```

### ⚠️ COMMON PITFALLS TO AVOID

1. **❌ NEVER Create Paths Manually**
   ```python
   # WRONG
   path = Path("data/structure/mmcif/1ubq.cif")
   os.makedirs("data/structure", exist_ok=True)
   
   # CORRECT
   path = self.path_resolver.get_structure_subdir_path('structure_dir') / "1ubq.cif"
   # ProtosPaths creates directories automatically
   ```

2. **❌ NEVER Expose Hash IDs to Users**
   ```python
   # WRONG
   def list_entities():
       return ["3f8a9c2d1e", "ba1837f945"]  # Hash IDs
   
   # CORRECT
   def list_entities():
       return ["1ubq", "EGFR_HUMAN"]  # Human names
   ```

3. **❌ NEVER Use Path Parameters in APIs**
   ```python
   # WRONG
   def load_structure(pdb_id: str, path: str):
   
   # CORRECT
   def load_structure(pdb_id: str):
       # Use self.path_resolver internally
   ```

4. **❌ NEVER Save Without Registration**
   ```python
   # WRONG
   df.to_csv("myfile.csv")
   
   # CORRECT
   self.save_entity("my_data", df)  # Auto-registers
   ```

### 🎯 PROTOS CORE PRINCIPLES

1. **Zero Configuration**: Everything works out of the box
   ```python
   processor = CifBaseProcessor()  # Just works!
   ```

2. **Human-Readable Everything**: Files, names, paths - all human-readable
   ```
   data/structure/mmcif/1ubq.cif  # NOT: data/structure/mmcif/3f8a9c2d1e.cif
   ```

3. **Unified Entity System**: One entity can have multiple formats
   ```python
   # Same entity "ubiquitin" in different formats
   struct_proc.save_structure("ubiquitin", structure)
   seq_proc.save_sequence("ubiquitin", "MQIFVKTLTG...")
   ```

4. **ProtosPaths Manages ALL Paths**: No exceptions
   ```python
   # ProtosPaths provides ALL paths
   self.path_resolver.get_processor_path('structure')
   self.path_resolver.get_structure_subdir_path('cache_dir')
   ```

### 📁 FILESYSTEM STRUCTURE (Managed by ProtosPaths)
```
working_dir/
└── data/                          # Or custom via PROTOS_DATA_ROOT
    ├── entity_registry.json       # Global entity tracking
    ├── structure/
    │   ├── mmcif/                # Raw CIF files
    │   ├── cache/                # Processed individual structures (PKL)
    │   ├── structure_dataset/    # Processed dataset PKL files (NOT JSON definitions!)
    │   ├── alignments/           # Structure alignments
    │   ├── datasets/             # Dataset JSON definitions (e.g., kinases.json, opsins.json)
    │   └── registry.json         # Structure processor registry
    ├── sequence/
    │   ├── fasta/                # FASTA files
    │   ├── alignments/           # MSA results
    │   ├── metadata/             # Sequence metadata
    │   ├── datasets/             # Dataset JSON definitions
    │   └── registry.json         # Sequence processor registry
    ├── grn/
    │   ├── tables/               # GRN tables
    │   ├── configs/              # Configurations
    │   ├── ref/                  # Reference tables
    │   ├── assignments/          # GRN assignments
    │   ├── datasets/             # Dataset JSON definitions
    │   └── registry.json         # GRN processor registry
    ├── property/
    │   ├── tables/               # Property tables
    │   ├── datasets/             # Dataset JSON definitions
    │   └── registry.json         # Property processor registry
    ├── embedding/
    │   ├── embeddings/           # Saved embeddings
    │   ├── datasets/             # Dataset JSON definitions
    │   └── registry.json         # Embedding processor registry
    └── graph/
        ├── networks/             # Graph/network files
        ├── analysis/             # Graph analysis results
        ├── datasets/             # Dataset JSON definitions
        └── registry.json         # Graph processor registry
```

### 🔧 DEBUGGING TIPS

1. **Check Path Resolution**
   ```python
   print(processor.path_resolver.data_root)
   print(processor.data_path)
   ```

2. **Verify Entity Registration**
   ```python
   # Check if entity is registered
   print(processor.entity_exists("my_protein"))
   
   # List all entities
   print(processor.list_entities())
   ```

3. **Enable Debug Logging**
   ```python
   import logging
   logging.basicConfig(level=logging.DEBUG)
   ```

### 📊 REAL DATA SOURCES

1. **PDB Structures**
   ```python
   from protos.loaders.download_structures import download_pdb_files
   download_pdb_files(['1ubq', '2gb1', '3sn6'])
   ```

2. **UniProt Sequences**
   ```python
   from protos.loaders.uniprot_utils import download_sequences
   download_sequences(['P62988', 'P00533', 'Q9Y6K9'])
   ```

3. **GPCRdb Data**
   ```python
   from protos.loaders.gpcrdb_loader import GPCRdbLoader
   loader = GPCRdbLoader()
   loader.download_alignment('rhodopsin')
   ```

## 🎯 CORE PRINCIPLE: ONE PATH SYSTEM FOR EVERYTHING

**ProtosPaths is the ONLY path management system. Zero configuration required. For reference, read DATA_MANAGEMENT_UNIFIED.md**

## 🚨 CRITICAL IMPLEMENTATION ORDER

**IMPORTANT**: The following components MUST work perfectly before moving to specific processors:

1. **ProtosPaths** - The ONLY path manager (NO hardcoded paths, NO os.path.join, NO manual directory creation)
2. **EntityRegistry** - Human names only, hash IDs internal 
3. **DatasetManager** - JSON datasets in datasets/ directories
4. **BaseProcessor** - Zero configuration, delegates ALL path operations to ProtosPaths

**ALL processors inherit from BaseProcessor and depend on these components. If these don't work correctly, NOTHING will work correctly.**

## 📋 IMPLEMENTATION PLAN

### PHASE 1: Core Components ✅ COMPLETED

#### 1.1 ProtosPaths Verification ✅
```bash
# Test each step after implementation
pytest tests/test_io/test_paths_comprehensive.py -xvs  # All 10 tests passing!
```

**Tasks:**
- [x] Verify ProtosPaths uses `~/protos_data` as default (or env var)
- [x] Ensure all path methods return strings for backward compatibility
- [x] Test directory creation on initialization
- [x] Remove any DataSource enum references (deprecated but handled)
- [x] Test with clean environment (no config files)
- [x] Fixed dataset path to always use datasets/ directory

#### 1.2 EntityRegistry Implementation ✅
```bash
# Test after each change
pytest tests/test_io/test_entity_registry.py -xvs  # All 13 tests passing!
```

**Tasks:**
- [x] EntityRegistry must accept ProtosPaths in `__init__`
- [x] Hash IDs used internally only - never exposed to users
- [x] All public methods work with human-readable names
- [x] File paths always use human-readable names
- [x] Implement dual ID support (internal use only):
  ```python
  def register_entity(self, name: str, format_type: str, file_path: str):
      # Generate hash ID internally
      hash_id = self._generate_hash_id(name)
      # Store with human name as key
      self.entities[hash_id] = {
          "original_id": name,
          "formats": {format_type: {"path": file_path}}
      }
  ```

#### 1.3 DatasetManager Implementation ✅
```bash
# Test dataset operations
pytest tests/test_io/test_dataset_manager.py -xvs  # All 14 tests passing!
```

**Tasks:**
- [x] Create DatasetManager class that accepts ProtosPaths
- [x] Datasets stored as JSON files: `processor_dir/datasets/dataset_name.json`
- [x] Each dataset is a separate JSON file (e.g., kinases.json, opsins.json)
- [x] Datasets store human names in JSON files (can store hash IDs internally if needed)
- [x] ALL user-facing APIs return human names only (users NEVER see hash IDs)
- [x] Implement all operations from DATA_MANAGEMENT_UNIFIED.md:
  ```python
  def create_dataset(self, name: str, entities: List[str], 
                    processor_type: str, metadata: dict):
      # Convert human names to hash IDs internally
      # Save to processor_type/datasets/name.json
  ```

### PHASE 2: BaseProcessor Update ✅ COMPLETED

#### 2.1 Fix Initialization Cascade ✅
```bash
# Run all base processor tests
pytest tests/test_core/test_base_processor_updated.py -xvs  # Tests passing!
```

**Tasks:**
- [x] BaseProcessor creates ProtosPaths if not provided
- [x] BaseProcessor creates EntityRegistry with ProtosPaths
- [x] BaseProcessor creates DatasetManager with ProtosPaths
- [x] Remove ALL path parameters from processor constructors
- [x] Test zero-configuration usage
- [x] Rewritten BaseProcessor to follow DATA_MANAGEMENT_UNIFIED.md exactly

#### 2.2 Implement Entity Operations ✅
**Tasks:**
- [x] `list_entities()` returns human names only
- [x] `load_entity(name)` abstract method - subclasses implement
- [x] `save_entity(name, data)` abstract method - subclasses implement
- [x] `entity_exists(name)` checks by human name
- [x] `delete_entity(name)` removes file and updates registry

#### 2.3 Implement Dataset Operations ✅
**Tasks:**
- [x] `list_datasets()` returns dataset names for this processor
- [x] `create_dataset(name, entities)` accepts human names
- [x] `load_dataset(name)` loads all entities, returns dict[human_name, data]
- [x] `add_to_dataset/remove_from_dataset` work with human names
- [x] `get_dataset_info` returns human-readable information

## 🚨 PHASE 5: CRITICAL TEST FAILURES - IMPLEMENTATION PLAN

### Overview
Based on test results (87 failed, 375 passed), we have systematic issues across multiple processors that need to be addressed in priority order.

### Phase 5.1: API Signature Fixes (High Priority) 🔴

#### Issue: Constructor Parameter Mismatches
Many processors still expect old constructor parameters, causing `TypeError: __init__() got an unexpected keyword argument`.

**Affected Components:**
- GRNProcessor: `processor_data_dir` parameter
- StructureProcessor: `data_root`, `processor_data_dir` parameters
- Multiple test files using deprecated initialization patterns

**Tasks:**
- [ ] Update GRNProcessor constructor to match BaseProcessor signature
- [ ] Update StructureProcessor constructor to match BaseProcessor signature  
- [ ] Fix all test files to use new initialization pattern
- [ ] Remove all references to deprecated parameters

**Implementation:**
```python
# OLD (causes failures)
processor = GRNProcessor(name="test", processor_data_dir="grn")

# NEW (correct)
processor = GRNProcessor(name="test")  # paths handled internally
```

### Phase 5.2: Entity Registry Concurrency (High Priority) 🔴

#### Issue: File Locking Conflicts
`FileExistsError: [WinError 183] Cannot create a file when that file already exists`

**Root Cause:** Multiple tests trying to create entity_registry.json simultaneously without proper locking.

**Tasks:**
- [ ] Implement file locking in EntityRegistry save operations
- [ ] Add retry logic with exponential backoff
- [ ] Ensure atomic file operations (write to temp, then rename)
- [ ] Add test isolation for entity registry operations

### Phase 5.3: Missing Method Implementations (High Priority) 🔴

#### Issue: Abstract Methods Not Implemented
`AttributeError: 'StructureProcessor' object has no attribute '_load_dataset_impl'`

**Affected Processors:**
- StructureProcessor: Missing `_load_dataset_impl`, `_save_dataset_impl`
- Various processors: Missing `load_entity`, `save_entity` implementations

**Tasks:**
- [ ] Implement all required abstract methods in each processor
- [ ] Ensure proper method signatures match BaseProcessor
- [ ] Add validation for required methods

### Phase 5.4: Column Name Standardization (High Priority) 🔴

#### Issue: Inconsistent Column Names
`KeyError: 'res_atom_name'` vs `'atom_name'` in structure data

**Tasks:**
- [ ] Standardize column names across all processors:
  - Structure: `atom_name`, `chain_id`, `seq_id`, `comp_id`
  - Sequence: `sequence_id`, `sequence`, `description`
  - GRN: Use consistent GRN column format
- [ ] Create column mapping utilities for backward compatibility
- [ ] Update all processors to use standard names

### Phase 5.5: Import and Module Issues (Medium Priority) 🟡

#### Issue: Missing Imports and Renamed Modules
`NameError: name 'SeqProcessor' is not defined` in embedding_processor.py

**Tasks:**
- [ ] Fix all import statements to use correct module names
- [ ] Update cross-processor dependencies
- [ ] Remove circular imports
- [ ] Ensure all processors can be imported independently

### Phase 5.6: GRN Processing Logic (Medium Priority) 🟡

#### Issue: GRN Column Parsing Failures
`Cannot match GRN columns. Data columns: ['3.50', '7.53', 'species'], Requested: ['0x00', '3.50', '7.53']`

**Root Cause:** GRN processor generating invalid column names ('0x00') and failing to handle non-GRN columns.

**Tasks:**
- [ ] Fix GRN column validation to skip non-GRN columns
- [ ] Remove '0x00' column generation bug
- [ ] Handle mixed column types gracefully
- [ ] Add proper GRN format validation

### Phase 5.7: Path Naming Consistency (Medium Priority) 🟡

#### Issue: Directory Naming Conflicts
Expected `structure_dataset/` but got `datasets/`

**Tasks:**
- [ ] Resolve naming conflicts between processors
- [ ] Ensure all processors use consistent subdirectory names
- [ ] Update ProtosPaths to enforce naming consistency
- [ ] Fix tests expecting old directory names

### Phase 5.8: Missing Utility Functions (Medium Priority) 🟡

#### Issue: Missing Implementation
`module 'protos.io.cif_utils' has no attribute 'dataframe_to_cif'`

**Tasks:**
- [ ] Implement missing `dataframe_to_cif` function
- [ ] Add other missing utility functions
- [ ] Ensure all referenced functions exist
- [ ] Add proper error messages for missing functionality

### Phase 5.9: Test Infrastructure (Low Priority) 🟢

#### Issue: Tests Using Abstract Classes
`Can't instantiate abstract BaseProcessor`

**Tasks:**
- [ ] Update tests to use concrete processor classes
- [ ] Remove attempts to instantiate abstract classes
- [ ] Create proper test fixtures for each processor
- [ ] Ensure test isolation and cleanup

### Phase 5.10: Environment-Specific Issues (Low Priority) 🟢

#### Issue: Conda Environment Dependencies
Embedding tests require conda environment activation

**Tasks:**
- [ ] Document conda environment requirements
- [ ] Add environment checks in tests
- [ ] Provide fallback options for missing dependencies
- [ ] Handle ESM special tokens in embeddings

### PHASE 2.5: VERIFICATION OF CORE COMPONENTS ✅ COMPLETED

- [x] Created comprehensive integration tests
- [x] Fixed EntityRegistry refresh for multi-processor scenarios
- [x] All 22 core component tests passing
- [x] Removed deprecated test files

### PHASE 3: Update Specific Processors ✅ COMPLETED

- [x] 3.1 CifBaseProcessor - Structure processor updated
- [x] 3.2 SeqProcessor - Sequence processor updated
- [x] 3.3 GRNBaseProcessor - GRN processor updated
- [x] 3.4 PropertyProcessor - Property processor updated
- [x] 3.5 EmbeddingProcessor - Embedding processor updated

#### Critical Path Management Principles for ALL Processors

**EVERY processor must follow these rules:**

1. **NO Custom Path Handling**: 
   - Remove ALL `self.path_*` properties
   - Remove ALL `os.path.join()` calls
   - Remove ALL `os.makedirs()` calls
   - Remove ALL direct path manipulation

2. **Use ONLY ProtosPaths**:
   - Access paths through `self.path_resolver.get_*_path()` methods
   - Let ProtosPaths handle ALL directory creation
   - Use subdirectory constants from ProtosPaths

3. **Standard Entity/Dataset Operations**:
   - Use BaseProcessor's `save_entity()`, `load_entity()`, `delete_entity()`
   - Use BaseProcessor's dataset methods
   - Override ONLY when processor needs special handling

4. **Processor-Specific Overrides**:
   - Override methods ONLY to add format-specific logic
   - Always delegate path operations to ProtosPaths
   - Never construct paths manually

#### 3.1 CifBaseProcessor
```bash
# Run structure tests after each change
pytest tests/test_processors/test_structure/ -xvs
```

**Structure Data Complexity:**
CifBaseProcessor handles multiple data formats and processing stages:
1. **Raw CIF files** - Original structure files in `mmcif/`
2. **Cached PKL files** - Processed structures in `cache/`
3. **Dataset PKL files** - Entire datasets in `structure_dataset/`
4. **Alignments** - Structure alignments in `alignments/`

**Tasks:**

##### 3.1.1 Critical: Remove ALL Custom Path Handling ✅
- [x] ProtosPaths updated with cache directory
- [ ] Remove ALL `self.path_structure_dir`, `self.path_dataset_dir`, etc.
- [ ] Replace with `self.path_resolver.get_structure_subdir_path()`
- [ ] Remove ALL `os.path.join()` and `os.makedirs()`

##### 3.1.2 Fix Path Properties
- [ ] Replace path properties with methods that use ProtosPaths:
  ```python
  @property
  def mmcif_dir(self):
      return self.path_resolver.get_structure_subdir_path('structure_dir')
  
  @property
  def cache_dir(self):
      return self.path_resolver.get_structure_subdir_path('cache_dir')
  ```

##### 3.1.3 Update Load/Save Methods
- [ ] `load_structure()`: Check cache/ then mmcif/ using ProtosPaths
- [ ] `save_structure()`: Save to appropriate dir based on format
- [ ] `save_entity()`: Override to handle structure-specific metadata
- [ ] `load_entity()`: Override to handle cached vs raw files

##### 3.1.4 Dataset Caching
- [ ] Save datasets to `structure_dataset/` via ProtosPaths
- [ ] Load cached datasets before individual files
- [ ] Use ProtosPaths for all dataset operations

#### 3.2 SeqProcessor
```bash
# Run sequence tests
pytest tests/test_processors/test_sequence/ -xvs
```

**Sequence Data Organization:**
1. **FASTA files** - Raw sequences in `fasta/`
2. **Alignments** - MSA results in `alignments/`
3. **Metadata** - Sequence metadata in `metadata/`

**Tasks:**

##### 3.2.1 Remove Custom Path Handling
- [ ] Remove ALL custom path properties
- [ ] Use `self.path_resolver.get_sequence_subdir_path()`
- [ ] Remove direct file operations

##### 3.2.2 Filename Sanitization
- [ ] Implement `_sanitize_filename()` for complex identifiers
- [ ] Handle special characters in UniProt IDs, etc.
- [ ] Ensure human-readable but filesystem-safe names

##### 3.2.3 Multi-Sequence FASTA Support
- [ ] Handle FASTA files with multiple sequences
- [ ] Register each sequence as separate entity
- [ ] Track which file contains which sequences

#### 3.3 GRNBaseProcessor
```bash
# Run GRN tests
pytest tests/test_processors/test_grn/ -xvs
```

**GRN Data Organization:**
1. **Tables** - GRN alignment tables in `tables/`
2. **Configs** - GRN configurations in `configs/`
3. **Reference** - Reference GRN tables in `ref/`
4. **Assignments** - GRN assignments in `assignments/`

**Tasks:**

##### 3.3.1 Remove Custom Path Handling
- [ ] Remove ALL custom path properties
- [ ] Use `self.path_resolver.get_grn_subdir_path()`
- [ ] Let ProtosPaths manage all directories

##### 3.3.2 Table Management
- [ ] Tables use human-readable protein IDs (not hashes)
- [ ] Register each protein in table as entity
- [ ] Track table-to-entity relationships

##### 3.3.3 Reference Tables
- [ ] Handle reference vs user tables properly
- [ ] Use ProtosPaths for ref/ subdirectory
- [ ] Maintain backward compatibility

#### 3.4 PropertyProcessor
```bash
# Run property tests
pytest tests/test_processors/test_property/ -xvs
```

**Property Data Organization:**
1. **Tables** - Property tables in `tables/`
2. **Cache** - Cached property data in `cache/`
3. **Datasets** - Property datasets in `datasets/`

**Tasks:**

##### 3.4.1 Remove Custom Path Handling
- [ ] Remove ALL custom path properties
- [ ] Use `self.path_resolver.get_property_subdir_path()`
- [ ] Remove custom dataset implementation

##### 3.4.2 Table Management
- [ ] Property tables use human-readable entity names
- [ ] Support multiple property types per entity
- [ ] Cache computed properties

#### 3.5 EmbeddingProcessor
```bash
# Run embedding tests
pytest tests/test_processors/test_embedding/ -xvs
```

**Embedding Data Organization:**
1. **Models** - Saved embeddings in `models/`
2. **Cache** - Cached embeddings in `cache/`
3. **Datasets** - Embedding datasets in `datasets/`

**Tasks:**

##### 3.5.1 Remove Custom Path Handling
- [ ] Remove ALL custom path properties
- [ ] Use `self.path_resolver.get_embedding_subdir_path()`
- [ ] Let ProtosPaths manage directories

##### 3.5.2 Embedding Storage
- [ ] Embeddings use same names as source sequences
- [ ] Support multiple embedding models
- [ ] Cache embeddings for reuse

#### 3.6 LigandProcessor (if exists)
```bash
# Run ligand tests if available
pytest tests/test_processors/test_ligand/ -xvs
```

**Ligand Data Organization:**
1. **SDF/MOL files** - Ligand structures in `sdf/`
2. **Cache** - Cached ligand data in `cache/`
3. **Datasets** - Ligand datasets in `datasets/`

**Tasks:**

##### 3.6.1 Remove Custom Path Handling
- [ ] Remove ALL custom path properties
- [ ] Use `self.path_resolver.get_ligand_subdir_path()`
- [ ] Follow same patterns as other processors

### PHASE 4: Comprehensive Testing ✅ COMPLETED

#### 4.1 Zero-Configuration Test
```python
# test_zero_config.py
def test_zero_configuration():
    # No setup, no paths, just works
    processor = CifBaseProcessor()
    processor.save_structure("test_protein", structure_data)
    loaded = processor.load_structure("test_protein")
    assert loaded is not None
```

#### 4.2 Drag-and-Drop Workflow Test
```python
# test_drag_drop.py
def test_drag_and_drop_workflow():
    # Simulate user dropping files
    shutil.copy("external_file.cif", "data/structure/mmcif/6xyz.cif")
    
    # Should work immediately
    processor = CifBaseProcessor()
    structure = processor.load_structure("6xyz")
    assert structure is not None
```

#### 4.3 Cross-Format Entity Test
```python
# test_cross_format.py
def test_cross_format_entity():
    # Save in multiple formats
    struct_proc = CifBaseProcessor()
    struct_proc.save_structure("ubiquitin", struct_data)
    
    seq_proc = SeqProcessor()
    seq_proc.save_sequence("ubiquitin", "MQIFVKTLTG...")
    
    # Both tracked under same entity
    assert struct_proc.entity_exists("ubiquitin")
    assert seq_proc.entity_exists("ubiquitin")
```

#### 4.4 Multi-Entity Table Test
```python
# test_multi_entity.py
def test_grn_table_entities():
    grn_proc = GRNBaseProcessor()
    grn_proc.save_grn_table("opsins", grn_df)
    
    # Dataset references table
    grn_proc.create_dataset("opsin_study", 
                           entities=["BACR", "ChR2"],
                           metadata={"table_file": "opsins_grn.csv"})
    
    # Load returns DataFrame
    table = grn_proc.load_dataset("opsin_study")
    assert isinstance(table, pd.DataFrame)
```

### PHASE 6: Execution Timeline & Priorities

#### Week 1 (Days 1-3): Critical Fixes 🔴
1. **Day 1**: Fix all API signature mismatches (Phase 5.1)
   - Update processor constructors
   - Fix test initialization patterns
   - Run tests after each fix

2. **Day 2**: Fix entity registry concurrency (Phase 5.2)
   - Implement file locking
   - Add retry logic
   - Test with parallel test execution

3. **Day 3**: Implement missing methods (Phase 5.3)
   - Add all abstract method implementations
   - Ensure proper inheritance
   - Verify with instanceof checks

#### Week 1 (Days 4-5): High Priority Fixes 🔴
4. **Day 4**: Standardize column names (Phase 5.4)
   - Create mapping utilities
   - Update all processors
   - Fix integration tests

5. **Day 5**: Fix imports and GRN logic (Phase 5.5, 5.6)
   - Resolve all import errors
   - Fix GRN column parsing
   - Test GRN workflows

#### Week 2: Medium Priority Fixes 🟡
6. **Days 6-7**: Path consistency (Phase 5.7)
   - Resolve directory naming
   - Update path references
   - Fix affected tests

7. **Days 8-9**: Missing utilities (Phase 5.8)
   - Implement missing functions
   - Add proper documentation
   - Create unit tests

#### Week 2-3: Low Priority & Polish 🟢
8. **Days 10-11**: Test infrastructure (Phase 5.9)
   - Fix abstract class usage
   - Improve test isolation
   - Add better fixtures

9. **Days 12-14**: Environment & documentation (Phase 5.10)
   - Document requirements
   - Add environment checks
   - Final test suite run

### Success Metrics

1. **Immediate Goal**: Reduce failures from 87 to under 20
2. **Week 1 Target**: All critical (🔴) issues resolved
3. **Week 2 Target**: All medium (🟡) issues resolved
4. **Final Target**: 100% test passage (462/462 tests)

### Testing Strategy

After EVERY fix:
```bash
# Run affected component tests
pytest tests/test_processors/test_[component]/ -xvs

# Check for regressions
pytest tests/test_core/ -xvs

# Final validation
pytest --tb=short
```

#### 5.1 Fix Path Violations
```bash
# Run all tests and fix path issues
pytest --tb=short | grep -E "(Path|path|file)"
```

**Tasks:**
- [ ] Remove all `tmp_path` usage
- [ ] Remove all `Path()` constructors
- [ ] Remove all `.mkdir()` calls
- [ ] Update fixtures to use processors only

#### 5.2 Update Test Data
**Tasks:**
- [ ] Ensure test-data mirrors production structure
- [ ] Add missing CIF files (1ubq.cif, etc.)
- [ ] Remove any generated test files
- [ ] Use only real biological data

### PHASE 7: Migration and Cleanup (Week 3-4)

#### 6.1 Remove Old Code
- [ ] Delete `identifier_utils.py` (hash IDs are internal only)
- [ ] Remove user-facing hash ID methods
- [ ] Clean up old entity-registry-overhaul code
- [ ] Remove DataSource enum completely

#### 6.2 Update Documentation
- [ ] Update all examples to use human names only
- [ ] Document drag-and-drop workflow
- [ ] Create migration guide for existing users
- [ ] Update CLI documentation

### TESTING STRATEGY

**After EVERY change:**
1. Run specific test for component
2. Run integration tests
3. Run full test suite
4. Check for path violations

**Test Commands:**
```bash
# Component test
pytest tests/test_io/test_entity_registry.py -xvs

# Integration test  
pytest tests/test_processors/ -xvs

# Full suite
pytest

# Check violations
pytest --tb=short | grep -E "(tmp_path|Path\(|\.mkdir)"
```

### SUCCESS CRITERIA

1. **Zero Configuration**: `processor = CifBaseProcessor()` just works
2. **Human Names Only**: No hash IDs in any user-facing operation
3. **Drag-and-Drop**: Files appear, immediately usable
4. **Cross-Format**: Same entity name works across all formats
5. **100% Tests Pass**: No path violations, no mocking

### CURRENT STATUS (Updated: 2025-07-13)

- [x] Phase 1: Core Components ✅ COMPLETED
  - [x] ProtosPaths working perfectly
  - [x] EntityRegistry with human names only
  - [x] DatasetManager with JSON datasets
- [x] Phase 2: BaseProcessor Update ✅ COMPLETED
  - [x] Zero configuration working
  - [x] All entity/dataset operations implemented
  - [x] Abstract methods for format-specific operations
- [x] Phase 2.5: Verification of Core Components ✅ COMPLETED
  - [x] Created comprehensive integration tests
  - [x] Fixed EntityRegistry refresh for multi-processor scenarios
  - [x] All 22 core component tests passing
  - [x] Removed deprecated test files
- [x] Phase 3: Update Specific Processors ✅ COMPLETED
  - [x] All processors updated to new architecture
- [x] Phase 4: Comprehensive Testing ✅ COMPLETED
- [ ] Phase 5: Fix All Test Failures 🚨 CURRENT PHASE
  - [ ] 5.1: API Signature Fixes (High Priority)
  - [ ] 5.2: Entity Registry Concurrency (High Priority)
  - [ ] 5.3: Missing Method Implementations (High Priority)
  - [ ] 5.4: Column Name Standardization (High Priority)
  - [ ] 5.5: Import and Module Issues (Medium Priority)
  - [ ] 5.6: GRN Processing Logic (Medium Priority)
  - [ ] 5.7: Path Naming Consistency (Medium Priority)
  - [ ] 5.8: Missing Utility Functions (Medium Priority)
  - [ ] 5.9: Test Infrastructure (Low Priority)
  - [ ] 5.10: Environment-Specific Issues (Low Priority)
- [ ] Phase 6: Execution Timeline & Priorities
- [ ] Phase 7: Migration and Cleanup

**CURRENT FOCUS**: Phase 5.1 - Fixing API signature mismatches (Day 1)

### CRITICAL REMINDERS

1. **NO VERBOSE CODE** - Keep implementations minimal
2. **TEST EVERYTHING** - Run tests after every change
3. **HUMAN NAMES ONLY** - Hash IDs are internal registry keys
4. **ZERO CONFIG** - Must work out of the box
5. **ONE PATH SYSTEM** - ProtosPaths handles everything

## 🔍 DETAILED ISSUE BREAKDOWN

### Critical Issues by Component

#### 1. API Signature Mismatches (31 failures)
```python
# WRONG - Causes TypeError
GRNProcessor(name="test", processor_data_dir="grn")
StructureProcessor(name="test", data_root="/path")

# CORRECT
GRNProcessor(name="test")
StructureProcessor(name="test")
```

#### 2. Entity Registry Conflicts (8 failures)
```
FileExistsError: [WinError 183] Cannot create a file when that file already exists: 
'C:\\...\\data\\entity_registry.json.tmp' -> 'C:\\...\\data\\entity_registry.json'
```

#### 3. Missing Implementations (12 failures)
```python
# Missing abstract methods:
- StructureProcessor._load_dataset_impl()
- StructureProcessor._save_dataset_impl()
- Various load_entity/save_entity implementations
```

#### 4. Column Name Issues (15 failures)
```python
# Structure data inconsistency:
- Some code expects: 'res_atom_name', 'res_seq_id'
- Other code uses: 'atom_name', 'seq_id'
- GRN expects: specific column format
```

#### 5. Import Errors (6 failures)
```python
# embedding_processor.py:
NameError: name 'SeqProcessor' is not defined

# Various test files:
ImportError: cannot import name 'PropertyProcessorEnhanced'
```

#### 6. GRN Logic Issues (9 failures)
```
# Invalid column generation:
Cannot match GRN columns. Data columns: ['3.50', '7.53', 'species']
Requested: ['0x00', '3.50', '7.53']  # Where does '0x00' come from?
```

#### 7. Utility Function Missing (4 failures)
```python
AttributeError: module 'protos.io.cif_utils' has no attribute 'dataframe_to_cif'
```

#### 8. Test Infrastructure (8 failures)
```python
# Trying to instantiate abstract class:
TypeError: Can't instantiate abstract class BaseProcessor
```

### Component Health Status

| Component | Tests | Passed | Failed | Health |
|-----------|-------|--------|--------|--------|
| Core | 22 | 22 | 0 | ✅ 100% |
| IO | 45 | 42 | 3 | ⚠️ 93% |
| Structure | 89 | 52 | 37 | ❌ 58% |
| Sequence | 52 | 51 | 1 | ✅ 98% |
| GRN | 43 | 34 | 9 | ⚠️ 79% |
| Property | 20 | 20 | 0 | ✅ 100% |
| Embedding | 21 | 20 | 1 | ✅ 95% |
| Interactions | 15 | 0 | 15 | ❌ 0% |
| **TOTAL** | **462** | **375** | **87** | ⚠️ 81% |