# PROTOS TODO

## 🎯 CORE PHILOSOPHY & CRITICAL PRINCIPLES

### ⚠️ CRITICAL: The Protos Promise
**Users work with names, never paths. Protos handles ALL file system complexity.**

**THIS IS NOT OPTIONAL - IT'S THE FUNDAMENTAL CONTRACT OF PROTOS**

### What This Means:
1. **Complete Abstraction** - Users NEVER see or manipulate file paths
2. **Name-Based Access** - Everything accessed by biological/dataset names  
3. **Universal Entity IDs** - One protein = one hash ID across all formats
4. **Transparent Operations** - Load/save/convert without format concerns
5. **Registry as Truth** - All lookups and discovery through registry

### ⚠️ WHY THIS MATTERS:
- **User Experience**: Scientists should think in biological terms, not file systems
- **Cross-Platform**: Path handling varies between Windows/Linux/Mac
- **Data Management**: Centralized control prevents data corruption/loss
- **Evolution**: System can reorganize files without breaking user code
- **Integration**: External tools can be swapped without user impact

### The Golden Rules:
```python
# ✅ CORRECT - Users work with names
processor.load_structure("1ABC")           # PDB ID
processor.load_sequence("P12345")          # UniProt ID
processor.load_grn_annotation("BR1_HUMAN") # Sequence name
processor.save_dataset("my_results", data) # Dataset name

# ❌ WRONG - Users never see paths or files
open("/path/to/file.csv")                 # NO direct file access
Path("data/structures/1ABC.cif")          # NO path construction
pd.read_csv("file.csv")                   # NO manual format handling
os.path.exists("some/path")               # NO file system checks
```

### Testing Principles:
0. **USE THE EXISTING MINICONDA ENV: 'protos'** - Activate it during development 
```bash
conda activate protos
```
1. **USE REAL TEST DATA** - Never mock biological data
2. **USE tests/test-data/** - Mirror production data structure  
3. **NO TEMPORARY FILES** - Use fixture data, not generated files
4. **TEST WITH REAL FILES** - CIF, FASTA, CSV from test-data
5. **FOLLOW PROTOS PHILOSOPHY** - Names only, never paths

### Current Violations:
- **170+ path violations in tests** - Breaking core abstraction
- **No universal entity system** - Can't track entities across formats
- **Incomplete registry** - Users can't discover available data
- **Direct file operations** - Bypassing processor abstraction

---

## 📁 PATH SYSTEM - ProtosPaths

**See [docs/PROTOSPATH.md](docs/PROTOSPATH.md) for complete documentation**

### Architecture Summary:
- **ProtosPaths** - Central path resolver (singleton pattern)
- **Global Configuration** - Set once via `ProtosPaths.set_data_root()`
- **Processor Integration** - BaseProcessor uses ProtosPaths automatically
- **Complete Abstraction** - Users never construct or see paths

### Key Components:
```python
# Global configuration (only in conftest.py for tests)
ProtosPaths.set_data_root("/path/to/data")

# Processors use it automatically
processor = CifBaseProcessor()  # No path parameters!
processor.load_structure("1ABC")  # Just names!
```

---

## 🔑 REGISTRY SYSTEM - Entity & Dataset Management

### Universal Entity ID System:
**One biological entity = One hash ID across ALL formats**

```python
# Same protein gets SAME entity ID everywhere
"P12345" → entity_id "a3f2d8c91b"
  ├── sequence format: /seq/a3f2d8c91b.fasta
  ├── structure format: /struct/a3f2d8c91b.cif  
  ├── GRN format: in table with entity_id column
  └── embedding format: /emb/a3f2d8c91b.pkl
```

### Registry Architecture:
1. **EntityRegistry** - Tracks all entities with hash IDs
   - Multi-format support (same entity, multiple formats)
   - Original ID → Hash ID lookups
   - Metadata per format
   
2. **Registry as Universal Translator**:
   ```python
   # User provides biological name
   processor.load_structure("1ABC")
   
   # Registry translates transparently
   entity_id = registry.resolve_identifier("1ABC")  # → "a3f2d8c91b"
   
   # System uses hash internally
   data = processor._load_entity_by_hash(entity_id)
   ```

3. **GRN Tables Special Case**:
   ```
   entity_id    | sequence_id | 1.50 | 2.50 | 3.50
   -------------|-------------|------|------|------
   a3f2d8c91b   | BR1_HUMAN   | L45  | V87  | I123
   b4e5f6a72c   | BR2_MOUSE   | L46  | V88  | I124
   ```

---

## 📋 TODOS BY PRIORITY

### 🔴 HIGH PRIORITY - Registry Implementation

#### Phase 1: Complete Multi-Format EntityRegistry ✅ COMPLETED
- [x] Basic EntityRegistry with hash IDs
- [x] Multi-format support per entity
- [x] Original ID lookups
- [x] Add resolve_identifier method:
  ```python
  def resolve_identifier(self, identifier: str) -> str:
      """Resolve any name to entity hash ID."""
      # Check if already hash
      # Lookup by original ID  
      # Generate new if needed
  ```

#### Phase 2: Update ALL Processors for Entity Support ✅ COMPLETED
- [✅] **CifBaseProcessor**:
  - [x] Basic entity registration
  - [x] load_structure uses resolve_identifier
  - [x] list_structures returns original IDs (not hashes!)
  - [x] Full multi-format awareness
  - [x] save_structure_as_entity method
  - [x] get_entity_id_for_pdb method
  
- [✅] **GRNBaseProcessor** (Special Case - Tables contain MULTIPLE entities):
  ```python
  def save_grn_table_with_entities(self, grn_df):
      # Add entity_id column as FIRST column
      entity_ids = [generate_entity_id(seq_id) for seq_id in grn_df.index]
      grn_df.insert(0, 'entity_id', entity_ids)
      
      # Register each row as an entity
      for seq_id, entity_id in zip(grn_df.index, entity_ids):
          register_entity(entity_id, "grn", seq_id, metadata={"in_table": table_name})
  ```
  
  **Important GRN Design Decisions**:
  - GRN tables contain MANY entities (one per row/sequence)
  - Each row gets an entity_id column based on sequence ID
  - Consider: Should single GRN entities be stored as:
    a) JSON files for individual annotated sequences?
    b) FASTA files (since GRN entity = annotated sequence)?
    c) Single-row CSV tables (current approach)?
  
- [✅] **SeqProcessor**:
  - [x] Parse FASTA headers for entity names
  - [x] Register each sequence as entity
  - [x] Support multi-sequence files
  - [x] load_sequence_entity method
  - [x] save_sequence_entity method
  - [x] list_sequence_entities method
  
- [✅] **EmbeddingProcessor**:
  - [x] Use same entity IDs as source sequences
  - [x] Register embeddings by model type
  - [x] _register_embedding_entities method
  - [x] load_embedding_entity method  
  - [x] list_embedding_entities method

#### Phase 3: Complete Processor Testing ✅ COMPLETED
- [x] Test resolve_identifier for all processors
- [x] Test list operations return names not hashes
- [x] Test multi-format entity scenarios
- [x] Test GRN tables with entity_id column
- [x] Test that all processors can load by both name and hash ID

### 🟡 MEDIUM PRIORITY - Cross-Format Operations

#### Phase 4: Cross-Format Workflows ✅ COMPLETED
- [x] Sequence → Structure (AlphaFold)
- [x] Structure → Sequence extraction  
- [x] Sequence → GRN assignment
- [x] Any format → Embeddings
- [x] Track conversion lineage in metadata

#### Phase 5: Migration Tools ✅ COMPLETED
- [x] Script to migrate existing data to entity system
- [x] Update all test data with proper entities
- [x] Documentation for users to migrate
- [x] Registration commands for manual file additions
- [x] Auto-registration for downloaded files
- [x] Scan function to detect unregistered files
- [x] Dataset creation from entity lists

### 🔴 HIGH PRIORITY - DOCUMENTATION FOR RELEASE

#### Comprehensive Documentation Task (NEW - CRITICAL FOR RELEASE)
**Goal**: Create complete documentation for pip-installable Protos release

##### Documentation Requirements:
1. **Code Documentation**:
   - [ ] Add concise docstrings to all files (module-level)
   - [ ] Document important functions/classes with clear descriptions
   - [ ] Review ALL comments - remove "you" language, use neutral matter-of-fact style
   - [ ] Ensure docstrings follow Google/NumPy style consistently

2. **Documentation Structure** (docs/):
   - [ ] **docs/protospath.md** - Path management system documentation
   - [ ] **docs/core_io.md** - Core I/O operations and philosophy
   - [ ] **docs/baseprocessor.md** - BaseProcessor architecture and usage
   - [ ] **docs/fileformats.md** - Supported file formats and handling
   - [ ] **docs/entities.md** - Entity system and cross-format tracking
   - [ ] **docs/registries.md** - Registry system and dataset management
   - [ ] **docs/processors/** - Individual processor documentation:
     - [ ] **cifbase_processor.md** - Structure processing
     - [ ] **grn_processor.md** - GRN annotation system
     - [ ] **seq_processor.md** - Sequence handling
     - [ ] **embedding_processor.md** - Embedding generation
     - [ ] **property_processor.md** - Property management
   - [ ] **docs/processing/** - Format-specific operations:
     - [ ] **alignments.md** - Sequence and structure alignment
     - [ ] **annotations.md** - GRN and other annotations
     - [ ] **conversions.md** - Cross-format conversions
   - [ ] **docs/cli/** - Command-line usage:
     - [ ] **installation.md** - Installation guide
     - [ ] **quickstart.md** - Getting started tutorial
     - [ ] **commands.md** - Complete CLI reference
     - [ ] **examples.md** - Real-world usage examples

3. **Review existing example**:
   - [ ] Check protos_review.py for documentation style guidance
   - [ ] Run tests to understand functionality for accurate documentation

4. **Documentation Standards**:
   - Use clear, neutral language (no "you" addressing)
   - Include code examples that follow Protos philosophy (names, not paths)
   - Link documentation to actual code where applicable
   - Focus on helping users understand usage and core philosophy

5. **Testing Documentation**:
   - [ ] Verify all code examples work
   - [ ] Ensure documentation matches current implementation
   - [ ] Test installation instructions on clean environment

### 🔴 CRITICAL - CORE PHILOSOPHY VIOLATIONS MUST BE FIXED FIRST

⚠️ **BREAKING THE CORE PROMISE: "Users work with names, never paths. Protos handles ALL file system complexity."**

#### Path System Violations Fixed ✅
- [x] Fix conftest.py to use single data root with ProtosPaths
- [x] Remove all tmp_path usage from test files (170+ violations) 
- [x] Eliminate manual directory creation (.mkdir() calls)

#### Remaining Path Issues (25 test failures) - HIGH PRIORITY
- [ ] **CRITICAL**: Fix attribute errors - processors missing path properties
  - `BaseProcessor` missing `path_structure_dir` attribute (downloads tests)
  - `GRNBaseProcessor` missing `path_grn_ref` attribute
  - Path operations returning strings instead of Path objects
- [ ] **CRITICAL**: Fix hardcoded path references in tests
  - Tests looking for files in `C:\Users\hidbe/protos_data\` instead of test-data
  - Missing imports of SeqProcessor in GRN assignment tests
  - Path concatenation errors (`TypeError: unsupported operand type(s) for /: 'str' and 'str'`)
- [ ] **CRITICAL**: Add missing test data files  
  - `1ubq.cif`, `1tqn.cif` referenced by tests but missing from test-data
  - Structure tests failing because files don't exist in test directory

#### Test System Issues (25 failures total)
- [ ] **Embedding Processor** (5 failures): Ankh model compatibility with PyTorch
- [ ] **Download System** (5 failures): Missing processor path attributes 
- [ ] **Structure Loading** (8 failures): Missing test data files (1ubq.cif, 1tqn.cif)
- [ ] **Path Operations** (6 failures): String/Path type mixing errors
- [ ] **Missing Imports** (1 failure): SeqProcessor not imported in test

### 🚨 IMMEDIATE ACTION PLAN (Fix in this order)

#### Step 1: Fix Missing Path Properties (HIGH PRIORITY)
1. Add missing `path_structure_dir` property to `BaseProcessor` class
2. Add missing `path_grn_ref` property to `GRNBaseProcessor` class  
3. Ensure all path properties return `Path` objects, not strings
4. Update processor constructors to use ProtosPaths correctly

#### Step 2: Fix Missing Imports (EASY WIN)
1. Add `SeqProcessor` import to `test_grn_assignment.py`
2. Verify all test files have required imports

#### Step 3: Add Missing Test Data Files (HIGH PRIORITY)
1. Add `1ubq.cif` to `tests/test-data/structure/mmcif/`
2. Add `1tqn.cif` to `tests/test-data/structure/mmcif/`
3. Verify structure files are valid and loadable

#### Step 4: Fix ProtosPaths Configuration Issues
1. Ensure conftest.py ProtosPaths configuration propagates to all processors
2. Debug why tests are looking in production path instead of test-data
3. Fix environment variable settings for test isolation

#### Step 5: Path Type Consistency  
1. Audit all processor path properties for consistent Path object returns
2. Fix string/Path concatenation errors
3. Update tests that assume string paths

### 🟡 MEDIUM PRIORITY - Path System Cleanup

#### Update Test Patterns (After core fixes)
- [ ] Replace direct Path() usage with processor methods
- [ ] Remove manual path construction (processor.data_path / "subdirectory")  
- [ ] Update tests to use processor.is_dataset_available() instead of file checks
- [ ] Fix environment variable manipulation for paths

### 🟢 LOW PRIORITY - Real Data Management

#### Implement Proper Real Data Loading
- [ ] Add loader functions for downloading test data during tests
- [ ] Implement init/cleanup CLI integration for reference data
- [ ] Create proper real data fixtures using ProtosPaths
- [ ] Remove direct file download code from test fixtures

### 🟢 LOW PRIORITY - Future Enhancements

#### Phase 6: Advanced Features
- [ ] Entity versioning support
- [ ] Entity relationship tracking
- [ ] Bulk operations optimization
- [ ] Export/import entity registry

#### Phase 7: Documentation & Examples
- [ ] Complete entity system user guide
- [ ] Cross-format workflow examples
- [ ] Migration guide from old system
- [ ] Best practices documentation

---

## 🐛 OTHER ISSUES TO FIX

### Test Issues:
- [ ] Fix 170+ path violations in tests
- [ ] Remove all Path() usage from tests
- [x] Fix GRN test data setup (column formats) - COMPLETED: Updated to use proper residue+position format
- [ ] Add proper test markers (@pytest.mark.slow, etc.)
- [ ] Fix embedding test timeouts (use small models)

### Test Results Summary (UPDATED - June 30, 2025 - AFTER PATH VIOLATION FIXES):

#### CURRENT TEST STATUS - PATH VIOLATIONS INTRODUCED NEW FAILURES:
- **Core tests**: 37/37 (100%) ✅  
- **IO tests**: 43/43 (100%) ✅
- **CLI tests**: 9/9 (100%) ✅
- **GRN processor tests**: 69/71 (97%, 2 failed due to missing path attributes) ⚠️
- **Sequence processor tests**: 40/41 (97%, 1 skipped) ✅
- **Embedding processor tests**: 45/50 (90%, 5 failed - Ankh model issues) ⚠️
- **Loader tests**: 46/52 (88%, 6 failed due to missing processor attributes) ⚠️
- **Structure processor tests**: 15/42 (36%, 27 failed - path issues + missing files) ❌
- **Workflow tests**: 3/7 (43%, 4 failed due to missing test files) ❌

#### TOTAL: 307/352 tests passing (87.2% success rate) - REGRESSION from path fixes

#### Critical Issues Introduced by Path Violation Fixes:

**1. Missing Processor Path Attributes (NEW - HIGH PRIORITY)**:
- `BaseProcessor` missing `path_structure_dir` property 
- `GRNBaseProcessor` missing `path_grn_ref` property
- Tests expect these attributes but they don't exist
- **Root cause**: Incomplete processor path property implementation

**2. Wrong Data Directory Usage (NEW - CRITICAL)**:
- Tests looking for files in `C:\Users\hidbe/protos_data\` (production path)
- Should be using `tests/test-data/` directory 
- ProtosPaths configuration not properly applied to all processors
- **Root cause**: Environment/path configuration not properly propagated

**3. Missing Test Data Files (NEW - HIGH PRIORITY)**:
- `1ubq.cif`, `1tqn.cif` referenced by multiple tests but don't exist
- Structure tests failing because expected files missing from test-data
- **Fix needed**: Add missing CIF files to test-data directory

**4. Path Type Errors (NEW - HIGH PRIORITY)**:
- `TypeError: unsupported operand type(s) for /: 'str' and 'str'`
- Path concatenation failing because properties return strings not Path objects
- **Root cause**: Inconsistent return types from path properties

**5. Missing Imports (NEW - EASY FIX)**:
- `NameError: name 'SeqProcessor' is not defined` in GRN tests
- Missing import statements in test files

#### Key Achievements:
1. **Complete Entity System**: All processors now support universal entity IDs
2. **Path System Stability**: ProtosPaths working reliably across all tests
3. **GRN System Robustness**: 71/71 tests passing, supports extended numbering
4. **Cross-Format Operations**: Sequence ↔ Structure workflows functional
5. **Registration System**: Entity registration and discovery working
6. **Data Management**: Clean separation of test vs reference data

#### Test Infrastructure Improvements:
1. Removed all example files - using only real test data
2. Fixed all path resolution issues
3. Updated validation to support extended GRN formats
4. Improved error handling in download workflows
5. Better test isolation and cleanup

#### Issues Identified from Notebook Review:
1. [x] GRN column warnings for metadata columns - FIXED
2. [x] Pandas SettingWithCopyWarning - FIXED
3. [x] Sequence loading returning None - FIXED
4. [x] Dataset loading 'path' error - FIXED
5. [x] Registration summary 'total_files' KeyError - FIXED
6. [x] Missing extract_sequence method - FIXED
7. [x] GRN validation allowing helix 0 incorrectly - FIXED

### CLI Updates:
- [ ] Update CLI embed.py to use new EmbeddingProcessor
- [ ] Ensure all CLI commands use entity system
- [ ] Add entity discovery commands

### Documentation:
- [ ] Update README.md with correct usage
- [ ] Document entity system
- [ ] Add workflow examples

---

## 📊 PROGRESS TRACKING

### Completed ✅:
- ProtosPaths architecture and documentation
- Basic EntityRegistry implementation  
- Multi-format entity support
- Hash-based entity ID generation
- BaseProcessor entity methods
- Phase 1: Complete Multi-Format EntityRegistry
- Phase 2: Update ALL Processors for Entity Support
- Phase 3: Complete Processor Testing
- Phase 4: Cross-Format Workflows
- Phase 5: Migration Tools and Data Management
- Added Entity System section to README
- Created comprehensive data management system
- Registration/download CLI commands
- Auto-registration for downloads

### In Progress 🔄:
- Investigate intermittent test failures (appears to be test isolation issue)
- Fix path violations in tests
- Update existing CLI commands to use entity system
- Advanced features (versioning, relationships)

### Blocked ❌:
- Cross-format operations (need processors done first)
- Migration tools (need final system design)

---

## 💡 IMPLEMENTATION NOTES

### Entity ID Generation:
```python
# ALWAYS hash the biological identifier, not filename
entity_id = generate_entity_id("P12345")     # UniProt
entity_id = generate_entity_id("1ABC")       # PDB
entity_id = generate_entity_id("BR1_HUMAN")  # Sequence name

# For FASTA with random names
>some_random_name
MKLSPADKTN...
entity_id = generate_entity_id("some_random_name")
```

### Processor Pattern:
```python
class AnyProcessor(BaseProcessor):
    def load_entity(self, identifier: str):
        # User provides name
        entity_id = self.registry.resolve_identifier(identifier)
        
        # Internal uses hash
        return self._load_by_hash(entity_id)
        
    def list_entities(self):
        # Return names, not hashes!
        entities = self.registry.list_entities(format_type=self.format)
        return [self.registry.get_original_id(e) for e in entities]
```

### Timeline:
- Week 1: Complete resolve_identifier and processor updates
- Week 2: Full testing suite  
- Week 3: Cross-format operations
- Week 4: Migration and documentation