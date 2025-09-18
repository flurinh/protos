# Entity System Refactoring: Detailed TODO Plan

## Quick Start Guide

### Environment Setup
```bash
# Activate the Protos conda environment
source /home/hidberf/miniconda/bin/activate protos  # Linux/WSL
# or
conda activate protos  # Windows

# Verify environment
python -c "import protos; print('Protos environment active')"
```

### Key Development Principles

1. **ProtosPaths Manages Everything**
   ```python
   # NEVER hardcode paths
   # ❌ WRONG
   path = Path("/home/user/data/structure/mmcif/1ubq.cif")
   
   # ✅ CORRECT - Let ProtosPaths handle it
   from protos.io.paths import ProtosPaths
   paths = ProtosPaths()  # Uses PROTOS_DATA_ROOT env var or default
   struct_path = Path(paths.get_processor_path("structure")) / "mmcif" / "1ubq.cif"
   ```

2. **Testing Pattern**
   ```python
   # Set up test environment (from conftest.py pattern)
   os.environ["PROTOS_DATA_ROOT"] = str(test_data_root)
   paths = ProtosPaths()  # Will use the env var
   
   # Ensure directories exist by calling getters
   paths.get_processor_path("structure")
   paths.get_global_registry_path()
   ```

3. **Current State**
   - Entity IDs are 10-character hashes (e.g., `3f8a9c2d1e`)
   - Migration tests in `migration_tests/` guide the refactoring
   - All public APIs already return human names (no changes needed there)

### Running Tests
```bash
# Run all migration tests
python -m pytest migration_tests/ -v

# Run specific test
python -m pytest migration_tests/test_entity_registry_core.py::TestEntityRegistryCore::test_entity_registration_creates_uuid -xvs

# Run with coverage
python -m pytest migration_tests/ --cov=protos.io.entity_registry
```

### Important Files
- `src/protos/io/entity_registry.py` - Core registry implementation
- `src/protos/io/paths/path_config.py` - ProtosPaths implementation
- `src/protos/core/base_processor.py` - Base processor with entity operations
- `migration_tests/` - Tests that guide the refactoring

---

## Overview
Implement a dual identifier system where:
- **Entity IDs**: UUIDs for internal tracking (never shown to users)
- **Entity Names**: Human-readable names for all user interactions
- **Relationships**: Track derivations, versions, and connections between entities

## Phase 0: Safe Registration and File Management ✅ COMPLETED (Jan 19, 2025)

### Implementation Summary
- ✅ Created InputManager for safe file registration workflow
- ✅ Implemented conflict resolution strategies (skip, version, replace)
- ✅ Added content hash-based duplicate detection
- ✅ Created RegistryHealthChecker for dataset integrity
- ✅ Centralized file format definitions in FormatRegistry
- ✅ Cleaned up IO module to remove duplicates and use centralized formats
- ✅ All 11 InputManager tests passing

### Data Acquisition Methods in Protos

1. **Reference Data** (Built-in):
   - Shipped with Protos installation
   - Initialized during `protos init-data`
   - Read-only, never modified by users
   - Examples: Standard GRN tables, reference structures

2. **Online Databank Loaders** (Automated):
   - Download from PDB, UniProt, ChEMBL, etc.
   - Automatic registration after download
   - Validates data integrity
   - Examples: `download_structures()`, `fetch_sequences()`

3. **User-Provided Data** (New safe workflow):
   - Users place files in designated input folder
   - Protos validates and registers them
   - Moves to appropriate data directories
   - **Replaces unsafe "drag-and-drop" approach**

### Critical: Safe User Data Registration Workflow

#### 0.1 Input Folder System ✅ COMPLETED
**File**: `src/protos/io/input_manager.py` 

- [x] **Implement input folder monitoring**
  ```python
  class InputManager:
      """Manages user data input and registration."""
      
      def __init__(self, paths: ProtosPaths):
          self.paths = paths
          self.input_dir = paths.get_input_directory()  # data/input/
          self.processed_dir = paths.get_processed_directory()  # data/input/processed/
          self.rejected_dir = paths.get_rejected_directory()  # data/input/rejected/
          
      def scan_input_folder(self) -> List[InputFile]:
          """Scan input folder for new files to process."""
          input_files = []
          
          # Group files by type
          for file_path in self.input_dir.iterdir():
              if file_path.is_file():
                  file_type = self._detect_file_type(file_path)
                  if file_type:
                      input_files.append(InputFile(
                          path=file_path,
                          file_type=file_type,
                          entity_name=self._extract_entity_name(file_path)
                      ))
          
          return input_files
      
      def process_input_files(self, 
                            conflict_strategy: ConflictResolutionStrategy = SKIP,
                            dry_run: bool = False) -> ProcessingReport:
          """Process all files in input folder."""
          report = ProcessingReport()
          
          for input_file in self.scan_input_folder():
              try:
                  # Get appropriate processor
                  processor = self._get_processor_for_type(input_file.file_type)
                  
                  # Validate file
                  validation = processor.validate_before_registration(input_file.path)
                  if not validation.is_valid:
                      self._move_to_rejected(input_file, validation.errors)
                      report.add_rejected(input_file, validation.errors)
                      continue
                  
                  # Register with safety checks
                  result = processor.register_file_safely(
                      input_file.path,
                      input_file.entity_name,
                      conflict_strategy=conflict_strategy,
                      dry_run=dry_run
                  )
                  
                  if not dry_run and result.success:
                      # Move to processed folder
                      self._move_to_processed(input_file)
                  
                  report.add_result(input_file, result)
                  
              except Exception as e:
                  report.add_error(input_file, str(e))
          
          return report
      
      def _move_to_processed(self, input_file: InputFile):
          """Move successfully processed file to processed folder."""
          timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
          processed_path = self.processed_dir / f"{timestamp}_{input_file.path.name}"
          shutil.move(str(input_file.path), str(processed_path))
      
      def _move_to_rejected(self, input_file: InputFile, errors: List[str]):
          """Move rejected file with error log."""
          rejected_path = self.rejected_dir / input_file.path.name
          error_log = self.rejected_dir / f"{input_file.path.stem}_errors.txt"
          
          shutil.move(str(input_file.path), str(rejected_path))
          error_log.write_text("\n".join(errors))
  ```

- [x] **Add CLI commands for input processing** ✅ COMPLETED
  ```python
  # In cli.py
  @click.command()
  @click.option('--dry-run', is_flag=True, help='Preview without making changes')
  @click.option('--strategy', type=click.Choice(['skip', 'version', 'replace']), 
                default='skip', help='Conflict resolution strategy')
  def register_input(dry_run, strategy):
      """Process and register files from input folder."""
      manager = InputManager(ProtosPaths())
      
      # Show what's in input folder
      input_files = manager.scan_input_folder()
      if not input_files:
          click.echo("No files found in input folder")
          return
      
      click.echo(f"Found {len(input_files)} files to process:")
      for f in input_files:
          click.echo(f"  - {f.path.name} ({f.file_type})")
      
      # Process with user confirmation
      if not dry_run:
          if not click.confirm("Proceed with registration?"):
              return
      
      # Process files
      report = manager.process_input_files(
          conflict_strategy=ConflictResolutionStrategy[strategy.upper()],
          dry_run=dry_run
      )
      
      # Show results
      report.display()
  ```

- [x] **Update ProtosPaths for input directories** ✅ COMPLETED (handled in InputManager)
  ```python
  # In path_config.py
  def get_input_directory(self) -> Path:
      """Get user input directory."""
      input_dir = self.data_root / "input"
      input_dir.mkdir(exist_ok=True)
      return input_dir
  
  def get_processed_directory(self) -> Path:
      """Get processed input directory."""
      processed = self.data_root / "input" / "processed"
      processed.mkdir(parents=True, exist_ok=True)
      return processed
  
  def get_rejected_directory(self) -> Path:
      """Get rejected input directory."""
      rejected = self.data_root / "input" / "rejected"
      rejected.mkdir(parents=True, exist_ok=True)
      return rejected
  ```

### Updated Data Directory Structure

```
data/
├── input/                    # User places new files here
│   ├── README.txt           # Instructions for users
│   ├── processed/           # Successfully registered files (archive)
│   └── rejected/            # Failed validation (with error logs)
├── reference/               # Read-only reference data from Protos
│   ├── grn/                # Standard GRN tables
│   └── structures/         # Example structures
├── structure/              # Registered structure data
│   ├── mmcif/              # User's registered structures
│   └── downloaded/         # From online databanks
├── sequence/               # Registered sequence data
│   ├── fasta/              # User's registered sequences
│   └── downloaded/         # From UniProt, etc.
└── entity_registry.json    # Central registry
```

### User Workflow for Adding Data

1. **Place files in input folder**:
   ```bash
   # User copies their files
   cp my_structure.cif ~/protos_data/input/
   cp protein_properties.csv ~/protos_data/input/
   ```

2. **Preview what will be registered**:
   ```bash
   protos register-input --dry-run
   # Shows:
   # Found 2 files to process:
   #   - my_structure.cif (structure)
   #   - protein_properties.csv (property)
   # 
   # Validation results:
   #   ✓ my_structure.cif - Valid PDB file
   #   ✓ protein_properties.csv - Valid property table
   ```

3. **Register with chosen strategy**:
   ```bash
   # Skip duplicates (default, safest)
   protos register-input
   
   # Or version conflicts
   protos register-input --strategy version
   
   # Or replace with backup
   protos register-input --strategy replace
   ```

4. **Review results**:
   ```bash
   # Success: files moved to processed/
   # Failed: files moved to rejected/ with error logs
   ```

#### 0.2 File Registration Safety ✅ COMPLETED
**File**: `src/protos/io/entity_registry.py`

- [x] **Implement safe registration workflow**
  ```python
  def register_file_safely(self, file_path: Path, entity_name: str, 
                          format_type: str, metadata: Dict = None) -> str:
      """
      Safely register a file with duplicate detection and validation.
      
      Steps:
      1. Validate file exists and is readable
      2. Compute content hash for duplicate detection
      3. Check if content already exists
      4. Check if name conflicts with existing entity
      5. Handle conflicts appropriately
      6. Register with transaction-like safety
      """
      # Validate file
      if not file_path.exists():
          raise FileNotFoundError(f"File {file_path} does not exist")
      
      # Compute content hash
      content_hash = self._compute_file_hash(file_path)
      
      # Check for duplicate content
      existing = self.find_by_content_hash(content_hash)
      if existing:
          # Same content already registered
          return self._handle_duplicate_content(existing, entity_name, format_type)
      
      # Check for name conflicts
      if self.entity_exists(entity_name):
          return self._handle_name_conflict(entity_name, file_path, format_type)
      
      # Safe registration with rollback
      return self._register_with_transaction(entity_name, file_path, format_type, metadata)
  ```

- [x] **Add content-based duplicate detection** ✅ COMPLETED
  ```python
  def _compute_file_hash(self, file_path: Path, chunk_size: int = 8192) -> str:
      """Compute SHA256 hash of file content."""
      hasher = hashlib.sha256()
      with open(file_path, 'rb') as f:
          while chunk := f.read(chunk_size):
              hasher.update(chunk)
      return hasher.hexdigest()
  
  def find_by_content_hash(self, content_hash: str) -> Optional[EntityInfo]:
      """Find entity with matching content hash."""
      for entity in self._registry.values():
          for format_data in entity['formats'].values():
              if format_data.get('content_hash') == content_hash:
                  return self._entity_to_info(entity)
      return None
  ```

#### 0.2 Conflict Resolution Strategies ✅ COMPLETED
**File**: `src/protos/io/conflict_resolver.py`

- [x] **Implement conflict resolution options**
  ```python
  class ConflictResolutionStrategy(Enum):
      SKIP = "skip"              # Skip registration, return existing
      VERSION = "version"        # Create versioned name (e.g., protein_v2)
      REPLACE = "replace"        # Replace existing (with backup)
      MERGE = "merge"            # Merge with existing (for tables)
      ASK = "ask"                # Prompt user for decision
  
  class ConflictResolver:
      def resolve_duplicate_content(self, existing: EntityInfo, 
                                   new_name: str) -> ResolutionResult:
          """Handle case where file content already exists."""
          if existing.entity_name == new_name:
              # Exact duplicate - safe to skip
              return ResolutionResult(action=SKIP, entity=existing)
          else:
              # Same content, different name - add as alias
              return ResolutionResult(action=ADD_ALIAS, entity=existing)
      
      def resolve_name_conflict(self, existing_name: str, 
                              new_file: Path) -> ResolutionResult:
          """Handle case where name already exists."""
          # Compare with existing file
          if self._files_identical(existing_file, new_file):
              return ResolutionResult(action=SKIP)
          else:
              # Different content, same name - version it
              versioned_name = self._generate_version_name(existing_name)
              return ResolutionResult(action=VERSION, new_name=versioned_name)
  ```

#### 0.3 Bulk Registration with Validation ✅ COMPLETED
**File**: `src/protos/core/base_processor.py`

- [x] **Add bulk registration methods**
  ```python
  def register_directory(self, directory: Path, 
                        strategy: ConflictResolutionStrategy = SKIP,
                        dry_run: bool = False) -> RegistrationReport:
      """
      Register all valid files in a directory.
      
      Returns report with:
      - Successfully registered files
      - Skipped files (duplicates)
      - Failed files (errors)
      - Conflicts resolved
      """
      report = RegistrationReport()
      
      # Scan directory for valid files
      valid_extensions = self._get_valid_extensions()
      for file_path in directory.glob(f"*{ext}" for ext in valid_extensions):
          try:
              # Extract entity name from filename
              entity_name = file_path.stem
              
              # Attempt registration
              result = self.register_file_safely(
                  file_path, entity_name, self.format_type, 
                  strategy=strategy, dry_run=dry_run
              )
              report.add_result(file_path, result)
              
          except Exception as e:
              report.add_error(file_path, str(e))
      
      return report
  ```

- [x] **Add validation before registration** ✅ COMPLETED
  ```python
  def validate_before_registration(self, file_path: Path) -> ValidationResult:
      """Validate file before attempting registration."""
      validations = []
      
      # Check file readability
      validations.append(self._check_readable(file_path))
      
      # Check file format
      validations.append(self._check_format(file_path))
      
      # Check for corruption
      validations.append(self._check_integrity(file_path))
      
      # Format-specific validation
      validations.append(self._validate_content(file_path))
      
      return ValidationResult(validations)
  ```

#### 0.4 Registry Health Checks and Dataset Integrity ✅ COMPLETED
**File**: `src/protos/io/registry_health.py`

- [x] **Detect unregistered files**
  ```python
  class RegistryHealthChecker:
      def __init__(self, entity_registry: EntityRegistry, paths: ProtosPaths):
          self.registry = entity_registry
          self.paths = paths
      
      def find_unregistered_files(self, processor_type: str) -> List[Path]:
          """Find files in data directories not in registry."""
          unregistered = []
          
          # Get processor-specific directories
          processor = self._get_processor(processor_type)
          data_dirs = processor.get_data_directories()
          
          # Scan for files
          for data_dir in data_dirs:
              if data_dir.exists():
                  for file_path in data_dir.glob(processor.file_pattern):
                      # Check if file is registered
                      relative_path = file_path.relative_to(self.paths.data_root)
                      if not self.registry.find_by_path(str(relative_path)):
                          unregistered.append(file_path)
          
          return unregistered
      
      def check_dataset_integrity(self, dataset_name: str, 
                                processor_type: str) -> DatasetIntegrityReport:
          """Check if all dataset entities exist and are accessible."""
          report = DatasetIntegrityReport()
          dataset = self.dataset_manager.load_dataset(dataset_name)
          
          for entity_info in dataset['entities']:
              entity_id = entity_info['entity_id']
              entity_name = entity_info['entity_name']
              
              # Check if entity still registered
              if not self.registry.get_entity_by_id(entity_id):
                  report.add_missing_registration(entity_name)
                  continue
              
              # Check if file exists
              entity_data = self.registry.get_entity_by_id(entity_id)
              file_path = Path(entity_data['formats'][processor_type]['file_path'])
              full_path = self.paths.data_root / file_path
              
              if not full_path.exists():
                  report.add_missing_file(entity_name, file_path)
              elif not full_path.is_file() or not os.access(full_path, os.R_OK):
                  report.add_inaccessible(entity_name, file_path)
              else:
                  report.add_valid(entity_name)
          
          return report
  ```

- [ ] **CLI commands for health checks** (TODO - needs implementation)
  ```python
  @click.command()
  @click.option('--processor', type=click.Choice(['structure', 'sequence', 'grn', 'property']),
                help='Check specific processor type')
  def check_registry(processor):
      """Check registry health and find unregistered files."""
      checker = RegistryHealthChecker(EntityRegistry(), ProtosPaths())
      
      # Find unregistered files
      if processor:
          unregistered = checker.find_unregistered_files(processor)
          if unregistered:
              click.echo(f"\nFound {len(unregistered)} unregistered {processor} files:")
              for path in unregistered[:10]:  # Show first 10
                  click.echo(f"  - {path.name}")
              if len(unregistered) > 10:
                  click.echo(f"  ... and {len(unregistered) - 10} more")
      
      # Check all datasets
      click.echo("\nDataset integrity check:")
      for dataset in dataset_manager.list_datasets():
          report = checker.check_dataset_integrity(dataset)
          report.display_summary()
  ```

#### 0.5 Dataset Merging (Creates New Datasets) - PENDING
**Important**: Merging ALWAYS creates new datasets, never modifies existing ones

- [ ] **Generic dataset merging** (TODO - needs implementation)
  ```python
  def merge_datasets(self, source_datasets: List[str], 
                    new_dataset_name: str,
                    merge_options: Dict = None) -> str:
      """
      Merge multiple datasets into a NEW dataset.
      Original datasets remain unchanged.
      """
      if self.dataset_exists(new_dataset_name):
          raise ValueError(f"Dataset '{new_dataset_name}' already exists")
      
      # Collect all unique entities
      all_entities = set()
      for dataset_name in source_datasets:
          dataset = self.load_dataset(dataset_name)
          all_entities.update(e['entity_name'] for e in dataset['entities'])
      
      # Create new dataset
      metadata = {
          'merged_from': source_datasets,
          'merge_date': datetime.now().isoformat(),
          'merge_options': merge_options or {}
      }
      
      return self.create_dataset(
          new_dataset_name,
          list(all_entities),
          metadata=metadata
      )
  ```

- [ ] **Property table merging (special case)**
  ```python
  def merge_property_tables(self, table_names: List[str],
                           new_table_name: str,
                           column_strategy: str = "union") -> str:
      """
      Merge property tables into NEW table.
      
      Column strategies:
      - "union": Include all columns from all tables
      - "intersection": Only common columns
      - "left": Keep columns from first table
      """
      if self.table_exists(new_table_name):
          raise ValueError(f"Table '{new_table_name}' already exists")
      
      tables = [self.load_table(name) for name in table_names]
      
      # Handle column conflicts
      if column_strategy == "union":
          # Combine all columns, NaN for missing values
          merged_df = pd.concat(tables, axis=0, join='outer')
      elif column_strategy == "intersection":
          # Only keep common columns
          common_cols = set(tables[0].columns)
          for table in tables[1:]:
              common_cols &= set(table.columns)
          merged_df = pd.concat([t[list(common_cols)] for t in tables])
      
      # Remove duplicate entities (keep first occurrence)
      merged_df = merged_df[~merged_df.index.duplicated(keep='first')]
      
      # Save as new table
      self.save_table(new_table_name, merged_df)
      
      # Create corresponding dataset
      return self.create_dataset(
          new_table_name,
          list(merged_df.index),
          metadata={'table': new_table_name, 'column_strategy': column_strategy}
      )
  ```

- [ ] **GRN table merging (special case)**
  ```python
  def merge_grn_tables(self, table_names: List[str],
                      new_table_name: str,
                      require_same_family: bool = True) -> str:
      """
      Merge GRN tables into NEW table.
      
      Special requirements:
      - Tables must use same GRN positions
      - Optionally require same protein family
      """
      if self.table_exists(new_table_name):
          raise ValueError(f"Table '{new_table_name}' already exists")
      
      tables = []
      for name in table_names:
          table = self.load_table(name)
          metadata = self.get_table_metadata(name)
          
          if require_same_family:
              # Check protein family consistency
              if 'protein_family' in metadata:
                  if not all(m.get('protein_family') == metadata['protein_family'] 
                           for m in table_metadata):
                      raise ValueError("Tables contain different protein families")
          
          tables.append(table)
      
      # Check GRN position consistency
      grn_positions = set(tables[0].columns)
      for table in tables[1:]:
          if set(table.columns) != grn_positions:
              raise ValueError("Tables have different GRN positions")
      
      # Merge (remove duplicates)
      merged_df = pd.concat(tables)
      merged_df = merged_df[~merged_df.index.duplicated(keep='first')]
      
      # Save as new table
      self.save_table(new_table_name, merged_df)
      
      return self.create_dataset(
          new_table_name,
          list(merged_df.index),
          metadata={'table': new_table_name, 'grn_positions': list(grn_positions)}
      )
  ```

### Data Storage Patterns Across Processors

#### File-Based Processors Storage Details

1. **StructureProcessor** (Dual Format):
   - **Raw**: `mmcif/{entity_name}.cif` - Original structure files
   - **Cached**: `cache/{entity_name}.pkl` - Processed DataFrames with annotations
   - **Priority**: Always prefer PKL if available (10-100x faster)
   - **Datasets**: Concatenated into single PKL in `structure_dataset/`
   - **Registry tracks both formats**: Must register CIF and PKL as same entity

2. **SequenceProcessor** (Single Format):
   - **Only Format**: `fasta/{entity_name}.fasta` - No caching needed
   - **Multi-sequence**: Single FASTA can contain multiple sequences
   - **Datasets**: Not concatenated, remain as individual FASTA files
   - **Alignments**: Stored separately in JSON format

3. **EmbeddingProcessor** (Flexible Format):
   - **Formats**: `.pt` (PyTorch) > `.npy` (NumPy) > `.pkl` (Pickle)
   - **Organization**: `{model_name}/{entity_name}.{ext}`
   - **Datasets**: Can be concatenated into tensor dictionaries
   - **Metadata**: Always paired with JSON metadata files

#### Table-Based Processors Storage Details

4. **GRNProcessor** (Table Storage):
   - **Tables**: `tables/{table_name}.csv` - Multiple entities per file
   - **Entity**: Row in table, not separate file
   - **Reference**: `ref/` directory for standard tables

5. **PropertyProcessor** (Table Storage):
   - **Tables**: `tables/{table_name}.csv` - Properties as columns
   - **Entity**: Row in table with properties
   - **No individual files**: All data in tables

### Critical Implementation Notes

- [ ] **Handle dual formats in StructureProcessor**
  ```python
  def load_entity(self, name: str) -> pd.DataFrame:
      # Check cache first (preferred)
      cache_path = self.path_cache_dir / f"{name}.pkl"
      if cache_path.exists():
          df = pd.read_pickle(cache_path)
          # Register both CIF and PKL under same entity
          self._ensure_dual_registration(name)
          return df
      
      # Fall back to CIF
      cif_path = self.path_structure_dir / f"{name}.cif"
      if cif_path.exists():
          df = self._parse_cif(cif_path)
          # Save to cache for next time
          df.to_pickle(cache_path)
          # Register both formats
          self._register_dual_formats(name, cif_path, cache_path)
          return df
  ```

- [ ] **Track format priority in registry**
  ```python
  "formats": {
    "structure": {
      "cif": {
        "file_path": "structure/mmcif/1ubq.cif",
        "priority": 2,
        "content_hash": "abc123..."
      },
      "pkl": {
        "file_path": "structure/cache/1ubq.pkl", 
        "priority": 1,  # Lower number = higher priority
        "content_hash": "def456...",
        "derived_from": "cif"
      }
    }
  }
  ```

## Phase 0.6: IO Module Cleanup ✅ COMPLETED (Jan 19, 2025)

### Consolidation Work Done
- ✅ Created `format_registry.py` with centralized file format definitions
- ✅ Defined ProcessorType and FormatCategory enums
- ✅ Created FileFormat dataclass with all format metadata
- ✅ Updated InputManager to use format registry instead of hardcoded mappings
- ✅ Updated BaseProcessor._get_valid_extensions() to use format registry
- ✅ Updated RegistryHealthChecker to use format registry
- ✅ Enhanced file_utils.py with better type hints and documentation
- ✅ Started removing deprecated generate_entity_id usage

### Key Improvements
1. **Single source of truth** for file formats and extensions
2. **Type-safe enums** for processor and format categories
3. **Automatic format detection** based on extension and content
4. **Consistent validation** using existing utilities (cif_utils, file_utils)
5. **Better separation of concerns** between simple utilities and complex handlers

---

## Phase 1: Core Infrastructure Changes ✅ COMPLETED (Previously)

### 1.1 Update Entity Registry Structure ✅
**File**: `src/protos/io/entity_registry.py`

- [ ] **Add UUID generation**
  ```python
  def _generate_entity_id(self) -> str:
      """Generate stable UUID for entity."""
      return str(uuid.uuid4())
  ```

- [ ] **Update entity structure**
  ```python
  # Current structure uses hash as key
  # New structure:
  {
    "entities": {
      "550e8400-e29b-41d4-a716-446655440000": {  # UUID key
        "entity_id": "550e8400-e29b-41d4-a716-446655440000",
        "entity_name": "1UBQ",  # Primary human name
        "canonical_name": "UBIQ_HUMAN",  # Normalized form
        "aliases": [...],
        "formats": {...},
        "relationships": {...}
      }
    }
  }
  ```

- [ ] **Add name-to-ID index for fast lookups**
  ```python
  {
    "name_index": {
      "1UBQ": "550e8400-e29b-41d4-a716-446655440000",
      "UBIQ_HUMAN": "550e8400-e29b-41d4-a716-446655440000",
      "P62988": "550e8400-e29b-41d4-a716-446655440000"
    }
  }
  ```

- [ ] **Implement relationship tracking**
  ```python
  "relationships": {
    "derived_from": ["entity_id1", "entity_id2"],  # Multi-parent
    "derived_entities": ["entity_id3", "entity_id4"],
    "version_of": "entity_id_previous",
    "versions": ["entity_id_next1", "entity_id_next2"]
  }
  ```

### 1.2 Create Alias Resolution System
**File**: `src/protos/io/alias_resolver.py` (new)

- [ ] **Implement AliasResolver class**
  - Case-insensitive matching
  - Organism suffix handling (e.g., _HUMAN, _MOUSE)
  - UniProt ID normalization
  - Chain notation handling
  - Priority-based resolution

- [ ] **Add to EntityRegistry**
  ```python
  def _resolve_name_to_id(self, name: str) -> Optional[str]:
      """Resolve any name/alias to entity ID."""
      # Direct lookup first
      if name in self.name_index:
          return self.name_index[name]
      
      # Try alias resolution
      return self.alias_resolver.resolve(name, self.name_index)
  ```

### 1.3 Update Registry Methods
**File**: `src/protos/io/entity_registry.py`

- [ ] **Modify register_entity**
  ```python
  def register_entity(self, name: str, format_type: str, 
                     file_path: str, metadata: Dict) -> str:
      """Register entity, return human name (not ID)."""
      entity_id = self._generate_entity_id()
      # Store with UUID key
      # Update name_index
      # Return human name to user
  ```

- [ ] **Add relationship methods**
  ```python
  def add_relationship(self, from_name: str, to_name: str, 
                      rel_type: str, metadata: Dict = None)
  def get_derived_entities(self, name: str) -> List[str]
  def get_parent_entities(self, name: str) -> List[str]
  def get_entity_lineage(self, name: str) -> Dict
  ```

- [ ] **Add ID/name translation helpers**
  ```python
  def get_entity_name(self, entity_id: str) -> Optional[str]
  def get_entity_id(self, name: str) -> Optional[str]
  ```

## Phase 2: Dataset Manager Updates

### 2.1 Store Entity IDs Internally
**File**: `src/protos/io/dataset_manager.py`

- [ ] **Update dataset structure**
  ```python
  {
    "dataset_id": "660e8400-e29b-41d4-a716-446655440001",  # Add UUID
    "name": "kinase_study",
    "processor_type": "structure",
    "entities": [  # Store IDs, not names
      {
        "entity_id": "550e8400-e29b-41d4-a716-446655440000",
        "entity_name": "EGFR",  # Cache for performance
        "added": "2024-01-15T10:30:00Z"
      }
    ]
  }
  ```

- [ ] **Update create_dataset**
  ```python
  def create_dataset(self, name: str, entity_names: List[str], 
                    metadata: Dict = None) -> str:
      # Generate dataset UUID
      dataset_id = str(uuid.uuid4())
      
      # Convert names to IDs
      entities = []
      for name in entity_names:
          entity_id = self.entity_registry.get_entity_id(name)
          if entity_id:
              entities.append({
                  "entity_id": entity_id,
                  "entity_name": name,  # Cache current name
                  "added": datetime.now().isoformat()
              })
  ```

- [ ] **Update all methods to translate IDs**
  - load_dataset: Return dict with human names as keys
  - add_to_dataset: Accept human names, store IDs
  - remove_from_dataset: Accept human names, remove by ID
  - get_dataset_entities: Return human names

## Phase 2.5: Processor Interface Standardization

### Naming Convention Issues to Address

1. **Inconsistent method names across processors**:
   - StructureProcessor: `load_structure()` 
   - SequenceProcessor: `load_sequence()`
   - GRNBaseProcessor: `load_grn_table()`
   - PropertyProcessor: `load_property_table()`
   - **Solution**: Standardize to consistent pattern

2. **Table-based vs File-based processors**:
   - **File-based**: Structure, Sequence, Embedding (entity = file)
   - **Table-based**: GRN, Property (entity = row in table)
   - **Solution**: Create two base classes or interfaces

### Proposed Standard Interface

#### For File-Based Processors:
```python
class FileBasedProcessor(BaseProcessor):
    # Standard entity operations
    def load(self, name: str) -> Any  # Replaces load_structure, load_sequence
    def save(self, name: str, data: Any, metadata: Dict = None)
    
    # Keep abstract methods from BaseProcessor
    def load_entity(self, name: str) -> Any  # Delegates to load()
    def save_entity(self, name: str, data: Any, metadata: Dict = None)  # Delegates to save()
    
    # Batch operations
    def load_multiple(self, names: List[str]) -> Dict[str, Any]
    def save_multiple(self, data_dict: Dict[str, Any])
```

#### For Table-Based Processors:
```python
class TableBasedProcessor(BaseProcessor):
    # Table operations (primary interface)
    def load_table(self, table_name: str) -> pd.DataFrame
    def save_table(self, table_name: str, df: pd.DataFrame, metadata: Dict = None)
    
    # Entity operations (row-based)
    def load_entity(self, name: str) -> pd.Series  # Returns row from table
    def save_entity(self, name: str, data: pd.Series, table_name: str)  # Updates row
    
    # Query operations
    def query_entities(self, table_name: str, filter_dict: Dict) -> pd.DataFrame
    def get_entity_from_table(self, table_name: str, entity_name: str) -> pd.Series
```

### Standardization Plan

- [ ] **Rename methods for consistency**:
  - `load_structure()` → `load()`
  - `save_structure()` → `save()`
  - `load_sequence()` → `load()`
  - `load_grn_table()` → `load_table()`
  - `load_property_table()` → `load_table()`

- [ ] **Add deprecation warnings**:
  ```python
  def load_structure(self, name: str):
      warnings.warn("load_structure is deprecated, use load()", DeprecationWarning)
      return self.load(name)
  ```

- [ ] **Update processor class names for consistency**:
  - StructureProcessor → StructureProcessor ✓ (already good)
  - SequenceProcessor → SequenceProcessor ✓ (already good)
  - GRNBaseProcessor → GRNProcessor (remove "Base")
  - PropertyProcessor → PropertyProcessor ✓ (already good)
  - EmbeddingProcessor → EmbeddingProcessor ✓ (already good)

## Phase 3: BaseProcessor Integration

### 3.1 Update BaseProcessor Methods
**File**: `src/protos/core/base_processor.py`

- [ ] **Ensure all entity operations hide IDs**
  ```python
  def list_entities(self) -> List[str]:
      """List entities - returns human names only."""
      # Get all entity IDs for this format
      entity_ids = self.entity_registry.list_entity_ids(self.processor_type)
      
      # Convert to human names
      names = []
      for entity_id in entity_ids:
          name = self.entity_registry.get_entity_name(entity_id)
          if name:
              names.append(name)
      return sorted(names)
  ```

- [ ] **Add relationship handling**
  ```python
  def create_derived_entity(self, parent_name: str, derived_name: str, 
                           data: Any, metadata: Dict = None):
      """Save entity and track derivation relationship."""
      # Save the entity
      self.save_entity(derived_name, data, metadata)
      
      # Add relationship
      self.entity_registry.add_relationship(
          from_name=parent_name,
          to_name=derived_name,
          rel_type="derived",
          metadata=metadata
      )
  ```

- [ ] **Add version management**
  ```python
  def save_entity_version(self, name: str, data: Any, 
                         version_tag: str, metadata: Dict = None):
      """Save new version of entity."""
      # Create versioned name
      versioned_name = f"{name}_{version_tag}"
      
      # Save entity
      self.save_entity(versioned_name, data, metadata)
      
      # Track version relationship
      self.entity_registry.add_relationship(
          from_name=name,
          to_name=versioned_name,
          rel_type="version",
          metadata={"version": version_tag}
      )
  ```

### 3.2 Update Processor Implementations

- [ ] **StructureProcessor** (`src/protos/processing/structure/structure_processor.py`)
  - Update load_entity to work with new registry
  - Update save_entity to register with UUID
  - Add extract_chain method with relationship tracking

- [ ] **SeqProcessor** (`src/protos/processing/sequence/sequence_processor.py`)
  - Similar updates as StructureProcessor
  - Track source when sequence extracted from structure

- [ ] **GRNBaseProcessor** (`src/protos/processing/grn/grn_processor.py`)
  - Handle multi-entity tables
  - Track which entities contribute to GRN alignment

- [ ] **PropertyProcessor** (`src/protos/processing/property/property_processor.py`)
  - Link properties to entity IDs internally
  - Return properties with human names

## Phase 4: Migration and Compatibility

### 4.1 Create Migration Tools
**File**: `src/protos/io/migration.py` (new)

- [ ] **Registry migration v1 to v2**
  ```python
  def migrate_registry_to_uuid():
      """Convert hash-based registry to UUID-based."""
      # Load old registry
      # Generate UUIDs for each entity
      # Create name_index
      # Save new format
  ```

- [ ] **Dataset migration**
  ```python
  def migrate_datasets_to_uuid():
      """Update datasets to use entity IDs."""
      # Convert entity lists to ID-based
      # Add dataset UUIDs
      # Update format
  ```

### 4.2 Backward Compatibility Layer

- [ ] **Add compatibility flag to EntityRegistry**
  ```python
  def __init__(self, paths: Optional[ProtosPaths] = None, 
               legacy_mode: bool = False):
      """Support legacy hash-based IDs if needed."""
  ```

- [ ] **Auto-detection of registry version**
  ```python
  def _detect_registry_version(self) -> str:
      """Detect if registry is v1 (hash) or v2 (UUID)."""
  ```

## Phase 5: Testing and Validation

### 5.1 Create Comprehensive Tests
**File**: `tests/test_entity_system_v2.py` (new)

- [ ] **Test UUID generation and stability**
- [ ] **Test alias resolution**
- [ ] **Test relationship tracking**
- [ ] **Test dataset operations with IDs**
- [ ] **Test migration tools**
- [ ] **Test backward compatibility**

### 5.2 Integration Tests

- [ ] **Test complete workflows**
  ```python
  # User saves structure
  proc.save_entity("my_kinase", structure)
  
  # Extract chain (creates relationship)
  chain_a = proc.extract_chain("my_kinase", "A")
  proc.save_entity("my_kinase_chain_A", chain_a)
  
  # Check relationships
  derived = registry.get_derived_entities("my_kinase")
  assert "my_kinase_chain_A" in derived
  ```

## Phase 6: Documentation Updates

### 6.1 Update User Documentation

- [ ] **Update all examples to show new patterns**
- [ ] **Document relationship tracking**
- [ ] **Document version management**
- [ ] **Add migration guide**

### 6.2 Update Developer Documentation

- [ ] **Document new registry structure**
- [ ] **Document UUID system**
- [ ] **Document alias resolution**
- [ ] **Document relationship queries**

## Implementation Order

1. **Core Infrastructure** (Phase 1)
   - EntityRegistry updates
   - Alias resolver
   - Basic UUID support

2. **Dataset and BaseProcessor** (Phases 2-3)
   - Dataset UUID support
   - BaseProcessor integration
   - Relationship methods

3. **Processor Updates** (Phase 3.2)
   - Update all processor implementations
   - Add relationship tracking

4. **Migration and Testing** (Phases 4-5)
   - Migration tools
   - Comprehensive testing
   - Bug fixes

5. **Documentation and Polish** (Phase 6)
   - Update all documentation
   - Final testing

## Key Design Decisions

1. **UUIDs over content hashes**: More stable, allows renaming
2. **Name index for performance**: Fast lookups without scanning
3. **Cached names in datasets**: Performance over perfect consistency
4. **Multi-parent relationships**: Real-world use cases need this
5. **Backward compatibility**: Smooth migration path

## Success Criteria

- [ ] Users never see entity IDs in any output
- [ ] All existing code continues to work
- [ ] Relationships are properly tracked
- [ ] Performance is maintained or improved
- [ ] Migration is smooth and reversible
- [ ] Documentation is comprehensive

## Example Usage After Implementation

```python
# User experience remains clean and simple
struct_proc = CifBaseProcessor()
seq_proc = SeqProcessor()

# Save structure
struct_proc.save_entity("therapeutic_target", structure_data)

# Extract sequence (automatic relationship)
sequence = struct_proc.extract_sequence("therapeutic_target")
seq_proc.save_entity("therapeutic_target", sequence)

# Create mutant version
mutant = create_mutation(structure_data, "A123V")
struct_proc.save_entity_version("therapeutic_target", mutant, "A123V")

# Query relationships
derived = struct_proc.get_derived_entities("therapeutic_target")
# Returns: ["therapeutic_target_A123V"]

# Create dataset (uses IDs internally)
struct_proc.create_dataset("drug_targets", [
    "therapeutic_target",
    "therapeutic_target_A123V",
    "reference_structure"
])

# Load dataset (returns human names)
targets = struct_proc.load_dataset("drug_targets")
for name, structure in targets.items():
    print(f"Processing {name}")  # Human-readable names
```

---

## Future Enhancements

### Rename Propagation System
Currently, after `rename_entity()`, datasets show the new name because `get_dataset_entities()` resolves IDs dynamically. However, the cached names in dataset JSON files become stale. We've added:
- `refresh_dataset_entities(dataset_name)` - Updates cached names for a single dataset
- `refresh_all_datasets()` - Updates all datasets for a processor

Consider implementing:
1. **Automatic propagation**: When `rename_entity()` is called, automatically refresh all affected datasets
2. **Lazy propagation**: Mark datasets as "stale" and refresh on next access
3. **Background job**: Periodic refresh of all dataset caches
4. **Event system**: Processors subscribe to entity rename events

This ensures consistency between the registry and dataset files while maintaining performance.