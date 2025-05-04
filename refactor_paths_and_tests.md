# UniprotDL Path System Refactoring Plan

## Core Principles

1. **File Format Driven Organization**: Directories should be organized by file format, not by data source
   - All FASTA files go in `sequence/fasta/` regardless of source (UniProt, local, etc.)
   - No data source-specific directories like `uniprot/`

2. **Standardized Path Management**: Use the ProtosPaths system consistently
   - Follow the existing standards in `path_config.py` and `path_constants.py`
   - Remove hardcoded paths and custom directory structures

## Current Issues

The `UniprotDL` class has several path-related issues:

1. **Source-specific Directory Structure**:
   - Uses `data/uniprot/fastas/` instead of the standard `data/sequence/fasta/`
   - Creates a custom path hierarchy that doesn't align with the project standards

2. **Hardcoded Database Path**:
   - Uses `data/Blast/db/` for database files
   - This path isn't integrated with the ProtosPaths system

3. **Manual Directory Management**:
   - Manual directory creation without leveraging ProtosPaths tools
   - Inconsistent verification of directory existence before file operations

4. **Test Environment Mismatch**:
   - Test fixtures create temporary directories that don't match the expected structure
   - Tests aren't properly validating path integration

## Detailed Analysis of Path Systems

### ProtosPaths Standard Structure
From `path_constants.py`, the standard sequence directories are:

```python
DEFAULT_SEQUENCE_SUBDIRS = {
    "fasta_dir": "fasta",          # All FASTA files go here
    "alignment_dir": "alignments", # All alignment files go here
    "metadata_dir": "metadata"     # All metadata files go here
}
```

The full path structure should be:
```
data_root/
└── sequence/
    ├── fasta/       # ALL sequence files regardless of source (UniProt, etc.)
    ├── alignments/  # ALL alignment files 
    └── metadata/    # ALL metadata files including dataset definitions
```

### Current UniprotDL Structure
Currently uses:
```
data/uniprot/         # Source-specific directory (problem\!)
├── datasets/         # Dataset definitions
├── table/            # Gene data tables
└── fastas/           # FASTA files from UniProt
```

And:
```
data/Blast/db/        # Database files (another custom path)
```

## Refactoring Plan

### 1. Update UniprotDL Class

#### A. Refactor Constructor to Use ProtosPaths

```python
def __init__(self,
             data_root=None,        # New: Data root parameter
             dataset='gprot_2585',  # Dataset name stays the same
             reload=False,
             limit=10):
    """
    Initialize UniprotDL with standardized paths from ProtosPaths.
    
    Args:
        data_root: Root directory for data. If None, uses default from ProtosPaths.
        dataset: Name of the dataset file (without extension).
        reload: Whether to reload data even if it exists.
        limit: Maximum number of sequences to process.
    """
    from protos.io.paths.path_config import ProtosPaths, DataSource, join_path
    self.paths = ProtosPaths(user_data_root=data_root, create_dirs=True)
    
    # Set dataset information
    self.dataset = dataset
    
    # Use standard path locations from ProtosPaths for all files
    self.sequence_dir = self.paths.get_processor_path('sequence', source=DataSource.USER)
    self.fasta_dir = self.paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    self.metadata_dir = self.paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    
    # Store paths for files
    self.dataset_file = join_path(self.metadata_dir, f"{dataset}.txt")
    self.gene_table_file = join_path(self.metadata_dir, f"gene_df_{dataset}.csv")
    
    # Note: We no longer need separate fastas, datasets, and table directories
    # Everything follows the file format organization principle
    
    # Initialize containers
    self.genes = []  # UniProt IDs from dataset
    self.data_list = []  # Intermediate storage
    self.gene_df = pd.DataFrame()  # Complete data
    
    # Settings
    self.reload = reload
    self.limit = limit
    
    # Load available dataset information
    self.available_datasets = self._get_available_datasets()
```

#### B. Simplified Dataset Loading Method

```python
def load_dataset(self):
    """
    Load UniProt IDs from a dataset file in the standard metadata directory.
    
    Returns:
        List of UniProt IDs
    """
    if not os.path.exists(self.dataset_file):
        raise FileNotFoundError(f"Dataset file not found: {self.dataset_file}")
        
    with open(self.dataset_file) as file:
        genes = [line.rstrip().split(' ') for line in file]
    self.genes = list(itertools.chain(*genes))
    return list(dict.fromkeys(self.genes))
```

#### C. Update FASTA Saving Method

```python
def save_uniprot_fasta(self, uniprot=None, mode='entry', ext='.fasta',
                      force_overwrite=True, validate=True):
    """
    Save sequence data to FASTA files in the standard fasta directory.
    
    Args:
        uniprot: UniProt ID or list of IDs (default: use all in gene_df)
        mode: 'entry' for individual files, 'database' for combined file
        ext: File extension to use (should be .fasta for consistency)
        force_overwrite: Whether to overwrite existing files
        validate: Whether to validate FASTA format
        
    Returns:
        List of saved file paths
    """
    # Determine which IDs to save
    uids = list(self.gene_df['uniprot']) if uniprot is None else (
        uniprot if isinstance(uniprot, list) else [uniprot])
    
    # Verify IDs exist in dataframe
    missing_uids = [uid for uid in uids if uid not in list(self.gene_df['uniprot'])]
    if missing_uids:
        raise ValueError(f"UniProt IDs not found in data: {missing_uids}")
    
    saved_files = []
    
    if mode == 'entry':
        # Save individual FASTA files in the standard fasta directory
        for uid in uids:
            row = self.gene_df[self.gene_df['uniprot'] == uid]
            
            # Extract sequence and info
            seq = str(row['seq'].values[0])
            info = str(row['info'].values[0]) if 'info' in row.columns else uid
            
            # Create sequence dictionary
            sequence_dict = {uid: seq}
            
            # Determine file path in standard location
            file_path = os.path.join(self.fasta_dir, f"{uid}{ext}")
            
            # Create directory if needed
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            
            # Write using standard fasta_utils function
            from protos.io.fasta_utils import write_fasta
            write_fasta(sequence_dict, file_path)
            saved_files.append(file_path)
            
    else:  # database mode
        # Save a combined FASTA file with all sequences
        db_path = os.path.join(self.metadata_dir, f"{self.dataset}_db{ext}")
        
        # Create combined dictionary
        sequences = {}
        for uid in uids:
            row = self.gene_df[self.gene_df['uniprot'] == uid]
            seq = str(row['seq'].values[0])
            sequences[uid] = seq
        
        # Create directory if needed
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        
        # Write using standard fasta_utils function
        from protos.io.fasta_utils import write_fasta
        write_fasta(sequences, db_path)
        saved_files.append(db_path)
    
    return saved_files
```

#### D. Simplified Table Management Methods

```python
def save_gene_table(self):
    """
    Save gene data table to the standard metadata directory.
    """
    if self.gene_df.empty:
        return False
        
    # Use standard metadata directory
    os.makedirs(os.path.dirname(self.gene_table_file), exist_ok=True)
    self.gene_df.to_csv(self.gene_table_file, index=False)
    return True
    
def load_gene_table(self):
    """
    Load gene data table from the standard metadata directory.
    """
    if not os.path.exists(self.gene_table_file):
        return None
        
    return pd.read_csv(self.gene_table_file)
```

#### E. Helper Method for Available Datasets

```python
def _get_available_datasets(self):
    """
    Get list of available datasets in the metadata directory.
    """
    if not os.path.exists(self.metadata_dir):
        return []
        
    # Look for .txt files in the metadata directory
    files = [f for f in os.listdir(self.metadata_dir) 
             if f.endswith('.txt') and os.path.isfile(os.path.join(self.metadata_dir, f))]
    
    # Remove extension to get dataset names
    return [os.path.splitext(f)[0] for f in files]
```

### 2. Update Test Fixtures

```python
@pytest.fixture
def test_data_root():
    """Create a temporary directory with ProtosPaths structure for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize ProtosPaths with the temp directory
        from protos.io.paths.path_config import ProtosPaths
        paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True)
        
        # Return the temp directory path
        yield temp_dir

@pytest.fixture
def test_dataset(test_data_root, test_uniprot_ids):
    """Create a test dataset in the proper metadata directory."""
    from protos.io.paths.path_config import ProtosPaths, DataSource, join_path
    
    # Create paths object
    paths = ProtosPaths(user_data_root=test_data_root)
    
    # Create dataset name
    dataset_name = "test_dataset"
    
    # Get metadata directory
    metadata_dir = paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    
    # Create dataset file with test IDs
    dataset_file = join_path(metadata_dir, f"{dataset_name}.txt")
    os.makedirs(os.path.dirname(dataset_file), exist_ok=True)
    
    with open(dataset_file, 'w') as f:
        f.write(" ".join(test_uniprot_ids[:2]))  # Use only 2 IDs for speed
    
    # Return paths object and dataset name
    yield paths, dataset_name
```

### 3. Updated Test Functions

```python
def test_uniprot_loader_initialization(test_dataset):
    """Test initializing the UniProt loader."""
    paths, dataset_name = test_dataset
    
    # Initialize the loader with the test data root
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name
    )
    
    # Verify paths were set up according to ProtosPaths standard
    assert loader.dataset == dataset_name
    assert loader.fasta_dir == paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert loader.metadata_dir == paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    assert os.path.exists(loader.dataset_file)
    
    # Verify dataset file exists
    assert os.path.exists(loader.dataset_file)
    
def test_uniprot_loader_save_fasta(test_dataset, monkeypatch):
    """Test saving UniProt sequences as FASTA files."""
    paths, dataset_name = test_dataset
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name
    )
    
    # Mock the download_gene_single_query to avoid network calls
    def mock_download(uniprot):
        return [uniprot, "SEQUENCE", f"{uniprot}_GENE", "SPECIES", "ORGANISM", "INFO", dataset_name]
    
    monkeypatch.setattr(loader, 'download_gene_single_query', mock_download)
    
    # Load dataset and download genes
    loader.load_dataset()
    test_id = loader.genes[0]
    
    # Add mock data
    loader.data_list.append(mock_download(test_id))
    loader.gene_df = pd.DataFrame(loader.data_list, columns=['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset'])
    
    # Save as FASTA files
    saved_files = loader.save_uniprot_fasta(uniprot=test_id, mode='entry')
    
    # Verify files were saved in the standard fasta directory
    assert len(saved_files) == 1
    saved_path = saved_files[0]
    
    # Verify file location is in the standard fasta directory
    expected_fasta_dir = paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert os.path.dirname(saved_path) == expected_fasta_dir
    
    # Verify filename includes the UniProt ID
    assert test_id in os.path.basename(saved_path)
    
    # Verify file content
    from protos.io.fasta_utils import read_fasta
    sequences = read_fasta(saved_path)
    assert test_id in sequences
    assert sequences[test_id] == "SEQUENCE"
```

## Implementation Strategy

1. **Phase 1**: Update the `UniprotDL` constructor and core path configuration
   - Adopt ProtosPaths for all path management
   - Eliminate custom directory structure

2. **Phase 2**: Update file operations to use standard locations
   - Save FASTA files to the standard fasta directory
   - Store metadata in the standard metadata directory

3. **Phase 3**: Update tests to use ProtosPaths
   - Create proper test fixtures that match the production structure
   - Verify path integration works correctly

4. **Phase 4**: Clean up and document
   - Remove any remaining hardcoded paths
   - Add clear documentation about the path structure

## Key Benefits

1. **Format-Based Organization**: Files are organized by their format, not their source
2. **Simplicity**: Single location for each file type
3. **Consistency**: Follows the project's established path structure
4. **Maintainability**: Easier to understand and extend
5. **Reliability**: Proper directory creation and verification

## Things to Avoid

1. **Source-Based Directories**: Don't create directories like `uniprot/` for data sources
2. **Hardcoded Paths**: Avoid any hardcoded paths like `data/Blast/db/`
3. **Custom Extensions**: Use standard extensions (.fasta) where possible
4. **Manual Path Construction**: Use path helpers instead of string concatenation

# Review of Current test_uniprot_loader.py

After analyzing the current test file, I've identified several issues that need to be addressed in the refactoring:

## Current Issues in the Test File

1. **Custom Directory Structure**: 
   - Creates a temporary directory with a custom structure (`uniprot/datasets`, `uniprot/table`, `uniprot/fastas`)
   - This structure doesn't match the format-based organization standard

2. **Path Construction**: 
   - Manually constructs paths using string concatenation and Path objects
   - Doesn't leverage the ProtosPaths system

3. **File Saving Pattern**:
   - The `test_uniprot_loader_save_to_standard_fasta` test creates a "standard" fasta directory but directly uses `write_fasta` function
   - Should test the integrated path handling in the `UniprotDL` class instead

4. **Testing Redundancy**:
   - Has both `test_uniprot_loader_save_fasta` and `test_uniprot_loader_save_to_standard_fasta` 
   - After refactoring, these should be unified since all files should go to the standard location

5. **Hard Dependencies on UniprotDL Structure**:
   - Tests depend on internal attributes like `loader.path_fastas` and `loader.path_database`
   - After refactoring, these paths would be derived from ProtosPaths

## Detailed Test Refactoring Plan

### 1. Create New Test Fixtures Using ProtosPaths

Replace the current fixtures with ones that use ProtosPaths:

```python
@pytest.fixture
def test_data_root():
    """Create a temporary directory with ProtosPaths structure for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize ProtosPaths with temp directory
        from protos.io.paths.path_config import ProtosPaths
        paths = ProtosPaths(user_data_root=temp_dir, create_dirs=True)
        
        yield temp_dir

@pytest.fixture
def prepared_test_environment(test_data_root, test_uniprot_ids):
    """Create a properly structured test environment with dataset file."""
    from protos.io.paths.path_config import ProtosPaths, DataSource, join_path
    
    # Create paths object
    paths = ProtosPaths(user_data_root=test_data_root)
    
    # Create dataset name
    dataset_name = "test_dataset"
    
    # Get metadata directory where dataset files should be stored
    metadata_dir = paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    
    # Create dataset file with test IDs
    dataset_file = join_path(metadata_dir, f"{dataset_name}.txt")
    
    # Make sure directory exists
    os.makedirs(os.path.dirname(dataset_file), exist_ok=True)
    
    # Write test IDs to dataset file
    with open(dataset_file, 'w') as f:
        f.write(" ".join(test_uniprot_ids[:2]))  # Use only 2 IDs for speed
    
    # Return paths object and dataset name
    yield paths, dataset_name
```

### 2. Update Basic Test Functions

```python
def test_uniprot_loader_initialization(prepared_test_environment):
    """Test initializing the UniProt loader with ProtosPaths."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader with data_root
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=2
    )
    
    # Verify initialization
    assert loader.dataset == dataset_name
    assert loader.limit == 2
    
    # Verify paths were set up correctly
    assert loader.fasta_dir == paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert loader.metadata_dir == paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    assert os.path.exists(loader.dataset_file)

def test_uniprot_loader_load_dataset(prepared_test_environment, test_uniprot_ids):
    """Test loading a dataset using the loader."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name
    )
    
    # Load the dataset
    loaded_ids = loader.load_dataset()
    
    # Verify loaded IDs
    assert len(loaded_ids) == 2
    assert loaded_ids[0] in test_uniprot_ids
    assert loaded_ids[1] in test_uniprot_ids
```

### 3. Create Mock for Network Tests

Since we want to test the path handling without relying on network calls:

```python
@pytest.fixture
def mock_uniprot_data():
    """Create mock data for testing without network calls."""
    def _create_mock_data(uniprot_id, sequence="TESTSEQUENCE"):
        return {
            "uniprot": uniprot_id,
            "seq": sequence,
            "gene": f"{uniprot_id}_GENE",
            "species": "TEST_SPECIES",
            "organism": "Test Organism",
            "info": f"Test info for {uniprot_id}",
            "dataset": "test_dataset"
        }
    return _create_mock_data
```

### 4. Update FASTA Saving Tests

```python
def test_uniprot_loader_save_fasta(prepared_test_environment, test_uniprot_ids, mock_uniprot_data, monkeypatch):
    """Test saving sequences as FASTA files in the standard location."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=1
    )
    
    # Mock download function to avoid network calls
    def mock_download(uniprot_id):
        data = mock_uniprot_data(uniprot_id)
        return [
            data["uniprot"], 
            data["seq"], 
            data["gene"], 
            data["species"], 
            data["organism"], 
            data["info"], 
            data["dataset"]
        ]
    
    monkeypatch.setattr(loader, "download_gene_single_query", mock_download)
    
    # Load dataset and mock download a gene
    loader.load_dataset()
    test_id = loader.genes[0]
    
    # Populate the data
    loader.data_list.append(mock_download(test_id))
    loader.gene_df = pd.DataFrame(loader.data_list, columns=['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset'])
    
    # Save as FASTA file (should go to standard fasta directory)
    saved_files = loader.save_uniprot_fasta(uniprot=test_id, mode='entry')
    
    # Verify file was created in the standard location
    assert len(saved_files) == 1
    
    # Verify correct location
    standard_fasta_dir = paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
    assert os.path.dirname(saved_files[0]) == standard_fasta_dir
    
    # Verify filename and content
    assert test_id in os.path.basename(saved_files[0])
    
    # Read the file
    loaded_sequences = read_fasta(saved_files[0])
    assert test_id in loaded_sequences
    assert loaded_sequences[test_id] == "TESTSEQUENCE"
```

### 5. Update Integration Test

```python
def test_uniprot_loader_integrated(prepared_test_environment, test_uniprot_ids, mock_uniprot_data, monkeypatch):
    """Integration test for the UniProt loader with standardized paths."""
    paths, dataset_name = prepared_test_environment
    
    # Initialize the loader
    loader = UniprotDL(
        data_root=paths.user_data_root,
        dataset=dataset_name,
        limit=2
    )
    
    # Mock download functions to avoid network calls
    def mock_download_single(uniprot_id):
        data = mock_uniprot_data(uniprot_id)
        return [
            data["uniprot"], 
            data["seq"], 
            data["gene"], 
            data["species"], 
            data["organism"], 
            data["info"], 
            data["dataset"]
        ]
    
    def mock_download_batch(*args, **kwargs):
        # Populate loader.data_list and loader.gene_df with test data
        for uid in loader.genes:
            loader.data_list.append(mock_download_single(uid))
        
        loader.gene_df = pd.DataFrame(
            loader.data_list, 
            columns=['uniprot', 'seq', 'gene', 'species', 'organism', 'info', 'dataset']
        )
        return loader.gene_df
    
    monkeypatch.setattr(loader, "download_gene_single_query", mock_download_single)
    monkeypatch.setattr(loader, "download_genes_single_query", mock_download_batch)
    
    # 1. Load dataset
    loader.load_dataset()
    assert len(loader.genes) == 2
    
    # 2. Download genes
    loader.download_genes_single_query(batchsize=2, save=True)
    assert len(loader.gene_df) == 2
    
    # 3. Save as FASTA in both modes
    entry_files = loader.save_uniprot_fasta(mode='entry')
    db_file = loader.save_uniprot_fasta(mode='database')
    
    # 4. Verify individual files
    for uid in loader.genes:
        # Find the corresponding file
        uid_file = [f for f in entry_files if uid in os.path.basename(f)]
        assert len(uid_file) == 1
        
        # Verify the file is in the standard fasta directory
        standard_fasta_dir = paths.get_sequence_subdir_path('fasta_dir', source=DataSource.USER)
        assert os.path.dirname(uid_file[0]) == standard_fasta_dir
        
        # Verify content
        loaded_sequences = read_fasta(uid_file[0])
        assert uid in loaded_sequences
    
    # 5. Verify database file
    assert len(db_file) == 1
    db_path = db_file[0]
    
    # Verify database file location (should be in metadata directory)
    metadata_dir = paths.get_sequence_subdir_path('metadata_dir', source=DataSource.USER)
    assert os.path.dirname(db_path) == metadata_dir
    
    # Verify database content
    db_sequences = read_fasta(db_path)
    assert len(db_sequences) == 2
    for uid in loader.genes:
        assert uid in db_sequences
```

## Transition Strategy

1. **Update UniprotDL First**: Implement the refactored UniprotDL class using ProtosPaths
2. **Create New Test Functions**: Create the new test fixtures and functions that use ProtosPaths
3. **Run Both Test Sets**: Keep the original tests temporarily and run both during transition
4. **Replace Original Tests**: Once the new tests pass, replace the original tests
5. **Final Cleanup**: Remove any references to the custom path structure

## Key Benefits of Test Refactoring

1. **Consistency**: Tests match the production code's use of ProtosPaths
2. **Isolation**: Tests use mocking to avoid network dependencies
3. **Clarity**: Path construction is handled by ProtosPaths, not manual string concatenation
4. **Standards**: All tests follow the format-driven organization principle

By refactoring both the UniprotDL class and its tests, we'll ensure that sequences from UniProt are properly stored in the standard sequence directory structure, and that all code follows the project's established path standards.

## Test Error Analysis and Resolution

After running the refactored tests, we found two types of errors:

### 1. Mock File Creation Issue in end_to_end_download_workflow

In `test_end_to_end_download_workflow`, the test fails with:
```
AssertionError: assert False
  where False = exists('C:\\Users\\hidbe\\AppData\\Local\\Temp\\tmpxpduyn9s\\structure\\mmcif\\1ABC.cif')
```

**Problem**: While we mocked the `retrieve_pdb_file` function to return a path to the downloaded file, we didn't actually create the file at that location. When the test tries to verify the file exists, it fails.

**Solution**: We need to not only mock the function calls but also create the expected output files:

```python
# Setup mock PDB download
mock_file_path = join_path(mmcif_dir, "1ABC.cif")
mock_retrieve.return_value = mock_file_path

# Create the mock file to satisfy the existence check
os.makedirs(os.path.dirname(mock_file_path), exist_ok=True)
with open(mock_file_path, 'w') as f:
    f.write("Mock CIF content")
```

### 2. Mock Data Sequence Mismatch in uniprot_download_and_save_to_standard_location

In `test_uniprot_download_and_save_to_standard_location`, the test fails with:
```
AssertionError: assert 'MRPSGTAGAALL...LRVAPQSSEFIGA' == 'MKPSLLHLLFLG...QETFEYMVSQGLL'
```

**Problem**: The actual sequence in the file doesn't match the expected sequence we specified in our mock. This is likely because:
1. The function is using real data from a network call despite our mock
2. The mock isn't being applied correctly
3. The sequence transformation in the save process is altering the data

**Solution**: We need to ensure our mock is properly intercepting the network call and that the sequence data passes through unchanged:

1. Check how the mock is applied and ensure it's intercepting at the right point
2. Add more logging to see what data is actually flowing through the system
3. Possibly simplify the test case to use very short, distinctive test sequences that make it easy to spot changes

### Key Test Fixes Needed

1. **Create Mock Files**: When using mocks that return file paths, ensure the files exist at those paths.

2. **Mock Network Calls Completely**: Ensure all network calls are properly mocked and test with distinctive data that makes problems obvious.

3. **Simplified Test Data**: Use shorter, clearly identifiable test data to make debugging easier.

4. **Consistent Path Usage**: Ensure all file operations use the same path helpers (`join_path`) consistently.

5. **Verify Mocks Are Applied**: Check that mocks are actually being called and intercepting the network requests.

## Updated Test Implementation Plan

1. Create helper methods for setting up mock files:
```python
def create_mock_file(path, content="Mock content"):
    """Create a mock file with the specified content."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)
    return path
```

2. For the UniProt download test, revise the mock implementation:
```python
# Use simple, distinctive test sequences
TEST_SEQUENCES = {
    "P00533": "SEQUENCE_FOR_P00533",
    "P01308": "SEQUENCE_FOR_P01308"
}

# Setup the mock with exact control
def mock_get_data(uid, **kwargs):
    if uid not in TEST_SEQUENCES:
        return pd.DataFrame()
    
    return pd.DataFrame({
        'Entry': [uid],
        'Entry Name': [f"{uid}_TEST"],
        'Sequence': [TEST_SEQUENCES[uid]],
        'Protein names': [f"Test protein {uid}"],
        'Organism': ['Test organism']
    })

# Ensure the mock is applied correctly
mock_get_uniprot.side_effect = mock_get_data
```

3. Add checks to verify mock function calls:
```python
# Verify mock was called with correct parameters
assert mock_get_uniprot.call_count >= 1
mock_get_uniprot.assert_any_call("P00533", reviewed=True)

# Verify the data in the dataframe directly before saving
assert loader.gene_df.iloc[0]['seq'] == TEST_SEQUENCES["P00533"]
```

By implementing these fixes, we can ensure the tests are robust and correctly verify the functionality without relying on network access.
