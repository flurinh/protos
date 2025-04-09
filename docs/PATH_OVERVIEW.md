# Protos Path Handling Architecture Overview

> **IMPORTANT BUG NOTICE**: There is a critical path resolution bug in `ProtosPaths._resolve_data_root()` where the `AUTO` source defaults to using `ref_data_root` instead of `user_data_root`. This causes dataset paths to incorrectly resolve to locations within the protos package directory instead of the user's project directory. See the [Known Issues](#known-issues) section for details and workarounds.

This document provides a detailed explanation of the path handling architecture in the Protos framework, focusing on the interactions between different components:

- **ProtosPaths**: Core path resolution system
- **BaseProcessor**: Abstract processor base class
- **CifBaseProcessor**: Structure-specific processor
- **DatasetManager**: Dataset management system
- **Registry System**: Data catalog and indexing

## 1. Path Resolution Flow

The path resolution in Protos follows this hierarchical flow:

```
ProtosPaths → BaseProcessor → StructBaseProcessor → DatasetManager → Registry
```

Each component builds on the path resolution capabilities of the previous one, with increasing specialization for specific data types.

## 2. ProtosPaths

`ProtosPaths` is the fundamental path management class in the system:

- **Location**: `protos/io/paths/path_config.py`
- **Purpose**: Centralizes all path resolution logic
- **Key Responsibilities**:
  - Defines standard directory structure
  - Resolves environment variables
  - Handles relative/absolute path conversion
  - Distinguishes between user and reference data

### Key ProtosPaths Methods:

```python
# Initialize with custom root directories
paths = ProtosPaths(
    user_data_root="/path/to/user/data",    # Where user data is stored
    ref_data_root="/path/to/reference/data", # Where reference data is stored
    create_dirs=True,                        # Create directories if missing
    validate=True                            # Validate directory structure
)

# Get processor-specific directory
structure_dir = paths.get_processor_path("structure", DataSource.USER)

# Get structure subdirectory
mmcif_dir = paths.get_structure_subdir_path("structure_dir", DataSource.USER)

# Get registry path
registry_path = paths.get_registry_path("structure", DataSource.USER)
```

## 3. BaseProcessor

`BaseProcessor` is the abstract base class for all processors in Protos:

- **Location**: `protos/core/base_processor.py`
- **Purpose**: Provides common functionality for data processing and path handling
- **Key Responsibilities**:
  - Initializes path resolution using ProtosPaths
  - Handles dataset operations through standardized API
  - Manages registry interaction
  - Provides methods for saving/loading data

### Initialization Flow:

1. **Path Resolution**:
   ```python
   # In BaseProcessor.__init__
   self.processor_type = self._get_processor_type()
   self.path_resolver = ProtosPaths(user_data_root=data_root)
   self.data_root = self.path_resolver.user_data_root
   self.data_path = self.path_resolver.get_processor_path(self.processor_type)
   ```

2. **Dataset Manager Setup**:
   ```python
   # In BaseProcessor.__init__
   self.dataset_manager = DatasetManager(
       processor_type=self._get_processor_type(),
       paths=self.path_resolver
   )
   ```

3. **Registry Loading**:
   ```python
   # In BaseProcessor._load_dataset_registry
   registry_path = self.path_resolver.get_registry_path(self._get_processor_type())
   ```

## 4. CifBaseProcessor

`CifBaseProcessor` specializes `BaseProcessor` for structural data:

- **Location**: `protos/processing/structure/struct_base_processor.py`
- **Purpose**: Handles structure-specific operations
- **Key Responsibilities**:
  - Defines structure-specific directories (mmcif, alignments, etc.)
  - Initializes structure-specific paths
  - Provides structure-specific dataset methods
  - Maps generic methods to structure concepts

### Initialization Flow:

1. **Base Initialization**:
   ```python
   # In CifBaseProcessor.__init__
   super().__init__(name=name, data_root=data_root, processor_data_dir=processor_data_dir)
   ```

2. **Structure-Specific Path Setup**:
   ```python
   # In CifBaseProcessor.__init__
   self.path_structure_dir = os.path.join(self.data_root, processor_data_dir, structure_dir)
   self.path_dataset_dir = os.path.join(self.data_root, processor_data_dir, dataset_dir)
   self.path_alignment_dir = os.path.join(self.data_root, processor_data_dir, alignments_dir)
   ```

## 5. DatasetManager

`DatasetManager` handles standardized dataset operations:

- **Location**: `protos/core/dataset_manager.py`
- **Purpose**: Centralizes dataset operations
- **Key Responsibilities**:
  - Creates and manages standard dataset format
  - Handles dataset loading/saving
  - Interacts with both processor-specific and global registries
  - Consolidates dataset operations across processor types

### Initialization Flow:

1. **Path Setup**:
   ```python
   # In DatasetManager.__init__
   self.paths = paths or ProtosPaths()
   self.global_registry = GlobalRegistry(self.paths)
   registry_path = self.paths.get_registry_path(processor_type)
   self.registry = DataRegistry(registry_path)
   self.dataset_dir = self._get_dataset_dir()
   ```

2. **Dataset Directory Resolution**:
   ```python
   # In DatasetManager._get_dataset_dir
   if self.processor_type == 'structure':
       return self.paths.get_structure_subdir_path('dataset_dir')
   elif self.processor_type == 'grn':
       return self.paths.get_grn_subdir_path('grn_dir')
   # ...
   ```

## 6. Registry System

The registry system catalogs data assets:

- **Locations**: 
  - `protos/io/data_access.py` (GlobalRegistry, DataRegistry)
  - Processor-specific registry.json files
  - Global registry.json file
- **Purpose**: Index and catalog datasets and other assets
- **Key Responsibilities**:
  - Maps dataset IDs to file paths
  - Stores metadata for datasets
  - Enables cross-processor data relationships
  - Supports dataset discovery

### Registry Flow:

1. **Dataset Registration**:
   ```python
   # Within DatasetManager.create_dataset
   self.registry.register_dataset(
       dataset_id=dataset_id,
       file_path=file_path,
       metadata={...}
   )
   
   self.global_registry.register_dataset(
       dataset_id=dataset_id, 
       file_path=file_path, 
       processor_type=self.processor_type,
       dataset_type=self.processor_type,
       metadata={...}
   )
   ```

2. **Dataset Path Resolution**:
   ```python
   # Within DatasetManager.load_dataset
   file_path = self.global_registry.get_dataset_path(dataset_id)
   
   if not file_path:
       # Try local registry as fallback
       file_path = self.registry.get_dataset_path(dataset_id)
   ```

## 7. Correct Usage Pattern

### In Application Code:

```python
# 1. Configure paths properly
project_dir = os.path.abspath("/path/to/project")
paths = ProtosPaths(user_data_root=project_dir, create_dirs=True)

# 2. Get processor-specific directories
structure_dir = paths.get_processor_path("structure", DataSource.USER)
mmcif_dir = os.path.join(structure_dir, "mmcif")

# 3. Initialize processor with proper paths
processor = CifBaseProcessor(
    name="my_processor",
    data_root=project_dir,
    processor_data_dir="structure"
)

# 4. Create and manage datasets through standardized API
dataset = processor.create_standard_dataset(
    dataset_id="my_dataset",
    name="My Dataset",
    description="Contains selected structures",
    content=structure_ids,
    metadata={"created_by": "my_application"}
)
```

## 8. Common Issues and Solutions

### Issue 1: Incorrect Path Resolution

**Problem**: Processor initializes with incorrect paths, leading to data being saved in unexpected locations.

**Solution**:
- Always use absolute paths for `data_root`
- Initialize `ProtosPaths` explicitly with the desired root directory
- Use `DataSource.USER` when writing data

### Issue 2: Registry Desynchronization

**Problem**: Local and global registries get out of sync, causing datasets to be invisible or inaccessible.

**Solution**:
- Always use the standardized dataset API (`create_standard_dataset`, etc.)
- Use `DatasetManager` directly for batch operations on datasets
- Avoid manually editing registry files

### Issue 3: Path Inconsistency Between Components

**Problem**: Different path resolution logic between components leads to fragmented data.

**Solution**:
- Share a single `ProtosPaths` instance across components
- Use component-specific path methods (`get_structure_subdir_path` vs generic path joining)
- Validate paths during initialization with `validate=True`

## 9. Best Practices

1. **Always use absolute paths** for `data_root` to avoid relative path confusion.
2. **Initialize `ProtosPaths` explicitly** rather than relying on defaults.
3. **Use `DataSource.USER` explicitly** for writing operations.
4. **Share path resolvers** between components rather than creating new ones.
5. **Use the standardized dataset API** rather than direct file operations.
6. **Validate paths** periodically with `paths.validate_directory_structure()`.
7. **Log resolved paths** at initialization time for debugging.

## 10. Known Issues and Roadmap for Improvement

The current path handling system has several critical issues that need to be addressed for reliable operation, especially if Protos is to be used with the Model Context Protocol (MCP) for LLM integration.

### Issue 1: AUTO Source Mode Uses ref_data_root Instead of user_data_root

**Problem:**

In `ProtosPaths._resolve_data_root()`, the default behavior when `DataSource.AUTO` is used is to return `ref_data_root` rather than `user_data_root`:

```python
def _resolve_data_root(self, source: DataSource) -> str:
    if source == DataSource.USER:
        return self.user_data_root
    elif source == DataSource.REFERENCE:
        return self.ref_data_root
    else:  # AUTO
        # Default to reference data for reading, but this should be refined
        # based on actual operation being performed
        return self.ref_data_root  # <-- BUG: Should use user_data_root in many cases
```

This causes `DatasetManager` instances to incorrectly default to using paths within the protos package directory (e.g., `C:\Users\hidbe\PycharmProjects\phd\protos\src\protos\data\structure\structure_dataset`) rather than the user's project directory.

**Impact:**

* Datasets get saved in unexpected locations
* Registry files get created in the wrong directory
* The DatasetManager's path differs from the processor's path

**Current Workarounds:**

1. **Option 1: Subclass ProtosPaths**

```python
class FixedProtosPaths(ProtosPaths):
    def _resolve_data_root(self, source: DataSource) -> str:
        if source == DataSource.USER:
            return self.user_data_root
        elif source == DataSource.REFERENCE:
            return self.ref_data_root
        else:  # AUTO
            return self.user_data_root  # Use user_data_root instead of ref_data_root
```

2. **Option 2: Always Specify DataSource.USER**

```python
# When getting processor paths
structure_dir = paths.get_processor_path("structure", DataSource.USER)

# When getting structure subdirectory paths
mmcif_dir = paths.get_structure_subdir_path("structure_dir", DataSource.USER)
```

3. **Option 3: Fix DatasetManager paths after initialization**

```python
# Fix dataset manager paths after processor initialization
processor.dataset_manager.dataset_dir = os.path.join(project_dir, "structure", "structure_dataset")
processor.dataset_manager.paths = your_correct_paths_instance
```

### Issue 2: DatasetManager Creates New ProtosPaths Instance

**Problem:**

When `BaseProcessor` initializes the `DatasetManager`, it correctly passes its `path_resolver`, but later in the chain, the `DatasetManager._get_dataset_dir()` method might still use its own path resolution logic, which can be affected by the AUTO source issue above.

**Current Workaround:**

After initializing a processor, explicitly correct the `DatasetManager.dataset_dir` path:

```python
# Create the processor with the correct data_root
processor = CifBaseProcessor(
    name="my_processor",
    data_root=project_dir
)

# Fix the DatasetManager paths
correct_dataset_dir = os.path.join(project_dir, "structure", "structure_dataset")
processor.dataset_manager.dataset_dir = correct_dataset_dir
```

## 11. Roadmap for Path System Improvements

For Protos to become more robust and usable, especially via MCP for LLM integration, the following improvements to the path system are planned:

### 1. Fix AUTO Source Mode

The immediate priority is to fix the default behavior of `DataSource.AUTO`:

```python
def _resolve_data_root(self, source: DataSource) -> str:
    if source == DataSource.USER:
        return self.user_data_root
    elif source == DataSource.REFERENCE:
        return self.ref_data_root
    else:  # AUTO
        # Use a more intelligent strategy: check user_data_root first, 
        # fall back to ref_data_root only for read operations if user path doesn't exist
        operation_type = self.current_operation_type  # New context tracking
        path_to_check = os.path.join(self.user_data_root, self.current_path_fragment)
        
        if operation_type == "write" or os.path.exists(path_to_check):
            return self.user_data_root
        else:
            return self.ref_data_root
```

### 2. Simplified Initialization API

Create a simplified initialization API that abstracts away the complexity of path handling:

```python
# Current verbose initialization
paths = ProtosPaths(
    user_data_root=str(data_dir.absolute()),
    ref_data_root=str(data_dir.absolute()),
    create_dirs=True
)

cp = CifBaseProcessor(
    name="structure_processor",
    data_root=str(data_dir.absolute()),
    processor_data_dir="structure"
)

# Planned simplified API
protos = ProtosAPI.init_project("/path/to/project")
cp = protos.get_processor("structure")
```

### 3. Operation Context Tracking

Add context tracking to make automatic source resolution more intelligent:

```python
# Internal implementation
class ProtosPaths:
    def __init__(self, ...):
        self.current_operation_type = None  # 'read' or 'write'
        self.current_path_fragment = None   # relative path being processed
    
    def with_operation_context(self, operation_type, path_fragment):
        self.current_operation_type = operation_type
        self.current_path_fragment = path_fragment
        return self
    
    # Usage in function
    def get_processor_path(self, processor_type, source=DataSource.AUTO):
        with self.with_operation_context('read', f"{processor_type}"):
            return self._resolve_path(...)
```

### 4. Stateless API Layer for MCP

Develop a stateless API layer on top of the current processors for MCP integration:

```python
@mcp_tool
def get_structure_sequence(pdb_id: str, chain_id: str, project_dir: str = None) -> str:
    """Get the amino acid sequence from a structure."""
    # Automatically handle path configuration
    project_dir = project_dir or os.environ.get("PROTOS_DATA_ROOT") or os.getcwd()
    protos = ProtosAPI.init_project(project_dir)
    
    # Use the processor through the API
    return protos.get_sequence(pdb_id, chain_id)
```

### 5. Consistent Error Handling

Standardize error handling for path-related operations:

```python
class ProtosPathError(Exception):
    """Base exception for path-related errors."""
    pass

class PathNotFoundError(ProtosPathError):
    """Raised when a path does not exist."""
    def __init__(self, path, operation):
        self.path = path
        self.operation = operation
        super().__init__(f"Path not found: {path} (during {operation})")
```

These improvements will significantly enhance the reliability and usability of the Protos path system, especially for integration with MCP and LLMs.