# Protos Data Management: The Complete Guide

## Table of Contents
1. [Core Principles](#core-principles)
2. [ProtosPaths: The Foundation](#protospaths-the-foundation)
3. [Entity, Registry, and Dataset System](#entity-registry-and-dataset-system)
4. [Processors: Format-Specific Data Handlers](#processors-format-specific-data-handlers)
5. [File Naming and Organization](#file-naming-and-organization)
6. [Implementation Guidelines](#implementation-guidelines)
7. [Common Workflows](#common-workflows)

---

## Core Principles

### 1. One Path System: ProtosPaths
**ProtosPaths is the ONLY path management system in Protos.** No component ever specifies or manages paths directly.

### 2. Zero Configuration
Protos works **out of the box** with zero setup. No configuration files, no environment variables, no path specifications needed.

### 3. Human-Readable Filesystem
All files use **human-readable names**. Hash IDs exist only as internal registry keys.

### 4. Unified Entity Tracking
Every biological object is an **entity** tracked across all formats through a central registry.

---

## ProtosPaths: The Foundation

ProtosPaths is the bedrock of Protos - a unified data management system that handles all path operations. It ensures consistent data organization, enables zero-configuration usage, and provides a single source of truth for data locations.

### The Fundamental Law

```python
# Users DON'T need to do this - processors handle it automatically!
# (Shown here only to illustrate the internal cascade)

# The initialization cascade that happens inside processors:
paths = ProtosPaths()              # 1. Path management
registry = EntityRegistry(paths)    # 2. Entity tracking  
dataset_manager = DatasetManager(paths)  # 3. Dataset management

# What users ACTUALLY do:
from protos.processing.structure import CifBaseProcessor
processor = CifBaseProcessor()     # Everything above happens automatically!
```

### What ProtosPaths Manages

```
working_dir/
└── data/                          # Default base (or user-specified)
    ├── entity_registry.json       # Central entity tracking
    ├── structure/                 # Structure data
    │   ├── mmcif/                # PDB/CIF files
    │   ├── cache/                # Processed individual structures (PKL)
    │   ├── structure_dataset/    # Processed dataset data (PKL files, NOT JSON!)
    │   ├── alignments/           # Structure alignments
    │   ├── datasets/             # Dataset JSON definitions
    │   │   ├── kinases.json      # Example dataset
    │   │   └── opsins.json       # Example dataset
    │   └── registry.json         # Structure processor registry
    ├── sequence/                 # Sequence data
    │   ├── fasta/               # FASTA files
    │   ├── alignments/          # All alignment results
    │   │   ├── pairwise/       # Pairwise alignments
    │   │   ├── multiple/       # Multiple sequence alignments
    │   │   └── mmseqs/         # MMseqs2 alignments
    │   ├── databases/           # MMseqs2 databases
    │   ├── metadata/            # Sequence metadata
    │   ├── datasets/            # Dataset JSON definitions
    │   └── registry.json        # Sequence processor registry
    ├── grn/                     # GRN data
    │   ├── tables/              # GRN annotation tables
    │   ├── reference/           # Reference GRN tables
    │   ├── configs/             # GRN configurations
    │   ├── assignments/         # GRN assignments
    │   ├── temp/                # Temporary files
    │   ├── datasets/            # Dataset JSON definitions
    │   └── registry.json        # GRN processor registry
    ├── property/                # Property data
    │   ├── tables/              # Property tables
    │   ├── datasets/            # Dataset JSON definitions
    │   └── registry.json        # Property processor registry
    ├── embedding/               # ML embeddings
    │   ├── embeddings/          # Saved embeddings
    │   ├── datasets/            # Dataset JSON definitions
    │   └── registry.json        # Embedding processor registry
    ├── ligand/                  # Ligand data
    │   ├── sdf/                 # SDF/MOL files
    │   ├── cache/               # Cached ligand data
    │   ├── datasets/            # Dataset JSON definitions
    │   └── registry.json        # Ligand processor registry
    ├── graph/                   # Graph/network data
    │   ├── networks/            # Graph/network files
    │   ├── analysis/            # Graph analysis results
    │   ├── datasets/            # Dataset JSON definitions
    │   └── registry.json        # Graph processor registry
    └── temp/                    # Temporary processing files
```

### ProtosPaths API

```python
class ProtosPaths:
    """Central path management for Protos."""
    
    def __init__(self, base_path: Optional[Path] = None):
        """
        Initialize paths.
        
        Args:
            base_path: Base data directory. Defaults to working_dir/data
        """
        self.base_path = base_path or Path.cwd() / "data"
        self.ensure_directories()
    
    # Core path methods
    def get_base_path(self) -> Path:
        """Get base data directory."""
        
    def get_processor_path(self, processor_type: str) -> Path:
        """Get processor-specific directory (e.g., data/structure)."""
        
    def get_registry_path(self, filename: str = "entity_registry.json") -> Path:
        """Get path to registry file."""
        
    def get_temp_path(self) -> Path:
        """Get temporary files directory."""
        
    def ensure_directories(self):
        """Create all required directories if they don't exist."""
```

### Implementation Rules

1. **Every component gets ProtosPaths**
   ```python
   def __init__(self, paths: Optional[ProtosPaths] = None):
       self.paths = paths or ProtosPaths()
   ```

2. **Never construct paths manually**
   ```python
   # ❌ WRONG
   filepath = Path("data/structure/mmcif/1ubq.cif")
   
   # ✅ CORRECT
   filepath = self.paths.get_processor_path("structure") / "mmcif" / "1ubq.cif"
   ```

3. **No path parameters in APIs**
   ```python
   # ❌ WRONG
   def load_structure(self, pdb_id: str, path: str):
   
   # ✅ CORRECT
   def load_structure(self, pdb_id: str):
       # Use self.paths internally
   ```

---

## Entity, Registry, and Dataset System

The management system consists of three interconnected components that work together to track and organize biological data.

### Entities: The Core Abstraction

An **entity** represents a single biological object (protein, complex, sequence) that can exist in multiple data formats.

#### Key Concepts

1. **Human-Readable Names**: What users see and use
   - PDB IDs: `1ubq`, `3SN6`, `7ZVL`
   - UniProt IDs: `P62988`, `EGFR_HUMAN`
   - Custom names: `my_protein_v2`, `kinase_complex_atp`

2. **Hash IDs**: Internal registry keys only
   - Generated once when entity is registered
   - Never shown to users
   - Never used in filenames
   - Only used for internal cross-referencing

#### Entity Structure in Registry

```json
{
  "entities": {
    "3f8a9c2d1e": {                      // Hash ID - internal only
      "original_id": "1ubq",              // Human-readable primary name
      "aliases": ["UBIQ_HUMAN", "P62988"], // Alternative names
      "formats": {
        "structure": {
          "file_path": "structure/mmcif/1ubq.cif",  // Human-readable path
          "metadata": {
            "resolution": 1.8,
            "method": "X-RAY",
            "chains": ["A"]
          }
        },
        "sequence": {
          "file_path": "sequence/fasta/1ubq.fasta",
          "metadata": {
            "length": 76,
            "organism": "Homo sapiens"
          }
        }
      },
      "created": "2024-01-15T10:30:00",
      "modified": "2024-01-20T14:22:00"
    }
  }
}
```

### Registry: The Central Ledger

The **EntityRegistry** maintains a complete inventory of all biological objects in the system.

#### Registry Responsibilities

1. **Entity Tracking**
   - Assigns hash IDs (internal use only)
   - Maps human names to hash IDs
   - Tracks aliases and alternative names
   - Records which formats exist for each entity

2. **File Path Mapping**
   - Maps entities to their file locations
   - All paths use human-readable names
   - Handles complex identifier sanitization

3. **Metadata Management**
   - Stores format-specific metadata
   - Tracks creation/modification times
   - Maintains relationships between entities

#### Registry Operations

```python
class EntityRegistry:
    def __init__(self, paths: Optional[ProtosPaths] = None):
        self.paths = paths or ProtosPaths()
        self.registry_file = self.paths.get_registry_path("entity_registry.json")
    
    # User-facing operations (work with human names)
    def register_entity(self, name: str, format_type: str, file_path: str, metadata: dict):
        """Register an entity - users provide human-readable name."""
        
    def find_entity(self, name: str) -> Optional[EntityInfo]:
        """Find entity by human-readable name or alias."""
        
    def list_entities(self) -> List[str]:
        """List all entities - returns human-readable names."""
    
    # Internal operations (use hash IDs)
    def _generate_hash_id(self, name: str) -> str:
        """Generate hash ID for internal use only."""
        
    def _resolve_to_hash(self, name: str) -> Optional[str]:
        """Convert human name to hash ID for internal lookups."""
```

### Datasets: Named Collections

A **dataset** is a named collection of related entities within a specific processor context.

#### Dataset Structure

```json
{
  "datasets": {
    "microbial_opsins": {
      "name": "Microbial Opsin Structures",
      "description": "Light-sensitive proteins from microorganisms",
      "content": ["3f8a9c2d1e", "ba1837f945", "a770526060"],  // Hash IDs for robustness
      "metadata": {
        "created": "2024-01-15",
        "organism": "various",
        "purpose": "structural comparison"
      }
    }
  }
}
```

#### Dataset Operations

```python
class DatasetManager:
    def __init__(self, paths: Optional[ProtosPaths] = None):
        self.paths = paths or ProtosPaths()
    
    # User-facing operations
    def create_dataset(self, name: str, entities: List[str], metadata: dict):
        """Create dataset - users provide human-readable entity names."""
        
    def load_dataset(self, name: str) -> Dataset:
        """Load dataset - returns human-readable information."""
        
    def list_dataset_contents(self, name: str) -> List[str]:
        """List entities in dataset - returns human-readable names."""
```

---

## Processors: Format-Specific Data Handlers

Processors are the **primary interface** for users. They handle specific biological data formats and automatically manage entities and datasets behind the scenes. Users typically ONLY interact with processors - not with ProtosPaths, EntityRegistry, or DatasetManager directly.

### Normal Usage Pattern

```python
# This is ALL a user needs to do:
from protos.processing.structure import CifBaseProcessor

# Just create a processor - everything else is automatic
processor = CifBaseProcessor()
processor.load_structure("1ubq")  # Works immediately!

# Advanced: Only if user wants non-default data location
from protos.io.paths import ProtosPaths
custom_paths = ProtosPaths(base_path="/my/custom/data")
processor = CifBaseProcessor(paths=custom_paths)
```

### BaseProcessor: The Foundation

```python
class BaseProcessor(ABC):
    """Abstract base class for all processors."""
    
    def __init__(self, name: str, paths: Optional[ProtosPaths] = None):
        self.name = name
        
        # Automatic initialization cascade:
        # 1. ProtosPaths (creates default if not provided)
        self.paths = paths or ProtosPaths()
        
        # 2. Processor-specific setup
        self.processor_type = self._get_processor_type()
        self.data_path = self.paths.get_processor_path(self.processor_type)
        
        # 3. Registry and Dataset Manager (created automatically)
        self.entity_registry = EntityRegistry(self.paths)
        self.dataset_manager = DatasetManager(self.paths)
        
        # User never needs to know about these internal components!
```

### How Processors Interact with Entities

Processors provide a clean interface for entity operations while the registry handles tracking behind the scenes.

#### Entity Operations

```python
class BaseProcessor(ABC):
    
    def list_entities(self) -> List[str]:
        """
        List all entities available in this format.
        Returns human-readable names, not hash IDs.
        
        Example:
            >>> processor.list_entities()
            ['1ubq', '2gb1', 'my_protein', 'EGFR_HUMAN']
        """
        # Registry returns hash IDs internally
        hash_ids = self.entity_registry.list_entities_by_format(self.processor_type)
        
        # Convert to human names for users
        return [self.entity_registry.get_human_name(hid) for hid in hash_ids]
    
    def load_entity(self, name: str) -> Any:
        """
        Load entity by human-readable name.
        
        Example:
            >>> structure = processor.load_entity("1ubq")
            >>> sequence = processor.load_entity("EGFR_HUMAN")
        """
        # Registry resolves human name to file path
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if not entity_info:
            # Try direct file access for unregistered files
            return self._try_load_unregistered(name)
        
        # Load from registered path
        return self._load_from_path(entity_info.file_path)
    
    def save_entity(self, name: str, data: Any, metadata: Optional[dict] = None):
        """
        Save entity with human-readable name.
        Automatically registers in entity system.
        
        Example:
            >>> processor.save_entity("my_protein", structure_data)
        """
        # Determine file path using human name
        file_path = self._get_entity_path(name)
        
        # Save the actual data
        self._save_to_path(file_path, data)
        
        # Register in entity system (hash ID created internally)
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=str(file_path),
            metadata=metadata or {}
        )
    
    def delete_entity(self, name: str):
        """
        Delete entity from this format.
        
        Example:
            >>> processor.delete_entity("old_structure")
        """
        # Find entity
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info:
            # Delete file
            Path(entity_info.file_path).unlink(missing_ok=True)
            # Update registry
            self.entity_registry.remove_format(name, self.processor_type)
    
    def entity_exists(self, name: str) -> bool:
        """
        Check if entity exists in this format.
        
        Example:
            >>> if processor.entity_exists("1ubq"):
            ...     structure = processor.load_entity("1ubq")
        """
        return self.entity_registry.has_entity(name, self.processor_type)
```

### How Processors Interact with Datasets

Datasets are collections of entities. Processors make it easy to work with groups of related data.

#### Dataset Operations

```python
class BaseProcessor(ABC):
    
    def list_datasets(self) -> List[str]:
        """
        List all available datasets for this processor.
        
        Example:
            >>> processor.list_datasets()
            ['training_set', 'test_structures', 'kinase_family']
        """
        return self.dataset_manager.list_datasets(self.processor_type)
    
    def create_dataset(self, dataset_name: str, entity_names: List[str], 
                      metadata: Optional[dict] = None):
        """
        Create a dataset from entity names.
        
        Example:
            >>> processor.create_dataset(
            ...     "kinase_study",
            ...     ["EGFR", "ERBB2", "ERBB3", "ABL1"],
            ...     metadata={"organism": "human", "family": "kinase"}
            ... )
        """
        # Convert human names to hash IDs (internal only)
        entity_ids = []
        for name in entity_names:
            entity_info = self.entity_registry.find_entity(name)
            if entity_info:
                entity_ids.append(entity_info.hash_id)
            else:
                print(f"Warning: Entity '{name}' not found")
        
        # Create dataset with hash IDs internally
        self.dataset_manager.create_dataset(
            name=dataset_name,
            entity_ids=entity_ids,
            processor_type=self.processor_type,
            metadata=metadata
        )
    
    def load_dataset(self, dataset_name: str) -> Dict[str, Any]:
        """
        Load all entities in a dataset.
        Returns dict mapping human names to data.
        
        Example:
            >>> structures = processor.load_dataset("kinase_study")
            >>> for name, structure in structures.items():
            ...     print(f"Processing {name}")
        """
        # Get dataset info
        dataset = self.dataset_manager.get_dataset(dataset_name, self.processor_type)
        if not dataset:
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        # Load each entity
        result = {}
        for entity_id in dataset.entity_ids:
            # Get human name for user-friendly output
            human_name = self.entity_registry.get_human_name(entity_id)
            try:
                data = self.load_entity(human_name)
                result[human_name] = data
            except Exception as e:
                print(f"Error loading {human_name}: {e}")
        
        return result
    
    def add_to_dataset(self, dataset_name: str, entity_names: List[str]):
        """
        Add entities to existing dataset.
        
        Example:
            >>> processor.add_to_dataset("kinase_study", ["SRC", "YES1"])
        """
        self.dataset_manager.add_entities_to_dataset(
            dataset_name, entity_names, self.processor_type
        )
    
    def remove_from_dataset(self, dataset_name: str, entity_names: List[str]):
        """
        Remove entities from dataset.
        
        Example:
            >>> processor.remove_from_dataset("kinase_study", ["ABL1"])
        """
        self.dataset_manager.remove_entities_from_dataset(
            dataset_name, entity_names, self.processor_type
        )
    
    def get_dataset_info(self, dataset_name: str) -> dict:
        """
        Get information about a dataset.
        
        Example:
            >>> info = processor.get_dataset_info("kinase_study")
            >>> print(f"Dataset contains {len(info['entities'])} structures")
            >>> print(f"Created on: {info['created']}")
        """
        dataset = self.dataset_manager.get_dataset(dataset_name, self.processor_type)
        if not dataset:
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        # Convert to user-friendly format
        return {
            'name': dataset.name,
            'entities': [self.entity_registry.get_human_name(eid) 
                        for eid in dataset.entity_ids],
            'metadata': dataset.metadata,
            'created': dataset.created,
            'modified': dataset.modified
        }
```

### Specialized Processors

#### CifBaseProcessor (Structure) - Special Caching Behavior

CifBaseProcessor has unique caching requirements due to the computational cost of parsing CIF files. Unlike other processors, it maintains a two-tier caching system:

1. **Individual Structure Cache** (`cache/`): Preprocessed PKL files for individual structures
2. **Dataset Cache** (`structure_dataset/`): Preprocessed PKL files for entire datasets

```python
class CifBaseProcessor(BaseProcessor):
    """
    Handles 3D protein structures with intelligent caching.
    
    Special behavior:
    - Individual structures can be cached as PKL files in cache/
    - Entire datasets can be saved as PKL files in structure_dataset/
    - Loading prioritizes PKL over CIF for performance
    """
    
    def __init__(self, name: str = "structures", paths: Optional[ProtosPaths] = None):
        super().__init__(name, paths)
    
    def load_structure(self, pdb_id: str, use_cache: bool = True) -> pd.DataFrame:
        """
        Load structure with caching support.
        
        Priority:
        1. Check cache/ for preprocessed PKL
        2. Load from mmcif/ and optionally save to cache
        """
        if use_cache:
            cache_path = self.data_path / "cache" / f"{pdb_id}.pkl"
            if cache_path.exists():
                return pd.read_pickle(cache_path)
        
        # Load from CIF
        cif_path = self.data_path / "mmcif" / f"{pdb_id}.cif"
        structure = self._parse_cif(cif_path)
        
        # Save to cache for next time
        if use_cache:
            structure.to_pickle(cache_path)
            
        return structure
    
    def save_data(self, dataset_id: str, data: pd.DataFrame = None, format: str = 'pkl'):
        """
        Save dataset to structure_dataset/.
        
        SPECIAL: Unlike other processors, saves the actual data as PKL,
        not just the dataset metadata JSON.
        """
        if format == 'pkl':
            # Save actual data to structure_dataset/
            dataset_path = self.data_path / "structure_dataset" / f"{dataset_id}.pkl"
            data.to_pickle(dataset_path)
            
        # Also create dataset JSON metadata
        super().create_dataset(dataset_id, list(data['pdb_id'].unique()))
    
    def load_data(self, dataset_id: str, format: str = 'pkl') -> pd.DataFrame:
        """
        Load dataset with fallback to individual structures.
        
        Priority:
        1. Check structure_dataset/ for preprocessed PKL
        2. Fall back to loading individual structures from cache/
        3. Fall back to loading individual structures from mmcif/
        """
        # Try preprocessed dataset
        if format == 'pkl':
            dataset_path = self.data_path / "structure_dataset" / f"{dataset_id}.pkl"
            if dataset_path.exists():
                return pd.read_pickle(dataset_path)
        
        # Fall back to loading individual structures
        dataset_info = self.get_dataset_info(dataset_id)
        structures = []
        for pdb_id in dataset_info['entities']:
            structure = self.load_structure(pdb_id)  # Uses cache
            structures.append(structure)
            
        return pd.concat(structures, ignore_index=True)
```

**Key Differences from Other Processors:**

1. **Dual PKL Storage**:
   - Other processors: Only dataset metadata (JSON) in `datasets/`
   - CifBaseProcessor: Both metadata (JSON) in `datasets/` AND data (PKL) in `structure_dataset/`

2. **Performance Optimization**:
   - Raw CIF parsing is computationally expensive
   - PKL files load 10-100x faster
   - Essential for large structure datasets

3. **Flexible Loading**:
   - Can load individual structures or entire datasets
   - Automatically uses cached version when available
   - Falls back gracefully through cache tiers

4. **Dataset as Single File**:
   - When working with many structures, a single PKL is more efficient
   - Reduces file I/O operations
   - Maintains consistent DataFrame structure

#### SeqProcessor (Sequence)
```python
class SeqProcessor(BaseProcessor):
    """Handles protein sequences."""
    
    def load_sequence(self, name: str) -> str:
        """Load sequence by name."""
        # Handle complex identifiers
        safe_name = self._sanitize_filename(name)
        filepath = self.data_path / "fasta" / f"{safe_name}.fasta"
```

#### GRNBaseProcessor (Generic Residue Numbering)

GRNBaseProcessor handles Generic Residue Numbering tables which provide standardized position numbering across protein families. It has a unique directory structure to support both user-generated tables and reference tables used for GRN assignment.

```python
class GRNBaseProcessor(BaseProcessor):
    """
    Handles GRN alignment tables and reference numbering systems.
    
    Directory structure:
    grn/
    ├── tables/          # User GRN tables (entities as rows)
    ├── ref/             # Reference GRN tables for assignment
    ├── configs/         # Configuration files for GRN assignment
    ├── datasets/        # Dataset metadata JSON files
    └── assignments/     # GRN assignment results
    """
    
    def load_entity(self, name: str) -> pd.Series:
        """
        Load a single GRN entity (row from a table).
        
        Returns the GRN positions for a specific protein.
        """
        # Search through tables for this entity
        # Returns Series with GRN positions as index
        
    def save_grn_table(self, table_name: str, grn_df: pd.DataFrame):
        """
        Save GRN table with human-readable protein IDs.
        
        Table format (saved to tables/):
        protein_id (index) | 1.50 | 2.50 | 3.50 | ...
        -------------------+------+------+------+----
        BACR_HALSA        |  A   |  L   |  V   | ...
        ChR2              |  A   |  L   |  I   | ...
        
        Dataset metadata (saved to datasets/):
        {
          "name": "opsin_grn_alignment",
          "entities": ["BACR_HALSA", "ChR2", ...],
          "metadata": {
            "grn_system": "Ballesteros-Weinstein",
            "protein_family": "microbial_rhodopsin"
          }
        }
        """
        # Save table to tables/
        table_path = self.data_path / "tables" / f"{table_name}.csv"
        grn_df.to_csv(table_path, index=True)  # Keep index for protein IDs
        
        # Create dataset metadata in datasets/
        dataset_info = {
            "name": table_name,
            "entities": list(grn_df.index),
            "table_file": f"../tables/{table_name}.csv",
            "metadata": {
                "grn_positions": list(grn_df.columns),
                "entity_count": len(grn_df),
                "created": datetime.now().isoformat()
            }
        }
        
        dataset_path = self.data_path / "datasets" / f"{table_name}.json"
        with open(dataset_path, 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        # Register each protein in the table
        for protein_id in grn_df.index:
            self.entity_registry.register_entity(
                name=protein_id,
                format_type="grn",
                file_path=str(table_path.relative_to(self.paths.data_root)),
                metadata={"table": table_name}
            )
    
    def load_reference_table(self, ref_name: str) -> pd.DataFrame:
        """
        Load reference GRN table for assignment.
        
        Reference tables are in ref/ and used for:
        - Template-based GRN assignment
        - Family-specific numbering systems
        - Standard reference alignments
        """
        ref_path = self.data_path / "ref" / f"{ref_name}.csv"
        return pd.read_csv(ref_path, index_col=0)
    
    def get_grn_config(self, config_name: str) -> dict:
        """
        Load GRN assignment configuration.
        
        Configs in configs/ define:
        - Motif patterns for position detection
        - Binding domain definitions
        - Assignment parameters
        """
        config_path = self.data_path / "configs" / f"{config_name}.json"
        with open(config_path, 'r') as f:
            return json.load(f)
```

**Key Features:**

1. **Table-Based Storage**: GRN data is stored in CSV tables where:
   - Each row is an entity (protein)
   - Each column is a GRN position (e.g., 1.50, 2.50)
   - Index contains human-readable protein names

2. **Reference Data Support**: 
   - `ref/`: Contains curated reference tables for GRN assignment
   - `configs/`: Contains JSON configuration for assignment algorithms
   - These are typically read-only and part of the Protos distribution

3. **Entity as Table Row**: Unlike file-based processors, GRN entities are rows in tables
   - `load_entity()` returns a Series (one row)
   - `save_entity()` adds/updates a row in a table

4. **Dataset Organization**:
   - `datasets/`: Contains JSON metadata about GRN tables
   - `tables/`: Contains actual CSV data files
   - Dataset JSON references the corresponding CSV in tables/

#### PropertyProcessor

PropertyProcessor manages experimental properties and metadata for entities across all format types. Like GRN, it uses a table-based approach but with even more flexibility in property assignment.

```python
class PropertyProcessor(BaseProcessor):
    """
    Handles experimental properties and metadata.
    
    Directory structure:
    property/
    ├── tables/          # Property data tables (CSV)
    ├── datasets/        # Dataset metadata (JSON only)
    ├── cache/           # Entity property cache
    ├── metadata/        # Property definitions
    └── assignments/     # Property assignment logs
    """
    
    def assign_property(self, entity_name: str, property_name: str, 
                       value: Any, dataset: str):
        """
        Assign a property to any entity.
        
        Works with entities from any processor:
        - Structure: "1UBQ", "2GB1"
        - Sequence: "BACR_HALSA", "sp|P02724|GLPA_ECOLI"
        - GRN: "ChR2_position_3.50"
        """
        
    def save_property_dataset(self, dataset_name: str, format: str = 'both'):
        """
        Save property dataset following strict separation:
        
        1. Dataset metadata → datasets/dataset_name.json
        2. Property data → tables/dataset_name.csv
        
        Example dataset JSON (datasets/opsin_properties.json):
        {
          "name": "opsin_properties",
          "entities": ["BACR_HALSA", "ChR2", "NpHR"],
          "data_file": "../tables/opsin_properties.csv",
          "metadata": {
            "properties": ["lambda_max", "photocycle", "pump_type"],
            "source": "literature_review",
            "created": "2024-01-15"
          }
        }
        
        Example data CSV (tables/opsin_properties.csv):
        entity_name  | lambda_max | photocycle | pump_type
        -------------+------------+------------+-----------
        BACR_HALSA   |    568     |   fast     |  proton
        ChR2         |    470     |   slow     |  channel
        """
```

**Key Features:**

1. **Strict File Separation**:
   - `datasets/`: Contains ONLY JSON metadata files
   - `tables/`: Contains ONLY CSV data files
   - JSON files reference their corresponding CSV

2. **Cross-Format Support**: Can assign properties to entities from any processor
   - Uses entity registry to resolve names
   - Maintains format-agnostic property storage

3. **Flexible Properties**: No fixed schema
   - Properties can be added dynamically
   - Supports any data type that can be serialized

4. **Human-Readable Throughout**:
   - Entity names in tables are human-readable
   - No hash IDs in any user-facing files

---

## File Naming and Organization

### Core Principle: Human-Readable Names Everywhere

All files in the Protos filesystem use human-readable names. Hash IDs are NEVER used in filenames or paths.

### File Naming Rules

1. **Use Original Biological Identifiers**
   ```
   structure/mmcif/1ubq.cif         ✓ CORRECT
   structure/mmcif/3f8a9c2d1e.cif   ✗ WRONG
   ```

2. **Sanitize Complex Identifiers**
   ```python
   def sanitize_filename(identifier: str) -> str:
       """Convert complex identifiers to safe filenames."""
       # sp|P02724|GLPA_ECOLI → sp_P02724_GLPA_ECOLI
       safe = identifier.replace('|', '_')
       safe = safe.replace('/', '_')
       safe = safe.replace(':', '_')
       # ... other unsafe characters
       return safe
   ```

3. **Consistent Naming Across Formats**
   ```
   structure/mmcif/1ubq.cif
   sequence/fasta/1ubq.fasta
   sequence/fasta/1ubq_A.fasta    # Chain A sequence
   ```

4. **Descriptive Names for Multi-Entity Files**
   ```
   grn/tables/microbial_opsins_grn.csv
   property/tables/kinase_inhibitor_properties.csv
   ```

### Directory Organization

```
data/
├── structure/
│   └── mmcif/
│       ├── 1ubq.cif              # PDB ID
│       ├── 3sn6.cif              # PDB ID
│       └── my_protein_v2.cif     # Custom name
├── sequence/
│   └── fasta/
│       ├── P62988.fasta          # UniProt ID
│       ├── EGFR_HUMAN.fasta      # Gene name
│       └── sp_P02724_GLPA_ECOLI.fasta  # Sanitized complex ID
├── grn/
│   └── tables/
│       ├── microbial_opsins_grn.csv     # Descriptive name
│       └── gpcr_family_grn.csv          # Descriptive name
└── property/
    └── tables/
        ├── opsin_properties.csv         # Descriptive name
        └── experimental_results_2024.csv # Descriptive with date
```

---

## Implementation Guidelines

### For New Components

```python
class NewProcessor(BaseProcessor):
    def __init__(self, name: str = "new_processor", paths: Optional[ProtosPaths] = None):
        # ALWAYS accept ProtosPaths
        super().__init__(name, paths)
    
    def save_data(self, entity_name: str, data: Any):
        # ALWAYS use human-readable names
        filepath = self.data_path / "subdir" / f"{entity_name}.ext"
        
        # Save file
        self._write_file(filepath, data)
        
        # ALWAYS register entity
        self.entity_registry.register_entity(
            name=entity_name,
            format_type=self.processor_type,
            file_path=str(filepath),
            metadata={...}
        )
```

### For File Operations

```python
# ✅ CORRECT - Use ProtosPaths
filepath = self.paths.get_processor_path("structure") / "mmcif" / f"{pdb_id}.cif"

# ❌ WRONG - Manual path construction
filepath = Path("data/structure/mmcif") / f"{pdb_id}.cif"

# ❌ WRONG - Hardcoded paths
filepath = f"/home/user/protos/data/structure/mmcif/{pdb_id}.cif"

# ❌ WRONG - Using hash IDs
filepath = self.data_path / f"{entity_hash_id}.cif"
```

### For User Interfaces

```python
# Users ALWAYS work with human names
processor.load_structure("1ubq")
processor.save_sequence("EGFR_HUMAN", sequence)
processor.create_dataset("kinase_study", ["EGFR", "ERBB2", "ERBB3"])

# Users NEVER see hash IDs
# ❌ WRONG
processor.load_entity("3f8a9c2d1e")
```

---

## Common Workflows

### 1. Drag-and-Drop Workflow

User drags PDB files into the mmcif folder:

```python
# Files appear as:
# data/structure/mmcif/6xyz.cif
# data/structure/mmcif/7abc.cif

# User code:
processor = CifBaseProcessor()

# Option 1: Direct load (works immediately)
structure = processor.load_structure("6xyz")

# Option 2: Create dataset
processor.create_dataset("my_structures", ["6xyz", "7abc"])
structures = processor.load_dataset("my_structures")
```

### 2. Cross-Format Workflow

```python
# Load structure
struct_proc = CifBaseProcessor()
structure = struct_proc.load_structure("1ubq")

# Extract and save sequence
sequence = struct_proc.extract_sequence(structure)
seq_proc = SeqProcessor()
seq_proc.save_sequence("1ubq", sequence)

# Both formats tracked under same entity
# Registry shows:
# "1ubq" -> {
#   "formats": ["structure", "sequence"],
#   "files": {
#     "structure": "structure/mmcif/1ubq.cif",
#     "sequence": "sequence/fasta/1ubq.fasta"
#   }
# }
```

### 3. Multi-Entity Table Workflow

```python
# Create GRN table
grn_data = {
    'protein_id': ['BACR_HALSA', 'ChR2', 'NpHR'],  # Human names
    '1.50': ['A', 'A', 'A'],
    '2.50': ['L', 'L', 'M'],
    # ...
}
grn_df = pd.DataFrame(grn_data)

grn_proc = GRNBaseProcessor()
grn_proc.save_grn_table("opsins_alignment", grn_df)

# Each protein automatically registered as entity
# File saved as: data/grn/tables/opsins_alignment_grn.csv
```

### 4. Property Import Workflow

```python
# Import experimental data
properties = pd.DataFrame({
    'entity_name': ['BACR_HALSA', 'ChR2', 'NpHR'],  # Human names
    'lambda_max': [568, 470, 590],
    'photocycle': ['fast', 'slow', 'fast']
})

prop_proc = PropertyProcessor()
prop_proc.import_properties("opsin_properties", properties)

# Query properties by entity
bacr_props = prop_proc.get_entity_properties("BACR_HALSA")
```

### 5. Complete Example: Entity and Dataset Operations

```python
# Normal user workflow - no setup needed!
from protos.processing.structure import CifBaseProcessor

# 1. Initialize processor (automatic setup)
processor = CifBaseProcessor()

# 2. List what's available
print("Available structures:", processor.list_entities())
print("Available datasets:", processor.list_datasets())

# 3. Load individual entities
structure = processor.load_entity("1ubq")
if processor.entity_exists("2gb1"):
    another = processor.load_entity("2gb1")

# 4. Save new entity (auto-registers)
processor.save_entity("my_complex", my_structure_data, 
                     metadata={"method": "cryo-EM", "resolution": 3.2})

# 5. Create a dataset
processor.create_dataset(
    "ubiquitin_study",
    ["1ubq", "2gb1", "my_complex"],
    metadata={"purpose": "comparative analysis"}
)

# 6. Work with datasets
structures = processor.load_dataset("ubiquitin_study")
for name, structure in structures.items():
    print(f"Processing {name}: {len(structure)} atoms")

# 7. Update dataset
processor.add_to_dataset("ubiquitin_study", ["1d3z", "1f9j"])

# 8. Get dataset info
info = processor.get_dataset_info("ubiquitin_study")
print(f"Dataset now contains: {info['entities']}")

# Everything above works with ZERO configuration!
# ProtosPaths, EntityRegistry, and DatasetManager all created automatically
```

### Key Points About Processor-Entity-Dataset Interaction

1. **Processors are the interface** - Users never need to import or create ProtosPaths, EntityRegistry, or DatasetManager

2. **Automatic cascade** - Creating a processor automatically creates:
   - ProtosPaths (with default data location)
   - EntityRegistry (for tracking entities)
   - DatasetManager (for managing datasets)

3. **Human names everywhere** - All processor methods accept and return human-readable names

4. **Hash IDs are invisible** - Used internally by registry but never exposed to users

5. **Registration is automatic** - Saving an entity automatically registers it

6. **Unregistered files work** - Can load files that aren't registered yet

7. **Cross-format aware** - Same entity can exist in multiple formats, tracked by registry

---

## Summary

The Protos data management system is built on three pillars:

1. **ProtosPaths**: The single, unified path management system
   - Zero configuration required
   - Works out of the box
   - Every component uses it

2. **Entity/Registry/Dataset**: The tracking and organization system
   - Entities use human-readable names everywhere
   - Hash IDs are internal registry keys only
   - Datasets are individual JSON files in `datasets/` subdirectories
   - Registry files track entity locations and metadata

3. **Processors**: Format-specific data handlers
   - All inherit from BaseProcessor
   - All use ProtosPaths
   - All register entities automatically
   - Each implements dataset loading according to its data format:
     - **Structure/Sequence/Embedding**: Return individual files as dict
     - **GRN/Property**: Return tables where each row is an entity

### Key Dataset Concepts

1. **Datasets are individual JSON files** in dedicated directories:
   - `structure/datasets/kinases.json`
   - `sequence/datasets/human_proteins.json`
   - `grn/datasets/opsin_alignment.json`
   - `property/datasets/experimental_data.json`
   - etc.

2. **Format-specific behavior**:
   - Single-entity formats (structure, sequence): Dataset = collection of file references
   - Multi-entity formats (GRN, property): Dataset = reference to table file

3. **Consistent interface**:
   - All processors use same dataset methods
   - Implementation varies by data format
   - Users get format-appropriate return values

This design ensures:
- **Simplicity**: Users work with familiar names
- **Consistency**: One path system, one registry
- **Discoverability**: Browse files and understand content
- **Robustness**: Internal tracking handles edge cases
- **Zero Setup**: Works immediately after installation