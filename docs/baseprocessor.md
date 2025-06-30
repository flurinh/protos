# BaseProcessor Architecture

BaseProcessor is the abstract base class that provides a common interface and shared functionality for all Protos data processors. It implements core patterns for data management, entity tracking, and cross-format integration.

## Overview

The BaseProcessor class establishes:
- Common interface for all processors
- Entity and dataset management
- Path resolution through ProtosPaths
- Registry integration
- Standardized data operations

## Class Hierarchy

```
BaseProcessor (Abstract)
    ├── CifBaseProcessor      # 3D protein structures
    ├── GRNBaseProcessor      # Generic Residue Numbering
    ├── SeqProcessor          # Sequence data
    ├── PropertyProcessor     # Experimental properties
    └── EmbeddingProcessor    # ML embeddings
```

## Core Interface

All processors implement these standard methods:

```python
class BaseProcessor(ABC):
    # Entity operations
    def list_entities(self) -> List[str]
    def load_entity(self, name: str) -> Any
    def save_entity(self, name: str, data: Any) -> str
    def delete_entity(self, name: str) -> bool
    
    # Dataset operations  
    def list_datasets(self) -> List[Dict]
    def create_dataset(self, dataset_id: str, **kwargs) -> bool
    def load_dataset(self, dataset_id: str) -> bool
    def delete_dataset(self, dataset_id: str) -> bool
    
    # Registry operations
    def get_entity_registry(self) -> EntityRegistry
    def resolve_identifier(self, identifier: str) -> str
```

## Initialization

BaseProcessor handles common initialization for all processors:

```python
from protos.core.base_processor import BaseProcessor

class CustomProcessor(BaseProcessor):
    def __init__(self, name: str = "custom", **kwargs):
        # Initialize base class
        super().__init__(
            name=name,
            processor_type="custom",  # Determines directory
            **kwargs
        )
        
        # Automatic setup:
        # - self.data_path set from ProtosPaths
        # - self.registry initialized
        # - self.logger configured
        # - Directories created if needed
```

## Key Features

### 1. Automatic Path Management

BaseProcessor integrates with ProtosPaths automatically:

```python
processor = AnyProcessor(name="analysis")

# These paths are automatically available:
processor.data_path          # Root data directory for processor
processor.path_registry      # Registry file path
processor.path_datasets      # Dataset directory

# Processor-specific paths added by subclasses:
struct_processor.path_mmcif  # Structure-specific
seq_processor.path_fasta     # Sequence-specific
```

### 2. Entity Management

Standardized entity operations across all processors:

```python
# List available entities (returns names, not hash IDs)
entities = processor.list_entities()
# ["1ubq", "2gb1", "P12345", "MY_PROTEIN"]

# Load entity by name (automatic ID resolution)
data = processor.load_entity("1ubq")

# Save new entity (automatic registration)
processor.save_entity("NEW_PROTEIN", data)
# - Generates deterministic entity ID
# - Saves data in appropriate format
# - Registers in entity registry
```

### 3. Dataset Management

Datasets provide logical grouping of entities:

```python
# Create dataset
processor.create_dataset(
    dataset_id="kinase_study",
    name="Human kinase structures",
    description="Structures for kinase inhibitor study",
    content=["1atp", "2src", "3erk", "4bcr"],  # Entity names
    metadata={
        "organism": "Homo sapiens",
        "family": "protein kinase",
        "date": "2024-01-15"
    }
)

# Load dataset (loads all entities)
processor.load_dataset("kinase_study")

# Iterate over dataset
for entity_name in processor.iter_dataset("kinase_study"):
    data = processor.load_entity(entity_name)
    # Process each entity
```

### 4. Registry Integration

Automatic registry management for all operations:

```python
# Registry tracks all entities
registry = processor.get_entity_registry()

# Resolve any identifier to entity ID
entity_id = processor.resolve_identifier("1ubq")  # PDB ID
entity_id = processor.resolve_identifier("P12345")  # UniProt ID
entity_id = processor.resolve_identifier("3f8a9c2d1e")  # Already an ID

# Query registry
entities = registry.get_entities_by_format("structure")
metadata = registry.get_entity_metadata(entity_id)
```

### 5. Logging System

Built-in logging for debugging and monitoring:

```python
class CustomProcessor(BaseProcessor):
    def process_data(self, data):
        self.logger.info(f"Processing {len(data)} items")
        
        try:
            result = self._analyze(data)
            self.logger.debug(f"Analysis complete: {result}")
        except Exception as e:
            self.logger.error(f"Processing failed: {e}")
            raise
```

## Abstract Methods

Subclasses must implement format-specific operations:

```python
class CustomProcessor(BaseProcessor):
    @abstractmethod
    def _load_entity_data(self, entity_id: str, file_path: Path) -> Any:
        """Load data from file."""
        pass
    
    @abstractmethod
    def _save_entity_data(self, data: Any, file_path: Path) -> bool:
        """Save data to file."""
        pass
    
    @abstractmethod
    def _get_entity_metadata(self, data: Any) -> Dict:
        """Extract metadata from data."""
        pass
```

## Common Patterns

### 1. Processor Initialization

```python
# Standard initialization pattern
processor = SpecificProcessor(
    name="my_analysis",      # Unique processor instance name
    preload=False,           # Don't load data on init
    verbose=True             # Enable detailed logging
)

# Data path automatically set from ProtosPaths
print(processor.data_path)  # /path/to/protos_data/processor_type
```

### 2. Batch Operations

```python
# Process multiple entities efficiently
entity_names = ["1ubq", "2gb1", "1crn", "3nir"]

results = {}
for name in entity_names:
    try:
        data = processor.load_entity(name)
        result = processor.analyze(data)
        results[name] = result
    except Exception as e:
        processor.logger.warning(f"Failed to process {name}: {e}")

# Save results as new dataset
processor.save_analysis_results(results, "batch_analysis")
```

### 3. Cross-Format Integration

```python
# BaseProcessor enables cross-format workflows
struct_proc = CifBaseProcessor()
seq_proc = SeqProcessor()

# Load structure
structure = struct_proc.load_entity("1ubq")

# Extract sequence (format-specific method)
sequence = struct_proc.extract_sequence(structure)

# Save in sequence processor (same entity name)
seq_proc.save_entity("1ubq", sequence)

# Both processors now track the same entity
```

### 4. Metadata Management

```python
# Rich metadata support
processor.save_entity(
    name="MY_EXPERIMENT",
    data=experimental_data,
    metadata={
        "date": "2024-01-15",
        "method": "X-ray crystallography",
        "resolution": 1.8,
        "conditions": {
            "pH": 7.4,
            "temperature": 298
        }
    }
)

# Retrieve metadata
metadata = processor.get_entity_metadata("MY_EXPERIMENT")
```

## Extending BaseProcessor

### Creating a Custom Processor

```python
from protos.core.base_processor import BaseProcessor
from pathlib import Path
from typing import Any, Dict, List
import json

class JsonProcessor(BaseProcessor):
    """Example processor for JSON data."""
    
    def __init__(self, name: str = "json_data", **kwargs):
        super().__init__(name, processor_type="json", **kwargs)
        self.path_json = self.data_path / "files"
        self.path_json.mkdir(exist_ok=True)
    
    def _load_entity_data(self, entity_id: str, file_path: Path) -> Any:
        """Load JSON data from file."""
        with open(file_path, 'r') as f:
            return json.load(f)
    
    def _save_entity_data(self, data: Any, file_path: Path) -> bool:
        """Save data as JSON."""
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
        return True
    
    def _get_entity_metadata(self, data: Any) -> Dict:
        """Extract metadata from JSON data."""
        return {
            "keys": list(data.keys()) if isinstance(data, dict) else None,
            "type": type(data).__name__,
            "size": len(str(data))
        }
    
    def _get_entity_file_path(self, entity_id: str) -> Path:
        """Get file path for entity."""
        return self.path_json / f"{entity_id}.json"
```

### Adding Processor-Specific Methods

```python
class JsonProcessor(BaseProcessor):
    # ... base implementation ...
    
    def query_json(self, entity_name: str, json_path: str) -> Any:
        """Query JSON data using path notation."""
        data = self.load_entity(entity_name)
        
        # Navigate JSON path (e.g., "results.proteins[0].name")
        current = data
        for part in json_path.split('.'):
            if '[' in part and ']' in part:
                key, index = part.split('[')
                index = int(index.rstrip(']'))
                current = current[key][index]
            else:
                current = current[part]
        
        return current
    
    def merge_entities(self, entity_names: List[str], 
                      output_name: str) -> str:
        """Merge multiple JSON entities."""
        merged_data = {}
        
        for name in entity_names:
            data = self.load_entity(name)
            merged_data[name] = data
        
        # Save merged result
        return self.save_entity(output_name, merged_data)
```

## Best Practices

### 1. Use Base Class Methods

Leverage inherited functionality rather than reimplementing:

```python
# ✅ GOOD - Use base class method
entities = self.list_entities()

# ❌ BAD - Reimplementing listing
import os
entities = [f for f in os.listdir(self.path_json) if f.endswith('.json')]
```

### 2. Consistent Error Handling

Follow the base class error patterns:

```python
def load_entity(self, name: str) -> Any:
    try:
        entity_id = self.resolve_identifier(name)
        file_path = self._get_entity_file_path(entity_id)
        
        if not file_path.exists():
            self.logger.warning(f"Entity {name} not found")
            return None
            
        return self._load_entity_data(entity_id, file_path)
        
    except Exception as e:
        self.logger.error(f"Failed to load {name}: {e}")
        raise
```

### 3. Metadata Standards

Include rich metadata for tracking:

```python
def _get_entity_metadata(self, data: Any) -> Dict:
    return {
        # Standard fields
        "format": self.processor_type,
        "processor": self.name,
        "timestamp": datetime.now().isoformat(),
        
        # Format-specific fields
        "specific_field": data.get('important_value'),
        "statistics": self._calculate_stats(data)
    }
```

### 4. Type Hints and Documentation

Use type hints and docstrings consistently:

```python
def process_entity(
    self,
    entity_name: str,
    parameters: Dict[str, Any],
    output_name: Optional[str] = None
) -> str:
    """
    Process an entity with given parameters.
    
    Args:
        entity_name: Name of entity to process
        parameters: Processing parameters
        output_name: Optional name for output entity
        
    Returns:
        Name of created output entity
        
    Raises:
        ValueError: If parameters are invalid
        FileNotFoundError: If entity doesn't exist
    """
    # Implementation...
```

## Advanced Features

### 1. Lazy Loading

BaseProcessor supports lazy loading patterns:

```python
class LargeDataProcessor(BaseProcessor):
    def __init__(self, name: str = "large_data", preload: bool = False):
        super().__init__(name, processor_type="large", preload=preload)
        self._data_cache = {}
    
    def load_entity(self, name: str, lazy: bool = True):
        if lazy and name in self._data_cache:
            return self._data_cache[name]
        
        data = super().load_entity(name)
        
        if lazy:
            self._data_cache[name] = data
        
        return data
```

### 2. Validation Hooks

Add validation at multiple points:

```python
def save_entity(self, name: str, data: Any, **kwargs) -> str:
    # Pre-save validation
    if not self._validate_data(data):
        raise ValueError(f"Invalid data for {name}")
    
    # Call parent save
    entity_id = super().save_entity(name, data, **kwargs)
    
    # Post-save verification
    self._verify_saved_entity(entity_id)
    
    return entity_id
```

### 3. Event Hooks

Implement event notifications:

```python
class EventProcessor(BaseProcessor):
    def __init__(self, name: str = "events"):
        super().__init__(name, processor_type="event")
        self._event_handlers = {}
    
    def on_entity_saved(self, handler):
        """Register save event handler."""
        self._event_handlers['save'] = handler
    
    def save_entity(self, name: str, data: Any) -> str:
        entity_id = super().save_entity(name, data)
        
        # Trigger event
        if 'save' in self._event_handlers:
            self._event_handlers['save'](name, entity_id)
        
        return entity_id
```

## Summary

BaseProcessor provides:
- Unified interface across all data types
- Automatic path and registry management
- Entity and dataset abstractions
- Cross-format integration support
- Extensible architecture for custom processors

By inheriting from BaseProcessor, new processors automatically gain access to the full Protos infrastructure while maintaining format-specific flexibility.