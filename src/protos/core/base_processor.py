"""
Base processor class for the protos framework.

This implementation follows DATA_MANAGEMENT_UNIFIED.md principles:
- Zero configuration required
- Human-readable names for all operations
- ProtosPaths handles ALL path management
- Hash IDs are internal only
"""

import os
import logging
import pickle
import json
import pandas as pd
import numpy as np
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple, Set
from abc import ABC, abstractmethod

# Import path management utilities
from protos.io.paths import ProtosPaths

# Import our new implementations
from protos.io.entity_registry import EntityRegistry
from protos.io.dataset_manager import DatasetManager


class BaseProcessor(ABC):
    """
    Base class for all processors in the protos framework.
    
    Key principles from DATA_MANAGEMENT_UNIFIED.md:
    - Zero configuration - works out of the box
    - Human-readable names for all user-facing operations
    - ProtosPaths handles ALL path operations
    - Automatic entity registration
    - Hash IDs are internal only - never exposed to users
    """
    
    def __init__(self, 
                 name: str = "processor",
                 paths: Optional[ProtosPaths] = None):
        """
        Initialize processor with automatic setup.
        
        Args:
            name: Processor instance name (default: "processor")
            paths: Optional ProtosPaths instance (created if not provided)
            
        Example:
            >>> processor = StructureProcessor()  # Zero configuration!
        """
        self.name = name
        self.paths = paths or ProtosPaths()
        self.processor_type = self._get_processor_type()
        self.data_path = Path(self.paths.get_processor_path(self.processor_type))
        
        # Create registries with ProtosPaths
        self.entity_registry = EntityRegistry(paths=self.paths)
        self.dataset_manager = DatasetManager(
            processor_type=self.processor_type,
            paths=self.paths,
            entity_registry=self.entity_registry
        )
        
        # Initialize data storage
        self.data = None
        self.metadata = {
            "processor_type": self.__class__.__name__,
            "name": name,
            "created_at": datetime.now().isoformat(),
            "data_path": str(self.data_path),
        }
        
        # Set up logging
        self.logger = self._setup_logger()
        self.logger.info(f"Initialized {self.__class__.__name__} '{name}' at {self.data_path}")
    
    def _setup_logger(self) -> logging.Logger:
        """Set up a logger for this processor instance."""
        logger = logging.getLogger(f"{self.__class__.__name__}.{self.name}")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger
    
    def _get_processor_type(self) -> str:
        """
        Get processor type from class name.
        
        Returns:
            Processor type string (e.g., 'structure', 'grn')
        """
        class_name = self.__class__.__name__
        
        # Handle common patterns
        if 'Structure' in class_name or 'Cif' in class_name:
            return 'structure'
        elif 'GRN' in class_name:
            return 'grn'
        elif 'Sequence' in class_name or 'Seq' in class_name:
            return 'sequence'
        elif 'Property' in class_name:
            return 'property'
        elif 'Embedding' in class_name or 'EMB' in class_name:
            return 'embedding'
        elif 'Ligand' in class_name:
            return 'ligand'
        elif 'Graph' in class_name:
            return 'graph'
        else:
            # Default to lowercase class name without 'Processor'
            return class_name.replace('Processor', '').replace('Base', '').lower()
    
    # ========== Entity Operations (DATA_MANAGEMENT_UNIFIED.md) ==========
    
    def list_entities(self) -> List[str]:
        """
        List all entities of this processor type.
        Returns human-readable names only.
        
        Example:
            >>> processor.list_entities()
            ['1ubq', '3sn6', 'EGFR_HUMAN']
        """
        return self.entity_registry.list_entities(self.processor_type)
    
    @abstractmethod
    def load_entity(self, name: str) -> Any:
        """
        Load entity by human-readable name.
        Subclasses must implement format-specific loading.
        
        Example:
            >>> structure = processor.load_entity("1ubq")
            >>> sequence = processor.load_entity("EGFR_HUMAN")
        """
        pass
    
    @abstractmethod  
    def save_entity(self, name: str, data: Any, metadata: Optional[dict] = None):
        """
        Save entity with human-readable name.
        Subclasses must implement format-specific saving.
        
        Example:
            >>> processor.save_entity("my_protein", structure_data)
            >>> processor.save_entity("EGFR_HUMAN", sequence)
        """
        pass
    
    def entity_exists(self, name: str) -> bool:
        """
        Check if entity exists in this format.
        
        Example:
            >>> if processor.entity_exists("1ubq"):
            ...     structure = processor.load_entity("1ubq")
        """
        return self.entity_registry.entity_exists(name, self.processor_type)
    
    def delete_entity(self, name: str) -> bool:
        """
        Delete entity from this format.
        
        Example:
            >>> processor.delete_entity("old_structure")
            
        Returns:
            True if entity was deleted, False if not found
        """
        # Find entity
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info:
            # Delete file
            file_path = Path(entity_info.file_path)
            if not file_path.is_absolute():
                file_path = Path(self.paths.data_root) / file_path
            file_path.unlink(missing_ok=True)
            # Update registry
            self.entity_registry.remove_format(name, self.processor_type)
            return True
        return False
    
    # ========== Dataset Operations (DATA_MANAGEMENT_UNIFIED.md) ==========
    
    def list_datasets(self) -> List[str]:
        """
        List all available datasets for this processor.
        
        Example:
            >>> processor.list_datasets()
            ['training_set', 'test_structures', 'kinase_family']
        """
        return self.dataset_manager.list_datasets()
    
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
        return self.dataset_manager.create_dataset(dataset_name, entity_names, metadata)
    
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
        dataset = self.dataset_manager.load_dataset(dataset_name)
        
        # Load each entity
        result = {}
        for entity_name in dataset['entities']:
            try:
                data = self.load_entity(entity_name)
                if data is not None:
                    result[entity_name] = data
            except Exception as e:
                print(f"Error loading {entity_name}: {e}")
        
        return result
    
    def add_to_dataset(self, dataset_name: str, entity_names: List[str]):
        """
        Add entities to existing dataset.
        
        Example:
            >>> processor.add_to_dataset("kinase_study", ["SRC", "YES1"])
        """
        self.dataset_manager.add_to_dataset(dataset_name, entity_names)
    
    def remove_from_dataset(self, dataset_name: str, entity_names: List[str]):
        """
        Remove entities from dataset.
        
        Example:
            >>> processor.remove_from_dataset("kinase_study", ["ABL1"])
        """
        self.dataset_manager.remove_from_dataset(dataset_name, entity_names)
    
    def get_dataset_info(self, dataset_name: str) -> dict:
        """
        Get information about a dataset.
        
        Example:
            >>> info = processor.get_dataset_info("kinase_study")
            >>> print(f"Dataset contains {len(info['entities'])} structures")
            >>> print(f"Created on: {info['created']}")
        """
        return self.dataset_manager.get_dataset_info(dataset_name)
    
    # ========== Helper Methods ==========
    
    def get_subdirectory_path(self, subdir_type: str) -> Path:
        """
        Get path to a processor-specific subdirectory.
        
        Args:
            subdir_type: Subdirectory type (e.g., 'cache_dir', 'temp_dir')
            
        Returns:
            Path to subdirectory
        """
        # Delegate to processor-specific methods in ProtosPaths
        method_name = f"get_{self.processor_type}_subdir_path"
        if hasattr(self.paths, method_name):
            method = getattr(self.paths, method_name)
            return Path(method(subdir_type))
        else:
            # Fallback to processor directory
            return self.data_path / subdir_type
    
    def _sanitize_filename(self, name: str) -> str:
        """
        Sanitize a filename to be filesystem-safe while keeping it readable.
        
        Args:
            name: Original name
            
        Returns:
            Sanitized filename
        """
        # Replace problematic characters
        replacements = {
            '/': '_',
            '\\': '_',
            ':': '_',
            '*': '_',
            '?': '_',
            '"': '_',
            '<': '_',
            '>': '_',
            '|': '_',
            ' ': '_'
        }
        
        safe_name = name
        for char, replacement in replacements.items():
            safe_name = safe_name.replace(char, replacement)
        
        return safe_name
    
    # ========== Data Operations (to be overridden by subclasses) ==========
    
    def save_data(self, name: str, data: Any, file_format: str = 'pkl'):
        """
        Save data in specified format. 
        Subclasses should override for format-specific saving.
        
        Args:
            name: Dataset name
            data: Data to save
            file_format: Format to save in
        """
        file_path = self.data_path / f"{name}.{file_format}"
        
        if file_format == 'pkl':
            with open(file_path, 'wb') as f:
                pickle.dump(data, f)
        elif file_format == 'json':
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
        elif file_format == 'csv' and isinstance(data, pd.DataFrame):
            data.to_csv(file_path)
        elif file_format == 'npy' and isinstance(data, np.ndarray):
            np.save(file_path, data)
        else:
            raise ValueError(f"Unsupported format: {file_format}")
    
    def load_data(self, name: str, file_format: str = 'pkl') -> Any:
        """
        Load data in specified format.
        Subclasses should override for format-specific loading.
        
        Args:
            name: Dataset name
            file_format: Format to load from
            
        Returns:
            Loaded data
        """
        file_path = self.data_path / f"{name}.{file_format}"
        
        if not file_path.exists():
            return None
        
        if file_format == 'pkl':
            with open(file_path, 'rb') as f:
                return pickle.load(f)
        elif file_format == 'json':
            with open(file_path, 'r') as f:
                return json.load(f)
        elif file_format == 'csv':
            return pd.read_csv(file_path)
        elif file_format == 'npy':
            return np.load(file_path)
        else:
            raise ValueError(f"Unsupported format: {file_format}")
    
    def is_dataset_available(self, name: str) -> bool:
        """
        Check if a dataset exists.
        
        Args:
            name: Dataset name
            
        Returns:
            True if dataset exists
        """
        return self.dataset_manager.dataset_exists(name)