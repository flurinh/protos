"""
Data access module for the Protos framework.

This module provides functionality for accessing and managing datasets,
including registry management and data IO operations.
"""

import os
import json
import pickle
import logging
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple, Set
from pathlib import Path
from datetime import datetime
from abc import ABC, abstractmethod

from .paths import (
    ProtosPaths, 
    DataSource,
    ensure_directory
)

# Configure logger
logger = logging.getLogger(__name__)

class Dataset:
    """
    Standard dataset class for the Protos framework.
    
    A dataset is a collection of related data items (e.g., structures, sequences)
    with associated metadata. This class provides a consistent interface for
    working with datasets across different processor types.
    """
    
    def __init__(self, 
                id: str,
                name: str,
                description: str,
                type: str,
                content: Union[List, Dict, Set],
                metadata: Optional[Dict[str, Any]] = None):
        """
        Initialize a dataset with metadata and content.
        
        Args:
            id: Unique identifier for the dataset
            name: Human-readable name
            description: Detailed description
            type: Processor type ('structure', 'grn', 'sequence', etc.)
            content: Dataset content (list of IDs, dictionary of values, etc.)
            metadata: Additional metadata
        """
        self.id = id
        self.name = name
        self.description = description
        self.type = type
        self.content = content
        self.metadata = metadata or {}
        self.creation_date = datetime.now().isoformat()
        self.last_modified = self.creation_date
        
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert dataset to dictionary for serialization.
        
        Returns:
            Dictionary representation of the dataset
        """
        return {
            "id": self.id,
            "name": self.name,
            "description": self.description,
            "type": self.type,
            "content": self.content,
            "metadata": self.metadata,
            "creation_date": self.creation_date,
            "last_modified": self.last_modified
        }
        
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Dataset':
        """
        Create dataset from dictionary.
        
        Args:
            data: Dictionary containing dataset information
            
        Returns:
            Dataset instance
        """
        dataset = cls(
            id=data["id"],
            name=data["name"],
            description=data["description"],
            type=data["type"],
            content=data["content"],
            metadata=data.get("metadata", {})
        )
        
        # Restore timestamps if available
        if "creation_date" in data:
            dataset.creation_date = data["creation_date"]
        if "last_modified" in data:
            dataset.last_modified = data["last_modified"]
            
        return dataset
    
    def save(self, file_path: str) -> None:
        """
        Save dataset to file.
        
        Args:
            file_path: Path to save the dataset
        """
        # Update modification timestamp
        self.last_modified = datetime.now().isoformat()
        
        # Ensure directory exists
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        
        # Save to file
        with open(file_path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
    
    @classmethod
    def load(cls, file_path: str) -> 'Dataset':
        """
        Load dataset from file.
        
        Args:
            file_path: Path to the dataset file
            
        Returns:
            Dataset instance
        """
        with open(file_path, 'r') as f:
            data = json.load(f)
            
        return cls.from_dict(data)
    
    def add_item(self, item: Any) -> None:
        """
        Add an item to the dataset content.
        
        Args:
            item: Item to add
        """
        if isinstance(self.content, list):
            if item not in self.content:
                self.content.append(item)
        elif isinstance(self.content, dict):
            if isinstance(item, tuple) and len(item) == 2:
                key, value = item
                self.content[key] = value
            else:
                raise ValueError("For dictionary content, item must be a (key, value) tuple")
        elif isinstance(self.content, set):
            self.content.add(item)
        else:
            raise TypeError(f"Unsupported content type: {type(self.content)}")
        
        # Update modification timestamp
        self.last_modified = datetime.now().isoformat()
    
    def remove_item(self, item: Any) -> bool:
        """
        Remove an item from the dataset content.
        
        Args:
            item: Item to remove
            
        Returns:
            True if item was removed, False if not found
        """
        if isinstance(self.content, list):
            if item in self.content:
                self.content.remove(item)
                self.last_modified = datetime.now().isoformat()
                return True
        elif isinstance(self.content, dict):
            if item in self.content:
                del self.content[item]
                self.last_modified = datetime.now().isoformat()
                return True
        elif isinstance(self.content, set):
            if item in self.content:
                self.content.remove(item)
                self.last_modified = datetime.now().isoformat()
                return True
        
        return False
    
    def update_metadata(self, metadata: Dict[str, Any]) -> None:
        """
        Update dataset metadata.
        
        Args:
            metadata: New metadata to merge with existing
        """
        self.metadata.update(metadata)
        self.last_modified = datetime.now().isoformat()
    
    def __len__(self) -> int:
        """
        Get number of items in the dataset.
        
        Returns:
            Number of items
        """
        return len(self.content)
    
    def __contains__(self, item: Any) -> bool:
        """
        Check if item is in the dataset.
        
        Args:
            item: Item to check
            
        Returns:
            True if item is in the dataset
        """
        if isinstance(self.content, dict):
            return item in self.content
        else:
            return item in self.content
    
    def __iter__(self):
        """
        Iterate over dataset items.
        
        Returns:
            Iterator over dataset content
        """
        return iter(self.content)
    
    def __str__(self) -> str:
        """
        Get string representation of dataset.
        
        Returns:
            String description
        """
        item_count = len(self.content)
        return f"Dataset(id={self.id}, name={self.name}, type={self.type}, items={item_count})"

class DataRegistry:
    """
    Registry for mapping dataset identifiers to file paths.
    
    This class manages a registry of datasets, making it easy to
    reference data by logical identifiers rather than file paths.
    """
    
    def __init__(self, registry_file: Optional[str] = None):
        """
        Initialize the data registry.
        
        Args:
            registry_file: Path to registry JSON file (default: data/registry.json)
        """
        self.registry_file = registry_file or os.path.join('data', 'registry.json')
        self.registry = self._load_registry()
    
    def _load_registry(self) -> Dict[str, Dict[str, Any]]:
        """Load registry from file or create if not exists."""
        if os.path.exists(self.registry_file):
            try:
                with open(self.registry_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"Error loading registry: {e}")
                return {}
        else:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.registry_file), exist_ok=True)
            return {}
    
    def _save_registry(self) -> None:
        """Save registry to file."""
        try:
            with open(self.registry_file, 'w') as f:
                json.dump(self.registry, f, indent=2)
        except Exception as e:
            logger.error(f"Error saving registry: {e}")
    
    def register_dataset(self, 
                        dataset_id: str, 
                        file_path: str, 
                        metadata: Optional[Dict[str, Any]] = None) -> None:
        """
        Register a dataset in the registry.
        
        Args:
            dataset_id: Unique identifier for the dataset
            file_path: Path to the dataset file
            metadata: Additional metadata for the dataset
        """
        self.registry[dataset_id] = {
            'path': file_path,
            'metadata': metadata or {},
            'timestamp': datetime.now().isoformat()
        }
        self._save_registry()
    
    def get_dataset_path(self, dataset_id: str) -> Optional[str]:
        """
        Get the file path for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            File path or None if not found
        """
        if dataset_id in self.registry:
            return self.registry[dataset_id]['path']
        return None
    
    def get_dataset_metadata(self, dataset_id: str) -> Optional[Dict[str, Any]]:
        """
        Get metadata for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            Metadata dictionary or None if not found
        """
        if dataset_id in self.registry:
            return self.registry[dataset_id].get('metadata', {})
        return None
    
    def list_datasets(self) -> List[str]:
        """
        List all registered datasets.
        
        Returns:
            List of dataset identifiers
        """
        return list(self.registry.keys())
    
    def remove_dataset(self, dataset_id: str) -> bool:
        """
        Remove a dataset from the registry.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            True if dataset was removed, False if not found
        """
        if dataset_id in self.registry:
            del self.registry[dataset_id]
            self._save_registry()
            return True
        return False
    
    def get_datasets_by_type(self, dataset_type: str) -> List[str]:
        """
        Get datasets of a specific type.
        
        Args:
            dataset_type: Type to filter by
            
        Returns:
            List of matching dataset identifiers
        """
        return [
            dataset_id for dataset_id, info in self.registry.items()
            if info.get('metadata', {}).get('type') == dataset_type
        ]
    
    def update_metadata(self, dataset_id: str, metadata: Dict[str, Any]) -> bool:
        """
        Update metadata for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            metadata: New metadata (merged with existing)
            
        Returns:
            True if successful, False if dataset not found
        """
        if dataset_id in self.registry:
            if 'metadata' not in self.registry[dataset_id]:
                self.registry[dataset_id]['metadata'] = {}
            
            self.registry[dataset_id]['metadata'].update(metadata)
            self._save_registry()
            return True
        return False


class GlobalRegistry:
    """
    Global registry that manages datasets across multiple processor types
    and handles both reference and user data.
    
    This registry provides a unified view of all datasets in the system,
    regardless of where they are physically stored.
    """
    
    def __init__(self, paths: Optional[ProtosPaths] = None):
        """
        Initialize the global registry.
        
        Args:
            paths: ProtosPaths instance for path resolution
        """
        self.paths = paths or ProtosPaths()
        self.registry_file = self.paths.get_global_registry_path()
        self.registry = self._load_registry()
        
        # Cache of processor-specific registries
        self._processor_registries = {}
        
    def _load_registry(self) -> Dict[str, Dict[str, Any]]:
        """Load registry from file or create if not exists."""
        if os.path.exists(self.registry_file):
            try:
                with open(self.registry_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.error(f"Error loading global registry: {e}")
                return {}
        else:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.registry_file), exist_ok=True)
            # Initialize with empty registry
            with open(self.registry_file, 'w') as f:
                json.dump({}, f, indent=2)
            return {}
    
    def _save_registry(self) -> None:
        """Save registry to file."""
        try:
            with open(self.registry_file, 'w') as f:
                json.dump(self.registry, f, indent=2)
        except Exception as e:
            logger.error(f"Error saving global registry: {e}")
    
    def _get_processor_registry(self, processor_type: str) -> DataRegistry:
        """
        Get a processor-specific registry.
        
        Args:
            processor_type: Type of processor ('structure', 'grn', etc.)
            
        Returns:
            DataRegistry instance for the specified processor
        """
        if processor_type not in self._processor_registries:
            registry_file = self.paths.get_registry_path(processor_type, DataSource.USER)
            self._processor_registries[processor_type] = DataRegistry(registry_file)
        return self._processor_registries[processor_type]
    
    def register_dataset(self, 
                        dataset_id: str, 
                        file_path: str, 
                        processor_type: str,
                        dataset_type: Optional[str] = None,
                        source: DataSource = DataSource.USER,
                        metadata: Optional[Dict[str, Any]] = None) -> None:
        """
        Register a dataset in the global registry.
        
        Args:
            dataset_id: Unique identifier for the dataset
            file_path: Path to the dataset file
            processor_type: Type of processor that owns the dataset
            dataset_type: Type of dataset (e.g., 'structure', 'sequence')
            source: Source of the dataset (reference or user)
            metadata: Additional metadata for the dataset
        """
        metadata = metadata or {}
        metadata.update({
            'processor_type': processor_type,
            'dataset_type': dataset_type,
            'source': source.value
        })
        
        # Add to global registry
        self.registry[dataset_id] = {
            'path': file_path,
            'metadata': metadata,
            'timestamp': datetime.now().isoformat()
        }
        self._save_registry()
        
        # Add to processor-specific registry if it's user data
        if source == DataSource.USER:
            processor_registry = self._get_processor_registry(processor_type)
            processor_registry.register_dataset(dataset_id, file_path, metadata)
    
    def get_dataset_path(self, 
                        dataset_id: str, 
                        check_reference: bool = True) -> Optional[str]:
        """
        Get the file path for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            check_reference: Whether to check reference data if not found in user data
            
        Returns:
            File path or None if not found
        """
        # First check global registry
        if dataset_id in self.registry:
            return self.registry[dataset_id]['path']
        
        # Then check processor-specific registries if not found
        for processor_type in ['structure', 'grn', 'sequence', 'graph', 'property']:
            processor_registry = self._get_processor_registry(processor_type)
            path = processor_registry.get_dataset_path(dataset_id)
            if path is not None:
                # Add to global registry for future lookups
                metadata = processor_registry.get_dataset_metadata(dataset_id) or {}
                metadata['processor_type'] = processor_type
                metadata['source'] = DataSource.USER.value
                self.register_dataset(
                    dataset_id, path, processor_type, 
                    metadata.get('dataset_type'), DataSource.USER, metadata
                )
                return path
        
        # Not found
        return None
    
    def get_dataset_metadata(self, dataset_id: str) -> Optional[Dict[str, Any]]:
        """
        Get metadata for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            Metadata dictionary or None if not found
        """
        if dataset_id in self.registry:
            return self.registry[dataset_id].get('metadata', {})
        
        # Check processor-specific registries if not found
        for processor_type in ['structure', 'grn', 'sequence', 'graph', 'property']:
            processor_registry = self._get_processor_registry(processor_type)
            metadata = processor_registry.get_dataset_metadata(dataset_id)
            if metadata is not None:
                # Add to global registry for future lookups
                path = processor_registry.get_dataset_path(dataset_id)
                if path is not None:
                    metadata['processor_type'] = processor_type
                    metadata['source'] = DataSource.USER.value
                    self.register_dataset(
                        dataset_id, path, processor_type, 
                        metadata.get('dataset_type'), DataSource.USER, metadata
                    )
                return metadata
        
        return None
    
    def list_datasets(self, processor_type: Optional[str] = None) -> List[str]:
        """
        List all registered datasets.
        
        Args:
            processor_type: Optional processor type to filter by
            
        Returns:
            List of dataset identifiers
        """
        if processor_type is None:
            return list(self.registry.keys())
        else:
            return [
                dataset_id for dataset_id, info in self.registry.items()
                if info.get('metadata', {}).get('processor_type') == processor_type
            ]
    
    def remove_dataset(self, dataset_id: str) -> bool:
        """
        Remove a dataset from the registry.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            True if dataset was removed, False if not found
        """
        # Can only remove user datasets
        if dataset_id in self.registry:
            source = self.registry[dataset_id].get('metadata', {}).get('source')
            if source == DataSource.REFERENCE.value:
                logger.warning(f"Cannot remove reference dataset: {dataset_id}")
                return False
            
            processor_type = self.registry[dataset_id].get('metadata', {}).get('processor_type')
            if processor_type:
                processor_registry = self._get_processor_registry(processor_type)
                processor_registry.remove_dataset(dataset_id)
            
            del self.registry[dataset_id]
            self._save_registry()
            return True
        
        return False
    
    def get_datasets_by_type(self, dataset_type: str) -> List[str]:
        """
        Get datasets of a specific type.
        
        Args:
            dataset_type: Type to filter by
            
        Returns:
            List of matching dataset identifiers
        """
        return [
            dataset_id for dataset_id, info in self.registry.items()
            if info.get('metadata', {}).get('dataset_type') == dataset_type
        ]
    
    def update_metadata(self, dataset_id: str, metadata: Dict[str, Any]) -> bool:
        """
        Update metadata for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            metadata: New metadata (merged with existing)
            
        Returns:
            True if successful, False if dataset not found
        """
        if dataset_id in self.registry:
            # Cannot update certain metadata fields for reference datasets
            source = self.registry[dataset_id].get('metadata', {}).get('source')
            if source == DataSource.REFERENCE.value:
                protected_fields = ['source', 'processor_type', 'path']
                for field in protected_fields:
                    if field in metadata:
                        logger.warning(f"Cannot update {field} for reference dataset: {dataset_id}")
                        metadata.pop(field)
            
            if 'metadata' not in self.registry[dataset_id]:
                self.registry[dataset_id]['metadata'] = {}
            
            self.registry[dataset_id]['metadata'].update(metadata)
            self._save_registry()
            
            # Update processor-specific registry if it's user data
            if source == DataSource.USER.value:
                processor_type = self.registry[dataset_id].get('metadata', {}).get('processor_type')
                if processor_type:
                    processor_registry = self._get_processor_registry(processor_type)
                    processor_registry.update_metadata(dataset_id, metadata)
            
            return True
        
        return False
    
    def get_datasets_by_source(self, source: DataSource) -> List[str]:
        """
        Get datasets from a specific source.
        
        Args:
            source: Data source to filter by
            
        Returns:
            List of matching dataset identifiers
        """
        return [
            dataset_id for dataset_id, info in self.registry.items()
            if info.get('metadata', {}).get('source') == source.value
        ]
    
    def import_reference_data(self) -> int:
        """
        Import reference data into the registry.
        
        Scans the reference data directory for datasets and adds them
        to the global registry.
        
        Returns:
            Number of reference datasets imported
        """
        count = 0
        
        # Scan reference data directory for processor types
        ref_root = self.paths.ref_data_root
        for processor_type in os.listdir(ref_root):
            processor_path = os.path.join(ref_root, processor_type)
            if not os.path.isdir(processor_path):
                continue
            
            # Check for registry file
            registry_file = os.path.join(processor_path, 'registry.json')
            if os.path.exists(registry_file):
                try:
                    with open(registry_file, 'r') as f:
                        processor_registry = json.load(f)
                    
                    # Import datasets from this registry
                    for dataset_id, info in processor_registry.items():
                        file_path = info.get('path')
                        if file_path:
                            # Make path absolute if it's relative
                            if not os.path.isabs(file_path):
                                file_path = os.path.join(processor_path, file_path)
                            
                            metadata = info.get('metadata', {})
                            dataset_type = metadata.get('dataset_type', processor_type)
                            
                            # Add to global registry
                            self.register_dataset(
                                dataset_id, file_path, processor_type, dataset_type,
                                DataSource.REFERENCE, metadata
                            )
                            count += 1
                
                except Exception as e:
                    logger.error(f"Error importing reference data from {registry_file}: {e}")
        
        return count