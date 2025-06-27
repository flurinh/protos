"""
Data access module for the Protos framework.

This module provides functionality for accessing and managing datasets,
including registry management and data IO operations.
"""

import os
import json
import pickle
import logging
import hashlib
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple, Set
from pathlib import Path
from datetime import datetime
from abc import ABC, abstractmethod

from .paths import (
    ProtosPaths, 
    ensure_directory
)
from .paths.path_config import DataSource

# Configure logger
logger = logging.getLogger(__name__)


def generate_entity_id(content: str, prefix: Optional[str] = None) -> str:
    """
    Generate a standardized 10-character entity ID from content.
    
    For universal entity IDs across formats, do NOT use prefix.
    The same biological entity should have the same ID regardless of format.
    
    Args:
        content: String to hash (e.g., UniProt ID, PDB ID, sequence ID)
        prefix: DEPRECATED - kept for compatibility but ignored for universal IDs
        
    Returns:
        10-character hash ID
    """
    # For universal IDs, we ONLY hash the content, no prefix
    # This ensures P12345 gets the same ID whether it's a sequence, structure, or GRN
    
    # Generate SHA-256 hash
    hash_obj = hashlib.sha256(content.encode('utf-8'))
    
    # Take first 10 characters of hex digest
    entity_id = hash_obj.hexdigest()[:10]
    
    return entity_id

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
        
        # Initialize entity registry
        entity_registry_path = os.path.join(os.path.dirname(self.registry_file), 'entity_registry.json')
        self.entity_registry = EntityRegistry(entity_registry_path)
        
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
            registry_file = self.paths.get_registry_path(processor_type)
            self._processor_registries[processor_type] = DataRegistry(registry_file)
        return self._processor_registries[processor_type]
    
    def register_dataset(self, 
                        dataset_id: str, 
                        file_path: str, 
                        processor_type: str,
                        dataset_type: Optional[str] = None,
                        source: str = 'user',
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
            'source': source
        })
        
        # Add to global registry
        self.registry[dataset_id] = {
            'path': file_path,
            'metadata': metadata,
            'timestamp': datetime.now().isoformat()
        }
        self._save_registry()
        
        # Add to processor-specific registry if it's user data
        if source == 'user':
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
                metadata['source'] = 'user'
                self.register_dataset(
                    dataset_id, path, processor_type, 
                    metadata.get('dataset_type'), 'user', metadata
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
                    metadata['source'] = 'user'
                    self.register_dataset(
                        dataset_id, path, processor_type, 
                        metadata.get('dataset_type'), 'user', metadata
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
            if source == 'reference':
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
            if source == 'reference':
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
            if source == 'user':
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
                                'reference', metadata
                            )
                            count += 1
                
                except Exception as e:
                    logger.error(f"Error importing reference data from {registry_file}: {e}")
        
        return count


class EntityRegistry:
    """
    Registry for managing individual entities with hash-based IDs.
    
    This registry tracks individual data items (entities) across all processors,
    using consistent 10-character hash IDs for identification. It maintains
    relationships between entities and datasets, and maps original IDs to hash IDs.
    
    Supports multiple formats per entity - the same biological entity can have
    sequence, structure, GRN, and embedding formats all sharing the same entity ID.
    """
    
    def __init__(self, registry_file: Optional[str] = None):
        """
        Initialize the entity registry.
        
        Args:
            registry_file: Path to entity registry JSON file
        """
        self.registry_file = registry_file or os.path.join('data', 'entity_registry.json')
        self.entities = self._load_registry()
        
        # Cache for reverse lookups (original_id -> entity_id)
        self._original_id_cache = self._build_original_id_cache()
        
        # Migrate old format if needed
        self._migrate_to_multi_format()
    
    def _load_registry(self) -> Dict[str, Dict[str, Any]]:
        """Load entity registry from file or create if not exists."""
        if os.path.exists(self.registry_file):
            try:
                with open(self.registry_file, 'r') as f:
                    data = json.load(f)
                    return data.get('entities', {})
            except Exception as e:
                logger.error(f"Error loading entity registry: {e}")
                return {}
        else:
            # Ensure directory exists
            os.makedirs(os.path.dirname(self.registry_file), exist_ok=True)
            return {}
    
    def _save_registry(self) -> None:
        """Save entity registry to file."""
        try:
            # Get all datasets that reference entities
            datasets = self._extract_datasets_info()
            
            registry_data = {
                'entities': self.entities,
                'datasets': datasets,
                'version': '1.0',
                'last_modified': datetime.now().isoformat()
            }
            
            with open(self.registry_file, 'w') as f:
                json.dump(registry_data, f, indent=2)
        except Exception as e:
            logger.error(f"Error saving entity registry: {e}")
    
    def _build_original_id_cache(self) -> Dict[str, str]:
        """Build cache for reverse lookups."""
        cache = {}
        for entity_id, info in self.entities.items():
            original_id = info.get('original_id')
            if original_id:
                # For multi-format entities, we don't use type prefix
                cache[original_id] = entity_id
        return cache
    
    def _migrate_to_multi_format(self) -> None:
        """Migrate old single-format entities to multi-format structure."""
        needs_save = False
        for entity_id, info in self.entities.items():
            # Check if this is old format (has 'type' field)
            if 'type' in info and 'formats' not in info:
                # Migrate to new format
                entity_type = info.pop('type')
                file_path = info.get('file_path')
                metadata = info.get('metadata', {})
                
                # Create formats structure
                info['formats'] = {
                    entity_type: {
                        'file_path': file_path,
                        'metadata': metadata,
                        'added_at': info.get('created', datetime.now().isoformat())
                    }
                }
                
                # Remove old fields
                if 'file_path' in info:
                    del info['file_path']
                if 'metadata' in info:
                    del info['metadata']
                    
                needs_save = True
        
        if needs_save:
            self._save_registry()
    
    def _extract_datasets_info(self) -> Dict[str, Dict[str, Any]]:
        """Extract dataset information from entity references."""
        datasets = {}
        
        # Collect all unique datasets from entities
        for entity_id, info in self.entities.items():
            entity_datasets = info.get('datasets', [])
            for dataset_id in entity_datasets:
                if dataset_id not in datasets:
                    datasets[dataset_id] = {
                        'entities': [],
                        'formats': set()
                    }
                datasets[dataset_id]['entities'].append(entity_id)
                # Track all formats in this dataset
                for fmt in info.get('formats', {}):
                    datasets[dataset_id]['formats'].add(fmt)
        
        # Convert sets to lists for JSON serialization
        for dataset_id in datasets:
            datasets[dataset_id]['formats'] = list(datasets[dataset_id]['formats'])
        
        return datasets
    
    def register_entity(self,
                       entity_id: str,
                       entity_type: str,
                       original_id: Optional[str] = None,
                       file_path: Optional[str] = None,
                       metadata: Optional[Dict[str, Any]] = None,
                       datasets: Optional[List[str]] = None) -> None:
        """
        Register an entity in the registry or add a new format to existing entity.
        
        Args:
            entity_id: Hash-based entity ID
            entity_type: Type of entity (structure, sequence, grn, embedding)
            original_id: Original identifier (e.g., PDB ID, sequence ID)
            file_path: Path to entity file
            metadata: Additional metadata
            datasets: List of dataset IDs this entity belongs to
        """
        # Initialize entity if new
        if entity_id not in self.entities:
            self.entities[entity_id] = {
                'original_id': original_id,
                'formats': {},
                'datasets': [],
                'created': datetime.now().isoformat(),
                'modified': datetime.now().isoformat()
            }
        
        # Add/update format information
        self.entities[entity_id]['formats'][entity_type] = {
            'file_path': file_path,
            'metadata': metadata or {},
            'added_at': datetime.now().isoformat()
        }
        
        # Update original_id if provided and not set
        if original_id and not self.entities[entity_id].get('original_id'):
            self.entities[entity_id]['original_id'] = original_id
        
        # Update datasets
        if datasets:
            for dataset_id in datasets:
                if dataset_id not in self.entities[entity_id]['datasets']:
                    self.entities[entity_id]['datasets'].append(dataset_id)
        
        # Update modified time
        self.entities[entity_id]['modified'] = datetime.now().isoformat()
        
        # Update cache
        if original_id:
            self._original_id_cache[original_id] = entity_id
        
        self._save_registry()
    
    def get_entity(self, entity_id: str) -> Optional[Dict[str, Any]]:
        """
        Get entity information by ID.
        
        Args:
            entity_id: Entity hash ID
            
        Returns:
            Entity information or None if not found
        """
        return self.entities.get(entity_id)
    
    def register_entity_format(self,
                             entity_id: str,
                             format_type: str,
                             file_path: str,
                             metadata: Optional[Dict[str, Any]] = None) -> bool:
        """
        Register a new format for an existing entity.
        
        Args:
            entity_id: Entity hash ID
            format_type: Format type (sequence, structure, grn, embedding)
            file_path: Path to file for this format
            metadata: Format-specific metadata
            
        Returns:
            True if successful
        """
        if entity_id not in self.entities:
            return False
            
        self.entities[entity_id]['formats'][format_type] = {
            'file_path': file_path,
            'metadata': metadata or {},
            'added_at': datetime.now().isoformat()
        }
        self.entities[entity_id]['modified'] = datetime.now().isoformat()
        self._save_registry()
        return True
    
    def get_entity_formats(self, entity_id: str) -> List[str]:
        """
        Get all available formats for an entity.
        
        Args:
            entity_id: Entity hash ID
            
        Returns:
            List of format types available
        """
        if entity_id not in self.entities:
            return []
        return list(self.entities[entity_id].get('formats', {}).keys())
    
    def get_entity_by_format(self, entity_id: str, format_type: str) -> Optional[Dict[str, Any]]:
        """
        Get entity info for a specific format.
        
        Args:
            entity_id: Entity hash ID
            format_type: Format type to retrieve
            
        Returns:
            Format-specific information or None
        """
        if entity_id not in self.entities:
            return None
        
        formats = self.entities[entity_id].get('formats', {})
        if format_type not in formats:
            return None
            
        return {
            'entity_id': entity_id,
            'original_id': self.entities[entity_id].get('original_id'),
            'format': format_type,
            **formats[format_type]
        }
    
    def resolve_identifier(self, identifier: str, format_type: Optional[str] = None) -> str:
        """
        Resolve any identifier to entity hash ID.
        
        This is the universal translator that allows users to work with
        biological names while the system uses hash IDs internally.
        
        Args:
            identifier: Any identifier (PDB ID, UniProt ID, sequence name, or hash)
            format_type: Optional format type for disambiguation
            
        Returns:
            Entity hash ID (creates new if not found)
            
        Examples:
            resolve_identifier("1ABC") → "a3f2d8c91b"
            resolve_identifier("P12345") → "b4e5f6a72c"
            resolve_identifier("a3f2d8c91b") → "a3f2d8c91b" (already a hash)
        """
        # Check if already a hash ID (10 chars, alphanumeric)
        if len(identifier) == 10 and identifier.isalnum():
            if self.entity_exists(identifier):
                # Valid existing hash ID
                if format_type is None:
                    return identifier
                # If format specified, verify entity has that format
                elif format_type in self.get_entity_formats(identifier):
                    return identifier
        
        # Try to find by original ID
        entity_id = self.find_entity_by_original_id(identifier, format_type)
        if entity_id:
            return entity_id
        
        # Not found - generate new hash ID from identifier
        # This allows automatic registration of new entities
        return generate_entity_id(identifier)
    
    def get_original_id(self, entity_id: str) -> Optional[str]:
        """
        Get the original identifier for an entity.
        
        Args:
            entity_id: Entity hash ID
            
        Returns:
            Original identifier or None if not found
        """
        if entity_id in self.entities:
            return self.entities[entity_id].get('original_id')
        return None
    
    def find_entity_by_original_id(self, original_id: str, format_type: Optional[str] = None) -> Optional[str]:
        """
        Find entity hash ID by original ID.
        
        Args:
            original_id: Original identifier (e.g., PDB ID)
            format_type: Optional format type for filtering
            
        Returns:
            Entity hash ID or None if not found
        """
        # Direct cache lookup (no type prefix in multi-format system)
        if original_id in self._original_id_cache:
            entity_id = self._original_id_cache[original_id]
            
            # If format_type specified, verify entity has that format
            if format_type:
                if format_type in self.entities[entity_id].get('formats', {}):
                    return entity_id
                else:
                    return None
            return entity_id
        
        # Fallback to searching all entities
        for entity_id, info in self.entities.items():
            if info.get('original_id') == original_id:
                if format_type is None or format_type in info.get('formats', {}):
                    return entity_id
        
        return None
    
    def list_entities(self, 
                     format_type: Optional[str] = None,
                     dataset: Optional[str] = None) -> List[str]:
        """
        List entity IDs, optionally filtered.
        
        Args:
            format_type: Filter by format type (structure, sequence, etc.)
            dataset: Filter by dataset membership
            
        Returns:
            List of entity IDs
        """
        entity_ids = []
        
        for entity_id, info in self.entities.items():
            # Filter by format if specified
            if format_type and format_type not in info.get('formats', {}):
                continue
            
            # Filter by dataset if specified
            if dataset and dataset not in info.get('datasets', []):
                continue
            
            entity_ids.append(entity_id)
        
        return entity_ids
    
    def add_entity_to_dataset(self, entity_id: str, dataset_id: str) -> bool:
        """
        Add an entity to a dataset.
        
        Args:
            entity_id: Entity hash ID
            dataset_id: Dataset ID
            
        Returns:
            True if successful
        """
        if entity_id not in self.entities:
            return False
        
        datasets = self.entities[entity_id].get('datasets', [])
        if dataset_id not in datasets:
            datasets.append(dataset_id)
            self.entities[entity_id]['datasets'] = datasets
            self.entities[entity_id]['modified'] = datetime.now().isoformat()
            self._save_registry()
        
        return True
    
    def remove_entity_from_dataset(self, entity_id: str, dataset_id: str) -> bool:
        """
        Remove an entity from a dataset.
        
        Args:
            entity_id: Entity hash ID
            dataset_id: Dataset ID
            
        Returns:
            True if successful
        """
        if entity_id not in self.entities:
            return False
        
        datasets = self.entities[entity_id].get('datasets', [])
        if dataset_id in datasets:
            datasets.remove(dataset_id)
            self.entities[entity_id]['datasets'] = datasets
            self.entities[entity_id]['modified'] = datetime.now().isoformat()
            self._save_registry()
            return True
        
        return False
    
    def get_dataset_entities(self, dataset_id: str) -> List[str]:
        """
        Get all entities in a dataset.
        
        Args:
            dataset_id: Dataset ID
            
        Returns:
            List of entity IDs in the dataset
        """
        return [
            entity_id for entity_id, info in self.entities.items()
            if dataset_id in info.get('datasets', [])
        ]
    
    def update_entity_metadata(self, entity_id: str, metadata: Dict[str, Any], format_type: Optional[str] = None) -> bool:
        """
        Update entity metadata.
        
        Args:
            entity_id: Entity hash ID
            metadata: New metadata to merge
            format_type: Optional format to update metadata for (updates all if None)
            
        Returns:
            True if successful
        """
        if entity_id not in self.entities:
            return False
        
        if format_type:
            # Update metadata for specific format
            if format_type in self.entities[entity_id].get('formats', {}):
                self.entities[entity_id]['formats'][format_type]['metadata'].update(metadata)
            else:
                return False
        else:
            # Update metadata for all formats
            for fmt in self.entities[entity_id].get('formats', {}):
                self.entities[entity_id]['formats'][fmt]['metadata'].update(metadata)
        
        self.entities[entity_id]['modified'] = datetime.now().isoformat()
        self._save_registry()
        
        return True
    
    def remove_entity(self, entity_id: str) -> bool:
        """
        Remove an entity from the registry.
        
        Args:
            entity_id: Entity hash ID
            
        Returns:
            True if entity was removed
        """
        if entity_id in self.entities:
            # Remove from cache
            original_id = self.entities[entity_id].get('original_id')
            if original_id:
                self._original_id_cache.pop(original_id, None)
            
            # Remove entity
            del self.entities[entity_id]
            self._save_registry()
            return True
        
        return False
    
    def entity_exists(self, entity_id: str) -> bool:
        """
        Check if an entity exists.
        
        Args:
            entity_id: Entity hash ID
            
        Returns:
            True if entity exists
        """
        return entity_id in self.entities
    
    def get_entity_stats(self) -> Dict[str, Any]:
        """
        Get statistics about entities in the registry.
        
        Returns:
            Dictionary with entity statistics
        """
        stats = {
            'total_entities': len(self.entities),
            'by_format': {},
            'by_dataset': {},
            'multi_format_entities': 0,
            'orphaned': 0  # Entities not in any dataset
        }
        
        # Count by format and multi-format entities
        for entity_id, info in self.entities.items():
            formats = info.get('formats', {})
            if len(formats) > 1:
                stats['multi_format_entities'] += 1
                
            for format_type in formats:
                stats['by_format'][format_type] = stats['by_format'].get(format_type, 0) + 1
            
            # Count by dataset
            datasets = info.get('datasets', [])
            if not datasets:
                stats['orphaned'] += 1
            for dataset_id in datasets:
                stats['by_dataset'][dataset_id] = stats['by_dataset'].get(dataset_id, 0) + 1
        
        return stats