"""
Dataset Manager for organizing collections of entities.

This module handles dataset creation, management, and operations,
working with EntityRegistry to track entities within datasets.
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Any, Set
from datetime import datetime

from .paths import ProtosPaths
from .entity_registry import EntityRegistry


class DatasetManager:
    """
    Manages datasets - collections of entities for a specific processor.
    
    Key principles:
    - Datasets are stored as individual JSON files in processor_dir/datasets/
    - Each dataset file contains entity names (human-readable)
    - Works with EntityRegistry to resolve entities
    """
    
    def __init__(self, processor_type: str, paths: Optional[ProtosPaths] = None,
                 entity_registry: Optional[EntityRegistry] = None):
        """
        Initialize DatasetManager with ProtosPaths and EntityRegistry.
        
        Args:
            processor_type: Type of processor (e.g., 'structure', 'sequence')
            paths: ProtosPaths instance for path management
            entity_registry: EntityRegistry instance for entity tracking
        """
        self.processor_type = processor_type
        self.paths = paths or ProtosPaths()
        self.entity_registry = entity_registry or EntityRegistry(self.paths)
        
        # Get datasets directory for this processor
        self.datasets_dir = Path(self.paths.get_dataset_path(processor_type, "dummy")[:-10])  # Remove dummy.json
        self.datasets_dir.mkdir(parents=True, exist_ok=True)
    
    def create_dataset(self, name: str, entities: List[str], 
                      metadata: Optional[Dict[str, Any]] = None) -> str:
        """
        Create a new dataset with given entities.
        
        Args:
            name: Dataset name (human-readable)
            entities: List of entity names (human-readable)
            metadata: Optional metadata for the dataset
            
        Returns:
            Dataset name
        """
        metadata = metadata or {}
        
        # Validate all entities exist
        missing_entities = []
        for entity in entities:
            if not self.entity_registry.entity_exists(entity):
                missing_entities.append(entity)
        
        if missing_entities:
            print(f"Warning: The following entities are not registered: {missing_entities}")
        
        # Create dataset structure
        dataset = {
            "name": name,
            "processor_type": self.processor_type,
            "entities": entities,  # Store human-readable names
            "metadata": metadata,
            "created": datetime.now().isoformat(),
            "modified": datetime.now().isoformat()
        }
        
        # Save to JSON file
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, name))
        with open(dataset_path, 'w') as f:
            json.dump(dataset, f, indent=2)
        
        return name
    
    def load_dataset(self, name: str) -> Dict[str, Any]:
        """
        Load dataset configuration.
        
        Args:
            name: Dataset name
            
        Returns:
            Dataset configuration dict
        """
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, name))
        
        if not dataset_path.exists():
            raise FileNotFoundError(f"Dataset '{name}' not found")
        
        with open(dataset_path, 'r') as f:
            return json.load(f)
    
    def list_datasets(self) -> List[str]:
        """
        List all datasets for this processor.
        
        Returns:
            List of dataset names
        """
        datasets = []
        
        if self.datasets_dir.exists():
            for file in self.datasets_dir.glob("*.json"):
                # Extract name without extension
                dataset_name = file.stem
                datasets.append(dataset_name)
        
        return sorted(datasets)
    
    def dataset_exists(self, name: str) -> bool:
        """
        Check if dataset exists.
        
        Args:
            name: Dataset name
            
        Returns:
            True if dataset exists
        """
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, name))
        return dataset_path.exists()
    
    def add_to_dataset(self, dataset_name: str, entities: List[str]):
        """
        Add entities to an existing dataset.
        
        Args:
            dataset_name: Dataset name
            entities: List of entity names to add
        """
        # Load existing dataset
        dataset = self.load_dataset(dataset_name)
        
        # Add new entities (avoid duplicates)
        existing = set(dataset['entities'])
        for entity in entities:
            if entity not in existing:
                dataset['entities'].append(entity)
        
        # Update modified time
        dataset['modified'] = datetime.now().isoformat()
        
        # Save updated dataset
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, dataset_name))
        with open(dataset_path, 'w') as f:
            json.dump(dataset, f, indent=2)
    
    def remove_from_dataset(self, dataset_name: str, entities: List[str]):
        """
        Remove entities from a dataset.
        
        Args:
            dataset_name: Dataset name
            entities: List of entity names to remove
        """
        # Load existing dataset
        dataset = self.load_dataset(dataset_name)
        
        # Remove entities
        entity_set = set(entities)
        dataset['entities'] = [e for e in dataset['entities'] if e not in entity_set]
        
        # Update modified time
        dataset['modified'] = datetime.now().isoformat()
        
        # Save updated dataset
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, dataset_name))
        with open(dataset_path, 'w') as f:
            json.dump(dataset, f, indent=2)
    
    def get_dataset_info(self, name: str) -> Dict[str, Any]:
        """
        Get detailed information about a dataset.
        
        Args:
            name: Dataset name
            
        Returns:
            Dataset information including entity details
        """
        dataset = self.load_dataset(name)
        
        # Add entity information
        entity_info = []
        for entity_name in dataset['entities']:
            info = self.entity_registry.find_entity(entity_name)
            if info:
                entity_info.append({
                    'name': entity_name,
                    'formats': self.entity_registry.get_entity_formats(entity_name)
                })
            else:
                entity_info.append({
                    'name': entity_name,
                    'formats': [],
                    'missing': True
                })
        
        return {
            'name': dataset['name'],
            'processor_type': dataset['processor_type'],
            'entity_count': len(dataset['entities']),
            'entities': entity_info,
            'metadata': dataset.get('metadata', {}),
            'created': dataset.get('created'),
            'modified': dataset.get('modified')
        }
    
    def delete_dataset(self, name: str):
        """
        Delete a dataset.
        
        Args:
            name: Dataset name
        """
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, name))
        
        if dataset_path.exists():
            os.remove(dataset_path)
    
    def get_dataset_entities(self, name: str) -> List[str]:
        """
        Get list of entities in a dataset.
        
        Args:
            name: Dataset name
            
        Returns:
            List of entity names
        """
        dataset = self.load_dataset(name)
        return dataset['entities']
    
    def update_metadata(self, name: str, metadata: Dict[str, Any]):
        """
        Update dataset metadata.
        
        Args:
            name: Dataset name
            metadata: New metadata to merge
        """
        dataset = self.load_dataset(name)
        
        # Update metadata
        dataset['metadata'].update(metadata)
        dataset['modified'] = datetime.now().isoformat()
        
        # Save updated dataset
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, name))
        with open(dataset_path, 'w') as f:
            json.dump(dataset, f, indent=2)
    
    def copy_dataset(self, source_name: str, target_name: str):
        """
        Create a copy of an existing dataset.
        
        Args:
            source_name: Source dataset name
            target_name: Target dataset name
        """
        # Load source dataset
        source = self.load_dataset(source_name)
        
        # Create new dataset with same entities
        self.create_dataset(
            target_name,
            source['entities'],
            source.get('metadata', {})
        )
    
    def merge_datasets(self, dataset_names: List[str], target_name: str):
        """
        Merge multiple datasets into a new one.
        
        Args:
            dataset_names: List of dataset names to merge
            target_name: Name for merged dataset
        """
        all_entities = set()
        merged_metadata = {}
        
        # Collect all entities and metadata
        for dataset_name in dataset_names:
            dataset = self.load_dataset(dataset_name)
            all_entities.update(dataset['entities'])
            
            # Merge metadata with source tracking
            for key, value in dataset.get('metadata', {}).items():
                if key not in merged_metadata:
                    merged_metadata[key] = value
        
        # Add merge information
        merged_metadata['merged_from'] = dataset_names
        merged_metadata['merge_date'] = datetime.now().isoformat()
        
        # Create merged dataset
        self.create_dataset(
            target_name,
            list(all_entities),
            merged_metadata
        )