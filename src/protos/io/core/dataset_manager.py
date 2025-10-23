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

from ..paths import get_protos_paths


class DatasetManager:
    """
    Manages datasets - collections of entities for a specific processor.
    
    Key principles:
    - Datasets are stored as individual JSON files in processor_dir/datasets/
    - Each dataset file contains entity names (human-readable)
    - Works with EntityRegistry to resolve entities
    """
    
    def __init__(self, processor_type: str):
        """
        Initialize DatasetManager - gets dependencies internally.
        
        Args:
            processor_type: Type of processor (e.g., 'structure', 'sequence')
        """
        self.processor_type = processor_type
        
        # Get singletons - no need to pass them
        self.paths = get_protos_paths()
        
        # Import here to avoid circular import
        from protos.io.core import get_registry
        self.entity_registry = get_registry()
        
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
        
        # Validate all entities exist and get their IDs
        entity_ids = []
        entity_names = []
        missing_entities = []
        
        for entity in entities:
            entity_info = self.entity_registry.find_entity(entity)
            if entity_info:
                entity_ids.append(entity_info.entity_id)
                entity_names.append(entity)
            else:
                missing_entities.append(entity)
        
        if missing_entities:
            import warnings
            warnings.warn(f"The following entities are not registered: {missing_entities}")
        
        # Create dataset structure with both IDs and names
        dataset = {
            "name": name,
            "processor_type": self.processor_type,
            "entities": entity_names,  # Keep for backward compatibility
            "entity_ids": entity_ids,  # New: store stable IDs
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
        
        # Handle both old and new dataset formats
        if 'entity_ids' in dataset:
            # New format: work with IDs
            existing_ids = set(dataset['entity_ids'])
            existing_names = set(dataset.get('entities', []))
            
            for entity in entities:
                entity_info = self.entity_registry.find_entity(entity)
                if entity_info and entity_info.entity_id not in existing_ids:
                    dataset['entity_ids'].append(entity_info.entity_id)
                    if entity not in existing_names:
                        dataset['entities'].append(entity)
        else:
            # Old format: work with names only
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
        
        # Get IDs of entities to remove
        ids_to_remove = set()
        for entity in entities:
            entity_info = self.entity_registry.find_entity(entity)
            if entity_info:
                ids_to_remove.add(entity_info.entity_id)
        
        if 'entity_ids' in dataset:
            # New format: remove by ID
            dataset['entity_ids'] = [eid for eid in dataset['entity_ids'] 
                                   if eid not in ids_to_remove]
            # Also update names for backward compatibility
            entity_set = set(entities)
            dataset['entities'] = [e for e in dataset['entities'] if e not in entity_set]
        else:
            # Old format: remove by name
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
        
        if 'entity_ids' in dataset:
            # New format: resolve IDs to current info
            for entity_id in dataset['entity_ids']:
                entity_data = self.entity_registry._registry.get(entity_id)
                if entity_data:
                    entity_name = entity_data['original_id']
                    entity_info.append({
                        'name': entity_name,
                        'formats': self.entity_registry.get_entity_formats(entity_name)
                    })
                else:
                    # Entity was deleted
                    # Try to get historic name from dataset
                    historic_idx = dataset['entity_ids'].index(entity_id)
                    historic_name = (dataset.get('entities', [])[historic_idx] 
                                   if historic_idx < len(dataset.get('entities', [])) 
                                   else 'Unknown')
                    entity_info.append({
                        'name': historic_name,
                        'formats': [],
                        'missing': True
                    })
        else:
            # Old format: use stored names
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
            'entity_count': len(dataset.get('entity_ids', dataset['entities'])),
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
            List of entity names (current names if entities exist, historic if deleted)
        """
        dataset = self.load_dataset(name)
        
        # If dataset has entity_ids, resolve to current names
        if 'entity_ids' in dataset:
            current_names = []
            for idx, entity_id in enumerate(dataset['entity_ids']):
                # Find entity by ID
                entity_data = self.entity_registry._registry.get(entity_id)
                if entity_data:
                    current_names.append(entity_data['original_id'])
                else:
                    # Entity was deleted - use historic name if available
                    if idx < len(dataset.get('entities', [])):
                        current_names.append(dataset['entities'][idx])
                    else:
                        # No historic name available
                        current_names.append(f"<deleted:{entity_id[:8]}>")
            return current_names
        else:
            # Backward compatibility: return stored names
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
    
    def refresh_dataset_entities(self, name: str):
        """
        Refresh cached entity names in a dataset from the registry.
        
        This updates the 'entities' array to match current entity names
        based on the entity IDs. Useful after entity renames.
        
        Args:
            name: Dataset name
        """
        dataset = self.load_dataset(name)
        
        # Only refresh if dataset has entity_ids
        if 'entity_ids' not in dataset:
            return  # Old format, nothing to refresh
            
        # Rebuild entities array from current registry state
        current_names = []
        for entity_id in dataset['entity_ids']:
            entity_data = self.entity_registry._registry.get(entity_id)
            if entity_data:
                current_names.append(entity_data['original_id'])
            else:
                # Keep historic name if entity was deleted
                # Find corresponding historic name
                try:
                    idx = dataset['entity_ids'].index(entity_id)
                    if idx < len(dataset.get('entities', [])):
                        current_names.append(dataset['entities'][idx])
                    else:
                        current_names.append(f"<deleted:{entity_id[:8]}>")
                except (ValueError, IndexError):
                    current_names.append(f"<deleted:{entity_id[:8]}>")
        
        # Update dataset
        dataset['entities'] = current_names
        dataset['modified'] = datetime.now().isoformat()
        
        # Save updated dataset
        dataset_path = Path(self.paths.get_dataset_path(self.processor_type, name))
        with open(dataset_path, 'w') as f:
            json.dump(dataset, f, indent=2)
    
    def refresh_all_datasets(self):
        """
        Refresh entity names in all datasets for this processor.
        
        Useful after bulk entity renames or registry updates.
        """
        for dataset_name in self.list_datasets():
            try:
                self.refresh_dataset_entities(dataset_name)
            except Exception as e:
                print(f"Warning: Could not refresh dataset '{dataset_name}': {e}")
