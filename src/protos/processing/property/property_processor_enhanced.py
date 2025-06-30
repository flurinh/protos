#!/usr/bin/env python3
"""
Enhanced PropertyProcessor for the Protos Framework

This processor manages properties assigned to entities across all format types.
Properties can be assigned to any entity (structure, sequence, grn, embedding, etc.)
and organized into named datasets for analysis and comparison.

Key Features:
- Entity-based property assignment using deterministic entity IDs
- Cross-format property support (same entity across multiple formats)
- Property dataset management and organization
- Advanced secondary selection (GRN columns, CIF atoms, etc.)
- Comprehensive metadata tracking and validation
- Integration with Protos entity registry system

Example Usage:
    # Initialize processor
    prop_proc = PropertyProcessor(name="opsin_properties")
    
    # Assign properties to entities
    prop_proc.assign_property("36c2c0da93", "lambda_max", 500.0, 
                             dataset="microbial_opsins")
    
    # Load property dataset
    prop_proc.load_property_dataset("microbial_opsins")
    
    # Query properties by entity
    properties = prop_proc.get_entity_properties("36c2c0da93")
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional, Union, Tuple

from protos.core.base_processor import BaseProcessor
from protos.io.data_access import generate_entity_id, EntityRegistry, GlobalRegistry
from protos.io.file_utils import save_json, load_json

logger = logging.getLogger(__name__)


class PropertyProcessor(BaseProcessor):
    """
    Enhanced Property Processor for assigning and managing entity properties.
    
    This processor handles property assignment to entities across all format types,
    with support for dataset organization and advanced selection criteria.
    """
    
    def __init__(self, name: str = "property_processor", 
                 processor_data_dir: str = "property", 
                 **kwargs):
        """
        Initialize the Property Processor.
        
        Args:
            name: Processor instance name
            processor_data_dir: Subdirectory for property data
            **kwargs: Additional arguments for BaseProcessor
        """
        super().__init__(
            name=name,
            processor_data_dir=processor_data_dir,
            **kwargs
        )
        
        # Set up data directories
        self.data_dirs = {
            'datasets': Path(self.data_path) / 'datasets',
            'metadata': Path(self.data_path) / 'metadata',
            'assignments': Path(self.data_path) / 'assignments',
            'cache': Path(self.data_path) / 'cache'
        }
        
        # Ensure directories exist
        for dir_path in self.data_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Initialize property tracking
        self.property_datasets = {}  # dataset_name -> property data
        self.entity_properties = {}  # entity_id -> {dataset: {prop: value}}
        self.property_metadata = {}  # property_name -> metadata
        self.dataset_metadata = {}   # dataset_name -> metadata
        
        # Secondary selection support (advanced feature)
        self.secondary_selections = {}  # entity_id -> selection criteria
        
        # Load existing data if available
        self._load_existing_data()
        
        logger.info(f"PropertyProcessor '{self.name}' initialized at {self.data_path}")
    
    def _load_existing_data(self):
        """Load existing property datasets and metadata."""
        try:
            # Load property datasets
            datasets_file = self.data_dirs['metadata'] / 'property_datasets.json'
            if datasets_file.exists():
                with open(datasets_file, 'r') as f:
                    dataset_info = json.load(f)
                    for dataset_name in dataset_info.get('datasets', []):
                        self._load_property_dataset_internal(dataset_name)
            
            # Load entity properties cache
            entity_cache = self.data_dirs['cache'] / 'entity_properties.json'
            if entity_cache.exists():
                with open(entity_cache, 'r') as f:
                    self.entity_properties = json.load(f)
            
            logger.info(f"Loaded {len(self.property_datasets)} property datasets")
            
        except Exception as e:
            logger.warning(f"Failed to load existing property data: {e}")
    
    def assign_property(self, 
                       entity_identifier: str,
                       property_name: str,
                       property_value: Any,
                       dataset_name: str,
                       metadata: Optional[Dict[str, Any]] = None,
                       secondary_selection: Optional[Dict[str, Any]] = None,
                       overwrite: bool = True) -> str:
        """
        Assign a property to an entity.
        
        Args:
            entity_identifier: Entity identifier (ID, PDB code, sequence name, etc.)
            property_name: Name of the property
            property_value: Value of the property
            dataset_name: Dataset to assign this property to
            metadata: Optional metadata about this property assignment
            secondary_selection: Advanced selection criteria (e.g., GRN column, atom ID)
            overwrite: Whether to overwrite existing property values
            
        Returns:
            Entity ID that was assigned the property
        """
        # Resolve entity identifier to entity ID
        entity_id = self._resolve_entity_identifier(entity_identifier)
        
        # Initialize dataset if it doesn't exist
        if dataset_name not in self.property_datasets:
            self._initialize_dataset(dataset_name)
        
        # Check for existing assignment
        if not overwrite and self._has_property(entity_id, property_name, dataset_name):
            existing_value = self.get_entity_property(entity_id, property_name, dataset_name)
            logger.warning(f"Property {property_name} already exists for entity {entity_id} "
                          f"with value {existing_value}. Use overwrite=True to replace.")
            return entity_id
        
        # Assign the property
        if entity_id not in self.entity_properties:
            self.entity_properties[entity_id] = {}
        
        if dataset_name not in self.entity_properties[entity_id]:
            self.entity_properties[entity_id][dataset_name] = {}
        
        # Store property value with metadata
        property_entry = {
            'value': property_value,
            'assigned_at': datetime.now().isoformat(),
            'metadata': metadata or {}
        }
        
        # Add secondary selection if provided
        if secondary_selection:
            property_entry['secondary_selection'] = secondary_selection
            self._store_secondary_selection(entity_id, secondary_selection)
        
        self.entity_properties[entity_id][dataset_name][property_name] = property_entry
        
        # Update dataset
        self._update_dataset(dataset_name, entity_id, property_name, property_entry)
        
        # Update property metadata
        self._update_property_metadata(property_name, property_value, metadata)
        
        # Save changes
        self._save_entity_properties()
        self.save_property_dataset(dataset_name)
        
        logger.info(f"Assigned property '{property_name}' = {property_value} "
                   f"to entity {entity_id} in dataset '{dataset_name}'")
        
        return entity_id
    
    def assign_properties_batch(self,
                               property_assignments: List[Dict[str, Any]],
                               dataset_name: str,
                               overwrite: bool = True) -> List[str]:
        """
        Assign multiple properties in batch for efficiency.
        
        Args:
            property_assignments: List of assignment dictionaries with keys:
                - entity_identifier: Entity identifier
                - property_name: Property name
                - property_value: Property value
                - metadata: Optional metadata
                - secondary_selection: Optional secondary selection
            dataset_name: Dataset name for all assignments
            overwrite: Whether to overwrite existing values
            
        Returns:
            List of entity IDs that were assigned properties
        """
        entity_ids = []
        
        # Initialize dataset
        if dataset_name not in self.property_datasets:
            self._initialize_dataset(dataset_name)
        
        for assignment in property_assignments:
            try:
                entity_id = self.assign_property(
                    entity_identifier=assignment['entity_identifier'],
                    property_name=assignment['property_name'],
                    property_value=assignment['property_value'],
                    dataset_name=dataset_name,
                    metadata=assignment.get('metadata'),
                    secondary_selection=assignment.get('secondary_selection'),
                    overwrite=overwrite
                )
                entity_ids.append(entity_id)
            except Exception as e:
                logger.error(f"Failed to assign property for {assignment['entity_identifier']}: {e}")
        
        logger.info(f"Batch assigned {len(entity_ids)} properties to dataset '{dataset_name}'")
        return entity_ids
    
    def get_entity_property(self,
                           entity_identifier: str,
                           property_name: str,
                           dataset_name: Optional[str] = None) -> Any:
        """
        Get a specific property value for an entity.
        
        Args:
            entity_identifier: Entity identifier
            property_name: Property name
            dataset_name: Specific dataset (searches all if None)
            
        Returns:
            Property value or None if not found
        """
        entity_id = self._resolve_entity_identifier(entity_identifier)
        
        if entity_id not in self.entity_properties:
            return None
        
        entity_props = self.entity_properties[entity_id]
        
        if dataset_name:
            # Search specific dataset
            if dataset_name in entity_props and property_name in entity_props[dataset_name]:
                return entity_props[dataset_name][property_name]['value']
            return None
        else:
            # Search all datasets
            for ds_name, ds_props in entity_props.items():
                if property_name in ds_props:
                    return ds_props[property_name]['value']
            return None
    
    def get_entity_properties(self,
                             entity_identifier: str,
                             dataset_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Get all properties for an entity.
        
        Args:
            entity_identifier: Entity identifier
            dataset_name: Specific dataset (all datasets if None)
            
        Returns:
            Dictionary of property_name -> value
        """
        entity_id = self._resolve_entity_identifier(entity_identifier)
        
        if entity_id not in self.entity_properties:
            return {}
        
        entity_props = self.entity_properties[entity_id]
        
        if dataset_name:
            # Return properties from specific dataset
            if dataset_name in entity_props:
                return {prop: data['value'] for prop, data in entity_props[dataset_name].items()}
            return {}
        else:
            # Return properties from all datasets
            all_props = {}
            for ds_name, ds_props in entity_props.items():
                for prop, data in ds_props.items():
                    # Prefix with dataset name if property exists in multiple datasets
                    prop_key = f"{ds_name}.{prop}" if prop in all_props else prop
                    all_props[prop_key] = data['value']
            return all_props
    
    def get_dataset_properties(self, dataset_name: str) -> pd.DataFrame:
        """
        Get all properties in a dataset as a DataFrame.
        
        Args:
            dataset_name: Dataset name
            
        Returns:
            DataFrame with entity_id as index and properties as columns
        """
        if dataset_name not in self.property_datasets:
            logger.warning(f"Dataset '{dataset_name}' not found")
            return pd.DataFrame()
        
        dataset = self.property_datasets[dataset_name]
        
        # Convert to DataFrame format
        rows = []
        for entity_id, properties in dataset.get('entities', {}).items():
            row = {'entity_id': entity_id}
            for prop_name, prop_data in properties.items():
                row[prop_name] = prop_data['value']
            rows.append(row)
        
        if not rows:
            return pd.DataFrame()
        
        df = pd.DataFrame(rows)
        df.set_index('entity_id', inplace=True)
        return df
    
    def filter_entities_by_property(self,
                                   dataset_name: str,
                                   property_filters: Dict[str, Any]) -> List[str]:
        """
        Filter entities by property values.
        
        Args:
            dataset_name: Dataset to filter
            property_filters: Dictionary of property_name -> filter_condition
                            Supports: value, {'gt': val}, {'lt': val}, {'in': [vals]}
            
        Returns:
            List of entity IDs matching the filter
        """
        df = self.get_dataset_properties(dataset_name)
        if df.empty:
            return []
        
        # Apply filters
        mask = pd.Series([True] * len(df), index=df.index)
        
        for prop_name, condition in property_filters.items():
            if prop_name not in df.columns:
                logger.warning(f"Property '{prop_name}' not found in dataset '{dataset_name}'")
                continue
            
            if isinstance(condition, dict):
                if 'gt' in condition:
                    mask &= df[prop_name] > condition['gt']
                elif 'lt' in condition:
                    mask &= df[prop_name] < condition['lt']
                elif 'in' in condition:
                    mask &= df[prop_name].isin(condition['in'])
                elif 'eq' in condition:
                    mask &= df[prop_name] == condition['eq']
            else:
                # Direct value comparison
                mask &= df[prop_name] == condition
        
        return df[mask].index.tolist()
    
    def create_property_dataset_from_file(self,
                                         file_path: str,
                                         dataset_name: str,
                                         entity_column: str = 'entity_id',
                                         format_type: Optional[str] = None) -> int:
        """
        Create a property dataset from a CSV/JSON file.
        
        Args:
            file_path: Path to data file (CSV or JSON)
            dataset_name: Name for the new dataset
            entity_column: Column containing entity identifiers
            format_type: Format type for entity resolution (auto-detect if None)
            
        Returns:
            Number of entities with properties assigned
        """
        file_path = Path(file_path)
        
        # Load data based on file type
        if file_path.suffix.lower() == '.csv':
            data = pd.read_csv(file_path)
        elif file_path.suffix.lower() == '.json':
            with open(file_path, 'r') as f:
                json_data = json.load(f)
            if isinstance(json_data, list):
                data = pd.DataFrame(json_data)
            else:
                data = pd.DataFrame([json_data])
        else:
            raise ValueError(f"Unsupported file format: {file_path.suffix}")
        
        if entity_column not in data.columns:
            raise ValueError(f"Entity column '{entity_column}' not found in data")
        
        # Initialize dataset
        self._initialize_dataset(dataset_name, {
            'source_file': str(file_path),
            'entity_column': entity_column,
            'format_type': format_type
        })
        
        # Assign properties for each entity
        assigned_count = 0
        property_columns = [col for col in data.columns if col != entity_column]
        
        for _, row in data.iterrows():
            entity_identifier = row[entity_column]
            
            try:
                entity_id = self._resolve_entity_identifier(entity_identifier, format_type)
                
                for prop_col in property_columns:
                    prop_value = row[prop_col]
                    if pd.notna(prop_value):  # Skip NaN values
                        self.assign_property(
                            entity_identifier=entity_id,
                            property_name=prop_col,
                            property_value=prop_value,
                            dataset_name=dataset_name,
                            metadata={'source_file': str(file_path)},
                            overwrite=True
                        )
                
                assigned_count += 1
                
            except Exception as e:
                logger.warning(f"Failed to assign properties for entity {entity_identifier}: {e}")
        
        logger.info(f"Created dataset '{dataset_name}' with {assigned_count} entities "
                   f"and {len(property_columns)} properties from {file_path}")
        
        return assigned_count
    
    def save_property_dataset(self, dataset_name: str, file_format: str = 'json'):
        """
        Save a property dataset to file.
        
        Args:
            dataset_name: Dataset to save
            file_format: Format to save ('json', 'csv', or 'both')
        """
        if dataset_name not in self.property_datasets:
            logger.warning(f"Dataset '{dataset_name}' not found")
            return
        
        dataset_dir = self.data_dirs['datasets'] / dataset_name
        dataset_dir.mkdir(exist_ok=True)
        
        dataset = self.property_datasets[dataset_name]
        
        if file_format in ['json', 'both']:
            # Save as JSON
            json_file = dataset_dir / f"{dataset_name}.json"
            with open(json_file, 'w') as f:
                json.dump(dataset, f, indent=2, default=str)
            logger.info(f"Saved dataset '{dataset_name}' to {json_file}")
        
        if file_format in ['csv', 'both']:
            # Save as CSV
            df = self.get_dataset_properties(dataset_name)
            if not df.empty:
                csv_file = dataset_dir / f"{dataset_name}.csv"
                df.to_csv(csv_file)
                logger.info(f"Saved dataset '{dataset_name}' to {csv_file}")
    
    def load_property_dataset(self, dataset_name: str) -> bool:
        """
        Load a property dataset from file.
        
        Args:
            dataset_name: Dataset name to load
            
        Returns:
            True if successfully loaded, False otherwise
        """
        return self._load_property_dataset_internal(dataset_name)
    
    def list_property_datasets(self) -> List[str]:
        """
        List all available property datasets.
        
        Returns:
            List of dataset names
        """
        # Combine loaded datasets and available dataset files
        loaded_datasets = set(self.property_datasets.keys())
        
        # Check for dataset files
        if self.data_dirs['datasets'].exists():
            for dataset_dir in self.data_dirs['datasets'].iterdir():
                if dataset_dir.is_dir():
                    loaded_datasets.add(dataset_dir.name)
        
        return sorted(list(loaded_datasets))
    
    def get_dataset_statistics(self, dataset_name: str) -> Dict[str, Any]:
        """
        Get statistics for a property dataset.
        
        Args:
            dataset_name: Dataset name
            
        Returns:
            Dictionary with dataset statistics
        """
        if dataset_name not in self.property_datasets:
            logger.warning(f"Dataset '{dataset_name}' not found")
            return {}
        
        dataset = self.property_datasets[dataset_name]
        df = self.get_dataset_properties(dataset_name)
        
        stats = {
            'dataset_name': dataset_name,
            'entity_count': len(dataset.get('entities', {})),
            'property_count': len(df.columns) if not df.empty else 0,
            'created_at': dataset.get('metadata', {}).get('created_at'),
            'modified_at': dataset.get('metadata', {}).get('modified_at'),
            'description': dataset.get('metadata', {}).get('description', ''),
        }
        
        if not df.empty:
            # Property statistics
            stats['properties'] = {}
            for col in df.columns:
                col_stats = {
                    'type': str(df[col].dtype),
                    'non_null_count': df[col].count(),
                    'null_count': df[col].isnull().sum()
                }
                
                if df[col].dtype in ['int64', 'float64']:
                    col_stats.update({
                        'mean': df[col].mean(),
                        'std': df[col].std(),
                        'min': df[col].min(),
                        'max': df[col].max()
                    })
                elif df[col].dtype == 'object':
                    col_stats.update({
                        'unique_count': df[col].nunique(),
                        'most_common': df[col].mode().iloc[0] if len(df[col].mode()) > 0 else None
                    })
                
                stats['properties'][col] = col_stats
        
        return stats
    
    # Advanced secondary selection methods
    def assign_grn_property(self,
                           entity_identifier: str,
                           grn_position: str,
                           property_name: str,
                           property_value: Any,
                           dataset_name: str) -> str:
        """
        Assign a property to a specific GRN position for an entity.
        
        Args:
            entity_identifier: Entity identifier
            grn_position: GRN position (e.g., "3.50", "7.45")
            property_name: Property name
            property_value: Property value
            dataset_name: Dataset name
            
        Returns:
            Entity ID
        """
        secondary_selection = {
            'type': 'grn_position',
            'grn_position': grn_position
        }
        
        return self.assign_property(
            entity_identifier=entity_identifier,
            property_name=f"{grn_position}.{property_name}",
            property_value=property_value,
            dataset_name=dataset_name,
            secondary_selection=secondary_selection
        )
    
    def assign_atom_property(self,
                            entity_identifier: str,
                            atom_selector: Dict[str, Any],
                            property_name: str,
                            property_value: Any,
                            dataset_name: str) -> str:
        """
        Assign a property to specific atoms in a structure.
        
        Args:
            entity_identifier: Entity identifier
            atom_selector: Atom selection criteria (e.g., {'chain': 'A', 'residue': 100})
            property_name: Property name
            property_value: Property value
            dataset_name: Dataset name
            
        Returns:
            Entity ID
        """
        secondary_selection = {
            'type': 'atom_selection',
            'atom_selector': atom_selector
        }
        
        selector_str = "_".join([f"{k}{v}" for k, v in atom_selector.items()])
        
        return self.assign_property(
            entity_identifier=entity_identifier,
            property_name=f"{selector_str}.{property_name}",
            property_value=property_value,
            dataset_name=dataset_name,
            secondary_selection=secondary_selection
        )
    
    # Internal helper methods
    def _resolve_entity_identifier(self, 
                                  identifier: str, 
                                  format_type: Optional[str] = None) -> str:
        """
        Resolve an entity identifier to an entity ID.
        
        Args:
            identifier: Entity identifier (could be entity_id, PDB code, etc.)
            format_type: Format type hint for resolution
            
        Returns:
            Resolved entity ID
        """
        # If already a valid entity ID (10 characters), return as-is
        if len(identifier) == 10 and identifier.isalnum():
            return identifier
        
        # Try to resolve through entity registry
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            
            if format_type:
                resolved_id = global_registry.entity_registry.resolve_identifier(
                    identifier, format_type=format_type
                )
                if resolved_id:
                    return resolved_id
            
            # Try all format types
            for fmt in ['structure', 'sequence', 'grn', 'embedding']:
                resolved_id = global_registry.entity_registry.resolve_identifier(
                    identifier, format_type=fmt
                )
                if resolved_id:
                    return resolved_id
        
        except Exception as e:
            logger.debug(f"Entity registry resolution failed: {e}")
        
        # Fall back to generating entity ID from identifier
        entity_id = generate_entity_id(identifier)
        logger.debug(f"Generated entity ID {entity_id} for identifier {identifier}")
        return entity_id
    
    def _initialize_dataset(self, dataset_name: str, metadata: Optional[Dict[str, Any]] = None):
        """Initialize a new property dataset."""
        base_metadata = {
            'created_at': datetime.now().isoformat(),
            'modified_at': datetime.now().isoformat(),
            'description': metadata.get('description', '') if metadata else ''
        }
        
        # Add additional metadata if provided
        if metadata:
            base_metadata.update(metadata)
            
        self.property_datasets[dataset_name] = {
            'metadata': base_metadata,
            'entities': {}
        }
        
        self.dataset_metadata[dataset_name] = self.property_datasets[dataset_name]['metadata']
        logger.info(f"Initialized property dataset '{dataset_name}'")
    
    def _update_dataset(self, dataset_name: str, entity_id: str, 
                       property_name: str, property_entry: Dict[str, Any]):
        """Update dataset with new property assignment."""
        if entity_id not in self.property_datasets[dataset_name]['entities']:
            self.property_datasets[dataset_name]['entities'][entity_id] = {}
        
        self.property_datasets[dataset_name]['entities'][entity_id][property_name] = property_entry
        self.property_datasets[dataset_name]['metadata']['modified_at'] = datetime.now().isoformat()
    
    def _has_property(self, entity_id: str, property_name: str, dataset_name: str) -> bool:
        """Check if entity has a specific property in a dataset."""
        return (entity_id in self.entity_properties and
                dataset_name in self.entity_properties[entity_id] and
                property_name in self.entity_properties[entity_id][dataset_name])
    
    def _update_property_metadata(self, property_name: str, property_value: Any, 
                                 metadata: Optional[Dict[str, Any]]):
        """Update metadata for a property."""
        if property_name not in self.property_metadata:
            self.property_metadata[property_name] = {
                'type': type(property_value).__name__,
                'first_seen': datetime.now().isoformat(),
                'usage_count': 0,
                'metadata': metadata or {}
            }
        
        self.property_metadata[property_name]['usage_count'] += 1
        self.property_metadata[property_name]['last_used'] = datetime.now().isoformat()
    
    def _store_secondary_selection(self, entity_id: str, selection: Dict[str, Any]):
        """Store secondary selection criteria for an entity."""
        if entity_id not in self.secondary_selections:
            self.secondary_selections[entity_id] = []
        
        self.secondary_selections[entity_id].append(selection)
    
    def _save_entity_properties(self):
        """Save entity properties cache to file."""
        cache_file = self.data_dirs['cache'] / 'entity_properties.json'
        with open(cache_file, 'w') as f:
            json.dump(self.entity_properties, f, indent=2, default=str)
    
    def _save_datasets_metadata(self):
        """Save datasets metadata."""
        metadata_file = self.data_dirs['metadata'] / 'property_datasets.json'
        metadata = {
            'datasets': list(self.property_datasets.keys()),
            'property_metadata': self.property_metadata,
            'updated_at': datetime.now().isoformat()
        }
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
    
    def _load_property_dataset_internal(self, dataset_name: str) -> bool:
        """Internal method to load a property dataset."""
        dataset_dir = self.data_dirs['datasets'] / dataset_name
        json_file = dataset_dir / f"{dataset_name}.json"
        
        if not json_file.exists():
            logger.warning(f"Dataset file {json_file} not found")
            return False
        
        try:
            with open(json_file, 'r') as f:
                dataset = json.load(f)
            
            self.property_datasets[dataset_name] = dataset
            
            # Rebuild entity_properties cache for this dataset
            for entity_id, properties in dataset.get('entities', {}).items():
                if entity_id not in self.entity_properties:
                    self.entity_properties[entity_id] = {}
                self.entity_properties[entity_id][dataset_name] = properties
            
            logger.info(f"Loaded property dataset '{dataset_name}' with "
                       f"{len(dataset.get('entities', {}))} entities")
            return True
            
        except Exception as e:
            logger.error(f"Failed to load dataset '{dataset_name}': {e}")
            return False
    
    # Entity integration methods
    def list_entities(self, dataset_name: Optional[str] = None) -> List[str]:
        """
        List entities that have properties assigned.
        
        Args:
            dataset_name: Specific dataset (all datasets if None)
            
        Returns:
            List of entity IDs
        """
        if dataset_name:
            if dataset_name in self.property_datasets:
                return list(self.property_datasets[dataset_name].get('entities', {}).keys())
            return []
        else:
            # All entities across all datasets
            all_entities = set()
            for dataset in self.property_datasets.values():
                all_entities.update(dataset.get('entities', {}).keys())
            return sorted(list(all_entities))
    
    def list_datasets(self) -> List[Dict[str, Any]]:
        """
        List property datasets with metadata.
        
        Returns:
            List of dataset information dictionaries
        """
        datasets = []
        for dataset_name in self.property_datasets:
            stats = self.get_dataset_statistics(dataset_name)
            datasets.append({
                'id': dataset_name,
                'entity_count': stats.get('entity_count', 0),
                'property_count': stats.get('property_count', 0),
                'description': stats.get('description', ''),
                'created_at': stats.get('created_at'),
                'type': 'property_dataset'
            })
        return datasets
    
    def get_entity_formats(self, entity_id: str) -> List[str]:
        """
        Get available formats for an entity through entity registry.
        
        Args:
            entity_id: Entity ID
            
        Returns:
            List of available format types
        """
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            
            entity_data = global_registry.entity_registry.get_entity(entity_id)
            if entity_data:
                return list(entity_data.get('formats', {}).keys())
            return []
            
        except Exception as e:
            logger.warning(f"Failed to get entity formats for {entity_id}: {e}")
            return []
    
    def __repr__(self) -> str:
        """String representation of the PropertyProcessor."""
        total_entities = len(self.entity_properties)
        total_datasets = len(self.property_datasets)
        
        info = [
            f"PropertyProcessor(name='{self.name}')",
            f"  Data path: {self.data_path}",
            f"  Datasets: {total_datasets}",
            f"  Entities with properties: {total_entities}"
        ]
        
        if self.property_datasets:
            info.append("  Property datasets:")
            for dataset_name in sorted(self.property_datasets.keys()):
                stats = self.get_dataset_statistics(dataset_name)
                info.append(f"    - {dataset_name}: {stats.get('entity_count', 0)} entities, "
                           f"{stats.get('property_count', 0)} properties")
        
        return "\n".join(info)