"""
Dataset management for the Protos framework.

This module provides functionality for working with standardized datasets
across different processor types. It integrates with the BaseProcessor class
to provide consistent dataset handling.
"""

import os
import json
import logging
from typing import Dict, List, Any, Optional, Union, Set
from datetime import datetime
from pathlib import Path

from protos.io.data_access import Dataset, DataRegistry, GlobalRegistry
from protos.io.paths import ProtosPaths, DataSource, ensure_directory

# Configure logger
logger = logging.getLogger(__name__)

class DatasetManager:
    """
    Manager for dataset operations across processors.
    
    This class provides methods for creating, loading, saving, and managing
    datasets across different processor types. It integrates with the path
    resolution system and registry management.
    """
    
    def __init__(self, processor_type: str, paths: Optional[ProtosPaths] = None):
        """
        Initialize the dataset manager.
        
        Args:
            processor_type: Type of processor ('structure', 'grn', etc.)
            paths: ProtosPaths instance for path resolution
        """
        self.processor_type = processor_type
        self.paths = paths or ProtosPaths()
        self.global_registry = GlobalRegistry(self.paths)
        
        # Initialize processor registry
        registry_path = self.paths.get_registry_path(processor_type)
        self.registry = DataRegistry(registry_path)
        
        # Dataset directory path
        self.dataset_dir = self._get_dataset_dir()
        print(f"Dataset Manager initialized at {self.dataset_dir}")
        ensure_directory(self.dataset_dir)
    
    def _get_dataset_dir(self) -> str:
        """
        Get the appropriate dataset directory for the processor type.
        
        Returns:
            Path to dataset directory
        """
        try:
            # Try processor-specific dataset directory first
            if self.processor_type == 'structure':
                return self.paths.get_structure_subdir_path('dataset_dir')
            elif self.processor_type == 'grn':
                return self.paths.get_grn_subdir_path('grn_dir')
            elif self.processor_type == 'sequence':
                return self.paths.get_sequence_subdir_path('metadata_dir')
            else:
                # Default to datasets subdirectory in processor directory
                processor_path = self.paths.get_processor_path(self.processor_type)
                return os.path.join(processor_path, 'datasets')
        except (ValueError, KeyError):
            # Fall back to basic dataset directory
            processor_path = self.paths.get_processor_path(self.processor_type)
            return os.path.join(processor_path, 'datasets')
    
    def create_dataset(self, 
                     dataset_id: str,
                     name: str,
                     description: str,
                     content: Union[List, Dict, Set],
                     metadata: Optional[Dict[str, Any]] = None) -> Dataset:
        """
        Create a new dataset.
        
        Args:
            dataset_id: Unique identifier for the dataset
            name: Human-readable name
            description: Detailed description
            content: Dataset content (list of IDs, dictionary, etc.)
            metadata: Additional metadata
            
        Returns:
            Created Dataset instance
        """
        # Create dataset object
        dataset = Dataset(
            id=dataset_id,
            name=name,
            description=description,
            type=self.processor_type,
            content=content,
            metadata=metadata or {}
        )
        
        # Save dataset to file
        file_path = os.path.join(self.dataset_dir, f"{dataset_id}.json")
        dataset.save(file_path)
        
        # Register in both local and global registries
        self.registry.register_dataset(
            dataset_id=dataset_id,
            file_path=file_path,
            metadata={
                "type": self.processor_type,
                "name": name,
                "description": description,
                "items": len(content),
                "created_at": dataset.creation_date
            }
        )
        
        self.global_registry.register_dataset(
            dataset_id=dataset_id, 
            file_path=file_path, 
            processor_type=self.processor_type,
            dataset_type=self.processor_type,
            metadata={
                "name": name,
                "description": description,
                "items": len(content),
                "created_at": dataset.creation_date
            }
        )
        
        logger.info(f"Created dataset '{dataset_id}' with {len(content)} items")
        return dataset
    
    def load_dataset(self, dataset_id: str) -> Optional[Dataset]:
        """
        Load a dataset by ID.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            Dataset instance or None if not found
        """
        # Try to get dataset path from registry
        file_path = self.global_registry.get_dataset_path(dataset_id)
        
        if not file_path:
            # Try local registry as fallback
            file_path = self.registry.get_dataset_path(dataset_id)
            
        if not file_path:
            # Try direct path construction if not in registry
            direct_path = os.path.join(self.dataset_dir, f"{dataset_id}.json")
            if os.path.exists(direct_path):
                file_path = direct_path
        
        if not file_path or not os.path.exists(file_path):
            logger.warning(f"Dataset '{dataset_id}' not found")
            return None
        
        # Load dataset from file
        try:
            dataset = Dataset.load(file_path)
            logger.info(f"Loaded dataset '{dataset_id}' with {len(dataset.content)} items")
            return dataset
        except Exception as e:
            logger.error(f"Error loading dataset '{dataset_id}': {e}")
            return None
    
    def save_dataset(self, dataset: Dataset) -> bool:
        """
        Save a dataset to file and update registry.
        
        Args:
            dataset: Dataset to save
            
        Returns:
            True if successful
        """
        # Ensure correct processor type
        if dataset.type != self.processor_type:
            logger.warning(
                f"Dataset '{dataset.id}' has type '{dataset.type}', expected '{self.processor_type}'"
            )
        
        # Save dataset to file
        file_path = os.path.join(self.dataset_dir, f"{dataset.id}.json")
        try:
            dataset.save(file_path)
            
            # Update registries
            self.registry.register_dataset(
                dataset_id=dataset.id,
                file_path=file_path,
                metadata={
                    "type": dataset.type,
                    "name": dataset.name,
                    "description": dataset.description,
                    "items": len(dataset.content),
                    "updated_at": dataset.last_modified
                }
            )
            
            self.global_registry.register_dataset(
                dataset_id=dataset.id, 
                file_path=file_path, 
                processor_type=dataset.type,
                dataset_type=dataset.type,
                metadata={
                    "name": dataset.name,
                    "description": dataset.description,
                    "items": len(dataset.content),
                    "updated_at": dataset.last_modified
                }
            )
            
            logger.info(f"Saved dataset '{dataset.id}' with {len(dataset.content)} items")
            return True
        except Exception as e:
            logger.error(f"Error saving dataset '{dataset.id}': {e}")
            return False
    
    def list_datasets(self) -> List[Dict[str, Any]]:
        """
        List available datasets for this processor type.
        
        Returns:
            List of dataset information dictionaries
        """
        # Get datasets from registry
        return self.registry.list_datasets()
    
    def delete_dataset(self, dataset_id: str) -> bool:
        """
        Delete a dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            True if deletion was successful
        """
        # Get file path from registry
        file_path = self.registry.get_dataset_path(dataset_id)
        
        if not file_path:
            # Try direct path construction
            file_path = os.path.join(self.dataset_dir, f"{dataset_id}.json")
        
        # Delete file if it exists
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                logger.info(f"Deleted dataset file: {file_path}")
            except Exception as e:
                logger.error(f"Error deleting dataset file: {e}")
                return False
        
        # Remove from registries
        self.registry.remove_dataset(dataset_id)
        self.global_registry.remove_dataset(dataset_id)
        
        return True
    
    def is_dataset_available(self, dataset_id: str) -> bool:
        """
        Check if a dataset is available.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            True if the dataset exists and is available
        """
        # Check registry first
        file_path = self.registry.get_dataset_path(dataset_id)
        
        if file_path and os.path.exists(file_path):
            return True
        
        # Try direct path
        direct_path = os.path.join(self.dataset_dir, f"{dataset_id}.json")
        return os.path.exists(direct_path)
    
    def get_dataset_info(self, dataset_id: str) -> Optional[Dict[str, Any]]:
        """
        Get information about a specific dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            Dataset information or None if not found
        """
        # Try local registry first, then global
        info = self.registry.get_dataset_metadata(dataset_id)
        if info:
            return {"id": dataset_id, **info}
        
        info = self.global_registry.get_dataset_metadata(dataset_id)
        if info:
            return {"id": dataset_id, **info}
        
        return None
    
    def update_dataset_metadata(self, dataset_id: str, metadata: Dict[str, Any]) -> bool:
        """
        Update metadata for a dataset.
        
        Args:
            dataset_id: Dataset identifier
            metadata: New metadata to merge with existing
            
        Returns:
            True if successful
        """
        # Load dataset
        dataset = self.load_dataset(dataset_id)
        if not dataset:
            return False
        
        # Update metadata and save
        dataset.update_metadata(metadata)
        return self.save_dataset(dataset)
    
    def merge_datasets(self, 
                      dataset_ids: List[str], 
                      new_dataset_id: str,
                      name: Optional[str] = None,
                      description: Optional[str] = None) -> Optional[Dataset]:
        """
        Merge multiple datasets into a new dataset.
        
        Args:
            dataset_ids: List of dataset IDs to merge
            new_dataset_id: ID for the merged dataset
            name: Optional name for the merged dataset
            description: Optional description for the merged dataset
            
        Returns:
            Merged dataset or None if error
        """
        # Load all datasets
        datasets = []
        for dataset_id in dataset_ids:
            dataset = self.load_dataset(dataset_id)
            if dataset:
                datasets.append(dataset)
            else:
                logger.warning(f"Dataset '{dataset_id}' not found, skipping")
        
        if not datasets:
            logger.error("No valid datasets to merge")
            return None
        
        # Determine content type from first dataset
        first_content = datasets[0].content
        content_type = type(first_content)
        
        # Prepare merged content
        if isinstance(first_content, list):
            merged_content = []
            for dataset in datasets:
                merged_content.extend(item for item in dataset.content if item not in merged_content)
        elif isinstance(first_content, dict):
            merged_content = {}
            for dataset in datasets:
                merged_content.update(dataset.content)
        elif isinstance(first_content, set):
            merged_content = set()
            for dataset in datasets:
                merged_content.update(dataset.content)
        else:
            logger.error(f"Unsupported content type: {content_type}")
            return None
        
        # Create new dataset
        name = name or f"Merged dataset ({', '.join(dataset_ids)})"
        description = description or f"Merged from datasets: {', '.join(dataset_ids)}"
        
        merged_dataset = self.create_dataset(
            dataset_id=new_dataset_id,
            name=name,
            description=description,
            content=merged_content,
            metadata={
                "source_datasets": dataset_ids,
                "merged_at": datetime.now().isoformat()
            }
        )
        
        return merged_dataset
    
    def convert_legacy_dataset(self, legacy_dataset_id: str, new_dataset_id: Optional[str] = None) -> Optional[Dataset]:
        """
        Convert a legacy dataset format to the new standardized format.
        
        Args:
            legacy_dataset_id: Legacy dataset identifier
            new_dataset_id: Optional new ID (defaults to legacy_dataset_id)
            
        Returns:
            Converted Dataset or None if conversion failed
        """
        # Try to load the legacy dataset using BaseProcessor methods
        from protos.core.base_processor import BaseProcessor
        
        # Create a temporary processor instance
        temp_processor = BaseProcessor(
            name="dataset_converter",
            processor_data_dir=self.processor_type
        )
        
        # Try to load legacy dataset
        try:
            legacy_data = temp_processor.load_data(legacy_dataset_id)
        except FileNotFoundError:
            logger.error(f"Legacy dataset '{legacy_dataset_id}' not found")
            return None
        
        # Determine content type and convert as needed
        content = legacy_data
        
        # For structure datasets, extract PDB IDs
        if self.processor_type == 'structure' and hasattr(legacy_data, 'get'):
            try:
                # Handle different legacy formats
                if isinstance(legacy_data, dict) and 'pdb_ids' in legacy_data:
                    content = legacy_data['pdb_ids']
                elif isinstance(legacy_data, dict) and all(isinstance(k, str) for k in legacy_data.keys()):
                    # Format where keys are dataset names and values are lists of IDs
                    if legacy_dataset_id in legacy_data:
                        content = legacy_data[legacy_dataset_id]
                    else:
                        # Use all values as content
                        content = []
                        for ids in legacy_data.values():
                            if isinstance(ids, list):
                                content.extend(ids)
            except Exception as e:
                logger.warning(f"Error processing structure dataset content: {e}")
        
        # Create new dataset
        new_id = new_dataset_id or legacy_dataset_id
        
        dataset = self.create_dataset(
            dataset_id=new_id,
            name=f"{legacy_dataset_id} (converted)",
            description=f"Converted from legacy dataset {legacy_dataset_id}",
            content=content,
            metadata={
                "original_dataset_id": legacy_dataset_id,
                "converted_at": datetime.now().isoformat(),
                "conversion_source": "legacy_format"
            }
        )
        
        return dataset