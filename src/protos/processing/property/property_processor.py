#!/usr/bin/env python3
"""
Fixed PropertyProcessor that treats properties as tables/DataFrames.

This version:
- Works with entire property tables (DataFrames) at once
- Each row is an entity
- Each column is a property
- Saves the entire table efficiently without per-cell logging
"""

import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union, Any
import json

from protos.core.base_processor import BaseProcessor

logger = logging.getLogger(__name__)


class PropertyProcessor(BaseProcessor):
    """
    Processor for managing properties as tables.
    
    Properties are stored as DataFrames where:
    - Each row represents an entity (indexed by entity ID)
    - Each column represents a property
    - Each table is a dataset
    """
    
    def __init__(self, name: str = "property_processor", paths=None):
        """Initialize PropertyProcessor."""
        super().__init__(name=name, paths=paths)
        self.processor_type = "property"
        
        # Store property tables as DataFrames
        self.property_tables: Dict[str, pd.DataFrame] = {}
        
        # Ensure directories exist
        self._ensure_directories()
    
    def _ensure_directories(self):
        """Ensure required directories exist."""
        # ProtosPaths handles directory creation automatically
        # We just need to use the proper subdirectory paths
        pass
    
    def create_property_table(self, 
                            dataset_name: str,
                            data: Union[pd.DataFrame, Dict[str, Dict[str, Any]]],
                            metadata: Optional[Dict[str, Any]] = None) -> pd.DataFrame:
        """
        Create a new property table.
        
        Args:
            dataset_name: Name for the dataset/table
            data: Either a DataFrame or dict of {entity_id: {prop: value}}
            metadata: Optional metadata for the dataset
            
        Returns:
            The created property table as DataFrame
        """
        # Convert to DataFrame if needed
        if isinstance(data, dict):
            df = pd.DataFrame.from_dict(data, orient='index')
        else:
            df = data.copy()
        
        # Ensure index is named properly
        if df.index.name is None:
            df.index.name = 'entity_id'
        
        # Store the table
        self.property_tables[dataset_name] = df
        
        # Save to disk
        self.save_property_table(dataset_name, metadata=metadata)
        
        logger.info(f"Created property table '{dataset_name}' with {len(df)} entities and {len(df.columns)} properties")
        
        return df
    
    def save_property_table(self, dataset_name: str, metadata: Optional[Dict[str, Any]] = None):
        """
        Save a property table to disk.
        
        Args:
            dataset_name: Name of the dataset to save
            metadata: Optional metadata to save with the dataset
        """
        if dataset_name not in self.property_tables:
            logger.warning(f"Dataset '{dataset_name}' not found")
            return
        
        df = self.property_tables[dataset_name]
        
        # Save CSV table
        tables_dir = self.get_subdirectory_path('tables_dir')
        csv_path = tables_dir / f"{dataset_name}.csv"
        df.to_csv(csv_path)
        
        # Save metadata
        if metadata is not None:
            datasets_dir = self.get_subdirectory_path('datasets_dir')
            metadata_path = datasets_dir / f"{dataset_name}.json"
            dataset_info = {
                'name': dataset_name,
                'metadata': metadata,
                'shape': list(df.shape),
                'properties': df.columns.tolist(),
                'entities': df.index.tolist(),
                'table_file': f"../tables/{dataset_name}.csv"
            }
            
            with open(metadata_path, 'w') as f:
                json.dump(dataset_info, f, indent=2)
        
        logger.info(f"Saved property table '{dataset_name}' to {csv_path}")
    
    def load_property_table(self, dataset_name: str) -> pd.DataFrame:
        """
        Load a property table from disk.
        
        Args:
            dataset_name: Name of the dataset to load
            
        Returns:
            The loaded property table as DataFrame
        """
        # Check if already loaded
        if dataset_name in self.property_tables:
            return self.property_tables[dataset_name]
        
        # Load from CSV
        tables_dir = self.get_subdirectory_path('tables_dir')
        csv_path = tables_dir / f"{dataset_name}.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"Property table '{dataset_name}' not found at {csv_path}")
        
        df = pd.read_csv(csv_path, index_col=0)
        self.property_tables[dataset_name] = df
        
        logger.info(f"Loaded property table '{dataset_name}' with {len(df)} entities and {len(df.columns)} properties")
        
        return df
    
    def add_property_column(self, dataset_name: str, property_name: str, 
                           values: Union[pd.Series, Dict[str, Any], Any]):
        """
        Add a new property column to an existing table.
        
        Args:
            dataset_name: Dataset to update
            property_name: Name of the new property
            values: Values for the property (Series, dict, or scalar)
        """
        if dataset_name not in self.property_tables:
            raise ValueError(f"Dataset '{dataset_name}' not found")
        
        df = self.property_tables[dataset_name]
        
        # Add the new column
        if isinstance(values, pd.Series):
            df[property_name] = values
        elif isinstance(values, dict):
            df[property_name] = pd.Series(values)
        else:
            # Scalar value - apply to all entities
            df[property_name] = values
        
        # Save updated table
        self.save_property_table(dataset_name)
        
        logger.info(f"Added property '{property_name}' to dataset '{dataset_name}'")
    
    def get_property_table(self, dataset_name: str) -> pd.DataFrame:
        """Get a property table."""
        if dataset_name not in self.property_tables:
            return self.load_property_table(dataset_name)
        return self.property_tables[dataset_name]
    
    def filter_by_property(self, dataset_name: str, property_name: str, 
                          condition: callable) -> pd.DataFrame:
        """
        Filter entities by property value.
        
        Args:
            dataset_name: Dataset to filter
            property_name: Property to filter by
            condition: Function that takes a value and returns True/False
            
        Returns:
            Filtered DataFrame
        """
        df = self.get_property_table(dataset_name)
        
        if property_name not in df.columns:
            raise ValueError(f"Property '{property_name}' not found in dataset")
        
        mask = df[property_name].apply(condition)
        return df[mask]
    
    def merge_property_tables(self, dataset_names: List[str], 
                             output_name: str, 
                             how: str = 'outer') -> pd.DataFrame:
        """
        Merge multiple property tables.
        
        Args:
            dataset_names: List of datasets to merge
            output_name: Name for the merged dataset
            how: Merge method ('outer', 'inner', 'left', 'right')
            
        Returns:
            Merged DataFrame
        """
        if not dataset_names:
            raise ValueError("No datasets provided to merge")
        
        # Load all tables
        tables = [self.get_property_table(name) for name in dataset_names]
        
        # Merge progressively
        merged = tables[0]
        for table in tables[1:]:
            merged = merged.join(table, how=how, rsuffix='_dup')
        
        # Remove duplicate columns
        merged = merged.loc[:, ~merged.columns.duplicated()]
        
        # Save as new dataset
        self.property_tables[output_name] = merged
        self.save_property_table(output_name)
        
        logger.info(f"Merged {len(dataset_names)} datasets into '{output_name}' "
                   f"with {len(merged)} entities and {len(merged.columns)} properties")
        
        return merged
    
    def get_entity_properties(self, entity_id: str, dataset_name: Optional[str] = None) -> Dict[str, Any]:
        """
        Get all properties for an entity.
        
        Args:
            entity_id: Entity ID
            dataset_name: Specific dataset to query (or all if None)
            
        Returns:
            Dictionary of property names to values
        """
        if dataset_name:
            df = self.get_property_table(dataset_name)
            if entity_id in df.index:
                return df.loc[entity_id].to_dict()
            return {}
        else:
            # Search all datasets
            all_props = {}
            for name in self.list_datasets():
                df = self.get_property_table(name)
                if entity_id in df.index:
                    all_props.update(df.loc[entity_id].to_dict())
            return all_props
    
    def list_datasets(self) -> List[str]:
        """List all available property datasets."""
        # From memory
        datasets = list(self.property_tables.keys())
        
        # From disk
        tables_dir = self.get_subdirectory_path('tables_dir')
        if tables_dir.exists():
            for csv_file in tables_dir.glob("*.csv"):
                dataset_name = csv_file.stem
                if dataset_name not in datasets:
                    datasets.append(dataset_name)
        
        return sorted(datasets)
    
    # Implement abstract methods from BaseProcessor
    def save_entity(self, name: str, data: Any, metadata: Optional[dict] = None):
        """Save entity data."""
        # For PropertyProcessor, entities are rows in tables
        logger.warning("PropertyProcessor manages entities as rows in tables. "
                      "Use create_property_table or add_property_column instead.")
    
    def load_entity(self, name: str) -> Optional[Any]:
        """Load entity data."""
        # For PropertyProcessor, entities are rows in tables
        logger.warning("PropertyProcessor manages entities as rows in tables. "
                      "Use get_property_table to access entity properties.")
        return None
    
    # Compatibility methods for old API
    def assign_property(self, entity_identifier: str, property_name: str, 
                       property_value: Any, dataset_name: str = "default", **kwargs):
        """
        Assign a single property value (compatibility method).
        
        This method provides compatibility with the old API but uses
        the new table-based approach internally.
        """
        # Get or create the dataset
        if dataset_name in self.property_tables:
            df = self.property_tables[dataset_name]
        else:
            # Create empty DataFrame with proper structure
            df = pd.DataFrame(columns=[], index=pd.Index([], name='entity_id'))
            self.property_tables[dataset_name] = df
        
        # Ensure the property column exists
        if property_name not in df.columns:
            df[property_name] = pd.Series(dtype='object')
        
        # Add the entity if it doesn't exist
        if entity_identifier not in df.index:
            # Add new row with NaN values
            df.loc[entity_identifier] = pd.NA
        
        # Set the property value
        df.loc[entity_identifier, property_name] = property_value
        
        # Save
        self.save_property_table(dataset_name)
        
        return entity_identifier
    
    def load_property_dataset(self, dataset_name: str) -> pd.DataFrame:
        """Load a property dataset (alias for load_property_table)."""
        return self.load_property_table(dataset_name)