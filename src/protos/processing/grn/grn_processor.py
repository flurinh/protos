"""
GRN Processor with BaseProcessor integration.

This module provides a GRNProcessor class that extends BaseProcessor
to provide standardized data management capabilities for GRN tables.
"""

# Import all GRN utilities from consolidated module
from protos.processing.grn.grn_utils import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string,
    sort_grns,
    get_grn_interval,
    get_seq,
    get_annot_seq,
    map_grn_to_color,
    GRNConfigManager,
    sort_grns_str
)
import os
import re
import json
import logging
import numpy as np
import pandas as pd
import plotly.graph_objs as go
from pathlib import Path
from datetime import datetime
from protos.core.base_processor import BaseProcessor
from protos.io.data_access import generate_entity_id
from typing import Dict, List, Optional, Any, Union


class GRNProcessor(BaseProcessor):
    """
    Processor for Generic Residue Numbering (GRN) data.
    
    Handles loading, saving, and processing of GRN tables which map
    sequence positions to standardized numbering schemes.
    
    Key changes from legacy version:
    - NO path parameters in constructor
    - ALL paths managed by ProtosPaths
    - Implements abstract methods from BaseProcessor
    - Uses human-readable names for all operations
    """
    
    def __init__(self,
                 name="grn_processor",
                 paths=None,
                 dataset=None,
                 preload=True):
        """
        Initialize the GRN processor.
        
        Args:
            name: Processor instance name
            paths: ProtosPaths instance (created if not provided)
            dataset: Dataset ID to load (or list of datasets to merge)
            preload: Whether to load the dataset on initialization
        """
        # Initialize BaseProcessor with ProtosPaths
        super().__init__(name=name, paths=paths)
        
        # Initialize GRN-specific attributes
        self.dataset = None
        self.ids = []
        self.grns = []
        self.features = pd.DataFrame()
        self.maps = {}
        self.map = pd.DataFrame()
        
        # Always use dot notation (x notation is deprecated)
        self.using_dot_notation = True
        self.dot_to_x = {}  # Maps dot notation (3.50) to x notation (3x50) - DEPRECATED
        self.x_to_dot = {}  # Maps x notation (3x50) to dot notation (3.50) - DEPRECATED
        
        # Track entity IDs for GRN tables
        self._grn_entity_map = {}  # Maps table names to entity IDs
        
        # Set dataset and preload if specified
        if dataset is not None:
            if isinstance(dataset, list):
                self.dataset = 'merged_' + '_'.join(dataset)
                if preload:
                    self.load_and_merge_grn_tables(datasets=dataset)
            else:
                self.dataset = dataset
                if preload:
                    self.load_grn_table(dataset_id=dataset)
                    self.map = pd.DataFrame(columns=self.grns)
        else:
            self.dataset = None
            self.data = pd.DataFrame()
    
    # Path properties using ProtosPaths
    
    @property
    def path_grn_dir(self):
        """Get path to GRN tables directory."""
        # Use ProtosPaths to get the tables directory
        return Path(self.paths.get_subdir_path('grn', 'table_dir'))
    
    @property
    def path_ref_dir(self):
        """Get path to reference GRN tables directory."""
        # Use ProtosPaths to get the reference directory
        return Path(self.paths.get_subdir_path('grn', 'ref_dir'))
    
    @property
    def path_config_dir(self):
        """Get path to GRN configs directory."""
        # Use ProtosPaths to get the configs directory
        return Path(self.paths.get_subdir_path('grn', 'configs_dir'))
    
    @property
    def path_temp_dir(self):
        """Get path to GRN temp directory."""
        # Use ProtosPaths to get the temp directory
        return Path(self.paths.get_subdir_path('grn', 'temp_dir'))
    
    def load_reference_table(self, ref_name: str) -> pd.DataFrame:
        """
        Load a reference GRN table from the reference/ directory.
        
        Reference tables are used for GRN assignment and contain
        curated alignments for specific protein families.
        
        Args:
            ref_name: Name of reference table (without .csv extension)
            
        Returns:
            DataFrame with reference GRN alignment
        """
        ref_path = self.path_ref_dir / f"{ref_name}.csv"
        if not ref_path.exists():
            raise FileNotFoundError(f"Reference table not found: {ref_path}")
        
        self.logger.info(f"Loading reference table from {ref_path}")
        df = pd.read_csv(ref_path, index_col=0)
        # Normalize column names to standard format (e.g., '1.5' -> '1.50')
        df.columns = [normalize_grn_format(str(col)) for col in df.columns.tolist()]
        return df
    
    def save_reference_table(self, ref_name: str, reference_data: pd.DataFrame, 
                           metadata: Optional[dict] = None):
        """
        Save a reference GRN table to the reference/ directory.
        
        Reference tables are typically curated alignments used for
        GRN assignment to new sequences.
        
        Args:
            ref_name: Name for the reference table (without .csv extension)
            reference_data: DataFrame with reference GRN alignment
            metadata: Optional metadata about the reference
        """
        # Ensure reference directory exists
        ref_path = self.path_ref_dir / f"{ref_name}.csv"
        
        # Save the reference table
        reference_data.to_csv(ref_path, index=True, na_rep='-')
        
        # Save metadata if provided
        if metadata:
            meta_path = self.path_ref_dir / f"{ref_name}_metadata.json"
            metadata['created'] = datetime.now().isoformat()
            metadata['reference_name'] = ref_name
            with open(meta_path, 'w') as f:
                json.dump(metadata, f, indent=2)
                
        self.logger.info(f"Saved reference table to {ref_path}")
    
    def load_grn_config(self, config_name: str) -> dict:
        """
        Load a GRN configuration from the configs/ directory.
        
        Configurations define parameters for GRN assignment such as:
        - Motif patterns for position detection
        - Binding domain definitions
        - Assignment thresholds
        
        Args:
            config_name: Name of config file (without .json extension)
            
        Returns:
            Dictionary with configuration parameters
        """
        config_path = self.path_config_dir / f"{config_name}.json"
        if not config_path.exists():
            raise FileNotFoundError(f"GRN config not found: {config_path}")
        
        self.logger.info(f"Loading GRN config from {config_path}")
        with open(config_path, 'r') as f:
            return json.load(f)
    
    # Implement abstract methods from BaseProcessor
    
    def load_entity(self, name: str) -> Optional[pd.Series]:
        """
        Load a single GRN entity (row from a GRN table) by name.
        
        GRN entities are rows in tables, not individual files.
        This method searches all loaded tables for the entity.
        
        Args:
            name: Entity name (sequence ID)
            
        Returns:
            Series with GRN data or None if not found
        """
        # First check if we have data loaded
        if self.data is not None and not self.data.empty:
            if name in self.data.index:
                return self.data.loc[name]
        
        # Search in registered entities
        entity_info = self.entity_registry.find_entity(name, self.processor_type)
        if entity_info and 'table' in entity_info.metadata:
            # Load the table containing this entity
            table_name = entity_info.metadata['table']
            table_data = self.load_grn_table(table_name, return_only=True)
            if table_data is not None and name in table_data.index:
                return table_data.loc[name]
        
        # Search all available tables
        for dataset_info in self.list_datasets():
            if dataset_info.get('type') in ['grn_table', 'grn_reference']:
                table_id = dataset_info['id']
                try:
                    table_data = self.load_grn_table(table_id, return_only=True)
                    if table_data is not None and name in table_data.index:
                        # Register the entity for faster future lookups
                        self.entity_registry.register_entity(
                            name=name,
                            format_type=self.processor_type,
                            file_path=dataset_info['file_path'],
                            metadata={'table': table_id}
                        )
                        return table_data.loc[name]
                except:
                    continue
        
        return None
    
    def save_entity(self, name: str, data: pd.Series, metadata: Optional[dict] = None):
        """
        Save a single GRN entity.
        
        For GRN data, entities are typically saved as part of a table,
        not as individual files. This method adds/updates a row in the
        current GRN table.
        
        Args:
            name: Entity name (sequence ID)
            data: GRN data as a Series
            metadata: Optional metadata dict
        """
        if self.data is None or self.data.empty:
            # Initialize DataFrame with the entity data
            self.data = pd.DataFrame([data], index=[name])
            grn_cols = [col for col in data.index if self._is_grn_position(col)]
            self.grns = grn_cols
        else:
            # Add or update the entity in the current data
            self.data.loc[name] = data
        
        # Register the entity
        table_name = self.dataset or "current"
        # Get the actual table path from ProtosPaths
        table_path = Path(self.paths.get_subdir_path('grn', 'table_dir')) / f"{table_name}.csv"
        # Store relative path in registry
        relative_path = table_path.relative_to(Path(self.paths.data_root))
        
        self.entity_registry.register_entity(
            name=name,
            format_type=self.processor_type,
            file_path=str(relative_path),
            metadata={
                'table': table_name,
                'grn_positions': len([col for col in data.index if self._is_grn_position(col)]),
                **(metadata or {})
            }
        )
        
        # Update internal tracking
        if name not in self.ids:
            self.ids.append(name)
    
    def save_data(self, name: str, data: Any = None, file_format: str = 'csv', **kwargs):
        """
        Save GRN data with proper CSV handling.
        
        Overrides base method to ensure CSV files are saved with index=True
        for GRN tables (to preserve protein IDs).
        """
        if data is None:
            data = self.data
            
        file_path = Path(self.paths.get_subdir_path('grn', 'table_dir')) / f"{name}.{file_format}"
        # ProtosPaths already ensures directories exist
        
        if file_format == 'csv' and isinstance(data, pd.DataFrame):
            # GRN tables must save with index to preserve protein IDs
            data.to_csv(file_path, index=True)
        else:
            # Use parent's save_data for other formats
            return super().save_data(name, data, file_format)
            
        return str(file_path)
            
    def list_available_datasets(self):
        """
        List available GRN datasets in the data directory.
        
        Returns:
            List of dataset IDs
        """
        datasets = self.list_datasets()
        return [d["id"] for d in datasets if d.get("format") == "csv"]
        
    def save_grn_table(self, dataset_id=None, normalize_formats=True, **kwargs):
        """
        Save the current GRN table to a dataset.
        
        GRN tables use entity names as the index (row labels) and GRN positions as columns.
        No entity_id column is added - the index contains the entity names.
        
        Args:
            dataset_id: Dataset identifier (uses current dataset if None)
            normalize_formats: Whether to normalize GRN formats before saving
            **kwargs: Additional format-specific saving parameters
            
        Returns:
            Path to the saved file
            
        Raises:
            ValueError: If no dataset is specified and no current dataset
        """
        # Use provided dataset_id or current dataset
        if dataset_id is not None:
            self.dataset = dataset_id
        elif self.dataset is None:
            raise ValueError("No dataset specified.")
            
        # Create a copy of the data to avoid modifying the original
        if self.data is not None:
            # Only include GRN columns (from self.grns) for saving
            grn_columns = []
            metadata_columns = ['entity_id', 'family', 'species', 'name', 'grn_system', 'id']
            
            for col in self.data.columns:
                if col not in metadata_columns and self._is_grn_position(col):
                    grn_columns.append(col)
                    
            # If self.grns is populated, use it as the source of truth
            if self.grns:
                grn_columns = [col for col in self.grns if col in self.data.columns]
                
            data_to_save = self.data[grn_columns].copy()
            
            # Normalize GRN formats if requested
            if normalize_formats:
                # Track normalization
                normalized_count = 0
                invalid_grns = []
                
                # Get column names excluding the index
                original_columns = data_to_save.columns.tolist()
                normalized_columns = {}
                
                for col in original_columns:
                    # Skip metadata columns including entity_id
                    if col in ['entity_id', 'family', 'species', 'name', 'grn_system', 'id']:
                        continue
                        
                    # Validate and normalize if needed
                    is_valid, message = validate_grn_string(col)
                    
                    if not is_valid:
                        # Try to normalize
                        normalized_col = normalize_grn_format(col)
                        normalized_is_valid, _ = validate_grn_string(normalized_col)
                        
                        if normalized_is_valid and normalized_col != col:
                            self.logger.info(f"Normalized GRN format for saving: {col} -> {normalized_col}")
                            normalized_columns[col] = normalized_col
                            normalized_count += 1
                        else:
                            self.logger.warning(f"Saving with invalid GRN format: {col}")
                            invalid_grns.append(col)
                
                # Rename columns if any were normalized
                if normalized_columns:
                    data_to_save = data_to_save.rename(columns=normalized_columns)
                    
                # Log normalization summary
                if normalized_count > 0:
                    self.logger.info(f"Normalized {normalized_count} GRN formats during save operation")
                if invalid_grns:
                    self.logger.warning(f"Saving with {len(invalid_grns)} invalid GRN formats: {invalid_grns}")
                
            # GRN tables must preserve protein IDs as the index
            # No manipulation of index needed - it will be saved with index=True
            # Ensure index doesn't have a name to match standard GRN format
            data_to_save.index.name = None
        else:
            self.logger.warning("No data to save")
            data_to_save = pd.DataFrame()
            
        # Save to tables/ subdirectory using ProtosPaths
        table_path = Path(self.paths.get_subdir_path('grn', 'table_dir')) / f"{self.dataset}.csv"
        # Ensure index is preserved for protein IDs
        save_index = kwargs.pop('index', True)
        # Use '-' for missing values to maintain consistency with GRN conventions
        na_rep = kwargs.pop('na_rep', '-')
        data_to_save.to_csv(table_path, index=save_index, na_rep=na_rep, **kwargs)
        
        # Create dataset metadata in datasets/
        # Use relative path from datasets to tables directory
        dataset_info = {
            "name": self.dataset,
            "entities": list(data_to_save.index),
            "table_file": str(Path("../tables") / f"{self.dataset}.csv"),
            "metadata": {
                "grn_positions": [col for col in data_to_save.columns if self._is_grn_position(col)],
                "entity_count": len(data_to_save),
                "created": datetime.now().isoformat()
            }
        }
        
        dataset_path = Path(self.paths.get_subdir_path('grn', 'datasets_dir')) / f"{self.dataset}.json"
        with open(dataset_path, 'w') as f:
            json.dump(dataset_info, f, indent=2)
        
        self.logger.info(f"Saved GRN table '{self.dataset}' to {table_path}")
        self.logger.info(f"Saved dataset metadata to {dataset_path}")
        return str(table_path)
        
    def load_grn_table(self, dataset_id=None, low_memory=False, remove_duplicates=True,
                       normalize_formats=True, register_entities=True, return_only=False, **kwargs):
        """
        Load a GRN table from a dataset with entity support.
        
        Args:
            dataset_id: Dataset identifier (uses current dataset if None)
            low_memory: Whether to use pandas low_memory option
            remove_duplicates: Whether to remove duplicate protein IDs
            register_entities: Whether to register entities found in the table
            normalize_formats: Whether to normalize GRN formats (e.g., convert legacy loop formats)
            return_only: If True, return data without setting as current
            **kwargs: Additional format-specific loading parameters
            
        Returns:
            DataFrame containing the GRN table (if return_only=True)
            None (if return_only=False)
            
        Raises:
            FileNotFoundError: If the dataset doesn't exist
        """
        # Determine which dataset to load
        if dataset_id is not None:
            dataset_to_load = dataset_id
            if not return_only:
                self.dataset = dataset_id
        elif self.dataset is not None:
            dataset_to_load = self.dataset
        else:
            raise ValueError("No dataset specified.")
            
        # Try to load from tables/ directory first using ProtosPaths
        table_path = Path(self.paths.get_subdir_path('grn', 'table_dir')) / f"{dataset_to_load}.csv"
        
        if table_path.exists():
            self.logger.info(f"Loading GRN table from {table_path}")
            # Load treating '-' as NaN for consistency
            if 'na_values' not in kwargs:
                kwargs['na_values'] = ['-']
            # Read CSV and ensure column names are preserved as strings
            df = pd.read_csv(table_path, index_col=0, low_memory=low_memory, **kwargs)
            # Normalize column names to standard format (e.g., '1.5' -> '1.50')
            df.columns = [normalize_grn_format(str(col)) for col in df.columns.tolist()]
        else:
            # Try loading from dataset JSON
            dataset_path = Path(self.paths.get_subdir_path('grn', 'datasets_dir')) / f"{dataset_to_load}.json"
            if dataset_path.exists():
                with open(dataset_path, 'r') as f:
                    dataset_info = json.load(f)
                
                if 'table_file' in dataset_info:
                    # Resolve table file path
                    table_file = dataset_path.parent / dataset_info['table_file']
                    if table_file.exists():
                        self.logger.info(f"Loading GRN table from dataset reference: {table_file}")
                        df = pd.read_csv(table_file, index_col=0, low_memory=low_memory, **kwargs)
                        # Normalize column names to standard format (e.g., '1.5' -> '1.50')
                        df.columns = [normalize_grn_format(str(col)) for col in df.columns.tolist()]
                    else:
                        raise FileNotFoundError(f"Table file not found: {table_file}")
                else:
                    raise ValueError(f"Dataset '{dataset_to_load}' has no table_file reference")
            else:
                # No legacy fallback - all data must be in proper directories
                raise FileNotFoundError(f"GRN table '{dataset_to_load}' not found")
        
        # Process the loaded DataFrame to ensure it has the right format
        # Handle case where the protein_id column is present
        if "protein_id" in df.columns:
            df = df.set_index("protein_id")
        # Handle older GRN tables where the index might be called "uen" or "Unnamed: 0"
        elif "uen" in df.columns:
            df = df.set_index("uen")
        elif "Unnamed: 0" in df.columns:
            df = df.rename(columns={"Unnamed: 0": "protein_id"})
            df = df.set_index("protein_id")
            
        # Sort by index and fill NA values
        df = df.sort_index(ascending=True)
        df = df.fillna('-')
        
        # If return_only, work with a copy and don't modify instance state
        if return_only:
            work_df = df.copy()
        else:
            self.data = df
            
        # Remove duplicates if requested
        if remove_duplicates:
            if return_only:
                # For return_only, remove duplicates from the copy
                dup_mask = work_df.index.duplicated(keep='first')
                work_df = work_df[~dup_mask]
            else:
                self.remove_duplicate_ids()
            
        # Update processor state only if not return_only
        if not return_only:
            self.ids = self.data.index.tolist()
            self.grns = self.data.columns.tolist()
        else:
            # For return_only, work with the copy
            self.ids = work_df.index.tolist()
            self.grns = work_df.columns.tolist()
        
        # Register entities if requested and no entity_id column exists
        current_df = work_df if return_only else self.data
        if register_entities and 'entity_id' not in current_df.columns:
            try:
                from protos.io.data_access import GlobalRegistry
                global_registry = GlobalRegistry()
                
                # Register each sequence as an entity
                for seq_id in self.ids:
                    entity_id = generate_entity_id(str(seq_id))
                    global_registry.entity_registry.register_entity(
                        entity_id=entity_id,
                        entity_type="grn",
                        original_id=str(seq_id),
                        file_path=None,  # GRN entities are in tables
                        metadata={
                            "table": dataset_to_load,
                            "grn_positions": len([col for col in self.grns if col not in ['entity_id', 'family', 'species', 'name', 'grn_system', 'id']])
                        },
                        datasets=[dataset_to_load] if dataset_to_load else []
                    )
                self.logger.info(f"Registered {len(self.ids)} GRN entities")
            except Exception as e:
                self.logger.warning(f"Could not register GRN entities: {e}")
        
        # Validate and normalize GRN formats
        if normalize_formats:
            # Store original column names for reference
            original_columns = self.data.columns.tolist()
            normalized_columns = []
            
            # Track validation results
            invalid_grns = []
            normalized_grns = []
            
            for col in original_columns:
                # If it's a GRN, validate and normalize it
                try:
                    # Skip non-GRN columns that might be metadata
                    if not self._is_grn_position(col):
                        normalized_columns.append(col)
                        continue
                        
                    # Validate the GRN format
                    is_valid, message = validate_grn_string(col)
                    
                    if is_valid:
                        # Keep as is, it's already valid
                        normalized_columns.append(col)
                    else:
                        # Try to normalize
                        normalized_col = normalize_grn_format(col)
                        normalized_is_valid, _ = validate_grn_string(normalized_col)
                        
                        if normalized_is_valid:
                            self.logger.info(f"Normalized GRN format: {col} -> {normalized_col}")
                            normalized_columns.append(normalized_col)
                            normalized_grns.append((col, normalized_col))
                        else:
                            # If it's still not valid, use the original for now but track it
                            self.logger.warning(f"Invalid GRN format (keeping as is): {col}")
                            normalized_columns.append(col)
                            invalid_grns.append(col)
                except Exception as e:
                    self.logger.error(f"Error processing column {col}: {e}")
                    # Keep the original to avoid data loss
                    normalized_columns.append(col)
            
            # Rename columns if any were normalized
            if normalized_grns:
                # Create mapping dictionary for column renaming
                rename_map = {old: new for old, new in normalized_grns}
                # Rename columns in the DataFrame
                self.data = self.data.rename(columns=rename_map)
                
            # Log summary of normalization
            if invalid_grns:
                self.logger.warning(f"Found {len(invalid_grns)} invalid GRN formats: {invalid_grns}")
            if normalized_grns:
                self.logger.info(f"Normalized {len(normalized_grns)} GRN formats")
            
            # Update column list with normalized names
            self.grns = self.data.columns.tolist()
        
        # Sort GRNs in standard order
        self.grns = sort_grns_str(self.grns)
        
        # We need to preserve the original notation in the data
        # while also supporting both notations in our methods
        original_columns = self.data.columns.tolist()
        
        # Create bidirectional mappings between dot and x notations
        self.dot_to_x = {}
        self.x_to_dot = {}
        
        for col in original_columns:
            if 'n.' in col or 'c.' in col:
                # Special cases don't need conversion
                continue
                
            if 'x' in col:
                # Convert x notation to dot notation
                tm, pos = col.split('x')
                dot_notation = f"{tm}.{pos}"
                self.x_to_dot[col] = dot_notation
                self.dot_to_x[dot_notation] = col
            elif '.' in col:
                # Convert dot notation to x notation
                try:
                    # Handle loop notation (AB.CCC)
                    if len(col.split('.')[0]) == 2 and len(col.split('.')[1]) == 3:
                        # This is a loop GRN in the format AB.CCC - don't try to convert
                        continue
                    
                    # Handle standard GRN with dot notation
                    tm, pos = col.split('.')
                    x_notation = f"{tm}x{pos}"
                    self.dot_to_x[col] = x_notation
                    self.x_to_dot[x_notation] = col
                except Exception as e:
                    self.logger.error(f"Error processing dot notation {col}: {e}")
        
        # Always use dot notation (x notation is deprecated)
        # No need to detect notation anymore
        
        # Create normalized grns list - convert any x notation to dot notation
        normalized_grns = []
        for grn in self.grns:
            if 'x' in grn and grn[0].isdigit():
                # Convert x notation to dot notation
                dot_form = self.x_to_dot.get(grn)
                if dot_form and dot_form in original_columns:
                    normalized_grns.append(dot_form)
                else:
                    # Manual conversion if not in mapping
                    try:
                        tm, pos = grn.split('x')
                        normalized_grns.append(f"{tm}.{pos}")
                    except:
                        normalized_grns.append(grn)
            else:
                normalized_grns.append(grn)
        
        # Update the GRNs list with normalized forms
        self.grns = normalized_grns
        
        # Use the normalized GRNs to select columns from the data
        # First try direct match, then try flexible matching for dot notation
        columns_to_select = []
        for grn in self.grns:
            if grn in original_columns:
                columns_to_select.append(grn)
            else:
                # Try flexible matching for dot notation (1.50 matches 1.5, 1.05 matches 1.5)
                matched = False
                # Check if it's a loop region (2 digits before dot, 3 digits after)
                if '.' in grn and len(grn.split('.')) == 2:
                    parts = grn.split('.')
                    if len(parts[0]) == 2 and len(parts[1]) == 3:
                        # Loop region like 67.001 - use as is
                        if grn in original_columns:
                            columns_to_select.append(grn)
                            matched = True
                
                if not matched and '.' in grn and grn not in ['n.', 'c.']:
                    # Standard dot notation - try different formats
                    parts = grn.split('.')
                    if len(parts) == 2:
                        helix = parts[0]
                        position = parts[1].lstrip('0')  # Remove leading zeros
                        
                        # Try variations: 1.50 -> 1.5, 1.05 -> 1.5
                        variations = [
                            f"{helix}.{position}",  # 1.5
                            f"{helix}.{int(position):02d}",  # 1.50
                            f"{helix}.{int(position):01d}",  # 1.5 (redundant but safe)
                        ]
                        
                        for var in variations:
                            if var in original_columns:
                                columns_to_select.append(var)
                                matched = True
                                break
                
                if not matched:
                    # Column not found, this will raise an error
                    columns_to_select.append(grn)
        
        try:
            if return_only:
                work_df = work_df[columns_to_select]
                # Update self.grns to match actual column names
                self.grns = list(work_df.columns)
            else:
                self.data = self.data[columns_to_select]
                # Update self.grns to match actual column names
                self.grns = list(self.data.columns)
        except KeyError as e:
            self.logger.error(f"Failed to match column names: {e}")
            raise KeyError(f"Cannot match GRN columns. Data columns: {original_columns}, Requested: {self.grns}")
            
        # Log notation format
        self.logger.info("Using dot notation for GRN positions (x notation is deprecated)")
            
        if return_only:
            work_df = work_df.fillna('-')
        else:
            self.data = self.data.fillna('-')
        
        # Dataset metadata is saved by save_grn_table, no need to register here
        
        if return_only:
            return work_df
        else:
            return None
    
    def get_available_grn_tables(self):
        """
        Get available GRN tables (legacy method).
        
        Returns:
            List of dataset IDs
        """
        self.logger.warning("get_available_grn_tables() is deprecated, use list_available_datasets() instead")
        return self.list_available_datasets()
    
    def load_and_merge_grn_tables(self, datasets, remove_duplicates=True):
        """
        Load and merge multiple GRN tables.
        
        Args:
            datasets: List of dataset IDs to merge
            remove_duplicates: Whether to remove duplicate protein IDs
            
        Returns:
            DataFrame containing the merged GRN table
        """
        # Reset the current data and mappings
        self.data = pd.DataFrame()
        self.dot_to_x = {}
        self.x_to_dot = {}
        tables = []
        
        # First, load all the tables to see what notation formats we're working with
        dot_count = 0
        x_count = 0
        
        for dataset in datasets:
            # Load each dataset and get the data
            self.load_grn_table(dataset_id=dataset, remove_duplicates=remove_duplicates)
            grn_table = self.data.copy()  # Get the loaded data
            tables.append(grn_table)
            
            # Count notation formats
            dot_cols = [col for col in grn_table.columns if '.' in col and not ('n.' in col or 'c.' in col)]
            x_cols = [col for col in grn_table.columns if 'x' in col and not ('n.' in col or 'c.' in col)]
            
            dot_count += len(dot_cols)
            x_count += len(x_cols)
            
        if not tables:
            return pd.DataFrame()
            
        # Always use dot notation (x notation is deprecated)
        
        # Build a comprehensive bidirectional mapping 
        # for all GRN positions across all tables
        for table in tables:
            for col in table.columns:
                if 'n.' in col or 'c.' in col:
                    continue
                    
                if 'x' in col:
                    tm, pos = col.split('x')
                    dot_form = f"{tm}.{pos}"
                    self.x_to_dot[col] = dot_form
                    self.dot_to_x[dot_form] = col
                elif '.' in col:
                    tm, pos = col.split('.')
                    x_form = f"{tm}x{pos}"
                    self.dot_to_x[col] = x_form
                    self.x_to_dot[x_form] = col
        
        # Merge the tables
        merged_table = pd.concat(tables, axis=0)
        
        # Standardize column names to dot notation
        std_columns = []
        for col in merged_table.columns:
            if 'x' in col and col[0].isdigit():
                # Convert x notation to dot notation
                dot_form = self.x_to_dot.get(col)
                if dot_form:
                    std_columns.append(dot_form)
                else:
                    # Manual conversion
                    try:
                        tm, pos = col.split('x')
                        std_columns.append(f"{tm}.{pos}")
                    except:
                        std_columns.append(col)
            else:
                std_columns.append(col)
        
        # Rename columns to the standardized format
        merged_table.columns = std_columns
        
        # Sort columns by GRN position
        self.grns = sort_grns_str(merged_table.columns.tolist())
        
        # Reorder columns using the sorted GRNs
        valid_cols = [col for col in self.grns if col in merged_table.columns]
        self.data = merged_table[valid_cols]
        self.data = self.data.fillna('-')
        self.ids = self.data.index.tolist()
        self.grns = valid_cols
        
        # Update dataset info
        self.dataset = 'merged_' + '_'.join(datasets)
        
        # Save the merged dataset
        self.save_grn_table(self.dataset)
        
        return self.data
    
    def remove_duplicate_ids(self):
        """
        Remove duplicate protein IDs from the data.
        """
        self.logger.info("Removing duplicate IDs...")
        seen = set()
        duplics = [x for x in self.ids if x in seen or seen.add(x)]
        duplics = list(set(duplics))
        
        if len(duplics) > 0:
            self.logger.info(f"Found {len(duplics)} duplicate IDs, removing: {duplics}")
            singles = [x for x in self.ids if x not in duplics]
            df1 = self.data.loc[singles, :]
            
            # Fix for the TypeError: Create a list of Series for concatenation
            df2_rows = []
            for dupli in duplics:
                try:
                    # Get the first occurrence of each duplicate
                    row = self.data.loc[dupli, :].iloc[0]
                    df2_rows.append(row)
                except (IndexError, AttributeError):
                    # Handle case where the duplicate ID doesn't exist in the data
                    continue
            
            # Only concatenate if there are rows to add
            if df2_rows:
                df2 = pd.DataFrame(df2_rows, index=[dupli for dupli in duplics if dupli in self.data.index])
                self.data = pd.concat([df1, df2])
            else:
                self.data = df1
                
            self.ids = self.data.index.tolist()
    
    def get_seq_dict(self):
        """
        Get sequences from GRN table.
        
        Returns:
            Dictionary mapping protein IDs to sequences
        """
        # Use the data directly without reordering columns
        # This avoids notation conversion issues
        grn_table = self.data
        seqs = [get_seq(idx, grn_table) for idx in self.ids]
        return dict(zip(self.ids, seqs))
    
    def reset_data(self, reset_maps=False, reset_features=False):
        """
        Reset the processor data to the original dataset.
        
        Args:
            reset_maps: Whether to reset the maps dictionary
            reset_features: Whether to reset the features DataFrame
        """
        # Store current dataset ID
        current_dataset = self.dataset
        
        # Reset mappings to ensure clean reload
        self.dot_to_x = {}
        self.x_to_dot = {}
        
        # Reload the dataset with all the appropriate conversions
        self.load_grn_table(dataset_id=current_dataset)
        
        # Reset other data structures if requested
        if reset_features:
            self.features = pd.DataFrame(index=self.data.index.tolist())
        if reset_maps:
            self.maps = {}
            
        # Update IDs and GRNs lists (loaded by load_grn_table)
        self.ids = self.data.index.tolist()
        self.grns = self.data.columns.tolist()
    
    def apply_interval(self, grn_interval, apply_to_maps=True):
        """
        Limit data to specific GRN positions.
        
        Args:
            grn_interval: List of GRN positions to keep
            apply_to_maps: Whether to apply the interval to maps
        """
        # Convert the interval to dot notation if needed
        mapped_interval = []
        for grn in grn_interval:
            # Check if x notation needs conversion
            if 'x' in grn and grn[0].isdigit():
                # Convert x notation to dot notation
                dot_form = self.x_to_dot.get(grn)
                if dot_form and dot_form in self.data.columns:
                    mapped_interval.append(dot_form)
                else:
                    # Try manual conversion
                    try:
                        tm, pos = grn.split('x')
                        dot_grn = f"{tm}.{pos}"
                        if dot_grn in self.data.columns:
                            mapped_interval.append(dot_grn)
                        elif grn in self.data.columns:
                            # Fall back to original if conversion doesn't exist
                            mapped_interval.append(grn)
                    except:
                        if grn in self.data.columns:
                            mapped_interval.append(grn)
            elif grn in self.data.columns:
                # Already in the right notation or special notation (n., c.)
                mapped_interval.append(grn)
        
        if not mapped_interval:
            self.logger.warning(f"No valid GRNs in interval {grn_interval}")
            return
            
        self.data = self.data[mapped_interval]
        
        if apply_to_maps:
            self.apply_interval_to_map(mapped_interval)
            
        self.grns = self.data.columns.tolist()
    
    def filter_by_ids(self, ids_to_keep):
        """
        Filter data to include only the specified IDs.
        
        Args:
            ids_to_keep: List of protein IDs to keep
        """
        if self.data is None or len(self.data) == 0:
            return
            
        # Use a more explicit approach to avoid pandas dtype issues
        valid_ids = [id for id in self.ids if id in ids_to_keep]
        
        if valid_ids:
            self.data = self.data.loc[valid_ids]
            self.ids = valid_ids
        else:
            self.logger.warning(f"No valid IDs found in {ids_to_keep}")
    
    def save_temp_table(self, temp_name: str, temp_data: pd.DataFrame) -> Path:
        """
        Save a temporary GRN table to the temp/ directory.
        
        Temporary tables are used for intermediate results during
        processing and can be automatically cleaned up.
        
        Args:
            temp_name: Name for the temporary file
            temp_data: DataFrame to save
            
        Returns:
            Path to the saved temporary file
        """
        # Add timestamp to temp file name to avoid conflicts
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        temp_file = f"{temp_name}_{timestamp}.csv"
        temp_path = self.path_temp_dir / temp_file
        
        # Save the temp table
        temp_data.to_csv(temp_path, index=True, na_rep='-')
        self.logger.info(f"Saved temporary table to {temp_path}")
        
        return temp_path
    
    def clean_temp_files(self, older_than_hours: int = 24):
        """
        Clean up old temporary files.
        
        Args:
            older_than_hours: Remove files older than this many hours
        """
        import time
        
        temp_dir = self.path_temp_dir
        if not temp_dir.exists():
            return
            
        current_time = time.time()
        cutoff_time = current_time - (older_than_hours * 3600)
        
        removed_count = 0
        for temp_file in temp_dir.glob("*.csv"):
            if temp_file.stat().st_mtime < cutoff_time:
                temp_file.unlink()
                removed_count += 1
                
        if removed_count > 0:
            self.logger.info(f"Removed {removed_count} old temporary files")
    
    def get_grn_dict(self, reset_data=False, notation=None):
        """
        Get GRN dictionary mapping proteins to GRN positions.
        
        Args:
            reset_data: Whether to reset the data before getting the dictionary
            notation: Which notation format to use in the returned dictionary
                      None: use the format in the data (default)
                      'x': convert all positions to x notation (3x50)
                      'dot': convert all positions to dot notation (3.50)
            
        Returns:
            Dictionary mapping protein IDs to lists of GRN positions
        """
        if reset_data:
            self.reset_data()
            
        # Create the basic GRN dictionary with original column names
        grn_table = self.data
        residue_mask = grn_table.replace('-', pd.NA).notna()
        grn_dict = {uen: residue_mask.columns[residue_mask.loc[uen]].tolist() for uen in residue_mask.index}
        
        # If no specific notation requested, use whatever is in the data
        if notation is None:
            return grn_dict
            
        # Convert to requested notation format
        converted_dict = {}
        for uen, grn_list in grn_dict.items():
            converted_grns = []
            for grn in grn_list:
                # No conversion needed for n. and c. notations
                if 'n.' in grn or 'c.' in grn:
                    converted_grns.append(grn)
                    continue
                    
                # Convert to requested notation
                if notation == 'x' and '.' in grn:
                    # Convert dot to x
                    converted_grns.append(self.dot_to_x.get(grn, grn))
                elif notation == 'dot' and 'x' in grn:
                    # Convert x to dot
                    converted_grns.append(self.x_to_dot.get(grn, grn))
                else:
                    # Already in requested format
                    converted_grns.append(grn)
                    
            converted_dict[uen] = converted_grns
            
        return converted_dict
    
    def filter_data_by_occurances(self, threshold):
        """
        Filter data to include only GRN positions that occur in at least threshold proteins.
        
        Args:
            threshold: Minimum number of proteins that must have the GRN position
        """
        # Calculate the number of '-' entries in each column of the 'data' DataFrame
        non_existent_counts = (self.data == '-').sum()
        # Calculate the number of genes in each column
        gene_counts = len(self.data) - non_existent_counts
        # Get the column names where the gene count is greater than or equal to the threshold
        filtered_columns = gene_counts[gene_counts >= threshold].index
        # Update the 'data' DataFrame, the maps, and the 'grns' list using the 'apply_interval' function
        self.apply_interval(filtered_columns)
    
    def sort_columns(self):
        """
        Sort columns by GRN position.
        
        This method sorts the GRN columns according to the standard order:
        1. N-terminal regions first
        2. TM regions in numerical order
        3. Loop regions in order of TM numbers
        4. C-terminal regions last
        """
        # Get current columns
        current_columns = self.data.columns.tolist()
        
        # Skip if no columns
        if not current_columns:
            return
            
        try:
            # Convert to float values for sorting
            # First validate and normalize to ensure consistent handling
            normalized_columns = []
            for col in current_columns:
                # Skip non-GRN columns
                if col in ['family', 'species', 'name', 'grn_system', 'id']:
                    normalized_columns.append(col)
                    continue
                    
                # Validate and normalize
                is_valid, _ = validate_grn_string(col)
                if not is_valid:
                    normalized_col = normalize_grn_format(col)
                    normalized_is_valid, _ = validate_grn_string(normalized_col)
                    if normalized_is_valid:
                        normalized_columns.append(normalized_col)
                    else:
                        # Keep original if can't normalize
                        normalized_columns.append(col)
                else:
                    normalized_columns.append(col)
                    
            # Create a mapping from original to normalized
            normalization_map = {orig: norm for orig, norm in zip(current_columns, normalized_columns)}
            
            # Convert normalized columns to float values for sorting
            try:
                cols_unsorted = []
                metadata_cols = []
                column_types = {}  # To track special column types
                
                for i, col in enumerate(normalized_columns):
                    if col in ['family', 'species', 'name', 'grn_system', 'id']:
                        metadata_cols.append((i, col))
                        continue
                    
                    try:
                        float_val = parse_grn_str2float(col)
                        cols_unsorted.append(float_val)
                        
                        # Track column type
                        if 'n.' in col:
                            column_types[float_val] = 'n_term'
                        elif 'c.' in col:
                            column_types[float_val] = 'c_term'
                        elif '.' in col and len(col.split('.')[0]) == 2 and len(col.split('.')[1]) == 3:
                            column_types[float_val] = 'loop'
                        elif 'x' in col:
                            column_types[float_val] = 'standard'
                    except Exception as e:
                        self.logger.warning(f"Cannot convert {col} to float for sorting: {e}")
                        # Add to metadata to preserve
                        metadata_cols.append((i, col))
            except Exception as e:
                self.logger.error(f"Error converting GRNs to float for sorting: {e}")
                # Fall back to alphabetical sorting to avoid data loss
                self.data = self.data.reindex(columns=sorted(current_columns))
                self.grns = self.data.columns.tolist()
                return
                
            # Sort the float values
            sorted_float_values = sort_grns(cols_unsorted)
            
            # Prepare sorted column names list - handle both x and dot notation
            cols_sorted = []
            
            for float_val in sorted_float_values:
                col_type = column_types.get(float_val, 'standard')
                
                # Get the appropriate string representation
                if col_type == 'n_term' or col_type == 'c_term':
                    # N-terminal and C-terminal formats don't change
                    cols_sorted.append(parse_grn_float2str(float_val))
                elif col_type == 'loop':
                    # Loop format should remain as <closer helix><further helix>.<distance>
                    cols_sorted.append(parse_grn_float2str(float_val))
                else:
                    # Standard GRN - use the correct notation based on processor preference
                    # Always use dot notation for TM regions
                    cols_sorted.append(parse_grn_float2str(float_val, notation_type='dot'))
            
            # Add metadata columns at the beginning
            for idx, col in sorted(metadata_cols, key=lambda x: x[0]):
                cols_sorted.insert(idx, col)
            
            # Create a mapping between normalized sorted columns and original column names
            # This is needed because the actual column names might differ from our standard notation
            col_mapping = {}
            
            # First map normalized columns to original columns
            reverse_norm_map = {norm: orig for orig, norm in normalization_map.items()}
            
            # Map sorted columns back to original columns
            for sorted_col in cols_sorted:
                # First get the original column name if it was normalized
                original_col = reverse_norm_map.get(sorted_col, sorted_col)
                
                # If it's in the current columns, use it directly
                if original_col in current_columns:
                    col_mapping[sorted_col] = original_col
                # If it's a normalized version, find the corresponding original
                elif sorted_col in normalized_columns:
                    idx = normalized_columns.index(sorted_col)
                    col_mapping[sorted_col] = current_columns[idx]
                else:
                    # This shouldn't happen, but just in case
                    col_mapping[sorted_col] = sorted_col
            
            # Use the mapping to reorder columns - preserve any columns not mapped
            try:
                self.data = self.data.loc[:, [col_mapping.get(col, col) for col in cols_sorted]]
                self.grns = self.data.columns.tolist()
            except KeyError as e:
                self.logger.error(f"Error reordering columns: {e}")
                # Preserve data by not changing the order
                pass
                
        except Exception as e:
            self.logger.error(f"Error sorting columns: {e}")
            # Don't change column order if there's an error to preserve data
    
    # Maps methods
    def get_maps(self):
        """
        Get list of map names.
        
        Returns:
            List of map names
        """
        return list(self.maps.keys())
    
    def apply_interval_to_map(self, grn_interval):
        """
        Apply GRN interval to maps.
        
        Args:
            grn_interval: List of GRN positions to keep
        """
        maps = self.get_maps()
        
        # Convert interval to dot notation
        dot_interval = []
        
        for grn in grn_interval:
            # Handle x notation conversion
            if 'x' in grn and grn[0].isdigit():
                dot_form = self.x_to_dot.get(grn)
                if dot_form:
                    dot_interval.append(dot_form)
                else:
                    # Create dot notation if not in mapping
                    try:
                        tm, pos = grn.split('x')
                        dot_interval.append(f"{tm}.{pos}")
                    except:
                        dot_interval.append(grn)
            else:
                # Keep as is (includes n., c., and already dot notation)
                dot_interval.append(grn)
        
        # Apply to each map
        for map_key in maps:
            map_cols = self.maps[map_key].columns.tolist()
            
            # Filter to columns that exist in this map
            valid_interval = [col for col in dot_interval if col in map_cols]
            
            # Apply if we have valid columns
            if valid_interval:
                self.maps[map_key] = self.maps[map_key][valid_interval]
                self.logger.info(f"Applied interval to map '{map_key}': {len(valid_interval)} GRNs")
            else:
                self.logger.warning(f"No valid GRNs to apply to map '{map_key}'")
    
    # Entity-specific methods for GRN tables
    def load_grn_entity(self, identifier: str) -> Optional[pd.Series]:
        """
        Load a single GRN entity (one row from a GRN table).
        
        Args:
            identifier: Entity identifier (sequence ID or entity hash)
            
        Returns:
            Series with GRN annotations for the entity or None if not found
        """
        # Resolve identifier
        try:
            from protos.io.data_access import GlobalRegistry
            global_registry = GlobalRegistry()
            entity_id = global_registry.entity_registry.resolve_identifier(identifier, format_type="grn")
            
            # Get original ID 
            original_id = global_registry.entity_registry.get_original_id(entity_id)
            if not original_id:
                original_id = identifier
        except Exception as e:
            # Fallback to direct lookup
            self.logger.debug(f"Could not resolve entity ID: {e}")
            original_id = identifier
            entity_id = generate_entity_id(identifier)
        
        # Check if we have this entity in current data
        if self.data is not None and original_id in self.data.index:
            return self.data.loc[original_id]
        
        # If not loaded, try to find which table contains this entity
        try:
            global_registry = GlobalRegistry()
            entity_info = global_registry.entity_registry.get_entity(entity_id)
            if entity_info and 'grn' in entity_info.get('formats', {}):
                table_name = entity_info['formats']['grn']['metadata'].get('table')
                if table_name:
                    # Load the table containing this entity
                    self.load_grn_table(table_name)
                    if original_id in self.data.index:
                        return self.data.loc[original_id]
        except:
            pass
        
        self.logger.warning(f"GRN entity not found: {identifier}")
        return None
    
    def list_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all GRN entities (sequence IDs).
        
        Args:
            dataset: Optional dataset/table to filter by
            
        Returns:
            List of sequence IDs (not hash IDs!)
        """
        if dataset:
            # Load specific dataset first
            self.load_grn_table(dataset)
            return self.ids
        elif self.data is not None:
            # Return current loaded entities
            return self.ids
        else:
            # List all GRN entities from registry
            try:
                from protos.io.data_access import GlobalRegistry
                global_registry = GlobalRegistry()
                entity_ids = global_registry.entity_registry.list_entities(format_type="grn", dataset=dataset)
                
                # Convert to original IDs
                original_ids = []
                for entity_id in entity_ids:
                    original_id = global_registry.entity_registry.get_original_id(entity_id)
                    if original_id:
                        original_ids.append(original_id)
                return original_ids
            except:
                return []
    
    def list_grn_entities(self, dataset: Optional[str] = None) -> List[str]:
        """
        List all GRN entities (backward compatibility).
        
        Deprecated: Use list_entities() instead.
        
        Args:
            dataset: Optional dataset/table to filter by
            
        Returns:
            List of sequence IDs
        """
        return self.list_entities(dataset=dataset)
    
    def _is_grn_position(self, col: str) -> bool:
        """
        Check if a column name is a GRN position.
        
        Args:
            col: Column name to check
            
        Returns:
            True if column is a GRN position, False otherwise
        """
        # Skip metadata columns - expanded list to include common column names
        metadata_columns = [
            'family', 'species', 'name', 'grn_system', 'id', 'entity_id',
            'sequence_id', 'protein_id', 'uniprot_id', 'pdb_id', 'chain_id',
            'description', 'organism', 'gene_name', 'protein_name',
            'Unnamed: 0', 'index'  # Common pandas artifacts
        ]
        
        if col.lower() in [c.lower() for c in metadata_columns]:
            return False
            
        # Check for GRN format (e.g., "3.50", "7.53", "34.50", "n.50", "c.50")
        # Also allow x notation like "3x50", "7x53"
        # Extended to support 0.XX and 9.XX formats, and loop regions
        import re
        
        # Match N-terminal (n.1, n.10)
        n_term_pattern = re.compile(r'^n\.[0-9]+$')
        # Match C-terminal (c.1, c.10)
        c_term_pattern = re.compile(r'^c\.[0-9]+$')
        # Match standard GRN formats with x or dot notation (1.50, 3x50)
        standard_pattern = re.compile(r'^[0-9]+[x.][0-9]+$')
        
        return bool(n_term_pattern.match(col) or 
                   c_term_pattern.match(col) or 
                   standard_pattern.match(col))
    
    def list_datasets(self) -> List[Dict[str, Any]]:
        """
        List available GRN datasets (GRN tables).
        
        Each GRN table is a dataset containing multiple entity rows.
        
        Returns:
            List of dataset information dictionaries
        """
        datasets = []
        
        # List GRN tables in the tables directory
        tables_dir = Path(self.paths.get_subdir_path('grn', 'table_dir'))
        if tables_dir.exists():
            for table_file in tables_dir.glob("*.csv"):
                try:
                    # Quick scan to get row count and columns
                    df = pd.read_csv(table_file, nrows=1)
                    total_rows = sum(1 for _ in open(table_file)) - 1  # Subtract header
                    
                    # Get GRN columns (format X.XX)
                    grn_cols = [col for col in df.columns if self._is_grn_position(col)]
                    
                    datasets.append({
                        'id': f'tables/{table_file.stem}',
                        'type': 'grn_table',
                        'format': 'csv',
                        'entity_count': total_rows,
                        'grn_positions': len(grn_cols),
                        'file_path': str(table_file),
                        'file_size': table_file.stat().st_size
                    })
                except:
                    pass
        
        # List GRN tables in reference directory
        ref_dir = Path(self.paths.get_subdir_path('grn', 'ref_dir'))
        if ref_dir.exists():
            for table_file in ref_dir.glob("*.csv"):
                try:
                    df = pd.read_csv(table_file, nrows=1)
                    total_rows = sum(1 for _ in open(table_file)) - 1
                    grn_cols = [col for col in df.columns if self._is_grn_position(col)]
                    
                    datasets.append({
                        'id': f'ref/{table_file.stem}',
                        'type': 'grn_reference',
                        'format': 'csv',
                        'entity_count': total_rows,
                        'grn_positions': len(grn_cols),
                        'file_path': str(table_file),
                        'file_size': table_file.stat().st_size
                    })
                except:
                    pass
        
        # Also check dataset manager if available
        if self.dataset_manager:
            managed_datasets = self.dataset_manager.list_datasets()
            for ds in managed_datasets:
                if isinstance(ds, dict) and ds not in datasets:
                    datasets.append(ds)
        
        return datasets
    
    def get_entity_grn_positions(self, identifier: str) -> List[str]:
        """
        Get the GRN positions that have residues for a given entity.
        
        Args:
            identifier: Entity identifier (sequence ID or entity hash)
            
        Returns:
            List of GRN positions with residues (not '-')
        """
        entity_data = self.load_grn_entity(identifier)
        if entity_data is not None:
            # Return columns where value is not '-'
            return [col for col in entity_data.index if entity_data[col] != '-']
        return []