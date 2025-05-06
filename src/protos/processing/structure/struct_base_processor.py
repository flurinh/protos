"""
Structure Processor with BaseProcessor integration.

This module provides a CifBaseProcessor class that extends BaseProcessor
to provide standardized data management capabilities for structural data.
"""

import os
import time
import json
import warnings
import tempfile
import shutil
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Union, Set, Any, Tuple

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

from protos.io import cif_utils
from protos.core.base_processor import BaseProcessor
from protos.processing.structure.struct_utils import (
    load_structure as load_structure_util,
    STRUCT_COLUMNS,
    STRUCT_COLUMN_DTYPE,
    SORTED_STRUCT_COLUMNS,
    ALPHA_CARBON
)


class CifBaseProcessor(BaseProcessor):
    """
    Processor for structural data in mmCIF format.
    
    This class extends BaseProcessor to provide standardized data management
    capabilities for structural data, such as loading and saving structures,
    filtering by various criteria, and extracting specific features.
    """
    
    def __init__(
            self,
            name="cif_processor",
            data_root=None,
            processor_data_dir="structure",
            structure_dir="mmcif",
            dataset_dir="structure_dataset",
            alignments_dir="alignments",
            pdb_ids_file=None,
            limit=None,
            path=None,
            preload=False,
            remove_hetatm=False,
            allow_exception=False
    ):
        """
        Initialize the CIF processor.
        
        Args:
            name: Processor instance name
            data_root: Root directory for all data
            processor_data_dir: Subdirectory for structure data
            structure_dir: Subdirectory for mmCIF files
            dataset_dir: Subdirectory for dataset files
            alignments_dir: Subdirectory for alignment files
            pdb_ids_file: File containing PDB IDs to load
            limit: Maximum number of PDB files to process
            path: Legacy path parameter (deprecated, use data_root instead)
            preload: Whether to load pdb_ids on initialization
            remove_hetatm: Whether to remove HETATM records
            allow_exception: Whether to allow exceptions during processing
        """
        # Handle legacy path parameter
        if path is not None:
            import logging
            logging.warning("The 'path' parameter is deprecated, use 'data_root' and 'processor_data_dir' instead")
            # Extract data_root from path
            if os.path.isabs(path):
                data_root = os.path.dirname(path)
                processor_data_dir = os.path.basename(path)
            else:
                parts = path.rstrip('/').split('/')
                if len(parts) >= 2:
                    data_root = '/'.join(parts[:-1])
                    processor_data_dir = parts[-1]
                else:
                    processor_data_dir = path
        
        # Initialize BaseProcessor
        super().__init__(name=name, data_root=data_root, processor_data_dir=processor_data_dir)
        
        # Store the subdirectories
        self.structure_dir = structure_dir
        self.dataset_dir = dataset_dir
        self.alignments_dir = alignments_dir
        
        # Set up actual paths
        self.path_structure_dir = os.path.join(self.data_root, processor_data_dir, structure_dir)
        self.path_dataset_dir = os.path.join(self.data_root, processor_data_dir, dataset_dir)
        self.path_alignment_dir = os.path.join(self.data_root, processor_data_dir, alignments_dir)
        
        # Create directories if they don't exist
        os.makedirs(self.path_structure_dir, exist_ok=True)
        os.makedirs(self.path_dataset_dir, exist_ok=True)
        os.makedirs(self.path_alignment_dir, exist_ok=True)
        
        # Set up parameters
        self.pdb_ids_file = pdb_ids_file
        if pdb_ids_file is not None:
            self.path_pdb_ids_file = os.path.join(self.data_root, processor_data_dir, pdb_ids_file)
        else:
            self.path_pdb_ids_file = None
            
        # Configure processing options
        self.limit = limit
        self.remove_hetatm = remove_hetatm
        self.allow_exception = allow_exception
        
        # Set up data containers
        self.pdb_ids = []
        self.structure_filenames = []
        self.ref_seq = None
        self.dfl = []
        self.dfl_list = []
        self.chain_dict = {}
        self.chain_dict_ca_atom_ids = {}
        self.data = None
        self.structure_info = {}
        self.chain_ids = []
        
        # Setup temp directory for temporary files
        self.temp_dir = Path(os.path.join(self.data_root, processor_data_dir, 'temp_cif'))
        os.makedirs(self.temp_dir, exist_ok=True)
        
        # Initialize PDB IDs if requested
        if preload:
            self.initialize_pdb_ids()
    
    # Initialize and setup methods
    def initialize_pdb_ids(self):
        """Set up PDB IDs from file or directory listing."""
        if self.path_pdb_ids_file is not None:
            self.logger.info(f"Loading PDB ids from {self.path_pdb_ids_file}.")
            self.pdb_ids = self.get_pdb_ids_from_file()
        else:
            self.pdb_ids = self.get_available_pdb_files()
        
        # Update structure filenames based on PDB IDs
        if self.pdb_ids:
            self.structure_filenames = [os.path.join(self.path_structure_dir, pdb_id + '.cif')
                                      for pdb_id in self.pdb_ids]
        
        # Apply limit if specified
        if self.limit is None:
            self.limit = len(self.pdb_ids)
        
        if len(self.structure_filenames) > self.limit:
            self.structure_filenames = self.structure_filenames[:self.limit]
            self.pdb_ids = self.pdb_ids[:self.limit]
            
        self.logger.info(f"Initialized with {len(self.pdb_ids)} PDB IDs")
    
    # File and PDB ID management methods
    def get_pdb_ids_from_file(self, path_pdb_ids_file=None):
        """
        Get PDB IDs from a file.
        
        Args:
            path_pdb_ids_file: Path to file containing PDB IDs
            
        Returns:
            List of PDB IDs
            
        Raises:
            ValueError: If no file is specified
        """
        if path_pdb_ids_file is not None:
            self.path_pdb_ids_file = path_pdb_ids_file
        
        if self.path_pdb_ids_file is None:
            raise ValueError("No PDB IDs file specified")
        
        with open(self.path_pdb_ids_file, 'r') as f:
            pdb_ids = f.readlines()[0].split(',')
            pdb_ids = [pdb_id.strip() for pdb_id in pdb_ids]
        
        return pdb_ids
    
    def get_available_pdb_files(self):
        """
        Get PDB IDs based on available CIF files.
        
        Returns:
            List of PDB IDs from CIF files
        """
        all_files = os.listdir(self.path_structure_dir)
        # Filter the list to only include files with a '.cif' extension
        pdb_ids = [file.replace('.cif', '')
                   for file in all_files if file.lower().endswith('.cif')]
        return pdb_ids
    
    # Core structure loading methods
    def load_structure(self, pdb_id, apply_dtypes=True, debug=False):
        """
        Load a structure from a CIF file.
        
        Args:
            pdb_id: PDB ID to load
            apply_dtypes: Whether to apply proper data types to the loaded structure
            debug: Whether to print debug information during data type formatting
            
        Returns:
            DataFrame with structure data or None if loading fails
        """
        try:
            # Call the utility function to load the structure
            structure = load_structure_util(pdb_id, folder=self.path_structure_dir)
            
            # Ensure the structure is properly loaded
            if structure is None or len(structure) == 0:
                self.logger.warning(f"Failed to load structure for {pdb_id}")
                return None
            
            # Add the PDB ID to our list if not already present
            if pdb_id not in self.pdb_ids:
                self.pdb_ids.append(pdb_id)
            
            # Filter HETATM records if requested
            if self.remove_hetatm:
                structure = structure[structure["group"] == "ATOM"]
            
            # Register the structure data for tracking
            self._register_dataset(
                f"structure/{pdb_id}",
                {
                    "type": "structure",
                    "atom_count": len(structure),
                    "chains": structure["auth_chain_id"].unique().tolist(),
                }
            )
            
            # Add to our data frame list
            if pdb_id in [df['pdb_id'].iloc[0] for df in self.dfl if df is not None and not df.empty]:
                # Replace existing data
                for i, df in enumerate(self.dfl):
                    if df is not None and not df.empty and df['pdb_id'].iloc[0] == pdb_id:
                        self.dfl[i] = structure
                        break
            else:
                # Add new data
                self.dfl.append(structure)
            
            # Update the main data structure
            self.concat_data()
            
            # Apply data types if requested
            if apply_dtypes:
                self.format_data_types(debug=debug)
                
                if debug:
                    self.logger.info(f"Applied data types to structure {pdb_id}")
                    
                    # Print debugging info for critical columns
                    for col in ['x', 'y', 'z', 'auth_seq_id']:
                        if col in self.data.columns:
                            self.logger.info(f"Column {col} dtype: {self.data[col].dtype}")
                            self.logger.info(f"Column {col} sample: {self.data[col].head()}")
            
            # Return the loaded structure
            return structure
        except Exception as e:
            if self.allow_exception:
                self.logger.error(f"Error loading structure {pdb_id}: {str(e)}")
                return None
            else:
                raise
    
    def load_structures(self, pdb_ids=None, apply_dtypes=True, debug=False):
        """
        Load multiple structures.
        
        Args:
            pdb_ids: List of PDB IDs to load (defaults to self.pdb_ids)
            apply_dtypes: Whether to apply proper data types to loaded structures
            debug: Whether to print debug information during data type formatting
            
        Returns:
            Concatenated DataFrame of loaded structures
        """
        if pdb_ids is None:
            pdb_ids = self.pdb_ids
        
        # Clear existing data
        self.dfl = []
        
        # Load each structure with type conversion disabled initially
        # We'll do a single type conversion at the end for better performance
        for pdb_id in tqdm(pdb_ids, desc="Loading structures"):
            # Don't apply types for individual structures yet
            structure = self.load_structure(pdb_id, apply_dtypes=False)
            self.dfl.append(structure)
        
        # Update the data structure
        self.concat_data()
        
        # Apply data types once to the entire dataset if requested
        if apply_dtypes and self.data is not None and not self.data.empty:
            if debug:
                self.logger.info(f"Applying data types to {len(pdb_ids)} structures at once")
            
            # Apply data types to all structures
            self.format_data_types(debug=debug)
            
            if debug:
                self.logger.info(f"Data type formatting complete for {len(pdb_ids)} structures")
                # Print statistics about key columns
                for col in ['x', 'y', 'z', 'auth_seq_id']:
                    if col in self.data.columns:
                        self.logger.info(f"Column {col} dtype: {self.data[col].dtype}")
        
        return self.data
    
    # Data management methods
    def concat_data(self):
        """Concatenate all DataFrames in dfl into the data attribute."""
        if not self.dfl:
            return
        
        # Filter out None values and empty DataFrames
        valid_dfl = [df for df in self.dfl if df is not None and not df.empty]
        if not valid_dfl:
            return
        
        self.data = pd.concat(valid_dfl)
        self.reset_index()
    
    def reset_index(self):
        """Reset and restructure the data index."""
        if self.data is None or len(self.data) == 0:
            return
        
        self.data.reset_index(drop=True, inplace=True)
        pdb_id_mapping = {pdb_id: i for i, pdb_id in enumerate(self.data['pdb_id'].unique(), 1)}
        self.data['idx'] = self.data['pdb_id'].map(pdb_id_mapping)
        self.data.set_index(['idx', self.data.groupby('pdb_id').cumcount()], inplace=True)
        self.data.index.names = [None, None]
    
    def update_pdb_ids(self):
        """Update pdb_ids list from current data."""
        if self.data is not None:
            self.pdb_ids = self.data['pdb_id'].unique().tolist()
    
    def update_chain_data(self):
        """Update chain dictionaries after data changes."""
        self.update_pdb_chain_ids()
    
    def update_pdb_chain_ids(self):
        """Update chain_ids list from chain_dict."""
        # List of all chains in all structures (<pdb_id>_<chain_id>)
        self.chain_ids = list(self.chain_dict.keys())
        
    def format_data_types(self, debug=False):
        """
        Apply correct data types to all columns in the structure data.
        
        This function ensures that all columns have the correct data types
        according to the STRUCT_COLUMN_DTYPE mapping. It handles special cases
        like coordinates, sequence IDs, and atom IDs, ensuring that:
        - Coordinates (x, y, z) are float64
        - Sequence IDs (auth_seq_id, gen_seq_id) are int64
        - Chain IDs and other string fields are properly formatted
        
        Args:
            debug (bool): Whether to print debug information
            
        Returns:
            None (modifies self.data in place)
        """
        if self.data is None or self.data.empty:
            self.logger.warning("No data to format")
            return
            
        if debug:
            self.logger.info(f"Formatting data types for {len(self.pdb_ids)} structures")
            self.logger.info(f"Before formatting: {self.data.dtypes}")
            
        # Create a mapping of column names to desired data types
        column_types = {}
        
        # Start with the standard column types from STRUCT_COLUMN_DTYPE
        for col, dtype in STRUCT_COLUMN_DTYPE.items():
            if col in self.data.columns:
                column_types[col] = dtype
        
        # Add common additional column types that might be present
        additional_columns = {
            # Coordinate columns should be float
            'cartn_x': float, 'cartn_y': float, 'cartn_z': float,
            'model_x': float, 'model_y': float, 'model_z': float,
            
            # B-factor and occupancy should be float
            'b_factor': float, 'b_iso_or_equiv': float, 'occupancy': float,
            
            # Chain identifiers should be string
            'label_asym_id': str, 'auth_asym_id': str,
            'label_chain_id': str, 'auth_chain_id': str,
            
            # Atom identifiers should be string
            'label_atom_id': str, 'auth_atom_id': str,
            'type_symbol': str, 'element': str,
            
            # Residue identifiers should be string
            'label_comp_id': str, 'auth_comp_id': str,
            'ins_code': str, 'pdbx_PDB_ins_code': str,
            'alt_id': str, 'label_alt_id': str,
            
            # Model identifiers should be int
            'model_id': int, 'pdbx_PDB_model_num': int,
        }
        
        # Add additional columns if they exist in the data
        for col, dtype in additional_columns.items():
            if col in self.data.columns and col not in column_types:
                column_types[col] = dtype
        
        # Apply the data types, with proper error handling
        for col, dtype in column_types.items():
            try:
                if dtype == float:
                    # Convert to numeric first to handle potential string values
                    self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
                elif dtype == int:
                    # Handle potentially missing values by first converting to float
                    if not pd.api.types.is_numeric_dtype(self.data[col]):
                        self.data[col] = pd.to_numeric(self.data[col], errors='coerce')
                    
                    # Convert to integer where possible, keeping NaN as is
                    mask = ~pd.isna(self.data[col])
                    if mask.any():
                        # Create a new column with Int64 dtype (nullable integer)
                        int_values = pd.Series(index=self.data.index, dtype='Int64')
                        int_values.loc[mask] = self.data.loc[mask, col].astype('int64')
                        self.data[col] = int_values
                elif dtype == str:
                    # Convert to string while handling NaN values
                    self.data[col] = self.data[col].astype(str)
                    # Replace 'nan' strings with empty strings
                    self.data[col] = self.data[col].replace('nan', '')
                else:
                    # For other types, use normal astype
                    self.data[col] = self.data[col].astype(dtype)
                    
                if debug:
                    self.logger.info(f"Formatted column {col} to {dtype.__name__}")
                    
            except Exception as e:
                self.logger.warning(f"Failed to format column {col} to {dtype.__name__}: {str(e)}")
                if debug:
                    self.logger.info(f"Sample values: {self.data[col].head()}")
        
        # Ensure x, y, z coordinates are always float
        for coord in ['x', 'y', 'z']:
            if coord in self.data.columns:
                try:
                    self.data[coord] = pd.to_numeric(self.data[coord], errors='coerce')
                except Exception as e:
                    self.logger.warning(f"Failed to convert {coord} to float: {str(e)}")
        
        if debug:
            self.logger.info(f"After formatting: {self.data.dtypes}")
            self.logger.info("Data type formatting complete")
    
    # Persistence methods using BaseProcessor
    def save_data(self, dataset_id, data=None, **kwargs):
        """
        Save structure data to a dataset.
        
        Args:
            dataset_id: Dataset identifier
            data: Data to save (uses self.data if None)
            **kwargs: Additional parameters for BaseProcessor.save_data
            
        Returns:
            Path to saved file
        """
        if data is None:
            if self.data is None:
                raise ValueError("No data to save")
            data = self.data
        
        # Set default file format
        file_format = kwargs.pop("file_format", "pkl")
        
        # Save using BaseProcessor's save_data method
        path = super().save_data(dataset_id, data, file_format=file_format, **kwargs)
        
        # Register metadata
        if not dataset_id.startswith("structure/") and "/" not in dataset_id:
            # Standardize dataset ID format for structures
            dataset_id = f"structure/{dataset_id}"
        
        # Update metadata
        self._register_dataset(
            dataset_id, 
            {
                "type": "structure",
                "saved_at": datetime.now().isoformat(),
                "atom_count": len(data),
                "chains": data["auth_chain_id"].unique().tolist() if "auth_chain_id" in data.columns else [],
                "format": file_format,
            }
        )
        
        return path
        
    def load_data(self, dataset_id, apply_dtypes=True, debug=False, **kwargs):
        """
        Load structure data from a dataset.
        
        Args:
            dataset_id: Dataset identifier
            apply_dtypes: Whether to apply proper data types to loaded data
            debug: Whether to print debug information during data type formatting
            **kwargs: Additional parameters for BaseProcessor.load_data
            
        Returns:
            Loaded data
        """
        # Handle legacy dataset IDs by standardizing format
        if not dataset_id.startswith("structure/") and "/" not in dataset_id:
            # Try standardized format first
            std_dataset_id = f"structure/{dataset_id}"
            try:
                return super().load_data(std_dataset_id, **kwargs)
            except FileNotFoundError:
                # Fall back to original ID
                pass
        
        # Load using BaseProcessor's load_data method
        self.data = super().load_data(dataset_id, **kwargs)
        
        # Apply data types if requested
        if apply_dtypes and self.data is not None and not self.data.empty:
            if debug:
                self.logger.info(f"Applying data types to loaded dataset '{dataset_id}'")
                self.logger.info(f"Before formatting: {self.data.dtypes}")
            
            # Use the comprehensive format_data_types function
            self.format_data_types(debug=debug)
            
            if debug:
                self.logger.info(f"Data type formatting complete for dataset '{dataset_id}'")
                # Print info about key columns
                for col in ['x', 'y', 'z', 'auth_seq_id']:
                    if col in self.data.columns:
                        self.logger.info(f"Column {col} dtype: {self.data[col].dtype}")
        
        # Update state
        if self.data is not None:
            self.update_pdb_ids()
            self.dfl = [self.data[self.data['pdb_id'] == pdb_id] for pdb_id in self.pdb_ids]
        
        return self.data

    def filter_by_ids(self, pdb_ids: List[str]):
        """
        Filter data by PDB IDs.

        Args:
            pdb_ids: List of PDB IDs to keep
        """
        if self.data is None or len(self.data) == 0:
            return

        # Use a simpler approach to avoid pandas type issues
        filtered_rows = []
        pdb_ids_set = set(str(pid) for pid in pdb_ids)

        # Filter data
        for idx, row in self.data.iterrows():
            if str(row['pdb_id']) in pdb_ids_set:
                filtered_rows.append(row)

        # Create new DataFrame from filtered rows
        if filtered_rows:
            self.data = pd.DataFrame(filtered_rows)
        else:
            self.data = pd.DataFrame(columns=self.data.columns)

        # Update processor state
        self.pdb_ids = sorted(list(self.data['pdb_id'].unique()))

    def reset_data(self, preserve_ids: bool = False):
        """
        Reset the processor data.

        Args:
            preserve_ids: If True, keep the list of PDB IDs but reset all other data structures
        """
        # Store current PDB IDs if needed
        current_ids = self.pdb_ids.copy() if preserve_ids else []

        # Reset all data structures
        self.data = None  # Explicitly set to None
        self.dfl = []
        self.chain_dict = {}
        self.pdb_ids = current_ids  # Restore IDs if preserving

    # Dataset management methods
    def list_datasets(self):
        """
        List available structure datasets.
        
        Returns:
            List of dataset information dictionaries
        """
        # Use the standardized dataset API
        if self.dataset_manager is not None:
            return self.dataset_manager.list_datasets()
        else:
            # Fallback to legacy method for backward compatibility
            datasets = super().list_datasets()
            # Filter for structure datasets
            return [d for d in datasets if d.get("type") == "structure" or d.get("type") == "structure_dataset"]
        
    def load_dataset(self, dataset_name):
        """
        Load a predefined dataset of structures.
        
        Args:
            dataset_name: Name of the dataset
            
        Returns:
            Dataset object or None if loading failed
        """
        # Use the standardized dataset API
        if self.dataset_manager is not None:
            dataset = self.dataset_manager.load_dataset(dataset_name)
            if dataset is not None:
                # Extract PDB IDs from dataset content
                self.pdb_ids = dataset.content
                # Load the actual structures
                self.load_structures(self.pdb_ids)
                return dataset
        
        # Legacy fallback approach
        dataset_info = self.get_dataset_info(dataset_name)
        if dataset_info is not None and "pdb_ids" in dataset_info:
            self.pdb_ids = dataset_info["pdb_ids"]
            self.load_structures(self.pdb_ids)
            
            # Create a standardized dataset object if dataset manager is available
            if self.dataset_manager is not None:
                return self.dataset_manager.create_dataset(
                    dataset_id=dataset_name,
                    name=dataset_info.get("name", dataset_name),
                    description=dataset_info.get("description", f"Structure dataset with {len(self.pdb_ids)} PDB IDs"),
                    content=self.pdb_ids,
                    metadata={
                        "converted_from_legacy": True,
                        "original_metadata": dataset_info
                    }
                )
            return None
            
        # Very legacy approach: check datasets.json
        datasets_json_path = os.path.join(self.path_dataset_dir, 'datasets.json')
        try:
            with open(datasets_json_path, 'r') as f:
                datasets = json.load(f)
                
            if dataset_name not in datasets:
                self.logger.error(f"Dataset '{dataset_name}' not found")
                return None
                
            self.pdb_ids = datasets[dataset_name]
            self.load_structures(self.pdb_ids)
            
            # Create a standardized dataset object if dataset manager is available
            if self.dataset_manager is not None:
                return self.dataset_manager.create_dataset(
                    dataset_id=dataset_name,
                    name=dataset_name,
                    description=f"Structure dataset with {len(self.pdb_ids)} PDB IDs",
                    content=self.pdb_ids,
                    metadata={
                        "converted_from_legacy_json": True,
                        "source_file": datasets_json_path
                    }
                )
            return None
        except (FileNotFoundError, json.JSONDecodeError):
            self.logger.error(f"Could not load dataset '{dataset_name}'")
            return None
        
    def load_dataset(self, dataset_id, apply_dtypes=True, debug=False):
        """
        Load a dataset of structure PDB IDs and the corresponding structures.
        
        This method loads the dataset definition, then loads each structure in the dataset.
        Finally, it applies proper data type formatting to ensure all columns have 
        the correct types.
        
        Args:
            dataset_id: Dataset identifier
            apply_dtypes: Whether to apply proper data types to loaded structures
            debug: Whether to print debug information during data type formatting
            
        Returns:
            List of loaded PDB IDs
        """
        # Use the standardized dataset API if available
        if self.dataset_manager is not None:
            dataset = self.dataset_manager.load_dataset(dataset_id)
            if dataset is not None and hasattr(dataset, 'content'):
                pdb_ids = dataset.content
                self.logger.info(f"Loaded dataset '{dataset_id}' with {len(pdb_ids)} structures using dataset manager")
                
                # Load the structures
                self.load_structures(pdb_ids, apply_dtypes=apply_dtypes, debug=debug)
                return pdb_ids
        
        # Legacy approach: check datasets.json
        datasets_json_path = os.path.join(self.path_dataset_dir, 'datasets.json')
        if os.path.isfile(datasets_json_path):
            try:
                with open(datasets_json_path, 'r') as f:
                    datasets = json.load(f)
                
                if dataset_id in datasets:
                    pdb_ids = datasets[dataset_id]
                    self.logger.info(f"Loaded dataset '{dataset_id}' with {len(pdb_ids)} structures from datasets.json")
                    
                    # Load the structures
                    self.load_structures(pdb_ids, apply_dtypes=apply_dtypes, debug=debug)
                    return pdb_ids
                else:
                    self.logger.warning(f"Dataset '{dataset_id}' not found in datasets.json")
            except (FileNotFoundError, json.JSONDecodeError) as e:
                self.logger.error(f"Error loading dataset '{dataset_id}': {str(e)}")
        else:
            self.logger.warning(f"No datasets.json file found at {datasets_json_path}")
        
        # Check registry
        dataset_info = self.get_dataset_info(dataset_id)
        if dataset_info is not None and "pdb_ids" in dataset_info:
            pdb_ids = dataset_info["pdb_ids"]
            self.logger.info(f"Loaded dataset '{dataset_id}' with {len(pdb_ids)} structures from registry")
            
            # Load the structures
            self.load_structures(pdb_ids, apply_dtypes=apply_dtypes, debug=debug)
            return pdb_ids
        
        self.logger.error(f"Could not load dataset '{dataset_id}'")
        return []
            
    def create_dataset(self, 
                    dataset_id: str,
                    name: str,
                    description: str,
                    content: Optional[List[str]] = None,
                    metadata: Optional[Dict[str, Any]] = None) -> Optional[object]:
        """
        Create a new dataset.
        
        Args:
            dataset_id: Unique identifier for the dataset
            name: Human-readable name
            description: Detailed description
            content: List of PDB IDs (defaults to current pdb_ids)
            metadata: Additional metadata
            
        Returns:
            Created Dataset instance or None if not available
        """
        # Use current PDB IDs if content not provided
        if content is None:
            content = self.pdb_ids
            
        # Ensure we have PDB IDs
        if not content:
            self.logger.error("No PDB IDs provided or loaded to create dataset")
            return None
            
        # Use the standardized dataset API if available
        if self.dataset_manager is not None:
            return self.dataset_manager.create_dataset(
                dataset_id=dataset_id,
                name=name,
                description=description,
                content=content,
                metadata=metadata or {}
            )
        else:
            # Legacy fallback approach
            self._register_dataset(
                dataset_id,
                {
                    "type": "structure_dataset",
                    "pdb_ids": content,
                    "name": name,
                    "description": description,
                    "created_at": datetime.now().isoformat(),
                    "structure_count": len(content),
                    **(metadata or {})
                }
            )
            
            # Also update datasets.json for backward compatibility
            datasets_json_path = os.path.join(self.path_dataset_dir, 'datasets.json')
            datasets = {}
            if os.path.isfile(datasets_json_path):
                try:
                    with open(datasets_json_path, 'r') as file:
                        datasets = json.load(file)
                except json.JSONDecodeError:
                    self.logger.warning(f"Could not parse existing datasets file: {datasets_json_path}")
                    
            # Add or update the dataset
            datasets[dataset_id] = content
            
            # Save the updated datasets back to the file
            try:
                with open(datasets_json_path, 'w') as file:
                    json.dump(datasets, file, indent=4)
                self.logger.info(f"Dataset '{dataset_id}' saved successfully")
            except IOError as e:
                self.logger.error(f"Failed to save dataset '{dataset_id}': {str(e)}")
                
            return None  # No Dataset object in legacy mode
            
    def delete_dataset(self, dataset_id: str) -> bool:
        """
        Delete a dataset.
        
        Args:
            dataset_id: Dataset identifier
            
        Returns:
            True if deletion was successful
        """
        # Use the standardized dataset API if available
        if self.dataset_manager is not None:
            return self.dataset_manager.delete_dataset(dataset_id)
            
        # Legacy fallback approach
        if dataset_id in self.dataset_registry:
            # Remove from registry
            del self.dataset_registry[dataset_id]
            self._save_dataset_registry()
            
            # Also remove from datasets.json for backward compatibility
            datasets_json_path = os.path.join(self.path_dataset_dir, 'datasets.json')
            if os.path.isfile(datasets_json_path):
                try:
                    with open(datasets_json_path, 'r') as file:
                        datasets = json.load(file)
                    
                    if dataset_id in datasets:
                        del datasets[dataset_id]
                        
                        with open(datasets_json_path, 'w') as file:
                            json.dump(datasets, file, indent=4)
                        
                        self.logger.info(f"Dataset '{dataset_id}' deleted successfully")
                        return True
                except (json.JSONDecodeError, IOError) as e:
                    self.logger.error(f"Failed to update datasets.json when deleting '{dataset_id}': {str(e)}")
                    return False
            
            return True
        
        return False
    
    # Structure analysis methods
    def get_chains(self, pdb_id):
        """
        Get chain IDs for a structure.
        
        Args:
            pdb_id: PDB ID of the structure
            
        Returns:
            List of chain IDs
        """
        if self.data is None:
            return []
            
        return self.data[self.data['pdb_id'] == pdb_id]['auth_chain_id'].unique().tolist()
        
    def get_ca_coordinates(self, pdb_id, chain_id):
        """
        Get CA (alpha carbon) coordinates for a specific chain.
        
        Args:
            pdb_id: PDB ID of the structure
            chain_id: Chain ID
            
        Returns:
            Array of coordinates
            
        Raises:
            ValueError: If data is None or no CA atoms found
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        df = self.data[(self.data['pdb_id'] == pdb_id) &
                      (self.data['auth_chain_id'] == chain_id) &
                      (self.data['res_atom_name'] == ALPHA_CARBON)]
                      
        if df.empty:
            raise ValueError(f"No CA atoms found for {pdb_id} chain {chain_id}")
            
        return df[['x', 'y', 'z']].to_numpy().astype(np.float64)
        
    def get_sequence(self, pdb_id, chain_id='A'):
        """
        Get amino acid sequence for a chain.
        
        Args:
            pdb_id: PDB ID of the structure
            chain_id: Chain ID
            
        Returns:
            Amino acid sequence as string
            
        Raises:
            ValueError: If data is None
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        seq = self.data[(self.data['pdb_id'] == pdb_id) &
                        (self.data['auth_chain_id'] == chain_id) &
                        (self.data['res_atom_name'] == ALPHA_CARBON)]['res_name1l'].tolist()
                        
        return ''.join(seq)
        
    def get_seq_dict(self, load_file=False, version='v1', chain_type=None):
        """
        Get a dictionary mapping chain IDs to sequences.
        
        Args:
            load_file: Whether to load from a file
            version: Version of the data
            chain_type: Type of chain to filter for
            
        Returns:
            Dictionary mapping chain IDs to sequences
        """
        if load_file:
            return self.load_chain_dict_from_fasta(version)
            
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        chain_dict = {}
        
        # Filter data based on chain_type if specified
        if chain_type is not None and 'chain_type' in self.data.columns:
            backbone = self.data[(self.data['chain_type'] == chain_type) &
                                (self.data['res_atom_name'] == ALPHA_CARBON)]
        else:
            backbone = self.data[(self.data['res_atom_name'] == ALPHA_CARBON) &
                                (self.data['group'] == 'ATOM')]
                                
        for pdb_id in tqdm(self.pdb_ids, desc="Creating sequence dictionary"):
            chains_in_pdb = backbone[backbone['pdb_id'] == pdb_id]['auth_chain_id'].unique()
            for chain_id in chains_in_pdb:
                # Sort by sequence ID to ensure we get residues in order
                chain_data = backbone[(backbone['pdb_id'] == pdb_id) &
                                     (backbone['auth_chain_id'] == chain_id)].sort_values(by='gen_seq_id')
                                     
                seq = chain_data['res_name1l'].tolist()
                ca_atom_ids = chain_data['atom_id'].tolist()
                
                # Store atom IDs for reference
                self.chain_dict_ca_atom_ids[f"{pdb_id}_{chain_id}"] = ca_atom_ids
                
                # Store the sequence
                chain_dict[f"{pdb_id}_{chain_id}"] = ''.join(seq).replace(' ', 'X')
                
        self.chain_dict = chain_dict
        return chain_dict
        
    def load_chain_dict_from_fasta(self, version='v1'):
        """
        Load chain sequences from a FASTA file.
        
        Args:
            version: Version string for the file
            
        Returns:
            Dictionary mapping chain IDs to sequences
        """
        filename = os.path.join(self.path_dataset_dir, f'chain_dict_{version}.fas')
        if not os.path.exists(filename):
            self.logger.warning(f"Chain dictionary file {filename} not found")
            return {}
            
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Create a dictionary from the lines
        chain_dict = {lines[i][1:].strip(): lines[i + 1].strip().upper()
                     for i in range(0, len(lines), 2)}
                     
        self.chain_dict = chain_dict
        return chain_dict
        
    def save_chain_dict_to_fasta(self, version='v1'):
        """
        Save chain sequences to a FASTA file.
        
        Args:
            version: Version string for the file
        """
        if not self.chain_dict:
            self.logger.warning("No chain dictionary to save")
            return
            
        filename = os.path.join(self.path_dataset_dir, f'chain_dict_{version}.fas')
        with open(filename, 'w') as f:
            for key, value in self.chain_dict.items():
                f.write(f'>{key}\n{value}\n')
                
        self.logger.info(f"Saved chain dictionary to {filename}")
        
    def extract_binding_pocket(self, pdb_id, ligand='RET', distance=5.0):
        """
        Extract atoms in the binding pocket around a ligand.
        
        Args:
            pdb_id: PDB ID of the structure
            ligand: Ligand residue name
            distance: Maximum distance to include
            
        Returns:
            DataFrame with atoms in the binding pocket
            
        Raises:
            ValueError: If data is None
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        # Get ligand data
        ligand_data = self.data[(self.data['pdb_id'] == pdb_id) &
                               (self.data['res_name3l'] == ligand)]
                               
        if ligand_data.empty:
            self.logger.warning(f"No ligand '{ligand}' found in {pdb_id}")
            return pd.DataFrame()
            
        # Get chain with ligand
        chain_id = ligand_data['auth_chain_id'].iloc[0]
        
        # Get ligand coordinates
        ligand_coords = ligand_data[['x', 'y', 'z']].values.astype(float)
        
        # Get protein coordinates
        protein_data = self.data[(self.data['pdb_id'] == pdb_id) &
                                (self.data['auth_chain_id'] == chain_id) &
                                (self.data['group'] == 'ATOM')]
                                
        protein_coords = protein_data[['x', 'y', 'z']].values.astype(float)
        
        # Find atoms within distance of ligand
        binding_pocket_indices = []
        
        for ligand_coord in ligand_coords:
            # Calculate distances
            distances = np.sqrt(np.sum((protein_coords - ligand_coord) ** 2, axis=1))
            # Find atoms within distance
            close_indices = np.where(distances <= distance)[0]
            binding_pocket_indices.extend(close_indices)
            
        # Remove duplicates
        binding_pocket_indices = list(set(binding_pocket_indices))
        
        # Return the binding pocket atoms
        if binding_pocket_indices:
            return protein_data.iloc[binding_pocket_indices].copy()
        else:
            return pd.DataFrame()
    
    # Structure filtering methods
    def filter_by_pdb_ids(self, pdb_ids):
        """
        Filter data by PDB IDs.
        
        Args:
            pdb_ids: List of PDB IDs to keep
        """
        if self.data is None:
            return
            
        # Convert input to a set for faster operations
        pdb_ids_set = set(pdb_ids)
        
        # Check if all requested pdb_ids are present
        present_pdb_ids = set(self.data['pdb_id'].unique())
        missing_pdb_ids = pdb_ids_set - present_pdb_ids
        
        if missing_pdb_ids:
            self.logger.warning(f"The following PDB IDs are missing: {missing_pdb_ids}")
            
        # Filter pdb_ids
        self.pdb_ids = list(present_pdb_ids & pdb_ids_set)
        
        # Filter data
        self.data = self.data[self.data['pdb_id'].isin(self.pdb_ids)]
        
        # Update chains
        self.update_pdb_chain_ids()
        
    def get_chain(self, pdb_chain_id):
        """
        Get data for a specific chain.
        
        Args:
            pdb_chain_id: Combined PDB and chain ID (e.g., '1U19_A')
            
        Returns:
            DataFrame with chain data
            
        Raises:
            ValueError: If data is None or format is invalid
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        if '_' not in pdb_chain_id:
            raise ValueError(f"Invalid pdb_chain_id format: {pdb_chain_id}. Expected format: 'pdbid_chainid'")
            
        parts = pdb_chain_id.split('_')
        chain_id = parts[-1]
        pdb_id = '_'.join(parts[:-1])
        
        return self.data[(self.data['pdb_id'] == pdb_id) &
                        (self.data['auth_chain_id'] == chain_id)]
                        
    def get_backbone(self, pdb_chain_id):
        """
        Get backbone atoms for a chain.
        
        Args:
            pdb_chain_id: Combined PDB and chain ID
            
        Returns:
            DataFrame with backbone atoms
            
        Raises:
            ValueError: If data is None or format is invalid
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        if '_' not in pdb_chain_id:
            raise ValueError(f"Invalid pdb_chain_id format: {pdb_chain_id}. Expected format: 'pdbid_chainid'")
            
        parts = pdb_chain_id.split('_')
        chain_id = parts[-1]
        pdb_id = '_'.join(parts[:-1])
        
        return self.data[(self.data['pdb_id'] == pdb_id) &
                        (self.data['auth_chain_id'] == chain_id) &
                        (self.data['res_atom_name'] == ALPHA_CARBON)]
                        
    def apply_grn_interval(self, grn_interval):
        """
        Filter data to include only specified GRNs.
        
        Args:
            grn_interval: List of GRN values to include
            
        Raises:
            ValueError: If data is None
        """
        if self.data is None:
            raise ValueError("No data loaded. Load a structure first.")
            
        if 'grn' not in self.data.columns:
            self.logger.warning("Application of a GRN interval selection requires the 'grn' annotation")
            return
            
        self.data = self.data[self.data['grn'].isin(grn_interval)]
        self.update_pdb_ids()
    
    # Download methods
    @staticmethod
    def download_cif(pdb_id, save_dir=None):
        """
        Download a CIF file from the PDB.
        
        Args:
            pdb_id: PDB ID to download
            save_dir: Directory to save the file (default: uses standard path resolution)
            
        Returns:
            True if download successful, False otherwise
        """
        # Import here to avoid circular import
        from protos.io.paths import get_structure_path
        
        # Store original PDB ID for the URL
        original_pdb_id = pdb_id
        
        # Fix for scientific notation issue - keep original PDB ID even if it looks like a number
        pdb_id = pdb_id.lower()
        # Handle known problematic cases
        if pdb_id == '1.00E+12' or pdb_id == '1.0e+12' or pdb_id == '1e+12':
            pdb_id = '1e12'
            print(f"Fixed scientific notation for PDB ID: {original_pdb_id} -> {pdb_id}")
        
        # Use original PDB ID for the URL to ensure proper download
        url = f'https://files.rcsb.org/download/{original_pdb_id}.cif'
        
        try:
            # Get file path
            if save_dir is not None:
                # Ensure directory exists
                os.makedirs(save_dir, exist_ok=True)
                # Ensure we use the exact PDB ID for the filename, not a converted numeric representation
                file_path = os.path.join(save_dir, f'{original_pdb_id}.cif')
            else:
                # Use standard path resolution with the original PDB ID
                file_path = get_structure_path(original_pdb_id, create_if_missing=True)
            
            # Download the file
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            
            with open(file_path, 'wb') as f:
                f.write(response.content)
                
            return True
        except:
            print(f"HTTP error occurred while downloading {pdb_id}")
            return False
            
    def check_and_download_missing_cifs(self, dataset_name):
        """
        Check for missing CIF files and download them.
        
        Args:
            dataset_name: Name of the dataset to check
        """
        # Import here to avoid circular import
        from protos.io.paths import get_dataset_path
        
        # First try the standardized dataset API
        if self.dataset_manager is not None:
            dataset = self.dataset_manager.load_dataset(dataset_name)
            if dataset is not None and hasattr(dataset, 'content'):
                pdb_ids = dataset.content
            else:
                # Fall back to other methods
                dataset = None
        else:
            dataset = None
        
        # If standardized dataset didn't work, try legacy methods
        if dataset is None:
            # Try dataset registry first
            dataset_info = self.get_dataset_info(dataset_name)
            if dataset_info is not None and "pdb_ids" in dataset_info:
                pdb_ids = dataset_info["pdb_ids"]
            else:
                # Legacy approach: check datasets.json
                datasets_json_path = get_dataset_path(
                    "datasets", 
                    processor_type="structure", 
                    file_extension=".json"
                )
                
                try:
                    with open(datasets_json_path, 'r') as f:
                        datasets = json.load(f)
                    if dataset_name not in datasets:
                        self.logger.error(f"Dataset '{dataset_name}' not found")
                        return
                    pdb_ids = datasets[dataset_name]
                except (FileNotFoundError, json.JSONDecodeError):
                    self.logger.error(f"Could not load dataset '{dataset_name}'")
                    return
                
        # Get existing files
        existing_files = os.listdir(self.path_structure_dir)
        existing_pdb_ids = [file.replace('.cif', '').upper()
                        for file in existing_files if file.lower().endswith('.cif')]
                        
        # Find missing files
        missing_pdb_ids = [pdb_id for pdb_id in pdb_ids if pdb_id.upper() not in existing_pdb_ids]
        
        if not missing_pdb_ids:
            self.logger.info(f"All CIF files for dataset '{dataset_name}' are available")
            return
            
        self.logger.info(f"Downloading {len(missing_pdb_ids)} missing CIF files")
        
        # Download missing files
        for pdb_id in tqdm(missing_pdb_ids, desc=f"Downloading CIF files for {dataset_name}"):
            self.download_cif(pdb_id, save_dir=self.path_structure_dir)

    def to_cif(self, pdb_id: str) -> str:
        """Converts structure data for a given PDB ID to CIF format text."""
        if self.data is None or self.data.empty:
            raise ValueError("No data loaded in the processor.")

        # Filter data for the specific PDB ID
        structure_df = self.data[self.data['pdb_id'] == pdb_id]
        if structure_df.empty:
            raise ValueError(f"No data found for PDB ID '{pdb_id}'.")

        # Convert filtered DataFrame to CIF string
        try:
            cif_content = cif_utils.df_to_cif(structure_df, pdb_id=pdb_id)
            return cif_content
        except Exception as e:
            self.logger.error(f"Error converting DataFrame to CIF for {pdb_id}: {e}")
            raise

    def save_cif(self, pdb_id: str, output_path: str, versioned: bool = False,
                 force_overwrite: bool = False) -> str:
        """Saves structure data for a given PDB ID to a CIF file."""
        if self.data is None or self.data.empty:
            raise ValueError("No data loaded in the processor.")

        # Filter data for the specific PDB ID
        structure_df = self.data[self.data['pdb_id'] == pdb_id]
        if structure_df.empty:
            raise ValueError(f"No data found for PDB ID '{pdb_id}'.")

        # Ensure output directory exists
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

        # Write DataFrame to CIF file using cif_utils
        try:
            final_path = cif_utils.write_cif_file(
                file_path=output_path,
                df=structure_df,
                versioned=versioned,
                force_overwrite=force_overwrite
            )
            self.logger.info(f"Saved structure for {pdb_id} to {final_path}")
            return final_path
        except Exception as e:
            self.logger.error(f"Error writing CIF file for {pdb_id} to {output_path}: {e}")
            raise

    def save_temp_cif(self, pdb_id: str, suffix: Optional[str] = None) -> str:
        """Creates a temporary CIF file for a given PDB ID."""
        if self.data is None or self.data.empty:
            raise ValueError("No data loaded in the processor.")

        # Filter data for the specific PDB ID
        structure_df = self.data[self.data['pdb_id'] == pdb_id]
        if structure_df.empty:
            raise ValueError(f"No data found for PDB ID '{pdb_id}'.")

        # Ensure temp directory exists (should be created in __init__)
        os.makedirs(self.temp_dir, exist_ok=True)

        # Generate unique filename
        timestamp = int(time.time() * 1000)
        filename_suffix = f"_{suffix}" if suffix else ""
        temp_filename = f"{pdb_id}{filename_suffix}_{timestamp}.cif"
        temp_filepath = os.path.join(self.temp_dir, temp_filename)  # Use os.path.join for Path object

        # Write DataFrame to temp CIF file using cif_utils
        try:
            final_path = cif_utils.write_cif_file(
                file_path=str(temp_filepath),  # Convert Path to string if needed by write_cif_file
                df=structure_df,
                versioned=False,  # No versioning for temp files
                force_overwrite=True  # Okay to overwrite temp files if collision (unlikely)
            )
            self.logger.info(f"Saved temporary structure for {pdb_id} to {final_path}")
            return final_path
        except Exception as e:
            self.logger.error(f"Error writing temporary CIF file for {pdb_id}: {e}")
            raise

    def cleanup_temp_files(self, older_than_hours: int = 1) -> None:
        """Removes temporary CIF files older than the specified hours."""
        if not os.path.exists(self.temp_dir):
            self.logger.debug(f"Temporary directory {self.temp_dir} does not exist. Skipping cleanup.")
            return

        now = time.time()
        cutoff = now - (older_than_hours * 3600)
        cleaned_count = 0

        try:
            for filename in os.listdir(self.temp_dir):
                file_path = os.path.join(self.temp_dir, filename)
                if os.path.isfile(file_path) and filename.lower().endswith('.cif'):
                    try:
                        file_mod_time = os.path.getmtime(file_path)
                        if file_mod_time < cutoff:
                            os.remove(file_path)
                            self.logger.debug(f"Removed old temporary file: {file_path}")
                            cleaned_count += 1
                    except OSError as e:
                        self.logger.warning(f"Could not remove temporary file {file_path}: {e}")
        except OSError as e:
            self.logger.error(f"Error listing temporary directory {self.temp_dir}: {e}")

        self.logger.info(f"Cleaned up {cleaned_count} old temporary files from {self.temp_dir}")

    def export_structures(self, pdb_ids: List[str], output_dir: str, validate: bool = False) -> List[str]:
        """Exports multiple structures to CIF files in the specified directory."""
        if self.data is None or self.data.empty:
            raise ValueError("No data loaded in the processor.")

        os.makedirs(output_dir, exist_ok=True)
        exported_files = []
        skipped_ids = []

        for pdb_id in tqdm(pdb_ids, desc="Exporting Structures"):
            try:
                # Filter data for the specific PDB ID
                structure_df = self.data[self.data['pdb_id'] == pdb_id]
                if structure_df.empty:
                    self.logger.warning(f"No data found for PDB ID '{pdb_id}'. Skipping export.")
                    skipped_ids.append(pdb_id)
                    continue

                # Optional validation
                if validate:
                    validation_results = cif_utils.validate_cif_data(structure_df)
                    if not validation_results["is_valid"]:
                        self.logger.warning(
                            f"Validation failed for {pdb_id}. Issues: {validation_results['issues']}. Skipping export.")
                        skipped_ids.append(pdb_id)
                        continue  # Or decide to export anyway

                # Define output path
                output_path = os.path.join(output_dir, f"{pdb_id}.cif")

                # Write file using save_cif logic (or directly call write_cif_file)
                # Using save_cif might be redundant, calling write_cif_file is more direct
                final_path = cif_utils.write_cif_file(
                    file_path=output_path,
                    df=structure_df,
                    versioned=False,  # Typically no versioning during bulk export
                    force_overwrite=True  # Overwrite if file exists from previous run
                )
                exported_files.append(final_path)

            except Exception as e:
                self.logger.error(f"Failed to export structure {pdb_id}: {e}")
                skipped_ids.append(pdb_id)

        self.logger.info(f"Exported {len(exported_files)} structures to {output_dir}.")
        if skipped_ids:
            self.logger.warning(f"Skipped {len(skipped_ids)} structures: {skipped_ids}")

        return exported_files

    def filter_data_flexibly(self, filters: Dict[str, Any], inplace: bool = True) -> Optional[pd.DataFrame]:
        """
        Filters the internal data DataFrame based on flexible criteria.

        Args:
            filters: A dictionary where keys are column names (optionally with operators)
                     and values are the conditions.
                     Supported operators (append to column name with __):
                       - __eq: equals (default if no operator)
                       - __ne: not equals
                       - __gt: greater than
                       - __lt: less than
                       - __ge: greater than or equal to
                       - __le: less than or equal to
                       - __is_in: value must be in the provided list/set
                       - __not_in: value must not be in the provided list/set
                       - __startswith: string starts with value
                       - __endswith: string ends with value
                       - __contains: string contains value
                       - __isna: value is NaN/None/NaT
                       - __notna: value is not NaN/None/NaT
            inplace: If True, modifies the internal self.data DataFrame.
                     If False, returns a new filtered DataFrame without modifying internal state.

        Returns:
            If inplace=False, returns the filtered DataFrame.
            If inplace=True, returns None.

        Raises:
            ValueError: If a filter column does not exist or an operator is invalid.
            TypeError: If filter value type is incompatible with the operator/column.
        """
        if self.data is None or self.data.empty:
            self.logger.warning("No data loaded to filter.")
            return pd.DataFrame() if not inplace else None

        current_df = self.data.copy() # Work on a copy initially

        # Define operators and their corresponding pandas methods
        operator_map = {
            'eq': lambda series, val: series == val,
            'ne': lambda series, val: series != val,
            'gt': lambda series, val: pd.to_numeric(series, errors='coerce') > val,
            'lt': lambda series, val: pd.to_numeric(series, errors='coerce') < val,
            'ge': lambda series, val: pd.to_numeric(series, errors='coerce') >= val,
            'le': lambda series, val: pd.to_numeric(series, errors='coerce') <= val,
            'is_in': lambda series, val: series.isin(val),
            'not_in': lambda series, val: ~series.isin(val),
            'startswith': lambda series, val: series.astype(str).str.startswith(val, na=False),
            'endswith': lambda series, val: series.astype(str).str.endswith(val, na=False),
            'contains': lambda series, val: series.astype(str).str.contains(val, na=False, regex=False),
            'isna': lambda series, val: series.isna(), # Value is ignored for isna/notna
            'notna': lambda series, val: series.notna(), # Value is ignored for isna/notna
        }

        combined_mask = pd.Series(True, index=current_df.index) # Start with all True

        for filter_key, filter_value in filters.items():
            col_name = filter_key
            operator = 'eq' # Default operator

            # Parse column name and operator
            if '__' in filter_key:
                parts = filter_key.split('__', 1)
                col_name = parts[0]
                op_suffix = parts[1].lower()
                if op_suffix in operator_map:
                    operator = op_suffix
                else:
                    raise ValueError(f"Invalid filter operator: __{op_suffix}")

            # Check if column exists
            if col_name not in current_df.columns:
                raise ValueError(f"Filter column '{col_name}' not found in data.")

            # Get the pandas operation function
            op_func = operator_map[operator]

            # Apply the filter to create a mask for the current condition
            try:
                # Handle special cases for operators that don't use the value directly
                if operator in ['isna', 'notna']:
                    condition_mask = op_func(current_df[col_name], None)
                else:
                    condition_mask = op_func(current_df[col_name], filter_value)

                # Combine with the overall mask using AND logic
                combined_mask &= condition_mask

            except TypeError as e:
                raise TypeError(f"Type error applying filter '{filter_key}' with value '{filter_value}' "
                                f"on column '{col_name}' (dtype: {current_df[col_name].dtype}): {e}")
            except Exception as e:
                # Catch other potential errors like regex errors for 'contains' if regex=True
                 raise ValueError(f"Error applying filter '{filter_key}' with value '{filter_value}': {e}")

        # Apply the final mask
        filtered_df = current_df[combined_mask]

        if inplace:
            original_size = len(self.data)
            self.data = filtered_df
            # Important: Update other state if filtering occurred
            self.update_pdb_ids() # Refresh the list of PDB IDs present
            self.dfl = [self.data[self.data['pdb_id'] == pdb_id] for pdb_id in self.pdb_ids] # Rebuild dfl
            # Consider if chain_dict etc. also need refreshing based on usage
            self.update_chain_data()
            self.logger.info(f"Filtered data inplace. Size changed from {original_size} to {len(self.data)} rows.")
            return None
        else:
            return filtered_df

    def add_ligand(self,
                   target_pdb_id: str,
                   ligand_df: pd.DataFrame,
                   ligand_chain_id: str = 'L',
                   ligand_res_seq_id: int = 9001):
        """
        Adds ligand atom data to a specific existing structure in self.data.

        NOTE: This performs data merging, not geometric docking. The ligand_df
              coordinates are assumed to be in the desired final position relative
              to the target protein.

        Args:
            target_pdb_id: The PDB ID of the structure to add the ligand to.
            ligand_df: A DataFrame containing the ligand atoms, formatted similarly
                       to the main structure data (must include 'x', 'y', 'z',
                       'atom_name', 'res_name'). Coordinates should be pre-positioned.
            ligand_chain_id: The chain ID to assign to the ligand (default 'L').
                             Ensure this doesn't clash with existing protein chains
                             if separation is desired.
            ligand_res_seq_id: The residue sequence number to assign to the ligand
                               (default 9001). Ensure this is unique within the
                               target structure's ligand chain.

        Raises:
            ValueError: If target PDB ID is not loaded, ligand DataFrame is invalid,
                        or required columns are missing.
        """
        if self.data is None or self.data.empty:
            raise ValueError("Processor data (self.data) is not loaded.")

        if target_pdb_id not in self.pdb_ids:
            raise ValueError(f"Target PDB ID '{target_pdb_id}' not found in loaded structures.")

        # --- Validate Ligand DataFrame ---
        required_ligand_cols = ['x', 'y', 'z', 'atom_name', 'res_name']
        missing_cols = [col for col in required_ligand_cols if col not in ligand_df.columns]
        if missing_cols:
            raise ValueError(f"Ligand DataFrame is missing required columns: {missing_cols}")

        if ligand_df.empty:
             raise ValueError("Ligand DataFrame cannot be empty.")

        # Make a copy to avoid modifying the original ligand DF
        lig_df = ligand_df.copy()

        # --- Prepare Ligand DataFrame for Merging ---
        num_lig_atoms = len(lig_df)
        lig_res_name = lig_df['res_name'].iloc[0] # Assume single residue name for ligand

        # Assign standard columns
        lig_df['pdb_id'] = target_pdb_id
        lig_df['group'] = 'HETATM'
        lig_df['auth_chain_id'] = ligand_chain_id
        lig_df['auth_seq_id'] = ligand_res_seq_id
        lig_df['res_id'] = f"{lig_res_name}_{ligand_res_seq_id}_{ligand_chain_id}"
        lig_df['res_name3l'] = lig_res_name # Assume res_name is 3-letter code
        # Handle optional columns with defaults if missing
        if 'element' not in lig_df:
             lig_df['element'] = lig_df['atom_name'].str[0].str.upper() # Basic guess
        if 'occupancy' not in lig_df:
             lig_df['occupancy'] = 1.0
        if 'b_factor' not in lig_df:
             lig_df['b_factor'] = 20.0 # Default B-factor
        # Add other optional columns as needed with defaults (e.g., insertion=None)
        for col in ['insertion', 'alt_id', 'charge', 'model_num', 'auth_comp_id', 'auth_atom_name']:
             if col not in lig_df:
                 lig_df[col] = None # Or appropriate default, e.g., model_num=1

        # --- Renumber Atom IDs ---
        # Find max atom_id in the target protein structure *before* adding ligand
        target_protein_df = self.data[self.data['pdb_id'] == target_pdb_id]
        max_atom_id = target_protein_df['atom_id'].max()
        if pd.isna(max_atom_id): # Handle case where target might be empty initially
            max_atom_id = 0

        # Assign new sequential atom IDs to ligand atoms
        lig_df['atom_id'] = range(int(max_atom_id) + 1, int(max_atom_id) + 1 + num_lig_atoms)
        # Also add gen_seq_id if used (can just match atom_id here)
        if 'gen_seq_id' in self.data.columns:
             lig_df['gen_seq_id'] = lig_df['atom_id']

        # --- Merge DataFrames ---
        # Separate data for other PDB IDs
        other_data = self.data[self.data['pdb_id'] != target_pdb_id]

        # Concatenate the target protein and the prepared ligand data
        updated_target_df = pd.concat([target_protein_df, lig_df], ignore_index=True)

        # Combine with data from other PDB IDs
        self.data = pd.concat([other_data, updated_target_df], ignore_index=True)

        # --- Update Processor State ---
        self.reset_index() # Important to fix multi-index after concatenation
        self.dfl = [self.data[self.data['pdb_id'] == pdb_id] for pdb_id in self.pdb_ids] # Rebuild dfl
        # Consider if other state like chain_dict needs update
        self.update_chain_data()

        self.logger.info(f"Added ligand '{lig_res_name}' (Chain: {ligand_chain_id}, ResID: {ligand_res_seq_id}) "
                         f"with {num_lig_atoms} atoms to structure '{target_pdb_id}'.")
