"""
GRN Processor with BaseProcessor integration.

This module provides a GRNProcessor class that extends BaseProcessor
to provide standardized data management capabilities for GRN tables.
"""

# Import updated GRN utilities with improved loop handling
from protos.processing.schema.grn_utils_updated import (
    parse_grn_str2float,
    parse_grn_float2str,
    normalize_grn_format,
    validate_grn_string,
    sort_grns
)
# Keep necessary imports from the original module for backwards compatibility
from protos.processing.grn.grn_utils import (
    get_grn_interval,
    get_seq,
    get_annot_seq,
    map_grn_to_color,
    GRNConfigManager,
    sort_grns_str
)
import os
import pandas as pd
import plotly.graph_objs as go
from protos.core.base_processor import BaseProcessor


class GRNBaseProcessor(BaseProcessor):
    """
    Processor for Generic Residue Numbering (GRN) data.
    
    Handles loading, saving, and processing of GRN tables which map
    sequence positions to standardized numbering schemes.
    """
    
    def __init__(self,
                 name="grn_processor",
                 dataset=None,
                 data_root=None,
                 processor_data_dir="grn",
                 path=None,
                 preload=True):
        """
        Initialize the GRN processor.
        
        Args:
            name: Processor instance name
            dataset: Dataset ID to load (or list of datasets to merge)
            data_root: Root directory for all data
            processor_data_dir: Subdirectory for GRN data
            path: Legacy path parameter (deprecated, use data_root instead)
            preload: Whether to load the dataset on initialization
        """
        # Handle legacy path parameter
        if path is not None:
            import logging
            logging.warning("The 'path' parameter is deprecated, use 'data_root' and 'processor_data_dir' instead")
            # Extract data_root and processor_data_dir from path
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
        
        # Initialize GRN-specific attributes
        self.dataset = None
        self.ids = []
        self.grns = []
        self.features = pd.DataFrame()
        self.maps = {}
        self.map = pd.DataFrame()
        
        # Track notation format and mappings
        self.using_dot_notation = False
        self.dot_to_x = {}  # Maps dot notation (3.50) to x notation (3x50)
        self.x_to_dot = {}  # Maps x notation (3x50) to dot notation (3.50)
        
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
            data_to_save = self.data.copy()
            
            # Normalize GRN formats if requested
            if normalize_formats:
                # Track normalization
                normalized_count = 0
                invalid_grns = []
                
                # Get column names excluding the index
                original_columns = data_to_save.columns.tolist()
                normalized_columns = {}
                
                for col in original_columns:
                    # Skip metadata columns
                    if col in ['family', 'species', 'name', 'grn_system', 'id']:
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
                
            # Reset index if it's not named
            if data_to_save.index.name != "protein_id":
                data_to_save = data_to_save.reset_index()
                # Rename index column if it exists and is not named
                if "index" in data_to_save.columns:
                    data_to_save = data_to_save.rename(columns={"index": "protein_id"})
                elif self.data.index.name is None and "uen" not in data_to_save.columns:
                    # Handle older GRN tables where the index might be unnamed
                    first_col = data_to_save.columns[0]
                    data_to_save = data_to_save.rename(columns={first_col: "protein_id"})
        else:
            self.logger.warning("No data to save")
            data_to_save = pd.DataFrame()
            
        # Save using BaseProcessor's save_data method
        file_path = self.save_data(
            self.dataset, 
            data_to_save, 
            file_format="csv", 
            index=False,
            **kwargs
        )
        
        self.logger.info(f"Saved GRN table '{self.dataset}' to {file_path}")
        return file_path
        
    def load_grn_table(self, dataset_id=None, low_memory=False, remove_duplicates=True, normalize_formats=True, **kwargs):
        """
        Load a GRN table from a dataset.
        
        Args:
            dataset_id: Dataset identifier (uses current dataset if None)
            low_memory: Whether to use pandas low_memory option
            remove_duplicates: Whether to remove duplicate protein IDs
            normalize_formats: Whether to normalize GRN formats (e.g., convert legacy loop formats)
            **kwargs: Additional format-specific loading parameters
            
        Returns:
            DataFrame containing the GRN table
            
        Raises:
            FileNotFoundError: If the dataset doesn't exist
        """
        # Update dataset if provided
        if dataset_id is not None:
            self.dataset = dataset_id
        elif self.dataset is None:
            raise ValueError("No dataset specified.")
            
        # Load using BaseProcessor's load_data method
        self.logger.info(f"Loading GRN table '{self.dataset}'")
        df = self.load_data(
            self.dataset, 
            file_format="csv", 
            low_memory=low_memory, 
            **kwargs
        )
        
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
        self.data = df
        
        # Remove duplicates if requested
        if remove_duplicates:
            self.remove_duplicate_ids()
            
        # Update processor state
        self.ids = self.data.index.tolist()
        self.grns = self.data.columns.tolist()
        
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
                    if col in ['family', 'species', 'name', 'grn_system', 'id']:
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
        
        # Detect which notation is used in the data
        dot_cols = [col for col in original_columns if '.' in col and not ('n.' in col or 'c.' in col)]
        x_cols = [col for col in original_columns if 'x' in col and not ('n.' in col or 'c.' in col)]
        self.using_dot_notation = len(dot_cols) > len(x_cols)
        
        # Create normalized grns list - preserve the original notation for internal use
        # but keep track of alternate forms
        normalized_grns = []
        for grn in self.grns:
            if 'n.' in grn or 'c.' in grn:
                normalized_grns.append(grn)
            # Handle loop GRNs (AB.CCC format)
            elif '.' in grn and len(grn.split('.')[0]) == 2 and len(grn.split('.')[1]) == 3:
                normalized_grns.append(grn)
            elif self.using_dot_notation and 'x' in grn:
                # Convert x notation to dot notation if data uses dot
                dot_form = self.x_to_dot.get(grn)
                if dot_form and dot_form in original_columns:
                    normalized_grns.append(dot_form)
                else:
                    normalized_grns.append(grn)
            elif not self.using_dot_notation and '.' in grn:
                # Convert dot notation to x notation if data uses x
                x_form = self.dot_to_x.get(grn)
                if x_form and x_form in original_columns:
                    normalized_grns.append(x_form)
                else:
                    normalized_grns.append(grn)
            else:
                normalized_grns.append(grn)
        
        # Update the GRNs list with normalized forms
        self.grns = normalized_grns
        
        # Use the normalized GRNs to select columns from the data
        try:
            self.data = self.data[self.grns]
        except KeyError as e:
            self.logger.error(f"Failed to match column names: {e}")
            raise KeyError(f"Cannot match GRN columns. Data columns: {original_columns}, Requested: {self.grns}")
            
        # Log notation format
        self.logger.info(f"Using {'dot' if self.using_dot_notation else 'x'} notation for GRN positions")
            
        self.data.fillna('-', inplace=True)
        
        # Add to dataset registry with additional metadata
        self._register_dataset(
            self.dataset,
            {
                "type": "grn_table",
                "protein_count": len(self.ids),
                "grn_count": len(self.grns),
                "first_grn": self.grns[0] if self.grns else None,
                "last_grn": self.grns[-1] if self.grns else None,
                "invalid_grns": invalid_grns if 'invalid_grns' in locals() else []
            },
            None
        )
        
        return self.data
    
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
            # Load each dataset
            grn_table = self.load_grn_table(dataset_id=dataset, remove_duplicates=remove_duplicates)
            tables.append(grn_table)
            
            # Count notation formats
            dot_cols = [col for col in grn_table.columns if '.' in col and not ('n.' in col or 'c.' in col)]
            x_cols = [col for col in grn_table.columns if 'x' in col and not ('n.' in col or 'c.' in col)]
            
            dot_count += len(dot_cols)
            x_count += len(x_cols)
            
        if not tables:
            return pd.DataFrame()
            
        # Determine which notation to use for the merged table
        self.using_dot_notation = dot_count >= x_count
        
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
        
        # Standardize column names to the preferred notation
        std_columns = []
        for col in merged_table.columns:
            if 'n.' in col or 'c.' in col:
                std_columns.append(col)
            elif self.using_dot_notation and 'x' in col:
                std_columns.append(self.x_to_dot.get(col, col))
            elif not self.using_dot_notation and '.' in col:
                std_columns.append(self.dot_to_x.get(col, col))
            else:
                std_columns.append(col)
        
        # Rename columns to the standardized format
        merged_table.columns = std_columns
        
        # Sort columns by GRN position
        self.grns = sort_grns_str(merged_table.columns.tolist())
        
        # Handle notation conversion for sorted columns
        sorted_cols_in_data_notation = []
        for grn in self.grns:
            # Keep special notations as is
            if 'n.' in grn or 'c.' in grn:
                sorted_cols_in_data_notation.append(grn)
            # Convert to the data's preferred notation
            elif self.using_dot_notation and 'x' in grn:
                sorted_cols_in_data_notation.append(self.x_to_dot.get(grn, grn))
            elif not self.using_dot_notation and '.' in grn:
                sorted_cols_in_data_notation.append(self.dot_to_x.get(grn, grn))
            else:
                sorted_cols_in_data_notation.append(grn)
        
        # Reorder columns using the data's notation
        valid_cols = [col for col in sorted_cols_in_data_notation if col in merged_table.columns]
        self.data = merged_table[valid_cols]
        self.data.fillna('-', inplace=True)
        self.ids = self.data.index.tolist()
        self.grns = valid_cols
        
        # Update dataset info
        self.dataset = 'merged_' + '_'.join(datasets)
        
        # Register the merged dataset
        self._register_dataset(
            self.dataset,
            {
                "type": "grn_table",
                "merged_from": datasets,
                "protein_count": len(self.ids),
                "grn_count": len(self.grns),
                "notation": "dot" if self.using_dot_notation else "x",
                "first_grn": self.grns[0] if self.grns else None,
                "last_grn": self.grns[-1] if self.grns else None,
            },
            None
        )
        
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
        # Convert the interval to the notation used in the data
        mapped_interval = []
        for grn in grn_interval:
            # No conversion needed for n. and c. notations
            if 'n.' in grn or 'c.' in grn:
                if grn in self.data.columns:
                    mapped_interval.append(grn)
                continue
                
            # Handle TM region notation conversion
            if self.using_dot_notation and 'x' in grn:
                # Need to convert x notation to dot notation
                dot_form = self.x_to_dot.get(grn)
                if dot_form in self.data.columns:
                    mapped_interval.append(dot_form)
                elif grn in self.data.columns:
                    mapped_interval.append(grn)
            elif not self.using_dot_notation and '.' in grn:
                # Need to convert dot notation to x notation
                x_form = self.dot_to_x.get(grn)
                if x_form in self.data.columns:
                    mapped_interval.append(x_form)
                elif grn in self.data.columns:
                    mapped_interval.append(grn)
            elif grn in self.data.columns:
                # Already in the right notation
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
                    if self.using_dot_notation:
                        # Create dot notation for TM regions
                        x_notation = parse_grn_float2str(float_val)
                        tm, pos = x_notation.split('x')
                        dot_notation = f"{tm}.{pos}"
                        cols_sorted.append(dot_notation)
                    else:
                        # Use standard x notation
                        cols_sorted.append(parse_grn_float2str(float_val))
            
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
        
        # Build both dot and x notation versions of the interval
        dot_interval = []
        x_interval = []
        
        for grn in grn_interval:
            # Handle special notations
            if 'n.' in grn or 'c.' in grn:
                dot_interval.append(grn)
                x_interval.append(grn)
            # Handle TM regions
            elif 'x' in grn:
                x_interval.append(grn)
                dot_form = self.x_to_dot.get(grn)
                if dot_form:
                    dot_interval.append(dot_form)
                else:
                    # Create dot notation if not in mapping
                    tm, pos = grn.split('x')
                    dot_interval.append(f"{tm}.{pos}")
            elif '.' in grn:
                dot_interval.append(grn)
                x_form = self.dot_to_x.get(grn)
                if x_form:
                    x_interval.append(x_form)
                else:
                    # Create x notation if not in mapping
                    tm, pos = grn.split('.')
                    x_interval.append(f"{tm}x{pos}")
        
        # Apply to each map using the appropriate notation
        for map_key in maps:
            map_cols = self.maps[map_key].columns.tolist()
            
            # Determine which notation this map uses
            map_using_dot = len([col for col in map_cols 
                               if '.' in col and not ('n.' in col or 'c.' in col)]) > 0
            
            # Select the appropriate interval
            map_interval = dot_interval if map_using_dot else x_interval
            
            # Filter to columns that exist in this map
            valid_interval = [col for col in map_interval if col in map_cols]
            
            # Apply if we have valid columns
            if valid_interval:
                self.maps[map_key] = self.maps[map_key][valid_interval]
                self.logger.info(f"Applied interval to map '{map_key}': {len(valid_interval)} GRNs")
            else:
                self.logger.warning(f"No valid GRNs to apply to map '{map_key}'")
    
    # Additional methods for maps, feature population, etc. can be integrated similarly
    # by updating them to use the BaseProcessor functionality where appropriate