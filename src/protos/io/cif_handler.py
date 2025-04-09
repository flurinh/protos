"""
CIF file format handler for Protos.

This module provides a handler for reading and writing mmCIF/PDBx files.
"""

import os
import logging
import pandas as pd
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

# Import FormatHandler class definition to avoid circular import
# We only need the class definition, not the entire module
import sys
if 'protos.io.formats' in sys.modules:
    FormatHandler = sys.modules['protos.io.formats'].FormatHandler
else:
    # Define a minimal version if formats not imported yet
    class FormatHandler:
        def __init__(self, logger=None):
            import logging
            self.logger = logger or logging.getLogger(self.__class__.__name__)

        def read(self, file_path, **kwargs):
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"File not found: {file_path}")
            return self._read_impl(file_path, **kwargs)

        def write(self, file_path, data, **kwargs):
            os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)
            self._validate(data)
            self._write_impl(file_path, data, **kwargs)

        def _validate(self, data):
            pass

        def _read_impl(self, file_path, **kwargs):
            raise NotImplementedError("Subclasses must implement _read_impl")

        def _write_impl(self, file_path, data, **kwargs):
            raise NotImplementedError("Subclasses must implement _write_impl")

# Import cif_utils functions
from .cif_utils import (
    read_cif_file,
    write_cif_file,
    validate_cif_data,
    fix_cif_data,
    merge_structures,
    extract_bioassembly,
    REQUIRED_COLUMNS,
)

class CifHandler(FormatHandler):
    """
    Handler for mmCIF/PDBx files.

    This handler provides methods for reading, writing, and validating
    mmCIF/PDBx structure files using pandas DataFrames as the internal
    data representation.
    """

    def _read_impl(self, file_path: str, **kwargs) -> pd.DataFrame:
        """
        Read a mmCIF file into a DataFrame.

        This implementation uses the read_cif_file function from cif_utils
        to parse the CIF file into a standardized DataFrame format.

        Args:
            file_path: Path to the CIF file to read
            **kwargs: Additional arguments to pass to the parser

        Returns:
            DataFrame containing the atomic structure data

        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If parsing fails
        """
        try:
            # Ensure file path exists before trying to read
            if not os.path.exists(file_path):
                # Check for alternate path format
                alt_path = file_path.replace('/', os.path.sep)
                if os.path.exists(alt_path):
                    file_path = alt_path
                else:
                    alt_path = file_path.replace('\\', '/')
                    if os.path.exists(alt_path):
                        file_path = alt_path
                    else:
                        raise FileNotFoundError(f"File not found: {file_path}")

            return read_cif_file(file_path)
        except Exception as e:
            self.logger.error(f"Error reading CIF file {file_path}: {e}")
            raise

    def _write_impl(self, file_path: str, data: pd.DataFrame, **kwargs) -> str:
        """
        Write a DataFrame to a mmCIF file.

        This implementation uses the write_cif_file function from cif_utils
        to convert the DataFrame to a properly formatted CIF file.

        Args:
            file_path: Path to save the CIF file
            data: DataFrame with atomic structure data
            **kwargs: Additional arguments including:
                - versioned: Whether to add version numbering (default: False)
                - force_overwrite: Whether to overwrite existing files (default: False)

        Returns:
            Path to the written file

        Raises:
            ValueError: If required data is missing
            FileExistsError: If file exists and force_overwrite=False
        """
        try:
            # Extract relevant kwargs
            versioned = kwargs.get('versioned', False)
            force_overwrite = kwargs.get('force_overwrite', False)

            # Normalize path using os.path.normpath to ensure consistent path format
            normalized_path = os.path.normpath(file_path)

            # Use the cif_utils function to write the file
            written_path = write_cif_file(
                file_path=normalized_path,
                df=data,
                versioned=versioned,
                force_overwrite=force_overwrite
            )

            # Return the exact path of the written file (important for versioned files)
            return written_path

        except Exception as e:
            self.logger.error(f"Error writing CIF file {file_path}: {e}")
            raise

    def write_with_versioning(self, file_path: str, data: pd.DataFrame,
                             versioned: bool = True, force_overwrite: bool = False,
                             **kwargs) -> str:
        """
        Write a DataFrame to a mmCIF file with versioning support.

        Args:
            file_path: Path to the file
            data: DataFrame to write
            versioned: Whether to add version numbering (_v1, _v2, etc.)
            force_overwrite: Whether to allow overwriting existing files
            **kwargs: Additional parameters for _write_impl

        Returns:
            Path to the written file

        Raises:
            FileExistsError: If file exists and force_overwrite=False
        """
        self._validate(data)
        kwargs.update({
            'versioned': versioned,
            'force_overwrite': force_overwrite
        })
        return self._write_impl(file_path, data, **kwargs)
    
    def validate_data(self, data: pd.DataFrame) -> dict:
        """
        Validate structure data before writing to CIF.
        
        Args:
            data: DataFrame to validate
            
        Returns:
            Dictionary with validation results
        """
        return validate_cif_data(data)
    
    def fix_data(self, data: pd.DataFrame) -> pd.DataFrame:
        """
        Fix common issues in CIF data.
        
        Args:
            data: DataFrame with atomic structure data
            
        Returns:
            Corrected DataFrame
        """
        return fix_cif_data(data)
    
    def merge_structures(self, dfs: List[pd.DataFrame], 
                        chain_mapping: Optional[Dict[str, str]] = None) -> pd.DataFrame:
        """
        Merge multiple structure DataFrames into a single structure.
        
        Args:
            dfs: List of DataFrames to merge
            chain_mapping: Optional mapping to reassign chain IDs
            
        Returns:
            DataFrame with merged structure data
        """
        return merge_structures(dfs, chain_mapping)
    
    def extract_bioassembly(self, data: pd.DataFrame, 
                           assembly_id: int = 1,
                           include_chains: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Extract biological assembly from a structure DataFrame.
        
        Args:
            data: DataFrame with atomic structure data
            assembly_id: ID of the biological assembly to extract
            include_chains: Optional list of chains to include
            
        Returns:
            DataFrame with biological assembly structure
        """
        return extract_bioassembly(data, assembly_id, include_chains)
    
    def _validate(self, data: Any) -> None:
        """
        Validate that data is a DataFrame with required columns.
        
        Args:
            data: Data to validate
            
        Raises:
            ValueError: If data is not a DataFrame or missing required columns
        """
        if not isinstance(data, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame for CIF format")
        
        missing_cols = [col for col in REQUIRED_COLUMNS if col not in data.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns for CIF format: {missing_cols}")