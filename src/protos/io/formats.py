"""
Standardized file format handlers for Protos.

This module provides consistent interfaces for reading and writing
different file formats used in the Protos framework. All format handlers
follow a common interface and implement validation based on the
specifications defined in FILE_FORMATS.md.
"""

import os
import json
import pickle
import tempfile
import shutil
import uuid
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union, Tuple, BinaryIO, TextIO, ContextManager
from contextlib import contextmanager
import logging

class FormatHandler:
    """Base class for file format handlers."""
    
    def __init__(self, logger=None):
        """
        Initialize the format handler.
        
        Args:
            logger: Optional logger instance
        """
        self.logger = logger or logging.getLogger(self.__class__.__name__)
    
    def read(self, file_path: str, **kwargs) -> Any:
        """
        Read data from a file.
        
        Args:
            file_path: Path to the file
            **kwargs: Format-specific parameters
            
        Returns:
            Loaded data
            
        Raises:
            FileNotFoundError: If file doesn't exist
            ValueError: If file format is invalid
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        
        return self._read_impl(file_path, **kwargs)
    
    def write(self, file_path: str, data: Any, **kwargs) -> None:
        """
        Write data to a file.
        
        Args:
            file_path: Path to the file
            data: Data to write
            **kwargs: Format-specific parameters
            
        Raises:
            ValueError: If data is invalid for the format
        """
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(os.path.abspath(file_path)), exist_ok=True)
        
        # Validate data before writing
        self._validate(data)
        
        # Write data
        self._write_impl(file_path, data, **kwargs)
    
    def _validate(self, data: Any) -> None:
        """
        Validate data against format specifications.
        
        Args:
            data: Data to validate
            
        Raises:
            ValueError: If data is invalid
        """
        pass  # Base implementation does no validation
    
    def _read_impl(self, file_path: str, **kwargs) -> Any:
        """
        Implementation of file reading.
        
        Args:
            file_path: Path to the file
            **kwargs: Format-specific parameters
            
        Returns:
            Loaded data
        """
        raise NotImplementedError("Subclasses must implement _read_impl")
    
    def _write_impl(self, file_path: str, data: Any, **kwargs) -> None:
        """
        Implementation of file writing.
        
        Args:
            file_path: Path to the file
            data: Data to write
            **kwargs: Format-specific parameters
        """
        raise NotImplementedError("Subclasses must implement _write_impl")


class CSVHandler(FormatHandler):
    """Handler for CSV files."""
    
    def _read_impl(self, file_path: str, **kwargs) -> pd.DataFrame:
        """Read a CSV file into a DataFrame."""
        return pd.read_csv(file_path, **kwargs)
    
    def _write_impl(self, file_path: str, data: pd.DataFrame, **kwargs) -> None:
        """Write a DataFrame to a CSV file."""
        data.to_csv(file_path, **kwargs)
    
    def _validate(self, data: Any) -> None:
        """Validate that data is a DataFrame."""
        if not isinstance(data, pd.DataFrame):
            raise ValueError("Data must be a pandas DataFrame for CSV format")


class PickleHandler(FormatHandler):
    """Handler for pickle files."""
    
    def _read_impl(self, file_path: str, **kwargs) -> Any:
        """Read a pickle file."""
        with open(file_path, 'rb') as f:
            return pickle.load(f)
    
    def _write_impl(self, file_path: str, data: Any, **kwargs) -> None:
        """Write data to a pickle file."""
        protocol = kwargs.get('protocol', pickle.HIGHEST_PROTOCOL)
        with open(file_path, 'wb') as f:
            pickle.dump(data, f, protocol=protocol)


class JSONHandler(FormatHandler):
    """Handler for JSON files."""
    
    def _read_impl(self, file_path: str, **kwargs) -> Union[Dict, List]:
        """Read a JSON file."""
        with open(file_path, 'r') as f:
            return json.load(f)
    
    def _write_impl(self, file_path: str, data: Union[Dict, List], **kwargs) -> None:
        """Write data to a JSON file."""
        indent = kwargs.get('indent', 2)
        default = kwargs.get('default', str)
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=indent, default=default)
    
    def _validate(self, data: Any) -> None:
        """Validate that data is JSON-serializable."""
        try:
            json.dumps(data)
        except (TypeError, OverflowError):
            raise ValueError("Data is not JSON-serializable")


class FASTAHandler(FormatHandler):
    """Handler for FASTA files."""
    
    def _read_impl(self, file_path: str, **kwargs) -> Dict[str, str]:
        """
        Read a FASTA file.
        
        If the file doesn't start with '>', the entire content will be treated
        as a single sequence with key 'unnamed_sequence'.
        """
        sequences = {}
        current_id = None
        current_seq = []
        
        with open(file_path, 'r') as f:
            content = f.read().strip()
            
            # Handle files that don't start with '>' (user error)
            if not content.startswith('>'):
                self.logger.warning(f"FASTA file {file_path} doesn't start with '>'. Treating as single sequence.")
                # Remove whitespace and linebreaks
                sequence = ''.join(content.split())
                return {'unnamed_sequence': sequence}
            
            # Process normal FASTA format
            for line in content.split('\n'):
                line = line.strip()
                if not line:
                    continue
                
                if line.startswith('>'):
                    # Save previous sequence if any
                    if current_id is not None:
                        sequences[current_id] = ''.join(current_seq)
                    
                    # Start new sequence
                    # Extract ID (everything up to the first whitespace)
                    parts = line[1:].strip().split(None, 1)
                    current_id = parts[0]
                    current_seq = []
                else:
                    current_seq.append(line)
        
        # Save last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
        
        return sequences
    
    def _write_impl(self, file_path: str, data: Dict[str, str], **kwargs) -> None:
        """Write sequences to a FASTA file."""
        line_width = kwargs.get('line_width', 60)
        
        with open(file_path, 'w') as f:
            for seq_id, sequence in data.items():
                f.write(f">{seq_id}\n")
                # Write sequence in chunks
                for i in range(0, len(sequence), line_width):
                    f.write(f"{sequence[i:i+line_width]}\n")
    
    def _validate(self, data: Any) -> None:
        """Validate that data is a dictionary of sequences."""
        if not isinstance(data, dict):
            raise ValueError("FASTA data must be a dictionary mapping IDs to sequences")
        
        for seq_id, sequence in data.items():
            if not isinstance(seq_id, str):
                raise ValueError("FASTA sequence IDs must be strings")
            if not isinstance(sequence, str):
                raise ValueError("FASTA sequences must be strings")


class GRNTableHandler(CSVHandler):
    """Handler for GRN tables."""
    
    def _read_impl(self, file_path: str, **kwargs) -> pd.DataFrame:
        """Read a GRN table."""
        df = super()._read_impl(file_path, **kwargs)
        # Ensure protein_id is the first column
        if 'protein_id' not in df.columns:
            raise ValueError("GRN table must have a 'protein_id' column")
        
        # Set protein_id as index if specified
        if kwargs.get('set_index', True):
            df = df.set_index('protein_id')
        
        return df
    
    def _write_impl(self, file_path: str, data: pd.DataFrame, **kwargs) -> None:
        """Write a GRN table."""
        # Make a copy to avoid modifying the original
        df_to_write = data.copy()
        
        # Reset index if it's the protein_id to include it as a column
        if df_to_write.index.name == 'protein_id':
            df_to_write = df_to_write.reset_index()
        
        # Ensure protein_id is the first column
        if 'protein_id' in df_to_write.columns and list(df_to_write.columns)[0] != 'protein_id':
            cols = ['protein_id'] + [c for c in df_to_write.columns if c != 'protein_id']
            df_to_write = df_to_write[cols]
        
        # Ensure index=False is always set for GRN tables
        kwargs['index'] = False
        
        super()._write_impl(file_path, df_to_write, **kwargs)
    
    def _validate(self, data: Any, skip_cell_validation: bool = False) -> None:
        """
        Validate GRN table format.
        
        Args:
            data: The data to validate
            skip_cell_validation: Skip validation of individual cell values
        """
        super()._validate(data)
        
        # Check that protein_id exists
        if 'protein_id' not in data.columns and data.index.name != 'protein_id':
            raise ValueError("GRN table must have a 'protein_id' column or index")
        
        # Check that other columns are GRN positions
        non_id_cols = [c for c in data.columns if c != 'protein_id']
        for col in non_id_cols:
            if not isinstance(col, str):
                raise ValueError(f"GRN column name must be a string: {col}")
        
        # Check cell format if not empty (should be <amino acid><position> OR '-')
        if not skip_cell_validation:
            import re
            pattern = re.compile(r'^[A-Z][0-9]+$')
            
            for col in non_id_cols:
                for value in data[col]:
                    if pd.notna(value) and value != '':
                        if not isinstance(value, str):
                            raise ValueError(f"GRN value must be a string: {value}")
                        if not pattern.match(value):
                            raise ValueError(f"GRN value must be in format <amino acid><position>: {value}")


class PropertyTableHandler(CSVHandler):
    """Handler for property tables."""
    
    def _read_impl(self, file_path: str, **kwargs) -> pd.DataFrame:
        """Read a property table."""
        df = super()._read_impl(file_path, **kwargs)
        
        # Check for associated metadata file
        metadata_path = None
        
        # Try to find metadata file based on path
        base_path = os.path.splitext(file_path)[0]
        potential_metadata_paths = [
            f"{base_path}_metadata.pkl",
            f"{base_path.replace('_properties', '_metadata').replace('_identity', '_metadata')}.pkl"
        ]
        
        for path in potential_metadata_paths:
            if os.path.exists(path):
                metadata_path = path
                break
        
        # If metadata file exists, attach it to dataframe as an attribute
        if metadata_path and not kwargs.get('skip_metadata', False):
            try:
                with open(metadata_path, 'rb') as f:
                    metadata = pickle.load(f)
                    df.metadata = metadata
                    self.logger.debug(f"Loaded metadata from {metadata_path}")
            except Exception as e:
                self.logger.warning(f"Error loading metadata file {metadata_path}: {e}")
        
        return df
    
    def _write_impl(self, file_path: str, data: pd.DataFrame, **kwargs) -> None:
        """Write a property table."""
        super()._write_impl(file_path, data, **kwargs)
        
        # Write metadata if provided
        metadata = kwargs.get('metadata') or getattr(data, 'metadata', None)
        if metadata and not kwargs.get('skip_metadata', False):
            metadata_path = f"{os.path.splitext(file_path)[0]}_metadata.pkl"
            try:
                with open(metadata_path, 'wb') as f:
                    pickle.dump(metadata, f, protocol=pickle.HIGHEST_PROTOCOL)
                self.logger.debug(f"Saved metadata to {metadata_path}")
            except Exception as e:
                self.logger.warning(f"Error saving metadata to {metadata_path}: {e}")
    
    def _validate(self, data: Any) -> None:
        """Validate property table format."""
        super()._validate(data)
        
        # Check that protein_id exists
        if 'protein_id' not in data.columns:
            raise ValueError("Property table must have a 'protein_id' column")


class EmbeddingHandler(PickleHandler):
    """Handler for embedding files."""
    
    def _validate(self, data: Any) -> None:
        """Validate embedding format."""
        if not isinstance(data, dict):
            raise ValueError("Embedding data must be a dictionary")
        
        # Check that values are numpy arrays
        for key, value in data.items():
            if not isinstance(value, np.ndarray):
                raise ValueError(f"Embedding for {key} must be a numpy array")
        
        # Check that all embeddings have the same dimension
        dims = set(emb.shape[1] for emb in data.values() if len(emb.shape) > 1)
        if len(dims) > 1:
            raise ValueError(f"Inconsistent embedding dimensions: {dims}")


class GraphHandler(JSONHandler):
    """Handler for graph files."""
    
    def _validate(self, data: Any) -> None:
        """Validate graph format."""
        super()._validate(data)
        
        if not isinstance(data, dict):
            raise ValueError("Graph data must be a dictionary")
        
        for protein_id, graph_data in data.items():
            if not isinstance(graph_data, dict):
                raise ValueError(f"Graph data for {protein_id} must be a dictionary")
            
            if 'nodes' not in graph_data or 'edges' not in graph_data:
                raise ValueError(f"Graph data for {protein_id} must have 'nodes' and 'edges'")
            
            # Validate nodes
            for node in graph_data['nodes']:
                if 'id' not in node:
                    raise ValueError(f"Node in {protein_id} must have an 'id'")
            
            # Validate edges
            for edge in graph_data['edges']:
                if 'source' not in edge or 'target' not in edge:
                    raise ValueError(f"Edge in {protein_id} must have 'source' and 'target'")


class FormatRegistry:
    """Registry for file format handlers."""
    
    def __init__(self):
        """Initialize the format registry."""
        # Import CifHandler here to avoid circular import
        from protos.io.cif_handler import CifHandler
        
        self.handlers = {
            'csv': CSVHandler(),
            'pkl': PickleHandler(),
            'pickle': PickleHandler(),
            'json': JSONHandler(),
            'fasta': FASTAHandler(),
            'fas': FASTAHandler(),
            'grn': GRNTableHandler(),
            'emb': EmbeddingHandler(),
            'graph': GraphHandler(),
            'cif': CifHandler(),
            'mmcif': CifHandler(),
            'pdbx': CifHandler()
        }
    
    def get_handler(self, format_type: str) -> FormatHandler:
        """
        Get a handler for a specific format.
        
        Args:
            format_type: Format identifier
            
        Returns:
            Format handler
            
        Raises:
            ValueError: If format is not supported
        """
        if format_type not in self.handlers:
            raise ValueError(f"Unsupported format: {format_type}")
        
        return self.handlers[format_type]
    
    def register_handler(self, format_type: str, handler: FormatHandler) -> None:
        """
        Register a new format handler.
        
        Args:
            format_type: Format identifier
            handler: Format handler instance
        """
        self.handlers[format_type] = handler
    
    def infer_format(self, file_path: str) -> str:
        """
        Infer format from file extension.
        
        Args:
            file_path: Path to the file
            
        Returns:
            Format identifier
            
        Raises:
            ValueError: If format cannot be inferred
        """
        _, ext = os.path.splitext(file_path)
        if not ext:
            raise ValueError(f"Cannot infer format from file: {file_path}")
        
        ext = ext[1:].lower()  # Remove leading dot
        
        # Special case handling
        if ext == 'csv':
            # Check if it's a GRN table based on filename
            basename = os.path.basename(file_path).lower()
            if 'grn' in basename:
                return 'grn'
            return 'csv'
        
        if ext in self.handlers:
            return ext
        
        raise ValueError(f"Unsupported file extension: {ext}")
    
    def read(self, file_path: str, format_type: Optional[str] = None, **kwargs) -> Any:
        """
        Read data from a file.
        
        Args:
            file_path: Path to the file
            format_type: Optional format override
            **kwargs: Format-specific parameters
            
        Returns:
            Loaded data
        """
        if not format_type:
            format_type = self.infer_format(file_path)
        
        handler = self.get_handler(format_type)
        return handler.read(file_path, **kwargs)
    
    def write(self, file_path: str, data: Any, format_type: Optional[str] = None, **kwargs) -> None:
        """
        Write data to a file.
        
        Args:
            file_path: Path to the file
            data: Data to write
            format_type: Optional format override
            **kwargs: Format-specific parameters
        """
        if not format_type:
            format_type = self.infer_format(file_path)
        
        handler = self.get_handler(format_type)
        handler.write(file_path, data, **kwargs)


# Global instance for convenience
format_registry = FormatRegistry()


class TempFileManager:
    """Manager for temporary files with tracking and cleanup."""
    
    def __init__(self, base_dir: Optional[str] = None, prefix: str = "protos_"):
        """
        Initialize the temporary file manager.
        
        Args:
            base_dir: Base directory for temporary files (uses system temp dir if None)
            prefix: Prefix for temporary filenames
        """
        self.base_dir = base_dir
        self.prefix = prefix
        self.temp_files = set()
        self.temp_dirs = set()
        self.logger = logging.getLogger(self.__class__.__name__)
    
    @contextmanager
    def temp_file(self, suffix: str = "", content: Optional[Union[str, bytes]] = None, 
                  delete: bool = True) -> ContextManager[str]:
        """
        Create a temporary file and track it for cleanup.
        
        Args:
            suffix: Filename suffix (including extension)
            content: Optional content to write to the file
            delete: Whether to delete the file when context exits
            
        Yields:
            Path to the temporary file
        """
        if self.base_dir and not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir, exist_ok=True)
        
        # Create temporary file
        fd, path = tempfile.mkstemp(suffix=suffix, prefix=self.prefix, dir=self.base_dir)
        self.temp_files.add(path)
        
        try:
            # Close the file descriptor returned by mkstemp
            os.close(fd)
            
            # Write content if provided
            if content is not None:
                mode = 'wb' if isinstance(content, bytes) else 'w'
                with open(path, mode) as f:
                    f.write(content)
            
            yield path
            
        finally:
            # Clean up
            if delete and path in self.temp_files:
                try:
                    os.unlink(path)
                    self.temp_files.remove(path)
                except Exception as e:
                    self.logger.warning(f"Failed to delete temporary file {path}: {e}")
    
    @contextmanager
    def temp_dir(self, suffix: str = "", delete: bool = True) -> ContextManager[str]:
        """
        Create a temporary directory and track it for cleanup.
        
        Args:
            suffix: Directory name suffix
            delete: Whether to delete the directory when context exits
            
        Yields:
            Path to the temporary directory
        """
        # Create temporary directory
        path = tempfile.mkdtemp(suffix=suffix, prefix=self.prefix, dir=self.base_dir)
        self.temp_dirs.add(path)
        
        try:
            yield path
            
        finally:
            # Clean up
            if delete and path in self.temp_dirs:
                try:
                    shutil.rmtree(path)
                    self.temp_dirs.remove(path)
                except Exception as e:
                    self.logger.warning(f"Failed to delete temporary directory {path}: {e}")
    
    def create_temp_file(self, suffix: str = "", content: Optional[Union[str, bytes]] = None) -> str:
        """
        Create a temporary file without automatic cleanup.
        
        Args:
            suffix: Filename suffix (including extension)
            content: Optional content to write to the file
            
        Returns:
            Path to the temporary file
        """
        if self.base_dir and not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir, exist_ok=True)
        
        # Create temporary file
        fd, path = tempfile.mkstemp(suffix=suffix, prefix=self.prefix, dir=self.base_dir)
        self.temp_files.add(path)
        
        # Close the file descriptor returned by mkstemp
        os.close(fd)
        
        # Write content if provided
        if content is not None:
            mode = 'wb' if isinstance(content, bytes) else 'w'
            with open(path, mode) as f:
                f.write(content)
        
        return path
    
    def create_temp_dir(self, suffix: str = "") -> str:
        """
        Create a temporary directory without automatic cleanup.
        
        Args:
            suffix: Directory name suffix
            
        Returns:
            Path to the temporary directory
        """
        # Create temporary directory
        path = tempfile.mkdtemp(suffix=suffix, prefix=self.prefix, dir=self.base_dir)
        self.temp_dirs.add(path)
        
        return path
    
    def cleanup(self, path: Optional[str] = None) -> bool:
        """
        Clean up a specific temporary file/directory or all if none specified.
        
        Args:
            path: Path to clean up (all if None)
            
        Returns:
            True if cleanup was successful
        """
        if path:
            # Clean up specific path
            if path in self.temp_files:
                try:
                    os.unlink(path)
                    self.temp_files.remove(path)
                    return True
                except Exception as e:
                    self.logger.warning(f"Failed to delete temporary file {path}: {e}")
                    return False
            
            elif path in self.temp_dirs:
                try:
                    shutil.rmtree(path)
                    self.temp_dirs.remove(path)
                    return True
                except Exception as e:
                    self.logger.warning(f"Failed to delete temporary directory {path}: {e}")
                    return False
            
            return False
        
        else:
            # Clean up all
            success = True
            
            # Clean up files
            for path in list(self.temp_files):
                try:
                    if os.path.exists(path):
                        os.unlink(path)
                except Exception as e:
                    self.logger.warning(f"Failed to delete temporary file {path}: {e}")
                    success = False
            
            # Clean up directories
            for path in list(self.temp_dirs):
                try:
                    if os.path.exists(path):
                        shutil.rmtree(path)
                except Exception as e:
                    self.logger.warning(f"Failed to delete temporary directory {path}: {e}")
                    success = False
            
            self.temp_files.clear()
            self.temp_dirs.clear()
            
            return success
    
    def get_temp_files(self) -> List[str]:
        """
        Get list of tracked temporary files.
        
        Returns:
            List of temporary file paths
        """
        return list(self.temp_files)
    
    def get_temp_dirs(self) -> List[str]:
        """
        Get list of tracked temporary directories.
        
        Returns:
            List of temporary directory paths
        """
        return list(self.temp_dirs)
    
    def __del__(self):
        """Cleanup on garbage collection."""
        self.cleanup()


# Global instance for convenience
temp_manager = TempFileManager()


class ExternalToolHelper:
    """Helper for working with external tools that require temp files."""
    
    def __init__(self, tool_name: str, temp_manager: Optional[TempFileManager] = None):
        """
        Initialize the external tool helper.
        
        Args:
            tool_name: Name of the external tool
            temp_manager: Temp file manager to use (uses global if None)
        """
        self.tool_name = tool_name
        self.temp_manager = temp_manager or globals()['temp_manager']
        self.logger = logging.getLogger(f"ExternalToolHelper.{tool_name}")
    
    @contextmanager
    def tool_input(self, data: Any, format_type: str, suffix: Optional[str] = None, 
                   **kwargs) -> ContextManager[str]:
        """
        Prepare input for an external tool.
        
        Args:
            data: Data to write to temp file
            format_type: Format type for the data
            suffix: File suffix (default based on format)
            **kwargs: Additional format-specific parameters
            
        Yields:
            Path to the temporary input file
        """
        # Determine suffix if not provided
        if suffix is None:
            format_map = {
                'fasta': '.fa',
                'csv': '.csv',
                'json': '.json',
                'pkl': '.pkl',
                'pickle': '.pkl',
                'grn': '.csv',
            }
            suffix = format_map.get(format_type, f".{format_type}")
        
        # Create temp file and write data
        with self.temp_manager.temp_file(suffix=suffix) as temp_path:
            format_registry.write(temp_path, data, format_type=format_type, **kwargs)
            self.logger.debug(f"Created temp input file for {self.tool_name}: {temp_path}")
            yield temp_path
    
    @contextmanager
    def tool_output(self, suffix: str = "", expected_format: Optional[str] = None) -> ContextManager[Tuple[str, Any]]:
        """
        Prepare for output from an external tool.
        
        Args:
            suffix: File suffix for the output file
            expected_format: Format type to read after tool execution
            
        Yields:
            Tuple of (output file path, read function)
        """
        # Create temp file for output
        with self.temp_manager.temp_file(suffix=suffix) as temp_path:
            # Create a read function for the caller to use
            def read_output(**kwargs):
                if expected_format:
                    return format_registry.read(temp_path, format_type=expected_format, **kwargs)
                else:
                    # Try to infer format from suffix
                    try:
                        return format_registry.read(temp_path, **kwargs)
                    except Exception as e:
                        # Fall back to reading raw content
                        with open(temp_path, 'r') as f:
                            return f.read()
            
            yield temp_path, read_output
    
    @contextmanager
    def tool_workspace(self) -> ContextManager[str]:
        """
        Create a temporary workspace directory for the tool.
        
        Yields:
            Path to the temporary directory
        """
        with self.temp_manager.temp_dir(suffix=f"_{self.tool_name}") as temp_dir:
            self.logger.debug(f"Created temp workspace for {self.tool_name}: {temp_dir}")
            yield temp_dir
    
    def create_temp_copy(self, source_path: str, suffix: Optional[str] = None) -> str:
        """
        Create a temporary copy of a file.
        
        Args:
            source_path: Path to the source file
            suffix: Optional suffix for the temp file
            
        Returns:
            Path to the temporary copy
        """
        if suffix is None:
            _, suffix = os.path.splitext(source_path)
        
        temp_path = self.temp_manager.create_temp_file(suffix=suffix)
        shutil.copy2(source_path, temp_path)
        self.logger.debug(f"Created temp copy for {self.tool_name}: {temp_path}")
        
        return temp_path