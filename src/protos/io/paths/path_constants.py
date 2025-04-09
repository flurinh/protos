"""
Path constants for the Protos framework.

This module defines standard directory names and environment variables
used for path resolution throughout the Protos framework.
"""

import os
from typing import Dict, List
import importlib.resources as pkg_resources

# Environment variables for data paths
ENV_DATA_ROOT = "PROTOS_DATA_ROOT"  # For user data
ENV_REF_DATA_ROOT = "PROTOS_REF_DATA_ROOT"  # For reference data (optional override)

# Default data directories
DEFAULT_USER_DATA_ROOT = "data"  # Default location for user data (in working directory)

# Determine default reference data path within the package
try:
    # First try using importlib.resources to get package directory
    with pkg_resources.as_file(pkg_resources.files("protos") / "data") as ref_path:
        DEFAULT_REF_DATA_ROOT = str(ref_path)
except (ImportError, AttributeError, NotADirectoryError):
    # Fallback to environment variable or a relative path
    DEFAULT_REF_DATA_ROOT = os.environ.get(ENV_REF_DATA_ROOT, os.path.join(os.path.dirname(__file__), "../../../reference_data"))

# Default processor directories
DEFAULT_PROCESSOR_DIRS = {
    "structure": "structure",
    "grn": "grn",
    "sequence": "sequence",
    "graph": "graph",
    "property": "property",
    "embedding": "embedding",
    # Add test processors for test compatibility
    "test": "test",
    "test_processor": "test",
    "simple": "simple",
    "complex_processor_with_long_name": "complex_processor_with_long_name"
}

# Default structure subdirectories
DEFAULT_STRUCTURE_SUBDIRS = {
    "structure_dir": "mmcif",         # Directory for structure files
    "dataset_dir": "structure_dataset", # Directory for dataset files
    "alignments_dir": "alignments",   # Directory for alignment files
    "temp_dir": "temp_cif"           # Directory for temporary files
}

# Default test subdirectories (for test compatibility)
DEFAULT_TEST_SUBDIRS = {
    "dataset_dir": "datasets"         # Directory for test dataset files
}

# Default GRN subdirectories
DEFAULT_GRN_SUBDIRS = {
    "table_dir": "tables",           # Directory for GRN mapping tables (legacy)
    "grn_dir": "grn",                # Main GRN directory (current)
    "configs_dir": "configs",        # Directory for GRN configuration files
    "assignment_dir": "assignments"  # Directory for GRN assignments
}

# Default sequence subdirectories
DEFAULT_SEQUENCE_SUBDIRS = {
    "fasta_dir": "fasta",           # Directory for FASTA files
    "alignment_dir": "alignments",   # Directory for sequence alignments
    "metadata_dir": "metadata"      # Directory for sequence metadata
}

# Registry file names
DEFAULT_REGISTRY_FILENAME = "registry.json"
DEFAULT_GLOBAL_REGISTRY_FILENAME = "global_registry.json"

# Common path joining function that handles different OS path separators
def join_path(*args) -> str:
    """
    Join path components in a cross-platform way.
    
    Args:
        *args: Path components to join
        
    Returns:
        Joined path as a string
    """
    return os.path.normpath(os.path.join(*args))

# Standard file extensions
FILE_EXTENSIONS = {
    "structure": ".cif",
    "fasta": ".fasta",
    "dataset": ".json",
    "alignment": ".pkl",
    "table": ".csv",
    "embedding": ".npy"
}