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
    "ligand": "ligand"
}

# Default structure subdirectories
DEFAULT_STRUCTURE_SUBDIRS = {
    "structure_dir": "mmcif",         # Directory for structure files (raw CIF)
    "cache_dir": "cache",            # Directory for cached/processed structures
    "dataset_dir": "structure_dataset", # Directory for processed dataset PKL files
    "alignments_dir": "alignments",   # Directory for alignment files
    "temp_dir": "temp_cif",          # Directory for temporary files
    "datasets_dir": "datasets"       # Directory for dataset JSON definitions
}

# Default test subdirectories (for test compatibility)
DEFAULT_TEST_SUBDIRS = {
    "dataset_dir": "datasets"         # Directory for test dataset files
}

# Default GRN subdirectories
DEFAULT_GRN_SUBDIRS = {
    "table_dir": "tables",           # Directory for GRN mapping tables
    "ref_dir": "reference",          # Directory for reference GRN tables (changed from 'reference' to 'reference')
    "reference_dir": "reference",    # Alias for ref_dir
    "configs_dir": "configs",        # Directory for GRN configuration files
    "assignment_dir": "assignments", # Directory for GRN assignments
    "temp_dir": "temp",             # Directory for temporary files
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default sequence subdirectories
DEFAULT_SEQUENCE_SUBDIRS = {
    "fasta_dir": "fasta",           # Directory for FASTA files
    "alignment_dir": "alignments",   # Directory for sequence alignments
    "pairwise_alignment_dir": "alignments/pairwise",  # Directory for pairwise alignments
    "multiple_alignment_dir": "alignments/multiple",  # Directory for multiple alignments
    "mmseqs_alignment_dir": "alignments/mmseqs",      # Directory for MMseqs2 alignments
    "databases_dir": "databases",    # Directory for MMseqs2 databases
    "metadata_dir": "metadata",      # Directory for sequence metadata
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default property subdirectories
DEFAULT_PROPERTY_SUBDIRS = {
    "tables_dir": "tables",          # Directory for property tables
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default embedding subdirectories
DEFAULT_EMBEDDING_SUBDIRS = {
    "embeddings_dir": "embeddings",  # Directory for saved embeddings
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default ligand subdirectories
DEFAULT_LIGAND_SUBDIRS = {
    "sdf_dir": "sdf",               # Directory for SDF/MOL files
    "cache_dir": "cache",            # Directory for cached ligand data
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default graph subdirectories
DEFAULT_GRAPH_SUBDIRS = {
    "networks_dir": "networks",      # Directory for network files
    "analysis_dir": "analysis",      # Directory for graph analysis results
    "datasets_dir": "datasets"       # Directory for dataset definitions
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