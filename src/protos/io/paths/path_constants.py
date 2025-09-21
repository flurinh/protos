"""
Path constants for the Protos framework.

This module defines standard directory names and environment variables
used for path resolution throughout the Protos framework.
"""

import os
from typing import Dict, List
import importlib.resources as pkg_resources

# Environment variables for data paths
ENV_DATA_ROOT = "PROTOS_DATA_ROOT"  # For data directory override

# Default processor directories
DEFAULT_PROCESSOR_DIRS = {
    "structure": "structure",
    "grn": "grn",
    "sequence": "sequence",
    "graph": "graph",
    "property": "property",
    "embedding": "embedding",
    "ligand": "ligand",
    "input": "input",  # Added for InputManager
    "temp": "temp"     # Added for temporary files
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

# Default GRN subdirectories
DEFAULT_GRN_SUBDIRS = {
    "table_dir": "tables",           # Directory for GRN mapping tables
    "reference_dir": "reference",    # Directory for reference GRN tables
    "configs_dir": "configs",        # Directory for GRN configuration files
    "assignment_dir": "assignments", # Directory for GRN assignments
    "temp_dir": "temp",             # Directory for temporary files
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default sequence subdirectories
DEFAULT_SEQUENCE_SUBDIRS = {
    "entity_fasta_dir": "fasta/entities",        # Single-sequence FASTA files
    "dataset_fasta_dir": "fasta/datasets",       # Multi-sequence FASTA archives
    "alignment_dir": "alignments",               # Alignment outputs
    "pairwise_alignment_dir": "alignments/pairwise",
    "multiple_alignment_dir": "alignments/multiple",
    "mmseqs_alignment_dir": "alignments/mmseqs",
    "database_dir": "databases",                 # Future sequence databases
    "metadata_dir": "metadata",
    "datasets_dir": "datasets"
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
    "graphs_dir": "graphs",          # Directory for persisted graph entities
    "analysis_dir": "analysis",      # Directory for graph analysis results
    "datasets_dir": "datasets"       # Directory for dataset definitions
}

# Default input subdirectories (for InputManager)
DEFAULT_INPUT_SUBDIRS = {
    # Flat input directory - no subdirectories needed
    # Files are processed and moved to appropriate processor directories
}

# Default temp subdirectories
DEFAULT_TEMP_SUBDIRS = {
    # No subdirectories for temp - it's just a flat directory
}

# Registry file names
# DEPRECATED: DEFAULT_REGISTRY_FILENAME removed - processor-specific registries no longer used
DEFAULT_GLOBAL_REGISTRY_FILENAME = "global_registry.json"  # Used by unified EntityRegistry

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
