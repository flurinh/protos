"""
Centralized path management for Protos.

This module provides standardized path handling and resolution for all components
of the Protos framework, ensuring consistent data organization and access.
It handles both reference data (distributed with the package) and user data
(created at runtime).
"""

# Flag to indicate that the path module is available
_HAS_PATH_MODULE = True

# Explicitly make these values available for import
# This is needed for the BaseProcessor class to detect the path module

from .path_constants import (
    DEFAULT_USER_DATA_ROOT,
    DEFAULT_REF_DATA_ROOT,
    DEFAULT_PROCESSOR_DIRS,
    DEFAULT_STRUCTURE_SUBDIRS,
    DEFAULT_GRN_SUBDIRS,
    DEFAULT_SEQUENCE_SUBDIRS,
    ENV_DATA_ROOT,
    ENV_REF_DATA_ROOT,
    DEFAULT_REGISTRY_FILENAME,
    DEFAULT_GLOBAL_REGISTRY_FILENAME
)

from .path_config import (
    ProtosPaths,
    get_data_root,
    get_user_data_root,
    get_reference_data_root,
    ensure_directory,
    is_package_resource,
    DataSource,
    resolve_path,
    get_structure_path,
    get_grn_path,
    get_sequence_path,
    get_dataset_path
)

# Define the __all__ list to properly expose all public elements
# This is separate from the `from .module import x` imports
# We need to explicitly expose the _HAS_PATH_MODULE flag 
# for the BaseProcessor class to detect it
__all__ = [
    # Flags
    '_HAS_PATH_MODULE',
    
    # Constants
    'DEFAULT_USER_DATA_ROOT',
    'DEFAULT_REF_DATA_ROOT',
    'DEFAULT_PROCESSOR_DIRS',
    'DEFAULT_STRUCTURE_SUBDIRS',
    'DEFAULT_GRN_SUBDIRS',
    'DEFAULT_SEQUENCE_SUBDIRS',
    'ENV_DATA_ROOT',
    'ENV_REF_DATA_ROOT',
    'DEFAULT_REGISTRY_FILENAME',
    'DEFAULT_GLOBAL_REGISTRY_FILENAME',
    
    # Path Configuration
    'ProtosPaths',
    'get_data_root',
    'get_user_data_root',
    'get_reference_data_root',
    'ensure_directory',
    'is_package_resource',
    'DataSource',
    
    # Path Resolution
    'resolve_path',
    'get_structure_path',
    'get_grn_path',
    'get_sequence_path',
    'get_dataset_path'
]