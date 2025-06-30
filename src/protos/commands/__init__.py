"""
Protos command modules for data management and operations.
"""

from .register_data import (
    register_structure_file,
    register_sequence_file,
    scan_for_unregistered_files,
    bulk_register_files,
    clean_orphaned_entries,
    create_dataset_from_entities
)

from .download_with_registration import (
    download_and_register_structure,
    bulk_download_structures,
    download_sequence_database
)

__all__ = [
    'register_structure_file',
    'register_sequence_file',
    'scan_for_unregistered_files',
    'bulk_register_files',
    'clean_orphaned_entries',
    'create_dataset_from_entities',
    'download_and_register_structure',
    'bulk_download_structures',
    'download_sequence_database'
]