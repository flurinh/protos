"""
Base loader class for the protos framework.

This module provides a base class for all data loaders, handling:
- Data fetching from external sources (PDB, AlphaFold, etc.)
- Automatic entity registration
- Dataset creation for bulk downloads
- Integration with ProtosPaths, EntityRegistry, and DatasetManager

Key principles:
- Loaders handle data acquisition, processors handle data manipulation
- Zero configuration - works out of the box
- Human-readable names for all operations
- Automatic registration with entity tracking
"""

import logging
import os
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Union

# Import core components
from protos.io.paths import get_protos_paths


class BaseLoader(ABC):
    """
    Base class for all data loaders in the protos framework.
    
    Key responsibilities:
    - Fetch data from external sources
    - Register entities automatically
    - Create datasets for bulk operations
    - Handle path management via ProtosPaths
    
    Loaders are separate from processors:
    - Loaders: Fetch and register data
    - Processors: Manipulate and analyze data
    """
    
    # Subclasses should override this
    loader_type: Optional[str] = None
    
    def __init__(self, name: str = "loader"):
        """
        Initialize loader with automatic setup - NO PATH PARAMETERS!
        
        Args:
            name: Loader instance name (default: "loader")
            
        Example:
            >>> loader = StructureLoader()  # Zero configuration!
        """
        self.name = name
        
        # Get singletons - initialization happens automatically
        self.paths = get_protos_paths()
        
        # Import here to avoid circular import
        from protos.io.core import get_registry
        self.entity_registry = get_registry()
        
        # Determine loader type
        self.loader_type = self._get_loader_type()
        
        # Get processor path for this data type
        self.data_path = Path(self.paths.get_processor_path(self.loader_type))
        
        # Create dataset manager for this data type
        from protos.io.core.dataset_manager import DatasetManager
        self.dataset_manager = DatasetManager(self.loader_type)
        
        # Set up logging
        self.logger = self._setup_logger()
        self.logger.info(f"Initialized {self.__class__.__name__} '{name}' for {self.loader_type} data")
    
    def _setup_logger(self) -> logging.Logger:
        """Set up a logger for this loader instance."""
        logger = logging.getLogger(f"{self.__class__.__name__}.{self.name}")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter(
                '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
            )
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger
    
    def _get_loader_type(self) -> str:
        """
        Get loader type from class attribute or class name.
        
        Returns:
            Loader type string (e.g., 'structure', 'sequence')
        """
        # First check if subclass defined loader_type
        if self.loader_type is not None:
            return self.loader_type
            
        # Otherwise, infer from class name
        class_name = self.__class__.__name__
        
        # Handle common patterns
        if 'Structure' in class_name:
            return 'structure'
        elif 'Sequence' in class_name:
            return 'sequence'
        elif 'Property' in class_name:
            return 'property'
        elif 'Ligand' in class_name:
            return 'ligand'
        else:
            # Default to lowercase class name without 'Loader'
            return class_name.replace('Loader', '').replace('Base', '').lower()
    
    # ========== Core Loader Methods (must be implemented by subclasses) ==========
    
    @abstractmethod
    def fetch_entity(self, identifier: str, **kwargs) -> Optional[Path]:
        """
        Fetch a single entity from external source.
        
        Args:
            identifier: Entity identifier (PDB ID, UniProt ID, etc.)
            **kwargs: Additional source-specific parameters
            
        Returns:
            Path to downloaded file, or None if failed
            
        Example:
            >>> path = loader.fetch_entity("1ubq")
            >>> path = loader.fetch_entity("P00533", source="alphafold")
        """
        pass
    
    @abstractmethod
    def parse_identifier(self, identifier: str) -> Dict[str, Any]:
        """
        Parse and validate an identifier for this data type.
        
        Args:
            identifier: Raw identifier string
            
        Returns:
            Dict with parsed components (id, source, etc.)
            
        Example:
            >>> info = loader.parse_identifier("1ubq")
            >>> info = loader.parse_identifier("AF-P00533-F1")
        """
        pass
    
    # ========== Entity Registration Methods ==========
    
    def download_and_register(self, identifier: str, 
                            name: Optional[str] = None,
                            metadata: Optional[Dict[str, Any]] = None,
                            **kwargs) -> Optional[str]:
        """
        Download an entity and register it with the entity registry.
        
        Args:
            identifier: Entity identifier (PDB ID, UniProt ID, etc.)
            name: Human-readable name for registration (defaults to identifier)
            metadata: Additional metadata to store
            **kwargs: Additional source-specific parameters
            
        Returns:
            Registered entity name if successful, None if failed
            
        Example:
            >>> name = loader.download_and_register("1ubq")
            >>> name = loader.download_and_register("P00533", name="EGFR_HUMAN")
        """
        # Use identifier as name if not provided
        if name is None:
            name = identifier
            
        # Check if already registered
        if self.entity_registry.entity_exists(name, self.loader_type):
            self.logger.info(f"Entity '{name}' already registered for {self.loader_type}")
            return name
        
        try:
            # Fetch the entity
            file_path = self.fetch_entity(identifier, **kwargs)
            if file_path is None:
                self.logger.error(f"Failed to fetch entity: {identifier}")
                return None
            
            # Prepare metadata
            entity_metadata = metadata or {}
            entity_metadata.update({
                'source_id': identifier,
                'download_date': datetime.now().isoformat(),
                'loader': self.__class__.__name__
            })
            
            # Add parsed identifier info
            id_info = self.parse_identifier(identifier)
            entity_metadata.update(id_info)
            
            # Register with entity registry
            registered_name = self.entity_registry.register_entity(
                name=name,
                format_type=self.loader_type,
                file_path=str(file_path),
                metadata=entity_metadata
            )
            
            self.logger.info(f"Successfully registered entity '{registered_name}' from {identifier}")
            return registered_name
            
        except Exception as e:
            self.logger.error(f"Failed to download and register {identifier}: {e}")
            return None
    
    # ========== Bulk Operations ==========
    
    def download_batch(self, identifiers: List[str], 
                      dataset_name: Optional[str] = None,
                      create_dataset: bool = True,
                      **kwargs) -> Tuple[List[str], List[str]]:
        """
        Download multiple entities and optionally create a dataset.
        
        Args:
            identifiers: List of entity identifiers
            dataset_name: Name for the dataset (optional)
            create_dataset: Whether to create a dataset from successful downloads
            **kwargs: Additional source-specific parameters
            
        Returns:
            Tuple of (successful entity names, failed identifiers)
            
        Example:
            >>> success, failed = loader.download_batch(
            ...     ["1ubq", "2w9s", "3sn6"],
            ...     dataset_name="gpcr_structures"
            ... )
        """
        successful = []
        failed = []
        
        # Download each entity
        for identifier in identifiers:
            name = self.download_and_register(identifier, **kwargs)
            if name:
                successful.append(name)
            else:
                failed.append(identifier)
        
        # Create dataset if requested and we have successes
        if create_dataset and successful and dataset_name:
            try:
                self.dataset_manager.create_dataset(
                    name=dataset_name,
                    entities=successful,
                    metadata={
                        'created_by': self.__class__.__name__,
                        'source_identifiers': identifiers,
                        'download_date': datetime.now().isoformat(),
                        'success_count': len(successful),
                        'fail_count': len(failed)
                    }
                )
                self.logger.info(f"Created dataset '{dataset_name}' with {len(successful)} entities")
            except Exception as e:
                self.logger.warning(f"Failed to create dataset '{dataset_name}': {e}")
        
        # Log summary
        if failed:
            self.logger.warning(
                f"Downloaded {len(successful)}/{len(identifiers)} entities. "
                f"Failed: {failed}"
            )
        else:
            self.logger.info(f"Successfully downloaded all {len(successful)} entities")
        
        return successful, failed
    
    # ========== Helper Methods ==========
    
    def get_download_path(self, identifier: str, subdirectory: Optional[str] = None) -> Path:
        """
        Get the path where a downloaded file should be stored.
        
        Args:
            identifier: Entity identifier
            subdirectory: Optional subdirectory within the data path
            
        Returns:
            Path object for the file location
        """
        if subdirectory:
            base_path = self.data_path / subdirectory
        else:
            base_path = self.data_path
            
        base_path.mkdir(parents=True, exist_ok=True)
        return base_path
    
    def list_sources(self) -> List[str]:
        """
        List available data sources for this loader.
        
        Returns:
            List of source names
            
        Example:
            >>> loader.list_sources()
            ['rcsb', 'alphafold']
        """
        # Default implementation - subclasses can override
        return ['default']
    
    def validate_identifier(self, identifier: str) -> bool:
        """
        Check if an identifier is valid for this loader.
        
        Args:
            identifier: Entity identifier to validate
            
        Returns:
            True if valid, False otherwise
        """
        try:
            self.parse_identifier(identifier)
            return True
        except (ValueError, KeyError):
            return False
    
    def check_entity_exists(self, identifier: str, name: Optional[str] = None) -> bool:
        """
        Check if an entity is already registered.
        
        Args:
            identifier: Entity identifier
            name: Entity name (defaults to identifier)
            
        Returns:
            True if already registered, False otherwise
        """
        if name is None:
            name = identifier
        return self.entity_registry.entity_exists(name, self.loader_type)